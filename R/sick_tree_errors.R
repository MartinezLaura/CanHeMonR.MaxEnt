#' @title Sick tree error calculator
#' @description Calculate errors in automated detection of declining trees using visual inspection data as reference
#' @param r_pred_dir A directory where binary .tifs predicting presence as 1 and absence as 0 can be found for multiple tiles
#' @param tile Character vector. Names of tile(s) to run. 'ALL will run all tiles in r_pred_dir. Default is 'ALL'
#' @param thresh If the image data is non-binary, the value threhs can be set to split the
#' image values between presence (> thresh) and absence (<= thresh). Default is NA, in which case the image values are assumed binary.
#' @param vuln_classes A list of the classes you want to model.
#' The list can contain one or more vectors. Each vector represents a seperate vegetation class and response variable for the model
#' and the vector elements are synonyms used to describe that class.
#' The fist place in each vector will be used in the output name used to store the calibrated model, so it should not contain spaces.
#' The other places should appear as attributes in the field 'field_name' of pnts.
#' @param pnts SpatialPointsDataFrame of which one field contains the vuln_classes
#' @param radius The radius within which a presence point must be found for it to be considered 'correct'
#' @param field_name The field in pnts that contains the vuln_classes
#' @param abs_samp How many 'absence' pixels should be randomly selected from each tile to evaluate the absences? Default is 100.
#' @param parallel Should the code be run in parallel using the doParallel package? Default is FALSE.
#' @param nWorkers If running the ocde in parallel, how many workers should be used? Default is 4.
#' @param data_outp_dir The folder and filename prefix to save the sampled data to. No data is saved is data_outp_dir is NULL. Default is NULL.
#' @return A data frame with commission and ommission errors and sample sizes of presence and absence
#' @examples \dontrun{
#' }
#' @export
sick_tree_errors <- function(r_pred_dir, tile = 'ALL', thresh = NA, vuln_classes, pnts, radius, field_name, abs_samp = 100,
                                                   parallel = F, nWorkers = 4, data_outp_dir = NULL){



  if(is.factor( pnts@data[[field_name]])){
    pnts@data[[field_name]] <- droplevels(pnts@data[[field_name]])
  }

  #harvest all the tif files in the directories holding covariate/predictor images
  all_tifs <- list.files(r_pred_dir, recursive = F, full.names = T)
  all_tifs <- all_tifs[grepl('.tif',all_tifs)]
  #excluded the tif files in the unprojected folder
  all_tifs <- all_tifs[!grepl('raw', all_tifs)]
  all_tifs <- all_tifs[!grepl('.aux.xml', all_tifs)]

  #if you want to run all the tiles in a directory, harvest the available tilenames
  if (tile[1] == 'ALL'){
    #tile <- unlist(strsplit(gsub(pattern = "\\.tif$", "", basename(all_tifs)),"_"))
    tile <- unlist(gsub(pattern = "\\.tif$", "", basename(all_tifs)))
    # tile <- tile[1:length(tile) %% 3 == 0]
    tile <- unique(tile)
    # tile <- strsplit(tile,'-')
    # tile <- tile[1:length(tile) %% 2 == 0]
    # #only keep tiles that start  with 'pt'
    # tile <- tile[substr(tile,9,10) == 'pt']

    cat(length(tile),' tiles are considered\n')
  }

  tile_counter <- 0

  # a list to hold the outputs
  calval_dfs <- list()


  #set up the cluster for parallel processing
  if (parallel){
    try(parallel::stopCluster(cl), silent=T)
    # TO DO add a line that avoids allocating more workers than you have cores
    cl <- parallel::makeCluster(nWorkers)
    doParallel::registerDoParallel(cl)
  }
  #choose the appropriate operator for the foreach loop
  require(foreach)
  `%op%` <- if (parallel) `%dopar%` else `%do%`

  stime <- system.time({
    calval_dfs <- foreach::foreach(i = 1:length(tile), .combine = rbind.data.frame, .inorder=F, .multicombine=F, .errorhandling='remove') %op% {
      tile_i <- tile[i]
      pred_rs <- all_tifs[i]

      #an empty data frame to hold the data extracted for this tile
      tile_dat <- data.frame()
      ## you should only have one output tif per tile
      #check if you have any points in this tile
      #crop the calval to this tile
      pnts_tile <- raster::crop(pnts, raster::raster(pred_rs))

      #only proceed if you have training points in this tile
      if (!is.null(pnts_tile)){
        if (length(pnts_tile) >= 1 & any(pnts_tile@data$type == vuln_classes)){

          #read the predicted layers for this tile
          r_pred <- raster::raster(pred_rs)

          #if the image has not been made binary yet, do so
          if (!is.na(thresh)){
            r_pred <- r_pred > thresh
          }

          cat('Sampling data from ', basename( tile_i),'\n')

          #reproject the trainig pnts if necessary
          if (raster::projection(pnts) != raster::projection(r_pred)){
            pnts <- sp::spTransform(pnts, sp::CRS(raster::projection(r_pred)))
          }


          # extract the data for this tile for each class
          for (j in 1:length(vuln_classes)){
            print ('entro!')
            pres_train <- NULL

            class. <- vuln_classes[[j]]
            cat('sampling data for class ',class.[1],'\n')
            if(length(class.) > 1){
              cat('which also includes ',class.[-1],'\n')
            }

            # the sampling for reference points:
            pres_vis_tile <- pnts_tile[is.element(pnts_tile@data[[field_name]] , class.),]

            cat('For the class', class.[1],' this tile has ',length(pres_vis_tile),' presence points falling in it.\n')

            #get the predicted presence/absence for the reference presence points
            #extractions are done by counting the presences in the image within a certain radius of the point location
            if (length(pres_vis_tile) > 0){
              pres_dat <- data.frame(raster::extract(r_pred, pres_vis_tile, buffer = radius, fun = sum), pres_vis_tile$Id)
              colnames(pres_dat) <- c('pred','point_ID')
            }else{
              pres_dat <- data.frame()
            }
            #pres_dat[['obs']] <- 1

            #randomly sample pseudo-absence locations
            abs_loc <- dismo::randomPoints( r_pred, n = abs_samp, p = pres_vis_tile, warn=0 )

            #exclude pseude-absences that fall too close to presence locations
            dist_abs2pres <- sp::spDists(abs_loc, sp::coordinates(pnts_tile))
            mindist_abs2pres <- apply(dist_abs2pres, 1, min)
            abs_loc <- abs_loc[mindist_abs2pres > 2*radius,]

            #extract predicted values at the pseudo-absences
            #extractions are done by counting the presences in the image within a certain radius of the point location
            abs_dat <- data.frame(raster::extract(r_pred, abs_loc, buffer = radius, fun = sum),0)
            abs_dat <- stats::na.omit(abs_dat)
            colnames(abs_dat) <- c('pred','point_ID')
            #abs_dat[['obs']] <- 0


            if (nrow(abs_dat) == 0) {
              stop('could not get valid background point values; is there a layer with only NA values?')
            }
            if (nrow(abs_dat) < abs_samp/100) {
              stop('only got:', nrow(abs_dat), 'random background point values; is there a layer with many NA values?')
            }
            if (nrow(abs_dat) < abs_samp/10) {
              warning('only got:', nrow(abs_dat), 'random background point values; Small exent? Or is there a layer with many NA values?')
            }

            #join presence and absence data
            tile_dat_class <- rbind.data.frame(pres_dat, abs_dat)
            tile_dat_class <- cbind.data.frame(vis = c(rep(1,nrow(pres_dat)),rep(0,nrow(abs_dat))),tile_dat_class)
            #add the classname
            tile_dat_class$class <- rep(class., nrow(tile_dat_class) )
            #add the tilename
            tile_dat_class$tile <- tile_i

            #add the data for this class and tile, to the data for this tile
            tile_dat <- rbind.data.frame(tile_dat, tile_dat_class)
            #require(dismo)
            }
          }
        }

      tile_dat
    }
  })[3]
  cat("------------------------------------------\n")


  #+++++++++++++++++++++++++++++++++++++++++++++++
  # report performance statistics ----
  #+++++++++++++++++++++++++++++++++++++++++++++++

  if (parallel){
    cat('using \n',foreach::getDoParWorkers(),' parallel workers,\n')
  }else{
    cat('processing sequentially on a single worker \n')
  }
  cat('Estimated performance for ',length(tile),' tiles in ',round(stime/60),' minutes\n')

  #############################################
  # close the cluster set up forparallel processing
  if (parallel){
    parallel::stopCluster(cl)
  }


  #+++++++++++++++++++++++++++++++++++++++++++++++
  # save the extracted data ----
  #+++++++++++++++++++++++++++++++++++++++++++++++

  if (!is.null(data_outp_dir)){

    data_file <- paste0(data_outp_dir,'sicktree_performance_dfs.rdsdata')
    saveRDS(calval_dfs, file = data_file)
    cat('Wrote away ', data_file,'\n')
  }


  return(calval_dfs)
}
