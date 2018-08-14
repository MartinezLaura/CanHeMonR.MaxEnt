#' @title Treatment of tiles
#' @description Check if you want to read all the tiles from a folder or just the files needed via csv.
#' @param r_train_dir Path where the tiles are located
#' @param tile Or 'ALL' for reading all the files, or the path where the csv is located with the name of the rasters
#' @param text_train_dir Path where the textures are located
#' @param text Or 'ALL' for reading all the files, or the path where the csv is located with 
#' the sufix of the textures wanted
#' @return The list of files to process
#' @export
Tiles <- function(r_train_dir, tile, text_train_dir, text) {
  
  #if you want to run all the tiles in a directory, harvest the available tilenames
  if (tile[1] == 'ALL'){
    #harvest all the tif files in the directories holding covariate/predictor images
    all_tifs <- list.files(r_train_dir, recursive = F, full.names = T, pattern = ".tif", ignore.case = TRUE)
    all_tifs <- unique(all_tifs)
    all_tifs <- all_tifs[grepl('.tif',all_tifs)]
    #excluded the tif files in the unprojected folder
    all_tifs <- all_tifs[!grepl('orig_noPRJ', all_tifs)]
    #excluded the tif files ending like .aux.xml
    all_tifs <- all_tifs[!grepl('.aux.xml', all_tifs)]
    tile <- substr(basename(all_tifs),1,16)
    tile <- unique(tile)
    #only keep tiles that start  with 'pt'
    tile <- tile[substr(tile,1,2) == 'pt']
    cat(length(tile),' tiles are considered\n')
  }else{
    tifs <- read.csv(tile, header = FALSE)
    tile <- unlist(lapply(tifs$V1, function(x) paste0(r_train_dir,x)))
    tile <- tile[!grepl('.aux.xml', tile)]
    all_tifs <- unlist(tile)
    tile <- substr(basename(tile),1,16)
    tile <- unique(tile)
    tile = tile[substr(tile,1,2) == 'pt']
    cat(length(tile),' tiles are considered\n')
  }
  if (text[1] == 'ALL'){
    all_tifs <- append(all_tifs, list.files(text_train_dir, recursive = T, full.names = T, pattern = ".tif", ignore.case = TRUE))
    all_tifs <- unique(all_tifs)
  }else{
    texts <- read.csv(text, header = FALSE)
    texts <- lapply(texts$V1, function(x) list.files(text_train_dir, recursive = F, full.names = T, pattern = as.character(x)))
    texts <- texts[!grepl('.aux.xml', texts)]
    all_text <- unlist(texts)
    all_tifs <- append(all_tifs, all_text)
    all_tifs <- unique(all_tifs)
  }
  output<-list(all_tifs, tile)
  return(output)
}

#' @title Parallel enviroment
#' @description Set up the parallel enviroment if needed.
#' @param nWorkers Number of threads you want to use the maximum is Maxcores-1
#' @param data_outp_dir Path where you want to save the logfile
#' @param data_outp_name Name of the log file
#' @param parallel Boolean. True for parallel performance
#' @return The registration of parallel flag
#' @export
Par <- function(nWorkers, data_outp_dir, data_outp_name, parallel) {
  require(foreach)
  if (parallel){
    try(parallel::stopCluster(cl), silent=T)
    # Test that the cores assign not overlap the maximum cores available
    maxcl <- (parallel::detectCores(logical = TRUE)-1)
    if (nWorkers <= maxcl){
      cl <- parallel::makeCluster(nWorkers, outfile = paste0(data_outp_dir,"log-",data_outp_name,".txt"))
      doParallel::registerDoParallel(cl)
    }else{
      cl <- parallel::makeCluster(maxcl, outfile = paste0(data_outp_dir,"log-",data_outp_name,".txt"))
      doParallel::registerDoParallel(cl)
    }
  }else{
    cl <- NULL
  }
  
  #choose the appropriate operator for the foreach loop
  `%op%` <- if (parallel) `%dopar%` else `%do%`
  out <- list(`%op%`, cl)
  return(out)
}

#' @title Sample points of presence points and background points and save them in a serialized file
#' @description Sample points of presence points and background points and save them in a serialized file.
#' @param prob_tifs Boolean (FALSE) if not wanted, Directory if wanted. tifs of each raster to run with 0 vaue for the areas that we don't want to sample 
#' and 1 fore the ones that we want to sample. Default FALSE
#' @param tile_i Passed by the algorithm, name of the tile to run
#' @param all_tifs Passed by the algorithm, list of all the tiles to run
#' @param field_name Character. The field in AOI.filename that contains the vuln_classes
#' @param vuln_classes A list of the classes you want to model
#' The list can contain one or more vectors characters. Each element of the vector represents a seperate vegetation class and response variable for the model
#' and the vector elements are synonyms used to describe that class
#' The fist place in each vector will be used in the output name used to store the calibrated model, so it should not contain spaces.
#' The other places should appear as attributes in the field 'field_name' of Pols
#' @param Pols Passed by the algorithm, PolygonDataframe Object with all the visual points assesed
#' @param ninputs_tile Integer. Number of inputs that we have fore each tile, including the tile, for exemple number of textures
#' @param abs_samp Integer. How many 'absence' pixels should be randomly selected from eah tile to train the model. Default is 100.
#' @param randompt Boolean. if True random points will be added in the tiles that doesn't have any visual point. Default TRUE
#' @param tile_dat A dataframe with the points of the prevoius iteration, if the previous iteration is the initial one, an empy dataframe will be passed
#' @return A dataframe with the points done for the tile selected
#' @export
GetPoints <- function(tile_i, all_tifs,field_name, ninputs_tile, randompt, prob_tifs, Pols, vuln_classes, abs_samp, tile_dat){

  #grep all the textures needes for this tile
  pred_rs <- all_tifs[grepl(tile_i, all_tifs)]
  
  if (length(pred_rs) == ninputs_tile){
    #create the brick/stack of predictor layers for this tile
    pred_rs = normalizePath(pred_rs)
    r_train <- raster::stack(pred_rs)
    tile_i_multiversion_for_regexpr <- gsub("-", ".", tile_i)
    #adjust the name so the tile-specific label is removed, and names are consistent between tiles
    names(r_train) <- paste0('l',unlist(lapply(strsplit(names(r_train),tile_i_multiversion_for_regexpr,fixed=F),function(x){x[-1]})))
    #check if you have any points in this tile
    #crop the calval to this tile
    Pols_tile <- raster::crop(Pols, raster::raster(pred_rs[1]))
    
    if (!is.null(Pols_tile)){
      #only proceed if you have training points in this tile
      #pres_train_tile <- NULL
      # extract the data for this tile for each class
      for (j in 1:length(vuln_classes)){
        pres_dat <- data.frame()
        
        if (any(is.element(Pols_tile@data[[field_name]] , vuln_classes[[j]]))){
          #pres_train <- NULL
          class. <- vuln_classes[[j]]
          cat('Sampling data for class ',class.[1],'\n')
          # the sampling for points:
          pres_train_tile <- Pols_tile[is.element(Pols_tile@data[[field_name]] , class.),]
          #get covariates for presence locations
          pres_dat <- raster::extract(r_train, pres_train_tile, df = TRUE)
          #erase the ID column if exist
          pres_dat$ID <- NULL
          
          if(vuln_classes[[j]] == 'false_pos'){
            pres_dat <- cbind.data.frame(pres = c(rep(0,nrow(pres_dat))),pres_dat)
            pres_dat$class <- rep('Pb',nrow(pres_dat))
            #add the data for this class and tile, to the data for this tile
            tile_dat <- rbind.data.frame(tile_dat, pres_dat)
          }
          else{
            pres_dat <- cbind.data.frame(pres = c(rep(1,nrow(pres_dat))),pres_dat)
            pres_dat$class <- rep('Pb',nrow(pres_dat))
            #add the data for this class and tile, to the data for this tile
            tile_dat <- rbind.data.frame(tile_dat, pres_dat)
          }
        }
        else{
          cat('There are no points falling in this tile:', tile_i,'for the class:', vuln_classes[[j]], '\n')
        }
      }
    }
    if(randompt == TRUE && prob_tifs == FALSE){
      #get covariates for randomly sampled absence locations
      abs_dat <- data.frame()
      
      if (exists('pres_train_tile')){
        cat('Detecting absences for tiles with visual points\n')
        abs_loc <- dismo::randomPoints(r_train, n = abs_samp, p = pres_train_tile, warn=0)
        #exclude pseude-absences that fall too close (< 20 m) to presence locations
        dist_abs2pres <- sp::spDists(abs_loc, sp::coordinates(Pols_tile))
        mindist_abs2pres <- apply(dist_abs2pres, 1, min)
        abs_loc <- abs_loc[mindist_abs2pres > 20,]
      }else{
        cat('Detecting absences for tiles without visual points\n')
        abs_loc <- dismo::randomPoints(r_train, n = abs_samp, warn=0 )
      }
      abs_dat <- data.frame(raster::extract(r_train, abs_loc))
      abs_dat <- stats::na.omit(abs_dat)
      
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
      abs_dat <- cbind.data.frame(pres = c(rep(0,nrow(abs_dat))),abs_dat)
      abs_dat$class <- rep('Pb', nrow(abs_dat) )
      
      #add the data for this class and tile, to the data for this tile
      tile_dat <- rbind.data.frame(tile_dat, abs_dat)
    }
    if(prob_tifs == TRUE){
      #Pick the raster mask belonging to this tile with prob = 0 is NA and the others are probs
      prob_tif <-  prob_tifs[grepl(paste0(tile_i,'.tif'), prob_tifs)]
      prob_tif <- prob_tif[!grepl('.aux.xml', prob_tif)]
      prob_rast <- raster::raster(prob_tif)
      print(sum(raster::values(prob_rast), na.rm = TRUE))
      
      if (sum(raster::values(prob_rast), na.rm = TRUE) > 0){
        prob_rast[prob_rast<1] <- NA
        #get covariates for randomly sampled absence locations
        abs_dat <- data.frame()
        
        if (exists('pres_train_tile')){
          cat('Detecting absences for tiles with visual points\n')
          abs_loc <- dismo::randomPoints(r_train, n = abs_samp, p = pres_train_tile, warn=0, mask=prob_rast)
          #exclude pseude-absences that fall too close (< 20 m) to presence locations
          dist_abs2pres <- sp::spDists(abs_loc, sp::coordinates(Pols_tile))
          mindist_abs2pres <- apply(dist_abs2pres, 1, min)
          abs_loc <- abs_loc[mindist_abs2pres > 20,]
          
        }else{
          cat('Detecting absences for tiles without visual points\n')
          abs_loc <- dismo::randomPoints(r_train, n = abs_samp, warn=0, mask=prob_rast)
        }
        
        abs_dat <- data.frame(raster::extract(r_train, abs_loc))
        abs_dat <- stats::na.omit(abs_dat)
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
        abs_dat <- cbind.data.frame(pres = c(rep(0,nrow(abs_dat))),abs_dat)
        abs_dat$class <- rep('Pb', nrow(abs_dat) )
        
        #add the data for this class and tile, to the data for this tile
        tile_dat <- rbind.data.frame(tile_dat, abs_dat)
      }else{
        cat("We wont study this tile, doesnt have any ROI", tile_i,"\n")}
    }
  }else{
    cat("Not enoght layers for tile:", tile_i,"\n")
  }
  
  return(tile_dat)
}