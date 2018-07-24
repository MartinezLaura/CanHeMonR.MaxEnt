#' @title Sample training data for image classification from multiple image tiles
#' @description For each class in .shp polygon file, Sample training data for image classification from multiple image tiles using their raster bricks as predictors
#' @param r_train_dir A directory where .tifs for training can be found for multiple tiles
#' @param text_train_dir A directory where .tifs of the textures associated with r_train_dir
#' @param prob_tifs Boolean (FALSE) if not wanted, Directory if wanted. tifs of each raster to run with 0 vaue for the areas that we don't want to sample 
#' and 1 fore the ones that we want to sample. Default FALSE
#' @param tile Character vector or CVS file. Names of tile(s) to run a cvs file. 'ALL' will run all tiles in r_train_dir. Default is 'ALL'
#' @param text Character vector or CVS file. Names of text(s) to run in a cvs file. 'ALL' will run all tiles in text_train_dir. Default is 'ALL'
#' @param vuln_classes A list of the classes you want to model
#' The list can contain one or more vectors. Each vector represents a seperate vegetation class and response variable for the model
#' and the vector elements are synonyms used to describe that class
#' The fist place in each vector will be used in the output name used to store the calibrated model, so it should not contain spaces.
#' The other places should appear as attributes in the field 'field_name' of Pols
#' @param training_pol_filename Full path to the vector file (SpatialPointsDataFrame) of which one field contains the vuln.classes
#' @param field_name Character. The field in AOI.filename that contains the vuln_classes
#' @param ninputs_tile Integer. Number of inputs that we have fore each tile, including the tile, for exemple number of textures
#' @param data_outp_dir The folder where you want to save the sampled data
#' @param abs_samp Integer. How many 'absence' pixels should be randomly selected from eah tile to train the model. Default is 100.
#' @param parallel Boolean. Should the code be run in parallel using the doParallel package? Default is FALSE
#' @param nWorkers Integer. If running the ocde in parallel, how many workers should be used? Default is 4
#' @param data_outp_name Character. Name of the data to output
#' @param randompt Boolean. if True random points will be added in the tiles that doesn't have any visual point. Default TRUE
#' @param EOS if EOS true the for loop will be avoided if False will work with a for loop. Default FALSE
#' @return Saves a serialize object with the list with class-specific data frames of which the first column is the presence-absence response that can be used to train distribution model.
#' @examples \dontrun{
# tt <-  sample_points(r_train_dir <- '/EOS/projects/CANHEMON/data/ADS_CB_Buffer/',
#                               text_train_dir <- '/EOS/projects/CANHEMON/data/ADS_CB_Buffer/Textures/',
#                               tile = '/HDD/visual_interpretation/Buffer_list.csv',
#                               text = 'ALL', 
#                               prob_tifs = '/HDD/rasterbias/',
#                               vuln_classes = list(c('Pb')),
#                               training_pol_filename = '/HDD/visual_interpretation/visual_interpretation_ADS/ADS100_06032017_inspect_20171107.shp',
#                               field_name = 'type',
#                               ninputs_tile <- 27,
#                               data_outp_dir <- 'HDD/Move/',
#                               abs_samp = 10000,
#                               parallel = T,
#                               nWorkers = 16,
#                               data_outp_name = 'samp_Maxent_CB-Buffer',
#                               randompt = TRUE,
#                               EOS = FALSE)
#'                               }

#' @export
sample_points <- function(r_train_dir, text_train_dir, tile = 'ALL', text = 'ALL', prob_tifs = FALSE, vuln_classes, training_pol_filename, field_name, 
                          ninputs_tile, data_outp_dir, abs_samp = 100, parallel = F, nWorkers = 4, data_outp_name, randompt = TRUE, EOS = FALSE){
  require(maptools)
  require(raster)
  require(doParallel)
  require(dismo) 
  require(foreach)
  
  #Open the shapefile
  Pols <- raster::shapefile(training_pol_filename)
  
  #Tells if the dataframe is a factor, if yes erase the levels that are NULL. IN THIS CASE DATA IS NOT A FACTOR
  if(is.factor( Pols@data[[field_name]])){
    Pols@data[[field_name]] <- droplevels(Pols@data[[field_name]])
  }
  
  #only keep training points/polygons that fall in the vuln_classes to be considered
  Pols <- Pols[is.element(Pols@data[[field_name]] , unlist(vuln_classes)), ]
  
  all_tifs <- Tiles(r_train_dir, tile, text_train_dir, text)
  
  # a list to hold the outputs
  maxent_training_dfs <- list()
  tile_dat <- data.frame()
  #set up the cluster for parallel processing
  if(!EOS){
    `%op%` <- Par(nWorkers, data_outp_dir, data_outp_name, parallel) 
    
    stime <- system.time({maxent_training_dfs <- foreach(i = 1:length(tile), .combine = rbind.data.frame, .inorder=F) %op% {
      tile_i <- tile[i]
      cat(tile_i,'\n')
      
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
        else{
          #Pick the raster mask belonging to this tile with prob = 0 is NA and the othersw are probs
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
      
      #return the tile_dat at the end of each iteration
      tile_dat
    }
    cat("------------------------------------------\n")
    
    
    
    # report performance statistics 
    if (parallel){
      cat('using \n',foreach::getDoParWorkers(),' parallel workers,\n')
    }else{
      cat('processing sequentially on a single worker \n')
    }
    
    # close the cluster set up forparallel processing
    if (parallel){
      parallel::stopCluster(cl)
    }
    
    })
  }
  # save the extracted data ----
  if (!is.null(data_outp_dir)){
    data_file <- paste0(data_outp_name, '.rdsdata')
    data_file <- paste0(data_outp_dir,'/', data_file)
    saveRDS(maxent_training_dfs, file = data_file)
    cat('Wrote away ', data_file,'\n')
  }
  
  cat('Estimated in ',round(stime/60),' minutes\n')
  return()
}


