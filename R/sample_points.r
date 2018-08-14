#' @title Sample training data for image classification from multiple image tiles
#' @description For each class in .shp polygon file, Sample training data for image classification from multiple image tiles using their raster bricks as predictors.
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
#' @param EOS If EOS true the for loop will be avoided if False will work with a for loop. Default FALSE
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
Sample_points <- function(r_train_dir, text_train_dir, tile = 'ALL', text = 'ALL', 
                          prob_tifs = FALSE, vuln_classes, training_pol_filename, field_name, 
                          ninputs_tile, data_outp_dir, abs_samp = 100, parallel = F, nWorkers = 4, 
                          data_outp_name, randompt = TRUE, EOS = FALSE){
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
  
  out <- Tiles(r_train_dir, tile, text_train_dir, text)
  all_tifs <- unlist(out[1])
  tile <- unlist(out[2])
  
  # a list to hold the outputs
  maxent_training_dfs <- list()
  tile_dat <- data.frame()
  #set up the cluster for parallel processing
  if(!EOS){
    out <- Par(nWorkers, data_outp_dir, data_outp_name, parallel) 
    `%op%` <- unlist(out[[1]])
    cl <- out[[2]]
    stime <- system.time({maxent_training_dfs <- foreach(i = 1:length(tile), .combine = rbind.data.frame, .inorder=F) %op% {
      tile_i <- tile[i]
      cat(tile_i,'\n')

      tile_dat <- GetPoints(tile_i, all_tifs, field_name, ninputs_tile, randompt, prob_tifs, Pols, vuln_classes, abs_samp, tile_dat)
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
    
    # close the cluster set up for parallel processing
    if (parallel){
      parallel::stopCluster(cl)
    }
    
    })
  }
  if (EOS){
    tile_i = tile
    cat(tile_i,'\n')
    tile_dat <- GetPoints(tile_i, all_tifs, field_name, ninputs_tile, randompt, prob_tifs, Pols, vuln_classes, abs_samp, tile_dat)
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


