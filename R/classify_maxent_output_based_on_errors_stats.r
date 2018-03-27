#' @title Classify maxent output to binary presence absence maps
#' @description Classify maxent output to binary presence absence maps based on thresholds of
#' pixel-level probability, expected minimum crown size, and clump size
#' @param r_pred_dir A directory where binary .tifs predicting presence as 1 and absence as 0 can be found for multiple tiles
#' @param tile Character vector. Names of tile(s) to run. 'ALL will run all tiles in r_pred_dir. Default is 'ALL'
#' @param max_permittable_cutoff_final Cutoff-level for Maxent probability output. Typically produced by error_summaries.R
#' @param max_permittable_cutoff_npix_final Cutoff-level for the nr of pixels (clumpsize) of presence that a circles of radius radius
#' needs to contain to be maintained as presence. Typically produced by error_summaries.R
#' @param radius The radius of the disk in which clump size is counted
#' @param outp_dir the directory to write output images to.
#' @param parallel Should the code be run in parallel using the doParallel package? Default is FALSE.
#' @param nWorkers If running the ocde in parallel, how many workers should be used? Default is 4.
#' @return empty
#'
#' @export
classify_maxent_output_based_on_error_stats <- function(r_pred_dir,
                                                        tile = 'ALL',
                                                        max_permittable_cutoff_final  = max_permittable_cutoff_final,
                                                        max_permittable_cutoff_npix_final = max_permittable_cutoff_npix_final,
                                                        radius = 0.75,
                                                        outp_dir,
                                                        parallel = F, nWorkers = 4){



  #harvest all the tif files in the directories holding covariate/predictor images
  all_tifs <- list.files(r_pred_dir, recursive = F, full.names = T)
  all_tifs <- all_tifs[grepl('.tif',all_tifs)]
  #excluded the tif files in the unprojected folder
  all_tifs <- all_tifs[!grepl('orig_noPRJ', all_tifs)]
  #avoid crashing on .tif.aux.xml files
  all_tifs <- all_tifs[substr(all_tifs, nchar(all_tifs)-3,nchar(all_tifs)) == ".tif"]
  all_tifs <- all_tifs[!grepl('.aux.xml', all_tifs)]

  require(dplyr)

  for(i in (1:length(all_tifs))){
    prediction_tiff <- raster::raster(all_tifs[i])
    if (any((prediction_tiff[prediction_tiff> max_permittable_cutoff_final])>0)){
      #the implementation goes in two steps
      #1. apply the cutoff to the maxent output to make it binary
      #2. run a moving window through the image that removes pixel clusters smaller than the m
      #mw_focweight <- raster::focalWeight(prediction_tiff_binary, d = 0.75, type = 'circle')
      # mw_focweight <- raster::focalWeight(prediction_tiff, d = radius, type = 'circle')
      # mw_focweightNA <- mw_focweight
      # mw_focweightNA[mw_focweightNA == 0] <- NA
      # mw_focweightNA[mw_focweightNA != 0] <- 1
      # quantfun <- function(x, ...){(sum(na.omit(x)) > max_permittable_cutoff_npix_final)}
      #
      # raster::rasterOptions(progress = 'text')
      #
      outp_name <- gsub(".tif", '', paste0(outp_dir, "/",basename(all_tifs[i])))
      #
      # #this took 10 hours to run on my machine
      # filt_probs <- raster::focal(x = prediction_tiff[prediction_tiff> max_permittable_cutoff_final], w = mw_focweightNA, fun = quantfun,
      #                             filename = outp_name, pad = T,overwrite=T)
      # # filt_probs <- raster::focal(x = prediction_tiff > max_permittable_cutoff_final, w = mw_focweightNA, fun = quantfun,
      # #                             filename = outp_name, pad = T,overwrite=T)
      #prediction_tiff = prediction_tiff > max_permittable_cutoff_final
      prediction_tiff = prediction_tiff > max_permittable_cutoff_final
      shp <- poligonize(prediction_tiff,all_tifs[i])
      max_pix_area = units::set_units(((raster::res(prediction_tiff)^2)*max_permittable_cutoff_npix_final)[1],m^2)
      aux = shp %>% mutate(area = st_area(geometry)) %>%  filter(area>max_pix_area)
      outp_name <- gsub(".tif", '.shp', paste0(outp_dir, "/",basename(all_tifs[i])))
      st_write(aux, outp_name, delete_dsn = TRUE)
    }

  }
  return()
}

