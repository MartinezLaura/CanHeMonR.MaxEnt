#' @title Classify maxent output to binary presence absence maps
#' @description Classify maxent output to binary presence absence maps based on thresholds of
#' pixel-level probability, expected minimum crown size, and clump size
#' @param r_pred_dir A directory where binary .tifs predicting presence as 1 and absence as 0 can be found for multiple tiles
#' @param tile Character vector. Names of tile(s) to run. 'ALL will run all tiles in r_pred_dir. Default is 'ALL'
#' @param max_permittable_cutoff_final Cutoff-level for Maxent probability output. Typically produced by error_summaries.R
#' @param max_permittable_cutoff_npix_final Cutoff-level for the nr of pixels (clumpsize) of presence that a circles of radius radius
#' needs to contain to be maintained as presence. Typically produced by error_summaries.R
#' @param max_area The maximum are aof a polygon, all the polygons above this area will be filtered
#' @param res The resolution of the raster
#' @param radius The radius of the disk in which clump size is counted
#' @param outp_dir the directory to write output images to.
#' @return empty
#'
#' @export
classify_maxent_output_based_on_error_stats <- function(r_pred_dir,
                                                        tile = 'ALL',
                                                        max_permittable_cutoff_final  = max_permittable_cutoff_final,
                                                        max_permittable_cutoff_npix_final = max_permittable_cutoff_npix_final,
                                                        max_area,
                                                        res,
                                                        radius = 0.75,
                                                        outp_dir){
  
  
  if (tile[1] == 'ALL'){
    #harvest all the tif files in the directories holding covariate/predictor images
    all_tifs <- list.files(r_pred_dir, recursive = F, full.names = T, pattern = ".tif")
    all_tifs <- unique(all_tifs)
    all_tifs <- all_tifs[grepl('.tif',all_tifs)]
    #excluded the tif files in the unprojected folder
    all_tifs <- all_tifs[!grepl('orig_noPRJ', all_tifs)]
    #excluded the tif files ending like .aux.xml
    all_tifs <- all_tifs[!grepl('.aux.xml', all_tifs)]
    cat(length(all_tifs),' tiles are considered\n')
  }else{
    tifs <- read.csv(tile, header = FALSE)
    tile <- unlist(lapply(tifs$V1, function(x) paste0(r_pred_dir,x)))
    tile <- tile[!grepl('.aux.xml', tile)]
    all_tifs <- unlist(tile)
    cat(length(all_tifs),' tiles are considered\n')
  }
  
  require(dplyr)
  min_pix_area = units::set_units(((res^2)*max_permittable_cutoff_npix_final),m^2)
  max_pix_area = units::set_units((max_area),m^2)
  
  for(i in (1:length(all_tifs))){
    outp_name <- gsub(".tif", '.shp', paste0(outp_dir, "/",basename(all_tifs[i])))
    shp.exists <- file.exists(outp_name)
    
    if (!any(shp.exists)){
      prediction_tiff <- raster::raster(all_tifs[i])

      if (any((prediction_tiff[prediction_tiff> max_permittable_cutoff_final])>0)){     
        prediction_tiff = prediction_tiff > max_permittable_cutoff_final
        shp <- poligonize(prediction_tiff,all_tifs[i])

        aux = shp %>% mutate(area = st_area(geometry)) %>%  filter(area > min_pix_area) %>%  filter(area < max_pix_area)
        
        st_write(aux, outp_name, delete_dsn = TRUE)
      }else{
        print("Probabilities are lower than the threshold set. The shapefile won't be outputed")
      }
    }else{
      print("File already exists")
    }
    
  }
  return()
}

