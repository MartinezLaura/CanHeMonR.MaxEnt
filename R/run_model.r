#' @title Run a saved MaxEnt model in predictive mode on a tile of image data
#' @description Run a saved MaxEnt model on a the image data selected.
#' @param predictor_dir Path where predictor layers are held, rasters. If EOS = FALSE predictor_dir is a path
#' if EOS = True predictor_dir is the path plus the image to predict
#' @param text_train_dir Path where .tifs of the textures associated with r_train_dir. It is really important to avoid
#' errors on the execution to pass the same numer of textures per tile as in the MaxEnt model used
#' @param MaxEntmodel_dir Path where the MaxEnt model file is held
#' @param fname_MaxEntmodel_r Filename of the MaxEnt model saved in rdsdata format
#' @param output_dir Path to write the output to
#' @param rastername Character. Prefix to give the outputed raster image, for control versions
#' @param model_type Character. Type of model of maxent you want to use: raw, logistic or cloglog
#' @param EOS If EOS true the for loop will be avoided if False will work with a for loop. Default FALSE
#' @return A raster image for each tile with the probabilities or cummulative probabilities of presence for each class
#' @seealso Depends on: \code{\link[CanHeMonR.MaxEnt]{calibrate_model.r}}
#' @examples \dontrun{
#' Run_model(predictor_dir = "/H03_CANHEMON/Imagery/Portugal/ADS100/ortophotos_06032017/geotif/pt616000-4404000.tif",
#'                             text_train_dir <-'/home/martlur/Documents/TexturesAds/',
#'                             MaxEntmodel_dir = "/home/martlur/Documents/Dockers/docker6EOS/",
#'                             fname_MaxEntmodel_r = "samp10000_Pb.rdsdata",
#'                             output_dir = "/DATA/Results/Rcode/OutputRunSickTree",
#'                             rastername = "samp1000_",
#'                             model_type = 'cloglog',
#'                             loop = FALSE)
#' }
#' @export
Run_model <- function(predictor_dir,
                      text_train_dir,
                      MaxEntmodel_dir,
                      fname_MaxEntmodel_r,
                      output_dir,
                      rastername,
                      model_type,
                      EOS){
  
  require(raster)
  require(rJava)
  require(dismo)
  
  
  if (EOS != TRUE){
    # #harvest all the tif files in the directories holding covariate/predictor images
    # all_tifs <- list.files(predictor_dir, recursive = F, full.names = T, pattern = ".tif")
    # all_tifs <- unique(all_tifs)
    # all_tifs <- all_tifs[grepl('.tif',all_tifs)]
    # #excluded the tif files in the unprojected folder
    # all_tifs <- all_tifs[!grepl('orig_noPRJ', all_tifs)]
    # #excluded the tif files ending like .aux.xml
    # all_tifs <- all_tifs[!grepl('.sh', all_tifs)]
    # all_tifs <- all_tifs[!grepl('.aux.xml', all_tifs)]
    # cat(length(all_tifs),' tiles are considered\n')
    
    #load the model
    mod2 <- readRDS(paste0(MaxEntmodel_dir, fname_MaxEntmodel_r))
    for (i in 1:length(all_tifs)){
      tile_i <- all_tifs[i]
      raster::rasterOptions(progress="")
      
      # extract the name of the tile to save it
      rast = unlist(lapply(strsplit(tile_i,"/"),tail,1))
      tile <- strsplit(tail(rast, n=1),".", fixed = TRUE)[[1]][1]
      cat(rast,' processing\n')
      
      #Pick the textures to find the ones that belongs to this tile
      tifs <- list.files(text_train_dir, recursive = T, full.names = T, pattern = ".tif")
      tifs <- tifs[!grepl('.aux.xml', tifs)]
      pred_rs <- append(rast, tifs[grepl(tile_i, tifs)])
      raster_fnames <- unlist(pred_rs)
      raster_fnames[1]<-paste0(predictor_dir,'/',raster_fnames[1])
      r_pred <- raster::stack(raster_fnames)
      
      #adjust the layernames to make them conform the model syntax
      names(r_pred) <- paste0('l',substr(names(r_pred),17,nchar(names(r_pred))))
      
      # run the model and write the output away to a file
      if (model_type == 'raw'){
        px <- dismo::predict(r_pred, mod2,  progress = 'text', args="outputformat=raw")}
      
      if (model_type == 'logistic'){
        px <- dismo::predict(r_pred, mod2,  progress = 'text' ,args="outputformat=logistic")}
      
      if (model_type == 'cloglog'){
        px <- dismo::predict(r_pred, mod2,  progress = 'text' ,args="outputformat=cloglog")}
      
      fname_output_tif <- paste0(output_dir,rastername, "-",basename(tile), '.tif')
      raster::writeRaster(px,filename = fname_output_tif, overwrite = F)
      cat('Output saved to :',fname_output_tif,'\n')
      
      return()
    }
  }
  if(EOS){
    #load the model
    mod2 <- readRDS(paste0(MaxEntmodel_dir, fname_MaxEntmodel_r))
    
    #load the predictor layers
    rast = lapply(strsplit(predictor_dir,"/"),tail,1)
    tile = strsplit(unlist(rast),"[.]")[[1]][1]
    #Pick the textures to find the ones that belongs to this tile
    all_tifs <- list.files(text_train_dir, recursive = F, full.names = T, pattern = tile)
    #excluded the tif files in the unprojected folder
    all_tifs <- all_tifs[!grepl('orig_noPRJ', all_tifs)]
    
    pred_rs = rast
    pred_rs <- append(pred_rs, all_tifs)
    
    raster_fnames <- unlist(pred_rs)
    raster_fnames[1]<-paste(predictor_dir)
    r_pred <- raster::stack(raster_fnames)
    
    #adjust the layernames to make them conform the model syntax
    names(r_pred) <- paste0('l',substr(names(r_pred),17,nchar(names(r_pred))))
    
    # run the model and write the output away to a file
    if (model_type == 'raw'){
      px <- dismo::predict(r_pred, mod2,  progress = 'text', args="outputformat=raw")}
    
    if (model_type == 'logistic'){
      px <- dismo::predict(r_pred, mod2,  progress = 'text' ,args="outputformat=logistic")}
    
    if (model_type == 'cloglog'){
      px <- dismo::predict(r_pred, mod2,  progress = 'text' ,args="outputformat=cloglog")}
    
    fname_output_tif <- paste0(output_dir,rastername, "-",basename(tile), '.tif')
    raster::writeRaster(px,filename = fname_output_tif, overwrite = F)
    cat('Output saved to :',fname_output_tif,'\n')
    
    return()
  }
}
