#' @title Read Threshold
#' @description Read the thresholds calculated in: and create a dataframe with the information.
#' @param input_dir Full path where the serialized objects with the threshold informtion are \link[CanHeMonR.MaxEnt]{sick_tree_errors.r}
#' @return A dataframe contaning all the infor of the different threshold for the selected tiles
#' @export
ReadThresh <- function(input_dir) {
  ########## read all the confusion matrices.
  rdsfiles <- list.files(input_dir, full.names = T)
  rdsfiles <- rdsfiles[grep('.rdsdata', rdsfiles)]
  tt <- data.frame()
  for (i in rdsfiles){
    ttcut <- readRDS(i)
    ttcut[['cutoff']] <- as.numeric(unlist(strsplit(unlist(strsplit(basename(i),c('_')))[2],"sicktree")))
    tt <- rbind.data.frame(tt, ttcut)
  }
  return(tt)
}

#' @title Maximum cutoff
#' @description Selects, given a maximum error, the cut-off probability and the minimum amout of pixels that a tree have to have. 
#' This information is important to do the final polygonization of the raster images outputted by MaxEnt.
#' After that it binarize the rasters and transformates then into a polygon shapefile.
#' @param input_dir Full path where the serialized objects with the threshold informtion are
#' @param plots Boolean. If TRUE plots will be created. Default FALSE
#' Plot-1: %of trees correctly detected % of false positives
#' Plot-2: error in the tree dectection
#' Plot-3: Distribution of the number of pixels for the trees correctly decteted
#' @param stadisticspath Full path where the plots will be saved
#' @param prefix Prefix that you want to put to the plots in order to differenciate them. Default = NULl
#' @param max_error Integer. Maximum error you want to select 
#' @param cla Name of the attribute that you want to do the test on. Ex: 'Pb'
#' @param r_pred_dir Path where the rasters outputted by Maxent to binarize and poligonize are.
#' @param tile Path to a txtfile with the names of the tile you want to execute
#' If you want to run all the images of r_pred_dir, leave it empty. Default is 'ALL'
#' @param max_area Maximum area in cm. If the area detected on the image is above this value it will be erase
#' @param res Pixel resolution of the image. Ex: 0.30
#' @param outp_dir Path where you want to save the shapefiles of the binarized raster
#' @return Polygon shapefile
#' @export
MaxCutoff <- function(input_dir, plots, stadisticspath, prefix, max_error, cla, 
                      r_pred_dir, tile = 'ALL', max_area, res, outp_dir){
  
  tt <- ReadThresh(input_dir)
  
  #summarize the nr of trees inspected per threshold
  table(tt$vis,tt$class)/length(table(tt$cutoff))
  
  ####################
  #calculate errors from the data
  tt$class <- factor(tt$class)
  require(dplyr)
  #calculate omission and commission error
  errs <- tt %>% group_by(class, cutoff) %>% summarise(
    #false positive rate
    false_pos = 100*length(which((vis == 0) & (pred > 0)))/length(which(vis==0)),
    #false negative rate
    false_neg = 100*length(which((vis == 1) & (pred == 0)))/length(which(vis==1)),
    #smallest n_pixels detected in reference trees
    min_npix = min(pred[vis == 1]))
  
  ####################
  
  errs_long <- errs %>% tidyr::gather(key = error_type, value = error, - cutoff, -class, - min_npix)
  errs_long$correct <- 100 - errs_long$error
  errs_long$tree_type <- plyr::revalue(errs_long$error_type, c('false_pos' = 'mapped tree', 'false_neg' ='known tree'))
  #++++++++++++++++++++++++++++++
  max_permittable_cutoff <- errs_long %>% group_by(class) %>% filter( error_type == 'false_neg') %>% filter( error < max_error)
  max_permittable_cutoff_final <- max(max_permittable_cutoff$cutoff)
  cat('The final cut-off probability selected is: ',max_permittable_cutoff_final, ' with an error of: ', max_permittable_cutoff[max_permittable_cutoff$cutoff==0.7,  ]$error,'%','\n')
  
  # we might not go for the absolute smallest, so instead get the 2% or 5% quantile of the number of pixels...
  # the below table gives per class the minimum nr of pixels to keep (maximum cutoff) allowing 0%, 2% and 5% false negatives
  npix_ <- tt %>% group_by(class) %>% dplyr::filter((cutoff == max_permittable_cutoff_final) & (vis ==1) & (pred >= 1) ) %>% select(pred)
  
  max_permittable_cutoff_npix <- npix_ %>% group_by(class) %>% summarise(error_0 = min(pred, na.rm = T),
                                                                         error_2 = quantile(pred,0.02, na.rm = T),
                                                                         error_5 = quantile(pred,0.05, na.rm = T),
                                                                         error_10 = quantile(pred,0.10, na.rm = T),
                                                                         error_15 = quantile(pred,0.15, na.rm = T))
  
  #keep the 5% threshold for Pb trees
  max_permittable_cutoff_npix_final <- max_permittable_cutoff_npix %>% dplyr::filter(class==cla) %>% dplyr::select(error_5)
  max_permittable_cutoff_npix_final <- as.numeric(max_permittable_cutoff_npix_final)

  cat('The minimum amout of pixels that a crown has to have is: ',max_permittable_cutoff_npix_final, ' with an error of: 5%' ,'\n')  
  
  
  if(plots){
    require(ggplot2)
    #mapping the percent correct
    p <- ggplot(errs_long, aes(cutoff ,correct, colour = tree_type))
    p + geom_line()+xlab("Cut-off for probability [0,1]")+ylab("Trees correctly detected [%]") + facet_grid(class~.)+geom_vline(xintercept = max_permittable_cutoff_final, color= 'green')
    ggsave(paste0(prefix, '-TreesCorrectlyDetected.pdf'), plot = last_plot(),path = stadisticspath,  width = 40, height = 20, units = "cm")
    #mapping the errors
    p <- ggplot(errs_long, aes(cutoff, error, colour = error_type))
    p + geom_line()+xlab("Cut-off for probability [0,1]")+ylab("Error [%]") + facet_grid(class~.)+geom_vline(xintercept = max_permittable_cutoff_final, color= 'green')
    ggsave(paste0(prefix, '-TreesCorrectlyDetected-error.pdf'), plot = last_plot(),path = stadisticspath,  width = 40, height = 20, units = "cm")
    
    
    #++++++++++++++++++++++++++++++
    #at this cutoff, what is the smallest nr of pixels detected inside a true tree-disk?
    filter(errs_long, cutoff == max_permittable_cutoff_final)
    p2 <- ggplot(npix_, aes(x = pred)) + geom_histogram(binwidth = 5) + facet_grid(class~.) + xlab("nr of pixels above threshold in disk [r = 0.75 cm]")
    
    p2 + geom_vline(aes(xintercept = error_0), max_permittable_cutoff_npix, lwd = 1) +
      geom_vline(aes(xintercept = error_2), max_permittable_cutoff_npix, lwd = 1.5) +
      geom_vline(aes(xintercept = error_5), max_permittable_cutoff_npix, lwd = 2) +
      geom_vline(aes(xintercept = error_10), max_permittable_cutoff_npix, lwd = 2.5) +
      geom_vline(aes(xintercept = error_15), max_permittable_cutoff_npix, lwd = 3)+
      geom_vline(xintercept = max_permittable_cutoff_npix_final, color= 'green')
    #++++++++++++++++++++++++++++++
    #the error analysis tells you how to convert the Maxent continuous output to binary output.
    ggsave(paste0(prefix, '-MinimumNumberPixels.pdf'), plot = last_plot(),path = stadisticspath,  width = 40, height = 20, units = "cm")
    
  }
  
  classify_maxent_output_based_on_error_stats(r_pred_dir,
                                              tile = 'ALL',
                                              max_permittable_cutoff_final  = max_permittable_cutoff_final,
                                              max_permittable_cutoff_npix_final = max_permittable_cutoff_npix_final,
                                              max_area,
                                              res,
                                              outp_dir)
  
}

