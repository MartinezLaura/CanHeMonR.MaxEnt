training_pnt_filename <- '/HDD/visual_interpretation/visual_interpretation_ADS/ADS100_Pbpointsonly.shp'

prediction_path <- "/HDD/Maxent/ADS/PCA/rasters/"

conf_dir <- '/HDD/Maxent/ADS/PCA/shp/'
radius <- 3
tt<-data.frame()
#treat the txt file if needed
bu = read.table('/HDD/finalcode/CanHeMonR-MaxEnt/tilesPb.txt')
bu = unlist(bu$V1)
bu = as.vector(bu)
for (cutoff in c(0.1,0.2,0.3,0.4,0.5,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)){
  tt_cutoff <- CanHeMonR.MaxEnt::sick_tree_errors(r_pred_dir = prediction_path,
                          tile = bu,
                          thresh = cutoff,
                          vuln_classes = list(c('Pb')),
                          pnts = raster::shapefile(training_pnt_filename),
                          radius = radius,
                          field_name = 'type',
                          abs_samp = 100,
                          parallel = F,
                          nWorkers = 20,
                          data_outp_dir = paste0(conf_dir,'/thresh_',cutoff)
  )
  tt_cutoff['cutoff'] <- cutoff
  tt <- rbind.data.frame(tt,tt_cutoff)
}

########## read all the confusion matrices.
rdsfiles <- list.files(conf_dir, full.names = T)
rdsfiles <- rdsfiles[grep('.rdsdata', rdsfiles)]
tt <- data.frame()
for (i in rdsfiles){
  ttcut <- readRDS(i)
  ttcut[['cutoff']] <- as.numeric(unlist(strsplit(unlist(strsplit(basename(i),c('_')))[2],"sicktree")))
  tt <- rbind.data.frame(tt, ttcut)
}


#summarize the nr of trees inspected per threshold
table(tt$vis,tt$class)/length(table(tt$cutoff))

####################
#calculate errors from the data
#Hadley Wickham - group of packages is called 'tidyverse' includes dplyr ggplot2 . Google Cheatsheet data wrangling R

tt$class <- factor(tt$class)
require(dplyr)
detach(package:plyr)
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
detach(package:plyr)
#++++++++++++++++++++++++++++++
#which is the highest cutoff providing 0 false negatives
#%>% slice(which.max(cutoff)
max_permittable_cutoff <- errs_long %>% group_by(class) %>% filter( error_type == 'false_neg') %>% filter( error < 20)
max_permittable_cutoff
max_permittable_cutoff_final <- max(max_permittable_cutoff$cutoff)
max_permittable_cutoff_final

require(ggplot2)
#mapping the percent correct
p <- ggplot(errs_long, aes(cutoff ,correct, colour = tree_type))
p + geom_line()+xlab("Cut-off for probability [0,255]")+ylab("Trees correctly detected [%]") + facet_grid(class~.)+geom_vline(xintercept = max_permittable_cutoff_final, color= 'green')

#mapping the errors
p <- ggplot(errs_long, aes(cutoff, error, colour = error_type))
p + geom_line()+xlab("Cut-off for probability [0,255]")+ylab("Error [%]") + facet_grid(class~.)+geom_vline(xintercept = max_permittable_cutoff_final, color= 'green')
#max_permittable_cutoff_final = 35

#++++++++++++++++++++++++++++++
#at this cutoff, what is the smallest nr of pixels detected inside a true tree-disk?
filter(errs_long, cutoff == max_permittable_cutoff_final)
npix_ <- tt %>% group_by(class) %>% filter((cutoff == max_permittable_cutoff_final) & (vis ==1) & (pred >= 1) ) %>% select(pred)
p2 <- ggplot(npix_, aes(x = pred)) + geom_histogram(binwidth = 5) + facet_grid(class~.) + xlab("nr of pixels above threshold in disk [r = 0.75 cm]")
p2
# we might not go for the absolute smallest, so instead get the 2% or 5% quantile of the number of pixels...
# the below table gives per class the minimum nr of pixels to keep (maximum cutoff) allowing 0%, 2% and 5% false negatives
max_permittable_cutoff_npix <- npix_ %>% group_by(class) %>% summarise(error_0 = min(pred, na.rm = T),
                                                                       error_2 = quantile(pred,0.02, na.rm = T),
                                                                       error_5 = quantile(pred,0.05, na.rm = T),
                                                                       error_10 = quantile(pred,0.10, na.rm = T),
                                                                       error_15 = quantile(pred,0.15, na.rm = T))

#keep the 5% threshold for Pb trees
max_permittable_cutoff_npix_final <- max_permittable_cutoff_npix %>% dplyr::filter(class=='Pb') %>% dplyr::select(error_5)
max_permittable_cutoff_npix_final <- as.numeric(max_permittable_cutoff_npix_final)
max_permittable_cutoff_npix_final 

p2 + geom_vline(aes(xintercept = error_0), max_permittable_cutoff_npix, lwd = 1) +
  geom_vline(aes(xintercept = error_2), max_permittable_cutoff_npix, lwd = 1.5) +
  geom_vline(aes(xintercept = error_5), max_permittable_cutoff_npix, lwd = 2) +
  geom_vline(aes(xintercept = error_10), max_permittable_cutoff_npix, lwd = 2.5) +
  geom_vline(aes(xintercept = error_15), max_permittable_cutoff_npix, lwd = 3)+
  geom_vline(xintercept = max_permittable_cutoff_npix_final, color= 'green')
#++++++++++++++++++++++++++++++
#the error analysis tells you how to convert the Maxent continuous output to binary output.
#now do that conversion for all the tiles



r_pred_dir <- prediction_path
outputpath = '/home/martlur/Desktop/CANHEMON/outputs/docker6_latest/shp/'

CanHeMonR::classify_maxent_output_based_on_error_stats(r_pred_dir = prediction_path,
                                            tile = 'ALL',
                                            max_permittable_cutoff_final  = max_permittable_cutoff_final,
                                            max_permittable_cutoff_npix_final = max_permittable_cutoff_npix_final,
                                            radius = radius,
                                            outp_dir = outputpath,
                                            parallel = T, nWorkers = 20)
#this could be vectorized and uploaded to Marco's site


