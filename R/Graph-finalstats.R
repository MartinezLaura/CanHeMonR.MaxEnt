#' @title Outputdataset
#' @description Create a dataset of the visual inspected shapefile and save it.
#' @param shp_insp Path to the visual inspected shapefile
#' @param listrasters csv file with the silt of locations(rasternames) you want to study. 
#' @param prefix Prefix that you want to put to the plots in order to differenciate them. Default = NULL
#' @param savepath Path where you want to save the dataframe 
#' @return A dataframe with all the values of the validated shapefile. Used to create the stadistics of the model
#' @export
Createdataset <- function(shp_insp, listrasters, shp_poly_path, prefix, savepath){
  require(sf)
  require(dplyr)
  require(sp)
  require(ggplot2)
  require(reshape2)

  #Gather the shp belonging to this chunk
  shp <- list.files(shp_poly_path, pattern = "*.shp", recursive = F, full.names = F)
  shp <- shp[!grepl('.aux.xml', shp)]
  # shp = gsub(prefix, "", files)
  # names_shp = gsub(".shp", "", files)
  tifs <- read.csv(listrasters, header = FALSE)
  locs <- as.character(unlist(tifs$V1))
  
  shp_df <- sf::st_read(shp_insp)
  first = TRUE
  for (l in locs){
    filter = shp_df %>% select(location, type) %>% filter(location == l)
    if (nrow(filter) != 0){
      name_shp = gsub(".tif", ".shp", l)
      s = shp[grepl(name_shp, shp)]
      s = sf::st_read(paste0(shp_poly_path, s))
      s = s %>% select(area)
      s = sf::st_buffer(s, 0.9)
      inter = sf::st_intersection(s,filter)
      
      if (first){
        final = inter
        first = FALSE
      }else{
        final = rbind(final,inter)
      }
      
      
    }
    
  }
  saveRDS(final, paste0(savepath, prefix, "rdsdata"))
}

#' @title Validation plots
#' @description Creation of the plots related to the validation of the predictions of MaxEnt.
#' @param pathrds Path to the serialized object with the information of the visual inspected area
#' @param prefix Prefix that you want to put to the plots in order to differenciate them. Default = NULL
#' @param savepath Path where you want to save the dataframe 
#' @return A series of plots of Maxent results, true positives and False positives, Accurancy tratio, error, etc.
#' @export
Validationplots<- function(pathrds, validationweb = TRUE, savepath, prefix = NULL){
  final = readRDS(path)
  if (validationweb != TRUE){
    tmp = final %>% select(area, type)%>% filter(type != "oth")
    tmp$valid = 0
    tmp$valid[tmp$type=="t"]<-1
    tmp$valid[tmp$type=="t_n"]<-1
    tmp = tmp %>% select(area, valid, type)}
  else{
      tmp = final %>% select(area, is_valid) %>% filter(is_valid != "oth")
      tmp$valid = 0
      tmp$valid[tmp$is_valid=="t"]<-1
      tmp$valid[tmp$is_valid=="t_n"]<-1
      tmp = tmp %>% select(area, valid, is_valid)}
  tmp$area = round(tmp$area,2)
  graph1 = tmp %>% group_by(area) %>% summarise(f=sum(valid==0), t=sum(valid==1))
  graphhist = graph1 %>% mutate(cumulativesum_True = cumsum(t)/sum(t), cumulativesum_FalsePositive = cumsum(f)/sum(f), cumratio = cumulativesum_True/cumulativesum_FalsePositive)
  
  # Convert data to long-form with 'melt' from the reshape2 package.
  graphhist = melt(graphhist, id.vars=c("area"),
                   measure.vars=c("cumulativesum_FalsePositive", "cumulativesum_True"))
  
  
  #mapping the count of each
  ggplot(graph1, aes(x=area)) + geom_line(aes(y=cumsum(t), color=  "True Positives")) + 
    geom_line(aes(y=cumsum(f), color =  "False Positives"))+
    xlab("Area [m2]")+ylab("Cumulative number of retrived trees") 
  ggsave(paste0(savepath, prefix, "trees-Area.pdf"))
  
  #mapping the ratio of eachno cumulative
  ggplot(graph1, aes(x=area)) + 
    geom_line(aes(y=cumsum(f/sum(f)), color =  "False Positives"))+ 
    xlab("Area [m2]")+ylab("Cumulative ratio of retrived trees") +
    geom_smooth(aes(y=cumsum(f/sum(f))), method = "lm", formula = y ~ log(x), se = TRUE, color="black", size = 1, fill = "green") +
    
    ggsave(paste0(savepath, prefix,"ratio-area.pdf"))
  
  graph2 = graph1 %>% group_by(area) %>% summarise(FDR = (f/nrow(tmp)), PPV =(t/nrow(tmp)))
  
  
  #graph FP/TP ratio vs True positive
  ggplot(graph2, aes(x=area)) + 
    geom_line(aes(y=FDR, col ="False Discovery Rate (FP/Totalpop)"))+ geom_line(aes(y=PPV, col=  "Positive Predictive Value(TP/Totalpop)")) +
    xlab("Area [m2]")+ylab("Ratio")
  ggsave(paste0(savepath, prefix,"FPR-TPR-area.pdf"))

  graph4 <- graph1 %>% mutate(cumt = cumsum(t)/sum(t), cumf = cumsum(f)/sum(f))
  graph4$cumf[graph4$cumf== 0] <- 0.001 
  graph4$cumratio = cumsum(graph4$f)/cumsum(graph4$t)
  normalizer <- 1.0 / max(graph4$cumratio)
  
  ggplot(graph4,aes(y=cumratio, x=area)) +
    geom_line(data=graph4,
              aes(x = area,y = cumt / normalizer), 
              color = "Blue")+ ylab("Cumulative ratio FP/TP")+
    theme(axis.text.y = element_text(colour = "red"))+
    geom_line(color = "Red") +
    theme(axis.text.y.right = element_text(colour = "blue"))+
    scale_y_continuous(sec.axis = sec_axis(trans= ~.*normalizer,
                                           name = 'Cumulative ratio True positives')) + 
    labs(x = "Area [m2]")
  
  
  ggsave(paste0(savepath, prefix,"Two-axes-area.pdf"))
}