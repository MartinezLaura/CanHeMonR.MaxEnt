#' @title Intersect COS file
#' @description Intersect the shapefile resulting from the last step of Maxent with the COS file.
#' It also can be applied to any intersection between to shapefiles
#' @param outpath Path wwhere you want to save the resultin shapefiles
#' @param shppath Path where the Maxent shapefiles are
#' @param COSpath Full path plus name of the COS file
#' @param prefix Prefix of the Maxent shapefiles, if there is a prefix
#' @return Several shapefies with the points that intersect with the COS file
#' @export
COSoutput <- function(outpath, shppath, COSpath, prefix=NULL) {
  require(sf)
  require(dplyr)
  require(sp)
  start.time <- Sys.time()
  
  
  cos_shp = sf::st_read(COSpath, quiet = TRUE)
  
  all_shp <- list.files(shppath, recursive = F, full.names = T, pattern = ".shp")
  
  for (shps in (all_shp)){
    shps <- normalizePath(shps)
    outname <- paste0(outpath,basename(shps))
    f.exists <- file.exists(outname)
    if (any(f.exists)){
      print(paste0('File already exists:',outname))
    }else{
      loc = basename(shps)
      loc = sub(prefix, "", loc)
      loc = sub(".shp", ".tif", loc)
      print(loc)
      filter = cos_shp %>% select(location) %>% filter(location == loc)
      if (nrow(filter)>0){
        shp = sf::st_read(shps, quiet = TRUE)
        shp = shp %>% st_buffer(0)
        filter = filter %>% st_buffer(0)
        inter = sf::st_intersection(shp, filter)
        if (nrow(inter)>0){
          inter = as(inter,'Spatial')
          raster::shapefile(inter,outname)
        }
        else{
          print(paste0(shps,". All features have been filtered. No shapefile will be outputed"))
        }
      }else{
        print(paste0(shps,". All features have been filtered. No shapefile will be outputed"))
      }
    }
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(time.taken)
}

#' @title Intersect road CODfile
#' @description Intersect the shapefile with the COS and roads file to erase the points that intersects with the polygons.
#' @param outpath Path wwhere you want to save the resultin shapefiles.
#' @param shppath Path where the Maxent shapefiles are
#' @param RoadCOSpath Full path plus name of the road file from OSM
#' @param prefix Prefix of the Maxent shapefiles, if there is a prefix
#' @return Several shapefies with the points that don't intersect with the Road file
#' @export
Roadoutput <- function(outpath, shppath, RoadCOSpath, prefix=NULL) {
  require(sf)
  require(dplyr)
  require(sp)
  start.time <- Sys.time()
  
  
  cosroads_shp = sf::st_read(RoadCOSpath, quiet = TRUE)
  
  all_shp <- list.files(shppath, recursive = F, full.names = T, pattern = ".shp")
  
  for (shps in (all_shp)){
    shps <- normalizePath(shps)
    outname <- paste0(outpath,basename(shps))
    f.exists <- file.exists(outname)
    if (any(f.exists)){
      print(paste0('File already exists:',outname))
    }else{
      loc = basename(shps)
      loc = sub(prefix, "", loc)
      loc = sub(".shp", ".tif", loc)
      print(loc)
      filter = cosroads_shp %>% select(location) %>% filter(location == loc)
      if (nrow(filter)>0){
        shp = sf::st_read(shps, quiet = TRUE)
        shp = sf::st_centroid(shp)
        filter = sf::st_buffer(filter, 0)
        st_crs(shp) <- st_crs(filter)
        inter = sf::st_intersection(shp, filter)
        if (nrow(inter)>0){
          final = shp[!shp$geometry %in% inter$geometry,]
          if (nrow(final)>0){
            sf::st_write(final, outname)
          }else{
            print(paste0(shps,".All features have been filtered"))
          }
        }else{
          print(paste0(shps,".No filtering has been done. theres no intersecting polygons"))
          basetocopy = basename(outname)
          basetocopy = sub(".shp","", basetocopy)
          listtocopy = list.files(shppath, pattern = basetocopy)
          listtocopy = paste(shppath, listtocopy, sep = "")
          file.copy(listtocopy, outpath)
          
        }
      }else{
        print(paste0(shps,".No filtering has been done. theres no intersecting polygons"))
        basetocopy = basename(outname)
        basetocopy = sub(".shp","", basetocopy)
        listtocopy = list.files(shppath, pattern = basetocopy)
        listtocopy = paste(shppath, listtocopy, sep = "")
        file.copy(listtocopy, outpath)
      }
    }
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(time.taken)
}


#' @title Polygont to point
#' @description Creates a unique point shapefiles from several polygons shapefiles.
#' @param inpath Path where the polygons shapefiles are
#' @param outpath Path where you weant to save the resulting point shapefile
#' @param chunks Path and name of the shapefile with the chunks
#' @param name_chunk Name of the chunk that you want to process
#' @return A unique shapefile for the centroids of the polygons
#' @export
Pointshp <- function(inpath, outpath, chunks =  '/HDD/visual_interpretation/visual_interpretation_ADS/ADS100_Chunks.shp', 
                     name_chunk){
  start.time <- Sys.time()
  blocks_shp = sf::st_read(chunks, quiet = TRUE)
  all_shp <- list.files(inpath, recursive = F, full.names = T, pattern = ".shp")
  all_shp <- normalizePath(unique(all_shp))
  blockfil = blocks_shp %>% select(block,location) %>% filter(block == name_chunk)
  control = TRUE
  for (loc in (blockfil$location)){
    loc = as.character(loc)
    baseshp = unlist(strsplit(loc,"[.]"))[1][1]
    shp_path <- all_shp[grepl(baseshp,all_shp)]
    f.exists <- file.exists(shp_path)
    if (any(f.exists)){
      shp = sf::st_read(shp_path, quiet = TRUE)
      if (control == TRUE){
        shp_centroids = NULL
        #shp = filter(shp, area >= 2)
        shp_centroids = sf::st_centroid(shp)
        control = FALSE
      }else{
        #shp = filter(shp, area >= 2)
        shp_centroids = rbind(shp_centroids,sf::st_centroid(shp))
      }
    }
    
  }
  sf::write_sf(shp_centroids, paste0(outpath,'.shp'))
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(time.taken)
}


#' @title Point in polygon
#' @description Saves in a data frame all the intersections between a polygon shapefile and a point shapefile.
#' Used for the Crown delineation method
#' @param inpath Path where the polygons shapefiles are
#' @param outpath Path where you weant to save the resulting point shapefile
#' @param pointpath Shapefile where all the PB points are
#' @return A rdsdata dataframe called intersections
#' @export
Intersection <- function(inpath, outpath, pointpath){
  start.time <- Sys.time()
  point_shp = sf::st_read(pointpath)
  all_shp <- list.files(polypath, recursive = F, full.names = T, pattern = ".shp")
  all_shp <- unique(all_shp)
  intersections = data.frame()
  for (poly_path in (all_shp)){
    poly_shp = sf::st_read(poly_path, quiet = TRUE)
    aux = sf::st_intersection(poly_shp,point_shp)
    intersections = rbind.data.frame(intersections,aux)
  }
  saveRDS(intersections, paste0(outpath,'/intersection.rds'))
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(time.taken)
}


#' @title Buffered point method
#' @description Buffer out a point with an specific size and then compare the classes Pb and fals_pos to see if we can separate them in the post-process stage. 
#' It also creates the graphs to check this information.
#' A boxplot graph with all the quantiles and a graph to see how much fasle_pos can be extracted vs Pb retrieved.
#' @param pbtiles txt file with the tiles names that you want to analizse
#' @param visualpoints shapefile where all the PB points are
#' @param polyspath Path where the shapefiles to analise
#' @param rasterpath Path to the raster files associates with the shapefiles
#' @param prefix Prefix of the Maxent shapefiles, if there is a prefix
#' @param outfile Path where you want to save the stats and graphs generated
#' @return An rdsdata dataframe with the base info and the plots
#' @export
Buffer <- function(pbtiles, visualpoints, polyspath, rasterpath, prefix, outfile){  
  start.time <- Sys.time()
  tiles <- read.table(pbtiles)
  tiles <- unlist(tiles$V1)
  tiles <-  gsub(".tif", "", tiles)
  visual_shp = sf::st_read(visualpoints, quiet = TRUE)
  
  #initialize the data frames
  pb_min_dat =setNames(data.frame(matrix(ncol = 4, nrow = 0)))
  pb_max_dat = setNames(data.frame(matrix(ncol = 4, nrow = 0)))
  pb_fisrtq_dat = setNames(data.frame(matrix(ncol = 4, nrow = 0)))
  pb_med_dat = setNames(data.frame(matrix(ncol = 4, nrow = 0)))
  pb_thirdq_dat = setNames(data.frame(matrix(ncol = 4, nrow = 0)))
  
  fp_min_dat =setNames(data.frame(matrix(ncol = 4, nrow = 0)))
  fp_max_dat = setNames(data.frame(matrix(ncol = 4, nrow = 0)))
  fp_fisrtq_dat = setNames(data.frame(matrix(ncol = 4, nrow = 0)))
  fp_med_dat = setNames(data.frame(matrix(ncol = 4, nrow = 0)))
  fp_thirdq_dat = setNames(data.frame(matrix(ncol = 4, nrow = 0)))
  for (t in tiles){
    poly_shp = sf::st_read(paste0(polyspath,prefix,t,".shp"), quiet = TRUE)
    aux_pb = sf::st_intersection(poly_shp,visual_shp)
    if (nrow(aux_pb) > 0){
      buffer = sf::st_buffer(aux_pb, 8)
      trast = raster::stack(paste0(rasterpath,t,".tif"))
      counter = 1
      for(g in buffer$geometry){
        trans_sf = as.vector(sf::st_bbox(g))[c(1, 3, 2, 4)]
        values = raster::crop(trast,trans_sf)
        pdf(paste0(outfile,"truepos/", t,"-hist-",counter))
        raster::hist(values)
        dev.off()
        pdf(paste0(outfile,"truepos/", t,"-crop-",counter))
        raster::plotRGB(values)
        dev.off()
        aux_pb = summary(values)
        colnames(aux_pb) = c('r','g', 'b', 'n')
        pb_min_dat = rbind(pb_min_dat, aux_pb[1, ])
        pb_max_dat = rbind(pb_max_dat, aux_pb[5, ])
        pb_fisrtq_dat = rbind(pb_fisrtq_dat, aux_pb[2, ])
        pb_med_dat = rbind(pb_med_dat, aux_pb[3, ])
        pb_thirdq_dat = rbind(pb_thirdq_dat, aux_pb[4, ])
        #dev_dat = rbind(dev_dat, aux_pb)
        counter = counter + 1
      }
    }
    else{
      buffer = sf::st_buffer(poly_shp, 8)
      trast = raster::stack(paste0(rasterpath,t,".tif"))
      counter = 1
      buffer = buffer[sample(nrow(buffer), 5), ]
      for(g in buffer$geometry){
        trans_sf = as.vector(st_bbox(g))[c(1, 3, 2, 4)]
        values = raster::crop(trast,trans_sf)
        pdf(paste0(outfile, "falsepos/",t,"-hist-",counter))
        raster::hist(values)
        dev.off()
        pdf(paste0(outfile, "falsepos/",t,"-crop-",counter))
        raster::plotRGB(values)
        dev.off()
        aux_fp = summary(values)
        colnames(aux_pb) = c('r','g', 'b', 'n')
        fp_min_dat = rbind(fp_min_dat, aux_fp[1, ])
        fp_max_dat = rbind(fp_max_dat, aux_fp[5, ])
        fp_fisrtq_dat = rbind(fp_fisrtq_dat, aux_fp[2, ])
        fp_med_dat = rbind(fp_med_dat, aux_fp[3, ])
        fp_thirdq_dat = rbind(fp_thirdq_dat, aux_fp[4, ])
        #dev_dat = rbind(dev_dat, aux_pb)
        counter = counter + 1
      }
    }
  }
  
  #ploting info
  require(ggplot2)
  library(reshape2)
  require(ggpubr)
  
  #assign here the data frame you want to plot
  plots = c("_min_dat","_max_dat" ,"_fisrtq_dat","_med_dat" ,"_thirdq_dat")
  outfile <- "/HDD/testplots/"
  
  for (p in plots){
    df = eval(parse(text = paste0("fp",p)))
    names(df) = c('r','g', 'b', 'n')
    df$id = 1:nrow(df)
    final_data_fp <- melt(df, id='id')
    names(final_data_fp) <- c('buffer', 'band', 'value')
    final_data_fp$variables <- "False pos"
    
    #ggplot() + geom_boxplot(data = final_data_fp, aes(x = band, y = value, color = band, group = band), size = 0.5)
    #ggsave(paste0(outfile,"fp",p,".pdf"),plot = last_plot())
    
    df = eval(parse(text = paste0("pb",p)))
    names(df) = c('r','g', 'b', 'n')
    df$id = 1:nrow(df)
    final_data_tp <- melt(df, id='id')
    names(final_data_tp) <- c('buffer', 'band', 'value')
    final_data_tp$variables <- "True pos"
    final_data = merge(final_data_fp, final_data_tp, all = T)
    #ggplot() + geom_line(data = final_data_tp, aes(x = buffer, y = value, color = band, group = band), size = 0.5)
    #ggsave(paste0(outfile,"pb",p, ".pdf"),plot = last_plot())
    plot <- ggplot(final_data, aes(x = band, y = value)) + 
      geom_boxplot(aes(fill=variables), size = 0.5) + ggtitle(p)
    ggsave(paste0(outfile,p, ".pdf"),plot = plot)
    
    
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(time.taken)
}


#' @title Filter shapefile with th ADS 2012 evolution on NDVI
#' @description Filter shapefile with th ADS 2012 evolution on NDVI.
#' @param filtpath path tho the shapefiles that will act as filter
#' @param shppath shapefile that you want to filter
#' @param prefix Prefix of the Maxent shapefiles, if there is a prefix
#' @param outpath Path where you want to save the shp generated
#' @return a shapefile filtered
#' @export
ortho2012 <- function(filtpath, shppath, outpath, prefix){  
  require(sf)
  require(dplyr)
  require(sp)
  start.time <- Sys.time()
  
  all_shpfil <- list.files(filtpath, recursive = F, full.names = T, pattern = ".shp")
  basefil <- basename(all_shpfil)
  all_shp <- list.files(shppath, recursive = F, full.names = T, pattern = ".shp")
  
  for (shps in (all_shp)){
    shps <- normalizePath(shps)
    outname <- paste0(outpath,basename(shps))
    f.exists <- file.exists(outname)
    if (any(f.exists)){
      print(paste0('File already exists:',outname))
    }else{
      loc = basename(shps)
      loc = sub(prefix, "", loc)
      print(loc)
      filname = all_shpfil[grepl(loc, all_shpfil)]
      
      if (length(filname)>0){
        filshp <- sf::st_read(filname)
        shp = sf::st_read(shps, quiet = TRUE)
        shp = sf::st_centroid(shp)
        filshp = sf::st_buffer(filshp, 0)
        inter = sf::st_intersection(shp, filshp)
        if (nrow(inter)>0){
          final = shp[!shp$geometry %in% inter$geometry,]
          if (nrow(final)>0){
            sf::st_write(final, outname)
          }else{
            print(paste0(shps,".All features have been filtered"))
          }
        }else{
          print(paste0(shps,".No filtering has been done. theres no intersecting polygons"))
          file.copy(shps, outpath)
        }
      }else{
        print(paste0(shps,". No filtering has been done. There is no polygons to filter"))
        file.copy(shps, outpath)
      }
    }
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(time.taken)
}

#' @title NDVI by moving window
#' @description Calculates the NDVI given a moving window size.
#' @param rast raster to calculate the ndvi
#' @param sizemv Size desired for the moving window
#' @return the resulting NDVI raster
#' @export
NDVI <-function(rast, sizemv){
  return (raster::focal(x, w=matrix(1,sizemv,sizemv), fun=mean((rast[[4]]-rast[[1]]) / (rast[[4]]+rast[[1]]))))
  
}
