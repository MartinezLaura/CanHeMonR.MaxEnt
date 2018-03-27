#' @title poligonize Raster
#' @description Poligonize a raster with values between 0 and 1
#' @param rast to be poligonized
#' @param nameshp
#' @return The Layer object
#' @examples \dontrun{
#' }
#' @export
poligonize <- function(rast, nameshp){

  gdalformat = 'ESRI Shapefile'
  quiet=TRUE
  require(rgdal)
  require(sf)
  pypath <- Sys.which('gdal_polygonize.py')
  owd <- getwd()
  readpoly <- TRUE
  on.exit(setwd(owd))
  setwd(dirname(pypath))



  ##############
  #Creation of the shapefile to extract the radius to work with
  ##############
  raster::writeRaster(rast, {f <-paste0('/DATA/trash/',basename(nameshp))})
  outshape <- paste0('/DATA/trash/',basename(nameshp))
  outshape <- sub('.tif', '.shp', outshape)
  f.exists <- file.exists(outshape)
  if (any(f.exists)){
    cat('File already exists: %s',outshape)
  }
  rastpath <- normalizePath(f)
  system2('python', args=(sprintf('"%1$s" -mask "%2$s" "%2$s" -f "%3$s" "%4$s"',
                                    pypath, rastpath, gdalformat, outshape)))

  if (isTRUE(readpoly)) {
    shp <- sf::st_read(outshape)
    return(shp)
  }


}
