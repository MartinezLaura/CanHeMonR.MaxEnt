% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeometryOperations.R
\name{Pointshp}
\alias{Pointshp}
\title{Polygont to point}
\usage{
Pointshp(inpath, outpath,
  chunks = "/HDD/visual_interpretation/visual_interpretation_ADS/ADS100_Chunks.shp",
  name_chunk)
}
\arguments{
\item{inpath}{Path where the polygons shapefiles are}

\item{outpath}{Path where you weant to save the resulting point shapefile}

\item{chunks}{Path and name of the shapefile with the chunks}

\item{name_chunk}{Name of the chunk that you want to process}
}
\value{
A unique shapefile for the centroids of the polygons
}
\description{
Creates a unique point shapefiles from several polygons shapefiles.
}
