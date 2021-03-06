% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sick_tree_errors.r
\name{Sick_tree_errors}
\alias{Sick_tree_errors}
\title{Sick tree error calculator}
\usage{
Sick_tree_errors(r_pred_dir, tile = "ALL", prefix, vuln_classes,
  training_pnt_filename, radius = 2, field_name, abs_samp = 100,
  minthresh = 0, maxthresh = 1, stepthresh = 0.05, parallel = F,
  nWorkers = 4, data_outp_dir = NULL)
}
\arguments{
\item{r_pred_dir}{A directory where binary .tifs predicting presence as 1 and absence as 0 can be found for multiple tiles}

\item{tile}{path to a txtfile with the names of the tile you want to execute
If you want to run all the images of r_pred_dir, leave it empty. Default is 'ALL'}

\item{prefix}{Prefix that you want to add to the output file. EX: NorthPortugal NorthPortugal-> -sicktree_performance_dfs.rdsdata}

\item{vuln_classes}{A list of the classes you want to model.
The list can contain one or more vectors. Each vector represents a seperate vegetation class and response variable for the model
and the vector elements are synonyms used to describe that class.
The fist place in each vector will be used in the output name used to store the calibrated model, so it should not contain spaces.
The other places should appear as attributes in the field 'field_name' of pnts.}

\item{training_pnt_filename}{Full path to the SpatialPointsDataFrame of which one field contains the vuln_classes}

\item{radius}{The radius within which a presence point must be found for it to be considered 'correct'. Default is True}

\item{field_name}{Character. The field in AOI.filename that contains the vuln_classes}

\item{abs_samp}{How many 'absence' pixels should be randomly selected from each tile to evaluate the absences? Default is 100.}

\item{minthresh}{The minimum cut-off number you want to use, Default = 0}

\item{maxthresh}{The maximum cut-off number you want to use, Default = 1}

\item{stepthresh}{The step used to create de sequences of threshold from the minthresh to the max thresh. Default = 0.05}

\item{parallel}{Should the code be run in parallel using the doParallel package? Default is FALSE.}

\item{nWorkers}{If running the ocde in parallel, how many workers should be used? Default is 4.}

\item{data_outp_dir}{The folder to save the sampled data to. No data is saved is data_outp_dir is NULL. Default is NULL.}
}
\value{
A data frame with commission and ommission errors and sample sizes of presence and absence
}
\description{
Calculate errors in automated detection of declining trees using visual inspection data as reference.
}
\examples{
\dontrun{
}
}
