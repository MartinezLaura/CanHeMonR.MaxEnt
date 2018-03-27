#' @title Calibrate vegetation distribution models
#' @description For each class in .shp polygon file, calibrate a distribution model using a raster brick as predictors
#' @param vuln_classes A character vector of the classes you want to model. The should be presented in the column 'class' of training_df. Default 'ALL'
#' @param training_path Path to the rdsdata that contains the data.frame, with in the column 'pres' 1/0 to indicate presence absence, then covariate columns, and a colum 'class' groupin grows by
#' the land-cover class the data was sampled for. This df is typically generated by sample_points.r
#' @param model_outp_dir Path and filename prefix to save the model objects
#' @param name Character. Name that you want to give to the serialized object with the model
#' @param stadistics Boolean. If Ture stadistics of the model will be done and save. Take into account that it can take several time. Default FALSE 
#' @param myargs List. Arguments to pass to the maxent model, in the following format. Example: myargs <- c("noautofeature", "nohinge", "nothreshold", "noproduct","nolinear"). Default NULL
#' @param model_type Character. Type of model that we wnat to use to predict types: raw - logistic - cloglog
#' @param stadisticspath Path where you want to save the stadistics of the model. Default NULL
#' @return Serialized object with class-specific distribution models, using a data frame created from training points and covariate images
#' @seealso For more possible arguments for MaxEnt see: https://groups.google.com/forum/#!msg/Maxent/yRBlvZ1_9rQ/Fj8Two0lmHIJ
#' @seealso Depends on: \code{\link[CanHeMonR]{sick_tree_errors.r}}
#' @examples \dontrun{
#' tt <- calibrate_model(vuln_classes = list(c('Pb')), 
#'                       training_path = '/DATA/Results/RunSickTree_Output/ADS/Test/MaxentModel/MaxentCastelobranco_NA-buffer.rdsdata',
#'                       model_outp_dir = '/DATA/Results/RunSickTree_Output/ADS/Final_Iter1-2/MaxentModel/',
#'                       name='samp_NA_buffer_logistic', 
#'                       stadistics = TRUE,
#'                       myargs = c("noautofeature", "nohinge"),
#'                       model_type = 'cloglog',
#'                       stadisticspath = '/DATA/Results/RunSickTree_Output/ADS/Test/MaxentModel/Stadistics/')
#'                       }

#' @export
calibrate_model <- function(vuln_classes = 'ALL', training_path, model_outp_dir, name,
                            stadistics = FALSE, myargs = NULL, model_type, stadisticspath=NULL){
  
  require(maptools)
  require(raster)
  require(dismo)
  require(rJava)
  #option tospecify the maximum memory allocation pool for a Java Virtual Machine to avoid java.lang.OutOfMemoryError
  #It should be at most the total amount of phisical memory. Expressed in GB
  options(java.parameters = "-Xmx30g")
  
  training_df = readRDS(training_path)
  
  training_df = unique(training_df)
  
  #for each class that should be modelled
  if (vuln_classes == 'ALL'){
    vuln_classes <- unique(training_df$class)
  }
  
  for (class in vuln_classes){
    # get the data for this class
    class_rows <- which(training_df$class == class)
    class_resp <- training_df$pres[class_rows]
    class_pred <- within(training_df, rm(pres,class))[class_rows,]
    
    # calibrate the model
    dismo::maxent()
    if (is.null(myargs)){
      mod2 <- dismo::maxent(p = class_resp, x = class_pred)}
    else{
      mod2 <- dismo::maxent(p = class_resp, x = class_pred, args = myargs)
    }
    assign(paste0('mod.',class),mod2)

    # save the model
    model.file <- paste0(model_outp_dir,paste0(name,"-",class,".rdsdata"))
    saveRDS(mod2,file=model.file)
    
    cat('Wrote away ',model.file,'\n')
    
    if (stadistics ==TRUE){
      #Create de folds and select the 20% for training and the rest for test
      fold <- dismo::kfold(training_df, k=5)
      occtrain <- training_df[fold == 1, ]
      for (j in 2:5){
        occtest <- training_df[fold == j, ]
        #Create the trainning and test datasets for the validation
        trainclass_rows <- which(occtrain$class == class)
        trainclass_resp <- occtrain$pres[trainclass_rows]
        trainclass_pred <- within(occtrain, rm(pres,class))[trainclass_rows,]
        
        testclass_rows <- which(occtest$class == class)
        testclass_resp <- occtest$pres[testclass_rows]
        testclass_pred <- within(occtest, rm(pres,class))[testclass_rows,]
        
        if (is.null(myargs)){
          modValid <- dismo::maxent(p = trainclass_resp, x = trainclass_pred)}
        else{
          modValid <- dismo::maxent(p = trainclass_resp, x = trainclass_pred, args = myargs)}
        
        #plot showing importance of each variable
        pdf(paste0(stadisticspath, name, "-kfold", j, "-Variables.pdf"))
        plot(modValid)
        dev.off()
        
        # response curves
        pdf(paste0(stadisticspath, name, "-kfold", j, "-Response.pdf"))
        dismo::response(modValid)
        dev.off()
        
        #evaluate model by type:
        if (model_type == 'raw'){
          testp <- dismo::predict(modValid, subset(occtest, occtest$pres ==1),args="outputformat=raw")
          testa <- dismo::predict(modValid, subset(occtest, occtest$pres ==0),args="outputformat=raw")}
        if (model_type == 'logistic'){
          testp <- dismo::predict(modValid, subset(occtest, occtest$pres ==1),args="outputformat=logistic")
          testa <- dismo::predict(modValid, subset(occtest, occtest$pres ==0),args="outputformat=logistic")}
        
        if (model_type == 'cloglog'){
          testp <- dismo::predict(modValid, subset(occtest, occtest$pres ==1),args="outputformat=cloglog")
          testa <- dismo::predict(modValid, subset(occtest, occtest$pres ==0),args="outputformat=cloglog")}
        
        e3 <- dismo::evaluate(p= testp, a= testa)
        txtfilename = paste0(stadisticspath, name, "-kfold", j,"-Evaluation.rdsdata")
        saveRDS(e3, txtfilename)
        pdf(paste0(stadisticspath, name, "-kfold", j, "-AUC.pdf"))
        plot(e3,'ROC')
        dev.off()
        pdf(paste0(stadisticspath, name, "-kfold", j, "-kappa.pdf"))
        plot(e3,'kappa')
        dev.off()
        pdf(paste0(stadisticspath, name, "-kfold", j, "-TNR.pdf"))
        plot(e3,'TNR')
        dev.off()
        pdf(paste0(stadisticspath, name, "-kfold", j, "-TPR.pdf"))
        plot(e3,'TPR')
        dev.off()
      }
    }
  }
  return()
}