# Parallel computing preprocessing data
# the beginning of parallel computing work preparation

#=========================================================================================================
#' @title As a data preprocessing function, sets some global variables that are not visible to the user
#' @param Fpath : the path of file
#' @param fn : Name of the folder in which the soil data is stored
#' @param tname : Standard soil files, which can be used as sample files (under in the FolderName)
#' @param mc : Read file mode (using data sets, or read yourself)
#' @param icsv : Use df.input from the built-in dataset
#' @param itif: Use df.dem in the built-in data set
#' @param nblock : the number of blocks for data cutting
#' @param ncore : Computes the CPU's kernel in parallel(fill in according to the computer configuration)
#' @param Fc : the encoding of file
#'
#' @return NULL
#' @export Preprocess
#'
#' @importFrom utils read.csv
#' @importFrom stats sd
#' @importFrom raster res
#' @importFrom sp proj4string
#'
#' @examples
#'
#' #####################################################################
#' ##  Example code 1                                                 ##
#' ##  Select your own reading method, as shown below                 ##
#' #####################################################################
#' mydatas <- system.file("extdata", "all.input.csv", package = "ParallelDSM")
#' sampledatas <- system.file("extdata", "covariate", package = "ParallelDSM")
#' Preprocess(mydatas,sampledatas,"twi.tif")
#'
#' #####################################################################
#' ##  Example code 2 (It is highly recommended)                      ##
#' ##  If you want to use test cases, load the relevant data sets     ##
#' #####################################################################
#' #  Select the data set that comes with this package
#'
#' #data("df.input")
#' #data("df.dem")
#'
#' #####################################################################
#' ##  Use the data file references that come with this package       ##
#' #####################################################################
#' #sampledatas <- system.file("extdata", "covariate", package = "ParallelDSM")
#'
#' #####################################################################
#' ## Use preprocessor functions to process the data that is loaded in##
#' #####################################################################
#' #Preprocess(fn = sampledatas,mc = TRUE,icsv = df.input,itif = df.dem)
#'
#' ############################################################################
#' ## This function is the main function that performs parallel computations ##
#' ## The outpath field refers to the filename of the data output            ##
#' ## The mymodels field has three modes to choose from: QRF,RF and MLR      ##
#' ## ‘QRF’ stands for Random Forest Model Prediction Method                 ##
#' ## ‘RF’ stands for Machine Learning Model Prediction Method               ##
#' ## ‘MLR’ stands for Multiple Linear Regression Prediction Model           ##
#' ## 'from' and 'to' are reserved fields that can be left unused by the user##
#' ############################################################################
#'
#' #DsmParallel(outpath = "myoutputs",mymodels = "MLR",from=1,to=200)
#'
#'
#'
#'
#' @references{
#' Breiman, L. (2001). Random forests. Mach. Learn. 45, 5–32.
#' Meinshausen, N. (2006) "Quantile Regression Forests", Journal of Machine Learning Research 7,
#' 983-999 http://jmlr.csail.mit.edu/papers/v7/
#' Song, X.D., Ge, G.Q., Zhang, G.L. and Wu, H.Y. ParallelDSM: A R package for parallel soil mapping. Computers & Geosciences (to be available in 2021)
#' }
#'
dsm.env <- new.env()
dsm.env$df.input <- NULL
dsm.env$numColumn <- 0
dsm.env$meansx <- 0
dsm.env$sdsx <- 0
dsm.env$name_variable <- NULL
dsm.env$ids <- NULL
dsm.env$nblock <- 0
dsm.env$resolutions <- NULL
dsm.env$pro <- NULL
dsm.env$foldername <- NULL
dsm.env$nr <- 0
dsm.env$nc <- 0
dsm.env$ncore <- 0
dsm.env$sample.path <- NULL
dsm.env$df.dem <- NULL
dsm.env$max.change <- NULL
dsm.env$xtrain <- NULL
dsm.env$ytrain <- NULL
dsm.env$name.x.variable <- c()
dsm.env$qrf.variable <- NULL
dsm.env$rf.variable <- NULL
dsm.env$mlr.variable <- NULL
dsm.env$outputnames <- NULL
dsm.env$choicemodel <- "QRF"
Preprocess <- function(Fpath="",fn="",tname="",mc=NULL,icsv=NULL,itif=NULL,nblock=6,ncore=2,Fc=1){
  # Generate environment variables
  # Create a new environment variable to store global variables without exposing them to the public.

  # =========================================================

  # Import dependent packages
  # library(raster,warn.conflicts=F)
  # Sets the current path to the working path
  # setwd(getwd())
  # @Param : Fc =>  Set the file read format encoding
  # The file contains Chinese words to set GBK
  if(Fc == 1){
    Fc <- 'GBK'
  }else if(Fc == 0){
    Fc <- 'UTF-8'
  }else {
    Fc <- 'UTF-8'
  }

  if(is.null(mc) == TRUE){
    # @Param : df.input => Read the metadata file
    dsm.env$df.input<-read.csv(file = Fpath,sep=",",fileEncoding = Fc)
    #print(dsm.env$df.input)
  }else{
    dsm.env$df.input<-icsv
    #print(dsm.env$df.input)
  }

  # @Param : numCloumn => Number of columns of data
  dsm.env$numColumn <- length(names(dsm.env$df.input))

  #print(dsm.env$numColumn)
  #print("=======================")
  #print(dsm.env$df.input)

  # @Param : numColumn => The number of columns of data
  dsm.env$df.input <- dsm.env$df.input[,c(1:dsm.env$numColumn)]
  # Digital processing of data
  # Digitize the data format for a given number of columns
  for(item in 1:dsm.env$numColumn){
    dsm.env$df.input[[item]] <- as.numeric(dsm.env$df.input[[item]])
  }
  # The operation of averaging
  # The default first line here is a numeric variable
  # @Param : meansx => Mean value of data
  dsm.env$meansx <- apply(dsm.env$df.input[,c(2:dsm.env$numColumn)],2,mean,na.rm=T)
  # @Param : sdsx => Data standard deviation
  dsm.env$sdsx <- apply(dsm.env$df.input[,c(2:dsm.env$numColumn)],2,sd,na.rm=T)
  # Data is discretized and decentralized
  dsm.env$df.input[,c(2:dsm.env$numColumn)] <- scale(dsm.env$df.input[c(2:dsm.env$numColumn)])

  # Convert the data
  # @Param : name_variable => the name of variable
  dsm.env$name_variable <- names(dsm.env$df.input[1])
  # @Param : ids => As a subscript for no missing data
  index_array <- is.na(dsm.env$df.input[dsm.env$name_variable])
  dsm.env$ids <- which(index_array==FALSE)
  # Converts data into data frames
  dsm.env$df.nameVariable <- as.data.frame(dsm.env$df.input[dsm.env$ids,])
  dsm.env$df.nameVariable <- dsm.env$df.nameVariable[,c(1:dsm.env$numColumn)]
  # @Param : ids => Flags for special data processing
  dsm.env$ids <- which( dsm.env$df.nameVariable[dsm.env$name_variable] < 0.01 )
  # Reset data less than 0.01 to 0.01
  dsm.env$df.nameVariable[dsm.env$name_variable][dsm.env$ids] <- 0.01
  #================================================
  # @Param : df.all.sub => as a predictive variable
  dsm.env$df.all.sub <- NULL
  # judge if it's missing
  if(is.na(nblock) == FALSE){
    dsm.env$nblock <- nblock
  }else{
    # set the default value
    dsm.env$nblock <- 10
  }
  # judge if it's missing
  if(is.na(ncore) == FALSE){
    dsm.env$ncore <- ncore
  }else{
    # set the default value
    dsm.env$ncore <- 2
  }


  if(is.null(mc) == TRUE){
    # sample data(Standardized data)
    dsm.env$sample.path <- paste(fn,"/",tname,sep="")
    dsm.env$rmap_variable <- raster::raster(dsm.env$sample.path)

    # the datas for merge file.
    dsm.env$df.dem <- dsm.env$rmap_variable
    dsm.env$df.dem <- as(dsm.env$df.dem,"SpatialPointsDataFrame")
    dsm.env$df.dem <- as.data.frame(dsm.env$df.dem)
  }else{
    dsm.env$sample.path <- paste(fn,"/",tname,sep="")
    dsm.env$rmap_variable <- itif
    # the datas for merge file.

    #print(itif)
    #print("================================")
    #print(dsm.env$rmap_variable)

    dsm.env$df.dem <- as(itif,"SpatialPointsDataFrame")

    dsm.env$df.dem <- as.data.frame(dsm.env$df.dem)
  }



  # get information about data
  dsm.env$nr <- dsm.env$rmap_variable@nrows
  dsm.env$nc <- dsm.env$rmap_variable@ncols
  # calculation resolution
  dsm.env$resolutions <- res(dsm.env$rmap_variable)[1]
  # create projection
  dsm.env$pro <- proj4string(dsm.env$rmap_variable)
  # foldername
  dsm.env$foldername <- paste(fn,"/",sep="")
}
#=========================================================================================================
#  Parallelint function ====> Compute the function part in parallel
#=======================================================================================
#' @title Parallel computing initialization preparation
#'
#' @param mymodel : The models were selected, including QRF,RF and MLR.
#'
#' @return Represents whether the loading of the required variables and dependent packages is complete
#' @export
#'
#' @importFrom stats sd
#' @importFrom stats lm
#' @importFrom stats as.formula
#' @examples
#' \donttest{
#' Parallelinit(mymodel = "QRF")
#' }
#'
#'
Parallelinit <- function(mymodel) {
  # Parallel computation of the prepare function
  # Eliminate the dimensional
  dsm.env$max.change <- mean(dsm.env$df.nameVariable[[dsm.env$name_variable]]) + 3*sd(dsm.env$df.nameVariable[[dsm.env$name_variable]])
  dsm.env$ids <- which(dsm.env$df.nameVariable[[dsm.env$name_variable]] > dsm.env$max.change)
  dsm.env$df.nameVariable[dsm.env$ids,][dsm.env$name_variable] <- dsm.env$max.change
  # Get a set of variables
  mylens <- ncol(dsm.env$df.nameVariable)
  dsm.env$name.x.variable <- c()
  for (nums in 2:mylens){
    dsm.env$name.x.variable <- c(dsm.env$name.x.variable,names(dsm.env$df.nameVariable[nums]))
  }
  # Handle special value
  dsm.env$df.nameVariable$ln.variable <- log(dsm.env$df.nameVariable[[dsm.env$name_variable]])

  # Train a global prediction model
  dsm.env$df.input <- dsm.env$df.nameVariable
  # Select the model : QRF、RF、MLR

  dsm.env$xtrain <- dsm.env$df.input[,(names(dsm.env$df.input) %in% dsm.env$name.x.variable)]
  dsm.env$ytrain <- dsm.env$df.input$ln.variable
  dsm.env$qrf.variable <- quantregForest::quantregForest(x=dsm.env$xtrain, y=dsm.env$ytrain)

  if(mymodel == "MLR"){
    myformula <- paste("ln.variable","~");

    for (nums in 3:mylens-1){
      partformula <- paste(names(dsm.env$df.nameVariable[nums]),"+");
      myformula <- paste(myformula,partformula);
    }

    myformula <- paste(myformula,names(dsm.env$df.nameVariable[mylens]));
    fmla.ak05 <- as.formula(myformula);
    dsm.env$mlr.variable <- lm(fmla.ak05, data = dsm.env$df.input)
    print("=========");
    print(dsm.env$mlr.variable);
    print("=========");

  }else if(mymodel == "RF"){
    # get the formula of all
    myformula <- paste("ln.variable","~");

    for (nums in 3:mylens-1){
      partformula <- paste(names(dsm.env$df.nameVariable[nums]),"+");
      myformula <- paste(myformula,partformula);
    }

    myformula <- paste(myformula,names(dsm.env$df.nameVariable[mylens]));
    fmla.ak05 <- as.formula(myformula);
    dsm.env$rf.variable <- randomForest::randomForest(fmla.ak05, data = dsm.env$df.input, importance=TRUE)

    print("=========");
    print(dsm.env$rf.variable);
    print("=========");
  }else{
    dsm.env$xtrain <- dsm.env$df.input[,(names(dsm.env$df.input) %in% dsm.env$name.x.variable)]
    dsm.env$ytrain <- dsm.env$df.input$ln.variable
    dsm.env$qrf.variable <- quantregForest::quantregForest(x=dsm.env$xtrain, y=dsm.env$ytrain)

    print("=========");
    print(dsm.env$qrf.variable);
    print("=========");
  }
  #print("=========")
  #print(myformula)
  #print("=========")
  #print(fmla.ak05)
  #print("=========")
  #print(dsm.env$mlr.variable);
  #print("=========")
  #print(dsm.env$qrf.variable);
}
#===============================================================================================
# DsmParallel function =====> Main function
#===============================================================================================
#' @title DsmParallel computings
#' @param outpath : Output path of the result of the prediction file. The default is "output".
#' @param mymodels : The models were selected, including QRF,RF and MLR.
#' @param from : Which row to start cutting the matrix
#' @param to : Where does the last row of the cut matrix go
#'
#' @return NULL
#' @export DsmParallel
#'
#' @importFrom raster predict
#' @importFrom sp coordinates<-
#' @importFrom sp gridded<-
#' @importFrom rgdal writeGDAL
#'
#' @examples
#' \donttest{
#'
#' #####################################################################################
#' ## Review the documentation for the Preprocess function before using this function ##
#' #####################################################################################
#'
#' ############################################################################
#' ## This function is the main function that performs parallel computations ##
#' ## The outpath field refers to the filename of the data output            ##
#' ## The mymodels field has three modes to choose from: QRF,RF and MLR      ##
#' ## ‘QRF’ stands for Random Forest Model Prediction Method                 ##
#' ## ‘RF’ stands for Machine Learning Model Prediction Method               ##
#' ## ‘MLR’ stands for Multiple Linear Regression Prediction Model           ##
#' ## 'from' and 'to' are reserved fields that can be left unused by the user##
#' ############################################################################
#'
#' DsmParallel(outpath = "myoutputs",mymodels = "MLR",from=1,to=200)
#' }
#' @references{
#' Breiman, L. (2001). Random forests. Mach. Learn. 45, 5–32.
#' Meinshausen, N. (2006) "Quantile Regression Forests", Journal of Machine Learning Research 7,
#' 983-999 http://jmlr.csail.mit.edu/papers/v7/
#' Song, X.D., Ge, G.Q., Zhang, G.L. and Wu, H.Y. ParallelDSM: A R package for parallel soil mapping. Computers & Geosciences (to be available in 2021)
#' }
DsmParallel <- function(outpath,mymodels,from=NULL,to=NULL) {
  # Load the required functions
  Parallelinit(mymodel = mymodels)
  dsm.env$choicemodel <- mymodels
  dsm.env$outputnames <- outpath


  # Read / write between GDAL grid mapping and spatial objects
  # description :
  # =========
  # function reads or writes to the GDAL grid mapping.
  # If they can, they will set up the spatial reference system.
  # GDALinfo reports the size of the dataset and other parameters
  # Create2GDAL creates a GDAL dataset from the SpatialGridDataFrame object,
  # in particular being able to save to allow only replication and not creation
  # Build # GDAL driver format.   Cluster programming tools
  # =========

  # library(snowfall)
  requireNamespace("snowfall")


  #===================================================================================
  df.all.sub <- dsm.env$df.all.sub
  df.input <- dsm.env$df.input
  df.nameVariable <- dsm.env$df.nameVariable
  foldername <- dsm.env$foldername
  ids <- dsm.env$ids
  max.change <- dsm.env$max.change
  meansx <- dsm.env$meansx
  name_variable <- dsm.env$name.x.variable
  name.x.variable <- dsm.env$name.x.variable
  nblock <- dsm.env$nblock
  nc <- dsm.env$nc
  ncore <- dsm.env$ncore
  nr <- dsm.env$nr
  numColumn <- dsm.env$numColumn
  pro <- dsm.env$pro
  qrf.variable <- dsm.env$qrf.variable
  rf.variable <- dsm.env$rf.variable
  mlr.variable <- dsm.env$mlr.variable
  resolutions <- dsm.env$resolutions
  rmap_variable <- dsm.env$rmap_variable
  sample.path <- dsm.env$sample.path
  sdsx <- dsm.env$sdsx
  xtrain <- dsm.env$xtrain
  ytrain <- dsm.env$ytrain
  choicemodel <- dsm.env$choicemodel
  #===================================================================================
  ParallelComputingVariable <- function(idx) {
    warnings('off')

    # Parallel computations are performed for each predictive variable
    for(k in 1:length(name.x.variable)){
      # Interception of predicted values
      predictor.k <- GetPredictorSubset(name.x.variable[k], idx, nblock,foldername,nr,nc,resolutions,pro,from,to)
      # the mean of value
      meanx <- meansx[names(meansx)==name.x.variable[k]]
      # the sd of sdx
      sdx <- sdsx[names(sdsx)==name.x.variable[k]]
      # Eliminate the dimensional
      predictor.k[,1] <- (predictor.k[,1] - meanx)/sdx
      # The predictive variable is saved
      if(k==1) {
        df.all.sub <- predictor.k
      }else{
        s <- name.x.variable[k]
        df.all.sub[s] <- predictor.k[,1]}
    }
    # ====== Start parallel computing operations ======
    # The prediction of parallel computation is made according to the function of training prediction
    if(choicemodel == "QRF"){
      xtest <- df.all.sub[,(names(df.all.sub) %in% name.x.variable)]
      model.prediction <- predict(qrf.variable, xtest, what = c(0.05, 0.5, 0.95))

      # For variables that are not properly distributed, natural logarithm conversion is required,
      # and the predicted results require exponential function operation

      df.all.sub$variable.quantile05 <- exp(model.prediction[,1])
      df.all.sub$variable.quantile50 <- exp(model.prediction[,2])
      df.all.sub$variable.quantile95 <- exp(model.prediction[,3])


      # Build data box
      df.all2 <- as.data.frame(df.all.sub)
      # For the coordinate prediction of DF.ALL2
      coordinates(df.all2) <- c("x","y")
      # Grid dF.ALL2 / You can also see if the data is already grid
      gridded(df.all2) <- TRUE


      #output the idx_th block's predictions
      if(nblock == 1){
        # Determine if the file exists
        mydirs <- "outputall"
        if(!file.exists(mydirs)){
          dir.create(file.path(mydirs))
        }
        output.file.name1 <- paste("outputall/variable.quantile05_all.tif", sep = "")
        output.file.name2 <- paste("outputall/variable.quantile50_all.tif", sep = "")
        output.file.name3 <- paste("outputall/variable.quantile95_all.tif", sep = "")
      }else{
        # Determine if the file exists
        mydirs1 <- outpath
        if(!file.exists(mydirs1)){
          dir.create(file.path(mydirs1))
        }
        output.file.name1 <- paste(mydirs1,"/variable.quantile05_", idx, ".tif", sep = "")
        output.file.name2 <- paste(mydirs1,"/variable.quantile50_", idx, ".tif", sep = "")
        output.file.name3 <- paste(mydirs1,"/variable.quantile95_", idx, ".tif", sep = "")
      }

      writeGDAL(   dataset = df.all2["variable.quantile05"],  fname = output.file.name1,
                   drivername = "GTiff",  type = "Float32" )
      writeGDAL(   dataset = df.all2["variable.quantile50"],  fname = output.file.name2,
                   drivername = "GTiff",  type = "Float32" )
      writeGDAL(   dataset = df.all2["variable.quantile95"],  fname = output.file.name3,
                   drivername = "GTiff",  type = "Float32" )

    }else if(choicemodel == "RF"){
      xtest <- df.all.sub[,(names(df.all.sub) %in% name.x.variable)]
      model.prediction <- predict(rf.variable, newdata=xtest, type="class")

      df.all.sub$variable.quantileall <- exp(model.prediction)
      df.all2 <- as.data.frame(df.all.sub)
      coordinates(df.all2) <- c("x","y")
      gridded(df.all2) <- TRUE

      #output the idx_th block's predictions
      if(nblock == 1){
        # Determine if the file exists
        mydirs <- "outputall"
        if(!file.exists(mydirs)){
          dir.create(file.path(mydirs))
        }
        output.file.name <- paste("outputall/variable.quantile_rf_all.tif", sep = "")
      }else{
        # Determine if the file exists
        mydirs1 <- outpath
        if(!file.exists(mydirs1)){
          dir.create(file.path(mydirs1))
        }
        output.file.name <- paste(mydirs1,"/variable.quantile_rf_all_", idx, ".tif", sep = "")
      }

      writeGDAL(   dataset = df.all2["variable.quantileall"],  fname = output.file.name,
                   drivername = "GTiff",  type = "Float32" )



    }else if(choicemodel == "MLR"){
      xtest <- df.all.sub[,(names(df.all.sub) %in% name.x.variable)]
      #model.prediction <- predict(mlr.variable, newdata=xtest, interval="none")
      model.prediction <- predict(mlr.variable, newdata=xtest)




      df.all.sub$variable.quantileall <- exp(model.prediction)
      df.all2 <- as.data.frame(df.all.sub)
      coordinates(df.all2) <- c("x","y")
      gridded(df.all2) <- TRUE

      #output the idx_th block's predictions
      if(nblock == 1){
        # Determine if the file exists
        mydirs <- "outputall"
        if(!file.exists(mydirs)){
          dir.create(file.path(mydirs))
        }
        output.file.name <- paste("outputall/variable.quantile_mlr_all.tif", sep = "")
      }else{
        # Determine if the file exists
        mydirs1 <- outpath
        if(!file.exists(mydirs1)){
          dir.create(file.path(mydirs1))
        }
        output.file.name <- paste(mydirs1,"/variable.quantile_mlr_all_", idx, ".tif", sep = "")
      }

      writeGDAL(   dataset = df.all2["variable.quantileall"],  fname = output.file.name,
                   drivername = "GTiff",  type = "Float32" )
    }

    return (1)
  }
  #=====================================================================================================

  #===================================================================================
  # Cluster initialization setup kernel
  snowfall::sfInit(parallel=TRUE,cpus=dsm.env$ncore)

  mylibrary <- "
  snowfall::sfLibrary(snowfall)
  snowfall::sfLibrary(rgdal)
  snowfall::sfLibrary(raster)
  snowfall::sfLibrary(quantregForest)
  snowfall::sfLibrary(randomForest)
  snowfall::sfLibrary(stats)
  "
  eval(parse(text=mylibrary))

  # Loads the relevant dependency packages
  # Cluster operations using the Snowfall parallel computing function

  # Loading variables
  snowfall::sfExport("nblock","sample.path","rmap_variable", "nr", "nc",
                     "resolutions", "pro", "ncore", "name.x.variable",
                     "df.all.sub", "df.input", "meansx", "sdsx", "qrf.variable","rf.variable","mlr.variable","choicemodel","foldername","xtrain","df.nameVariable",
                     "ids","max.change","name_variable","numColumn","ytrain","from","to")

  # Start gets the current system time
  # and saves the run time by doing parallel operations on each partitioned block
  start <- Sys.time()
  rtest <-  snowfall::sfLapply(1:nblock, ParallelComputingVariable)
  print(Sys.time()-start)

  # End parallel returns resources such as memory
  snowfall::sfStop()


  if(dsm.env$nblock == 1){
    dsm.env$outputnames <- "outputall"
  }
  mystr_input <- paste(dsm.env$outputnames,"/",sep = "")
  f.i.d <- c(mystr_input)
  mystr_output <- paste(dsm.env$outputnames,"/",sep = "")
  f.o.d <- c(mystr_output)



  if(dsm.env$nblock != 1){

    if(dsm.env$choicemodel == "QRF"){

      f.iblock <- c("variable.quantile05_")
      mstrs <- paste(mystr_output,"variable.quantile05_all.tif",sep = "")
      f.suffix <- c(mstrs)
      MergingTiles(dsm.env$df.dem,f.i.d, f.iblock, dsm.env$nblock, f.o.d, f.suffix)

      f.iblock <- c("variable.quantile50_")
      mstrs <- paste(mystr_output,"variable.quantile50_all.tif",sep = "")
      f.suffix <- c(mstrs)
      MergingTiles(dsm.env$df.dem,f.i.d, f.iblock, dsm.env$nblock, f.o.d, f.suffix)

      f.iblock <- c("variable.quantile95_")
      mstrs <- paste(mystr_output,"variable.quantile95_all.tif",sep = "")
      f.suffix <- c(mstrs)
      MergingTiles(dsm.env$df.dem,f.i.d, f.iblock, dsm.env$nblock, f.o.d, f.suffix)

    }

    if(dsm.env$choicemodel == "MLR"){

      f.iblock <- c("variable.quantile_mlr_all_")
      mstrs <- paste(mystr_output,"variable.quantile_mlr_all.tif",sep = "")
      f.suffix <- c(mstrs)
      MergingTiles(dsm.env$df.dem,f.i.d, f.iblock, dsm.env$nblock, f.o.d, f.suffix)

    }

    if(dsm.env$choicemodel == "RF"){

      f.iblock <- c("variable.quantile_rf_all_")
      mstrs <- paste(mystr_output,"variable.quantile_rf_all.tif",sep = "")
      f.suffix <- c(mstrs)
      MergingTiles(dsm.env$df.dem,f.i.d, f.iblock, dsm.env$nblock, f.o.d, f.suffix)

    }


  }


}
#=====================================================================================================



