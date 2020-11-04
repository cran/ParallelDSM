# Parallel computing preprocessing data
# the beginning of parallel computing work preparation

#=========================================================================================================
#' @title As a data preprocessing function, sets some global variables that are not visible to the user
#' @param Fpath : the path of file
#' @param FEncoding : the encoding of file
#' @param nblock : the number of blocks for data cutting
#' @param ncore : Computes the CPU's kernel in parallel(fill in according to the computer configuration)
#' @param foldername : Name of the folder in which the soil data is stored
#' @param sample.name : Standard soil files, which can be used as sample files (under in the FolderName)
#'
#' @return NULL
#' @export Preprocessing
#'
#' @importFrom utils read.csv
#' @importFrom stats sd
#' @importFrom raster res
#' @importFrom sp proj4string
#'
#' @examples
#' mydatas <- system.file("extdata", "all.input.csv", package = "ParallelDSM")
#' sampledatas <- system.file("extdata", "covariate", package = "ParallelDSM")
#' Preprocessing(mydatas,1,2,2,sampledatas,"twi.tif")
#'
#'
#' @references{
#' Breiman, L. (2001). Random forests. Mach. Learn. 45, 5–32.
#' Meinshausen, N. (2006) "Quantile Regression Forests", Journal of Machine Learning Research 7,
#' 983-999 http://jmlr.csail.mit.edu/papers/v7/
#' Song, X.D., Ge, G.Q., Zhang, G.L. and Wu, H.Y. ParallelDSM: A R package for parallel soil mapping. Computers & Geosciences (to be available in 2021)
#' }
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
dsm.env$max.change <- NULL
dsm.env$xtrain <- NULL
dsm.env$ytrain <- NULL
dsm.env$name.x.variable <- c()
dsm.env$qrf.variable <- NULL
Preprocessing <- function(Fpath,FEncoding,nblock,ncore,foldername,sample.name) {
  # Generate environment variables
  # Create a new environment variable to store global variables without exposing them to the public.

  # =========================================================

  # Import dependent packages
  # library(raster,warn.conflicts=F)
  # Sets the current path to the working path
  # setwd(getwd())
  # @Param : FEncoding =>  Set the file read format encoding
  # The file contains Chinese words to set GBK
  if(FEncoding == 1){
    FEncoding <- 'GBK'
  }else if(FEncoding == 0){
    FEncoding <- 'UTF-8'
  }else {
    FEncoding <- 'UTF-8'
  }
  # @Param : df.input => Read the metadata file
  dsm.env$df.input<-read.csv(file = Fpath,sep=",",fileEncoding = FEncoding)
  # @Param : numCloumn => Number of columns of data
  dsm.env$numColumn <- length(names(dsm.env$df.input))
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
  # sample data(Standardized data)
  dsm.env$sample.path <- paste(foldername,"/",sample.name,sep="")
  dsm.env$rmap_variable <- raster::raster(dsm.env$sample.path)
  # get information about data
  dsm.env$nr <- dsm.env$rmap_variable@nrows
  dsm.env$nc <- dsm.env$rmap_variable@ncols
  # calculation resolution
  dsm.env$resolutions <- res(dsm.env$rmap_variable)[1]
  # create projection
  dsm.env$pro <- proj4string(dsm.env$rmap_variable)
  # foldername
  dsm.env$foldername <- paste(foldername,"/",sep="")
}
#=========================================================================================================
#  Parallelint function ====> Compute the function part in parallel
#=======================================================================================
#' @title Parallel computing initialization preparation
#'
#' @return Represents whether the loading of the required variables and dependent packages is complete
#' @export
#'
#' @importFrom stats sd
#'
#' @examples
#' \donttest{
#' Parallelinit()
#' }
#'
#'
Parallelinit <- function() {
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
  # The data backup
  dsm.env$df.input <- dsm.env$df.nameVariable
  dsm.env$xtrain <- dsm.env$df.input[,(names(dsm.env$df.input) %in% dsm.env$name.x.variable)]
  dsm.env$ytrain <- dsm.env$df.input$ln.variable
  dsm.env$qrf.variable <- quantregForest::quantregForest(x=dsm.env$xtrain, y=dsm.env$ytrain)
}
#========================================================================================
# dsmParallel function =====> Main function
#=====================================================================================================
#' @title dsmparallel computings
#' @param outpath : Output path of the result of the prediction file. The default is "output".
#'
#' @return NULL
#' @export dsmParallel
#'
#' @importFrom raster predict
#' @importFrom sp coordinates<-
#' @importFrom sp gridded<-
#' @importFrom rgdal writeGDAL
#'
#' @examples
#' \donttest{
#' dsmParallel(outpath = "myoutputs")
#' }
#' @references{
#' Breiman, L. (2001). Random forests. Mach. Learn. 45, 5–32.
#' Meinshausen, N. (2006) "Quantile Regression Forests", Journal of Machine Learning Research 7,
#' 983-999 http://jmlr.csail.mit.edu/papers/v7/
#' Song, X.D., Ge, G.Q., Zhang, G.L. and Wu, H.Y. ParallelDSM: A R package for parallel soil mapping. Computers & Geosciences (to be available in 2021)
#' }
dsmParallel <- function(outpath) {
  # Load the required functions
  Parallelinit()
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
  resolutions <- dsm.env$resolutions
  rmap_variable <- dsm.env$rmap_variable
  sample.path <- dsm.env$sample.path
  sdsx <- dsm.env$sdsx
  xtrain <- dsm.env$xtrain
  ytrain <- dsm.env$ytrain
  #===================================================================================
  ParallelComputingVariable <- function(idx) {
    warnings('off')
    # Parallel computations are performed for each predictive variable
    for(k in 1:length(name.x.variable)){
      # Interception of predicted values
      predictor.k <- GetPredictorSubset(name.x.variable[k], idx, nblock,foldername,nr,nc,resolutions,pro)
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
  snowfall::sfLibrary(stats)
  "
  eval(parse(text=mylibrary))

  # Loads the relevant dependency packages
  # Cluster operations using the Snowfall parallel computing function

  # Loading variables
  snowfall::sfExport("nblock","sample.path","rmap_variable", "nr", "nc",
                     "resolutions", "pro", "ncore", "name.x.variable",
                     "df.all.sub", "df.input", "meansx", "sdsx", "qrf.variable","foldername","xtrain","df.nameVariable",
                     "ids","max.change","name_variable","numColumn","ytrain")


  # Start gets the current system time
  # and saves the run time by doing parallel operations on each partitioned block
  start <- Sys.time()
  rtest <-  snowfall::sfLapply(1:nblock, ParallelComputingVariable)
  print(Sys.time()-start)

  # End parallel returns resources such as memory
  snowfall::sfStop()
}
#=====================================================================================================



