
#=====================================================================================================
#' @title A function that checks the parallel computation for missing data of MLR model.
#'
#' @param block  : The number of blocks for data cutting.
#' @param outputDirectory  : The directory of output files.
#'
#' @return NULL
#' @export Insepect_MLR
#'
#' @importFrom raster raster
#' @importFrom methods as
#'
#' @examples
#' \donttest{
#' Insepect_MLR(30, "./MlrOutput")
#' }
#'
Insepect_MLR <- function(block, outputDirectory) {
  name.x.temp <- NULL
  alldata <- NULL
  if(block == 1) {
    print("* Without parallel computing, the data is always intact *")
    return(-1)
  }

  for(idx in 1: block) {
    fileAddress <- paste(outputDirectory, "/variable.quantile_mlr_all_", idx, ".tif", sep = "")
    if(file.exists(fileAddress)) {
      print(fileAddress)
      fileRaster  <- raster(fileAddress, header=FALSE)
      dfRaster  <- as(fileRaster,"SpatialPointsDataFrame")
      allRaster <- data.frame(dfRaster)
      name.x.temp <- data.frame(allRaster)
      if(is.null(alldata)) {
        alldata = name.x.temp
      }else {
        names(alldata) = names(name.x.temp)
        alldata = rbind(alldata,name.x.temp)
      }
    }
  }
  myaddress1 <- paste(outputDirectory, "/variable.quantile_mlr_all", ".tif", sep = "")
  alldata1 <- raster(myaddress1, header = FALSE)
  alldata1 <- as(alldata1, "SpatialPointsDataFrame")
  alldata1 <- data.frame(alldata1)
  # ============== check =================
  # Deep random number, loop traversal comparison operation
  mynrows <- nrow(alldata)
  xnums <- sample(1: mynrows, floor(mynrows * 0.3))
  flag <- TRUE
  for(i in 1: length(xnums)) {
    for(j in 1: 3) {
      if(!is.na(alldata[xnums[i], j] != alldata1[xnums[i], j])) {
        if(alldata[xnums[i], j] != alldata1[xnums[i], j]) {
          flag <- FALSE
        }
      } else {
        flag <- FALSE
      }
    }
  }

  if(flag) {
    print("* Data integrity. *")
  } else {
    print("* Data is missing. The file path is as follows: *")
    print(fileAddress)
    print(myaddress1)
  }
}


#=====================================================================================================
#' @title A function that checks the parallel computation for missing data of RF model.
#'
#' @param block  : The number of blocks for data cutting.
#' @param outputDirectory  : The directory of output files.
#'
#' @return NULL
#' @export Insepect_RF
#'
#' @importFrom raster raster
#' @importFrom methods as
#'
#' @examples
#' \donttest{
#' Insepect_RF(30, "./RfOutput")
#' }
#'
Insepect_RF <- function(block, outputDirectory) {
  if(block == 1) {
    print("* Without parallel computing, the data is always intact *")
    return(-1)
  }
  name.x.temp <- NULL
  alldata <- NULL

  for(idx in 1: block) {
    fileAddress <- paste(outputDirectory, "/variable.quantile_rf_all", "_", idx, ".tif", sep = "")
    if(file.exists(fileAddress)) {
      print(fileAddress)
      fileRaster  <- raster(fileAddress, header=FALSE)
      dfRaster  <- as(fileRaster,"SpatialPointsDataFrame")
      allRaster <- data.frame(dfRaster)
      name.x.temp <- data.frame(allRaster)
      if(is.null(alldata)) {
        alldata = name.x.temp
      }else {
        names(alldata) = names(name.x.temp)
        alldata = rbind(alldata,name.x.temp)
      }
    }
  }
  myaddress1 <- paste(outputDirectory, "/variable.quantile_rf_all", ".tif", sep = "")
  alldata1 <- raster(myaddress1, header = FALSE)
  alldata1 <- as(alldata1, "SpatialPointsDataFrame")
  alldata1 <- data.frame(alldata1)
  # ============== check =================
  # Deep random number, loop traversal comparison operation
  mynrows <- nrow(alldata)
  xnums <- sample(1: mynrows, floor(mynrows * 0.3))
  flag <- TRUE
  for(i in 1: length(xnums)) {
    for(j in 1: 3) {
      if(!is.na(alldata[xnums[i], j] != alldata1[xnums[i], j])) {
        if(alldata[xnums[i], j] != alldata1[xnums[i], j]) {
          flag <- FALSE
        }
      } else {
        flag <- FALSE
      }
    }
  }

  if(flag) {
    print("* Data integrity. *")
  } else {
    print("* Data is missing. The file path is as follows: *")
    print(fileAddress)
    print(myaddress1)
  }
}


#=====================================================================================================
#' @title A function that checks the parallel computation for missing data of QRF model.
#'
#' @param block  : The number of blocks for data cutting.
#' @param outputDirectory  : The directory of output files.
#'
#' @return NULL
#' @export Insepect_QRF
#'
#' @importFrom raster raster
#' @importFrom methods as
#'
#' @examples
#' \donttest{
#' Insepect_QRF(30, "./QrfOutput")
#' }
#'
Insepect_QRF <- function(block, outputDirectory) {
  if(block == 1) {
    print("* Without parallel computing, the data is always intact *")
    return(-1)
  }
  # Init Variables
  name.x.temp <- NULL
  alldata <- NULL
  fragList <- c('05', '50', '95')
  for (i in 1: length(fragList)) {
    for(idx in 1: block) {
      fileAddress <- paste(outputDirectory, "/variable.quantile",fragList[i],"_", idx, ".tif", sep = "")
      print(fileAddress)
      if(file.exists(fileAddress)) {
        print(fileAddress)
        fileRaster  <- raster(fileAddress, header=FALSE)
        dfRaster  <- as(fileRaster,"SpatialPointsDataFrame")
        allRaster <- data.frame(dfRaster)
        name.x.temp <- data.frame(allRaster)
        if(is.null(alldata)) {
          alldata = name.x.temp
        }else {
          names(alldata) = names(name.x.temp)
          alldata = rbind(alldata,name.x.temp)
        }
      }
    } # for end
    # Read the total data set
    myaddress1 <- paste(outputDirectory, "/variable.quantile", fragList[i], "_all.tif", sep = "")
    alldata1 <- raster(myaddress1, header = FALSE)
    alldata1 <- as(alldata1, "SpatialPointsDataFrame")
    alldata1 <- data.frame(alldata1)
    # ============== check =================
    # Deep random number, loop traversal comparison operation
    mynrows <- nrow(alldata)
    xnums <- sample(1: mynrows, floor(mynrows * 0.3))
    flag <- TRUE
    for(i in 1: length(xnums)) {
      for(j in 1: 3) {
        if(!is.na(alldata[xnums[i], j] != alldata1[xnums[i], j])) {
          if(alldata[xnums[i], j] != alldata1[xnums[i], j]) {
            flag <- FALSE
          }
        } else {
          flag <- FALSE
        }
      }
    }

    if(flag) {
      print("* Data integrity. *")
    } else {
      print("* Data is missing. The file path is as follows: *")
      print(fileAddress)
      print(myaddress1)
    }
  } # for end
}


#=====================================================================================================
#' @title A function that checks the parallel computation for missing data
#'
#' @param model  : The models were selected, including QRF,RF and MLR.
#' @param block  : The number of blocks for data cutting.
#' @param outputDirectory  : The directory of output files.
#'
#' @return NULL
#' @export InsepectionVariable
#'
#' @importFrom raster raster
#' @importFrom methods as
#'
#' @examples
#' \donttest{
#' InsepectionVariable(model = "MLR", block = 30, outputDirectory = "MlrOutput")
#' }
#'
#'
#'
InsepectionVariable <- function(model = 'MLR', block, outputDirectory) {
  InspectFunc <- NULL
  Inspect_Func <- switch(
    model,
    MLR = Insepect_MLR,
    RF = Insepect_RF,
    QRF = Insepect_QRF
  )
  Inspect_Func(block, outputDirectory)
}
