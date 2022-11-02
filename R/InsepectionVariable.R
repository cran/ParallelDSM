
#=====================================================================================================
#' @title A function that checks the parallel computation for missing data
#'
#' @param model  : The models were selected, including QRF,RF and MLR.
#' @param block  : The number of blocks for data cutting.
#' @param outputDirectory  : The directory of output files.
#'
#' @return NULL
#' @export
#'
#' @importFrom raster raster
#' @importFrom methods as
#'
#' @examples
#' \donttest{
#' InsepectionVariable(model = "MLR", block = 30, outputDirectory = "mlrOutput")
#' }
#'
#'
#'
InsepectionVariable <- function(model = 'MLR', block, outputDirectory) {
  InspectFunc <- NULL
  # MLR method function
  Inspect_MLR <- function() {
    name.x.temp <- NULL
    alldata <- NULL
    if(block == 1) {
      print("* Without parallel computing, the data is always intact *")
      return(-1)
    }

    for(idx in 1: block) {
      fileAddress <- paste(outputDirectory, "/variable.quantile_mlr_all_", idx, ".tif", sep = "")
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

    if(!flag) {
      print("* Data integrity. *")
    } else {
      print("* Data is missing. The file path is as follows: *")
      print(fileAddress)
      print(myaddress1)
    }
  }

  # RF method function
  Inspect_RF <- function() {
    if(block == 1) {
      print("* Without parallel computing, the data is always intact *")
      return(-1)
    }
    name.x.temp <- NULL
    alldata <- NULL

    for(idx in 1: block) {
      fileAddress <- paste(outputDirectory, "/variable.quantile_rf_all", "_", idx, ".tif", sep = "")
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

    if(!flag) {
      print("* Data integrity. *")
    } else {
      print("* Data is missing. The file path is as follows: *")
      print(fileAddress)
      print(myaddress1)
    }
    print('RF')
    print(block)
    print(outputDirectory)
  }

  # QRF method function
  Inspect_QRF <- function() {
    if(block == 1) {
      print("* Without parallel computing, the data is always intact *")
      return(-1)
    }
    print('QRF')
    # Init Variables
    name.x.temp <- NULL
    alldata <- NULL
    fragList <- c('05', '50', '95')
    for (i in 1: length(fragList)) {
      for(idx in 1: block) {
        fileAddress <- paste(outputDirectory, "/variable.quantile",fragList[i],"_", idx, ".tif", sep = "")
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

      if(!flag) {
        print("* Data integrity. *")
      } else {
        print("* Data is missing. The file path is as follows: *")
        print(fileAddress)
        print(myaddress1)
      }
    } # for end
  }

  Inspect_Func <- switch(
    model,
    MLR = Inspect_MLR,
    RF = Inspect_RF,
    QRF = Inspect_QRF
  )
  Inspect_Func()
}
