
#=====================================================================================================
#' @title A function that checks the parallel computation for missing data
#'
#' @param myblock  : the number of blocks for data cutting
#'
#' @return NULL
#' @export
#'
#' @importFrom raster raster
#' @importFrom methods as
#'
#' @examples
#' \donttest{
#' iblock = 10
#' InsepectionVariable(myblock = iblock)
#' }
#'
#'
#'
InsepectionVariable <- function(myblock){
  if(myblock == 1){
    print("* Without parallel computing, the data is always intact *")
    return (-1)
  }
  # ====================================================
  # Stored as a temporary variable
  name.x.temp <- NULL
  alldata <- NULL
  alldata1 <- NULL
  mylist <- c("05","50","95")
  for (i in 1 : length(mylist)){
    # Create the merged data set
    for(idx in 1 : myblock){
      myaddress <- paste("output/variable.quantile",mylist[i],"_", idx, ".tif", sep = "")
      myraster  <- raster(myaddress,header=FALSE)
      dfraster  <- as(myraster,"SpatialPointsDataFrame")
      allraster <- data.frame(dfraster)
      #plot(allraster)
      name.x.temp <- data.frame(allraster)
      if(is.null(alldata)){
        # determines whether the value is null
        alldata = name.x.temp
      }else{
        names(alldata) = names(name.x.temp)
        alldata = rbind(alldata,name.x.temp)
      }
    }
    # Read the total data set
    myaddress1 <- paste("outputall/variable.quantile",mylist[i],"_all.tif",sep="")
    alldata1 <- raster(myaddress1,header=FALSE)
    alldata1 <- as(alldata1,"SpatialPointsDataFrame")
    alldata1 <- data.frame(alldata1)
    # ========================= check =========================
    # Deep random number, loop traversal comparison operation
    mynrows <- nrow(alldata)
    xnums <-sample(1:mynrows,floor(mynrows * 0.3))
    myflag <- TRUE
    for(i in 1:length(xnums)){
      for(j in 1:3){
        if(alldata[xnums[i],j] != alldata[xnums[i],j]){
          myflag <- FALSE
        }
      }
    }
    if(myflag){
      print("* Data integrity. *")
    }else{
      print("* Data is missing. The file path is as follows: *")
      print(myaddress)
      print(myaddress1)
    }
    # ========================= check =========================
  }
  # ====================================================
}
#=====================================================================================================
