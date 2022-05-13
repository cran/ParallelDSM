# The GetPredictorSubset function, which is suitable for cutting
# and integrating spatial data, is a helper function

# to get the a part of predictors: 1/nblock
# this func will output a df.all with all predictorss values and coordinates (x and y)
#=======================================================================================
#' @title calculation function for cutting spatial data (tool function,Not as an open function, only for function calls)
#'
#' @param predictor.name : the name of the predictor variable
#' @param iblock : sequence code of parallel computing
#' @param nblock : number of target blocks (integer)
#' @param fn : The passed value of a global variable
#' @param nr : The passed value of a global variable
#' @param nc : The passed value of a global variable
#' @param resolutions : The passed value of a global variable
#' @param pro : The passed value of a global variable
#' @param from : Which row to start cutting the matrix
#' @param to : Where does the last row of the cut matrix go
#'
#' @importFrom raster cellsFromExtent
#' @importFrom raster projection<-
#' @importFrom raster values<-
#' @importFrom methods as
#'
#' @return Parallel calculation of the cut part of the data box data
#' @export
#'
#' @examples
#' \donttest{
#' GetPredictorSubset("dem",4,10,"covariate",486,777,NULL,NULL,1,10)
#' }
#' @references{
#' Breiman, L. (2001). Random forests. Mach. Learn. 45, 5–32.
#' Meinshausen, N. (2006) "Quantile Regression Forests", Journal of Machine Learning Research 7,
#' 983-999 http://jmlr.csail.mit.edu/papers/v7/
#' }
GetPredictorSubset <- function(predictor.name,iblock,nblock,fn,nr,nc,resolutions,pro,from,to) {
  #The address of the target file that the file points to
  file.directory <- paste(fn,predictor.name, ".tif", sep = "")
  # Read the data and convert it into Raster data
  if(file.exists(file.directory))
  {
    raster.temp <- raster::raster(file.directory)
    # Control the amount of change

    # Determine whether the user passed in the relevant parameters
    if(is.null(from) == TRUE){
      if(iblock!=nblock) {
        nr.block <- nr%/%nblock
      } else  {
        nr.block <- nr - (nr%/%nblock)*(nblock-1)
      }
      row.offset <- (iblock-1)*nr%/%nblock
      e <- raster::extent(x=raster.temp,
                          r1= 0+row.offset,
                          r2= nr.block+row.offset,
                          c1=0,
                          c2=nc)
    }else{
      cha <- to - from + 1;

      # Cut the data in equal quantities
      if(iblock!=nblock) {
        # Using the divisor operator
        nr.block <- cha%/%nblock
      } else  {
        nr.block <- cha - (cha%/%nblock)*(nblock-1)
      }
      # get blocks by rows: this value for the first partition will be 0
      row.offset <- (iblock-1)*cha%/%nblock
      # e represents a restricted object (boundary object) extent(xmin,xmax,ymin,ymax)
      e <- raster::extent(x=raster.temp,
                          r1= from + row.offset,
                          r2= from + nr.block+row.offset,
                          c1=0,
                          c2=nc)
    }

    # Get a grid object
    rowcol <- raster::rowColFromCell(raster.temp, cellsFromExtent(raster.temp,e))
    # rowColFromCell(r, c(5,15)) c(x,y)->cellsFromExtent(r,bb) bb->e->boundary object
    # Is a grid* the value of the block (rectangular area) of the object
    v <- raster::getValuesBlock(raster.temp,
                        row=rowcol[1,1],
                        nrows=(rowcol[nrow(rowcol),1] - rowcol[1,1]),
                        col=rowcol[1,2],
                        ncols=(rowcol[nrow(rowcol),2] - rowcol[1,2]))
    x <- raster::raster(ncol=(rowcol[nrow(rowcol),2] - rowcol[1,2]),
                nrow=(rowcol[nrow(rowcol),1] - rowcol[1,1]),
                xmn=e@xmin,
                xmx=e@xmax-resolutions,
                ymn=e@ymin+resolutions,
                ymx=e@ymax)
    # raster object's set of coordinate values
    projection(x) <- pro
    values(x) <- v
    # Force the object to a given class 'spatial point data frame'
    df <- as(x,"SpatialPointsDataFrame")
    # Then convert the raster object into a data frame
    df <- as.data.frame(df)
    # Then replace the first line of df with the predicted name
    names(df)[1] <- predictor.name
    # Return DataFrame Data
    return (df)
  }
  else
  {
    return(-1)
  }
}
#=======================================================================================
