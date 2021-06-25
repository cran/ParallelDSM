
#=====================================================================================================
#' @title A function that combines the results of parallel cutting into a single file
#'
#' @param df_dem  : The predicted source file before merging
#' @param f.i.d  : Enter the absolute path to the file
#' @param f.iblock  : The filename prefix of the resulting result
#' @param n.block  : The number of blocks cut is calculated in parallel
#' @param f.o.d  : The absolute output path of the file
#' @param f.suffix  : The suffix for the output of the file
#'
#' @return 1
#' @export MergingTiles
#'
#' @examples
#' \donttest{
#' # you must have a file, which is name "myres"
#' # Merging files, for example:
#' # f.input.directory <- c("e:/test/")
#' # f.input.iblock <- c("sics030_")
#' # n.block <- 100
#' # f.output.directory <- c("e:/test/myoutput")
#' # f.output.suffix <- c("sics030_together.tif")
#' # Naming rules: file.name.directory + file.name.iblock + ".tif"
#'
#' rmap_dem <- raster("E:/test/dem.tif")
#' spdf_dem <- as(rmap_dem,"SpatialPointsDataFrame")
#' df_dem <- as.data.frame(spdf_dem)
#'
#' # mergeing results together
#' n.block <- 100
#' f.i.d <- c("E:/test/mapping/")
#' f.o.d <- c("E:/test/mapping_merge/")
#' f.iblock <- c("mlr.ak05.")
#' f.suffix <- c("mlr.ak05.tif")
#' MergingTiles(df_dem, f.i.d, f.iblock, n.block, f.o.d, f.suffix)
#'}
#'
MergingTiles <- function(df_dem, f.i.d, f.iblock, n.block, f.o.d, f.suffix)
{

  df.output <- df_dem; names(df.output)[1] <- c("res");  df.output$res <- -999
  pixel.from <- 0;  pixel.to <- 0

  # Naming rules for the ith block
  f.output.suffix <- paste(f.o.d, f.suffix, sep = "")
  f.prefix <- paste(f.i.d, f.iblock, sep = "")

  # Output relative path names
  print("=============================** ParallelDSM information **=============================")
  print("Soil mapping is successfully performed, and the predicted results are stored as follows:")
  myprefix <- getwd()
  allAddress <- paste(myprefix,'/',substring(f.prefix,1,nchar(f.prefix) - 1),'.tif',sep = "")
  print(allAddress)

  for(i in 1:n.block)
  {
    in.file.name <- paste(f.prefix, i, ".tif", sep = "")

    #Reading a tile
    rmap_i <- raster(in.file.name); spdf_i <- as(rmap_i,"SpatialPointsDataFrame")
    df.i <- as.data.frame(spdf_i)

    #Update the row numbers of from and to
    if(i==1)  {
      pixel.from <- 1
      pixel.to <- nrow(df.i)
    }  else  if(i>1)  {
      pixel.from <- pixel.to+1;  pixel.to <- pixel.from+nrow(df.i)-1
    }

    df.output[ c(pixel.from:pixel.to), 1] <- df.i[,1]
  }

  df.output <- as.data.frame(df.output); coordinates(df.output) <- c("x","y"); gridded(df.output) <- TRUE

  writeGDAL(dataset = df.output["res"], fname = f.suffix, drivername = "GTiff", type = "Float32")
  return (1)
}
#=====================================================================================================






