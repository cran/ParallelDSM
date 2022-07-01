
#=====================================================================================================
#' @title Black box test function to test whether R package was installed successfully
#'
#' @param myflag  : Mark the successful installation of the R package
#'
#' @return NULL
#' @export
#'
#' @examples
#' flag <- ""
#' smalltesttoy(flag)
#'
#'
smalltesttoy <- function(myflag){
  if(!is.null(myflag)){
    print("* * * *R package has been installed successfully, through the black box testing, welcome to use the R package: ParallelDSM!! * * * *")
  }else{
    print("* * * *Input error or black box test fails, please retest or check whether the package is installed successfully!!* * * *")
  }
}
#=====================================================================================================
