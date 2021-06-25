
#===================================================================================================================
#' @title For the gap between the predicted value and expected value of the model, the model validates the function
#'
#' @param pred  : Value predicted by the model
#' @param actual  : The real value
#'
#' @return Vector-value after model accuracy verification
#'
#' @importFrom("stats", "cor")
#'
#' @examples
#' test.pred <- c(2,4,5,7,2,4)
#' test.obs <- c(1,2,3,4,5,6)
#' myres <- CVfunction(test.pred,test.real)
#' print(myres)
#'
#'
CVfunction <- function(pred, actual)
{
  e <- pred - actual
  ae <- abs(e)
  e2 <- e^2
  spred <- (pred - mean(pred))^2
  sactual <- (actual - mean(actual))^2
  spa <- (pred - mean(pred)) * (actual - mean(actual))
  df_cv <- data.frame(e,ae,e2,spred,sactual,spa)
  means <- apply(df_cv, MARGIN=2,FUN=mean)

  ME <- means[c(1)]
  MAE <- means[c(2)]
  RMSE <- sqrt(means[c(3)])

  sx <- means[c(4)]
  sy <- means[c(5)]
  sxy <- means[c(6)]
  RPD <- sd(actual) / RMSE
  pearcor <- cor(pred,actual)
  lin <- 2*sxy/(sx + sy + (mean(pred) - mean(actual))^2)# Lins's Concordance Correlation Coefficient

  tempd <- data.frame(pred,actual)
  r2.lm = lm(actual ~ pred, data=tempd)
  r2 <- summary(r2.lm)$r.squared[1]

  va <- c(ME, MAE, RMSE, RPD, pearcor, r2, lin)
  NAME <- c("ME", "MAE", "RMSE", "RPD","Pearson Correlation","Coefficient of Determination (R2)","Lins's Concordance")
  cv <- data.frame(NAME, va)
  return (cv)
}
#=====================================================================================================
