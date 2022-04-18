#' This function returns the systematic error of the model
#'
#' @param O observed values
#' @param P simulated values
#'
#' @return systematic error RMSEs
#' @export
#'
rmsesCoucheney=function(O,P){
  Nb.O=length(O)
  sum.error.egs<-0
  resreg <- stats::lm(P~O)
  Preg <- stats::fitted.values(resreg)
  for(k in 1:length(P)){
    error.egs<-(Preg[k]-O[k])^2
    sum.error.egs<-sum.error.egs+error.egs
  }
  sum.error.egs=unname(sum.error.egs)
  return(((1/Nb.O)*sum.error.egs)^0.5)
}


#' This function returns the unsystematic error of the model
#'
#' @param O observed values
#' @param P simulated values
#'
#' @return unsystematic error RMSEu
#' @export
#'
rmseuCoucheney=function(O,P){
  Nb.O=length(O)
  sum.error.egu <- 0
  resreg <- stats::lm(P~O)
  Preg <- stats::fitted.values(resreg)
  for(k in 1:length(P)){
    error.egu <- (P[k]-Preg[k])^2
    sum.error.egu <- sum.error.egu+error.egu
  }
  sum.error.egu=unname(sum.error.egu)
  return(((1/Nb.O)*sum.error.egu)^0.5)
}


#' This function returns the relative mean square error in %
#'
#' @param O observed values
#' @param P simulated values
#'
#' @return Relative Mean Square Error (RMSE)
#' @export
#'
rrmseCoucheney=function(O,P){
  O.mean=mean(O)
  return(100*Metrics::rmse(O,P)/O.mean)
}


#' This function returns the percentage of the systematic error of the model
#'
#' @param rmse_ Relative Mean Square Error
#' @param rmses_ systematic error of the model
#'
#' @return percentage of the systematic error
#' @export
#'
pRMSEs=function(rmse_, rmses_)
{
  rmses_**2 / rmse_**2
}


#' This function returns the percentage of the systematic error of the model
#'
#' @param rmse_ Relative Mean Square Error
#' @param rmseu_ unsystematic error of the model
#'
#' @return percentage of the unsystematic error
#' @export
#'
pRMSEu=function(rmse_, rmseu_)
{
  rmseu_**2 / rmse_**2
}


#' This function returns the modelling Efficiency
#'
#' @param O observed values
#' @param P simulated values
#'
#' @return Modelling Efficiency (EF)
#' @export
#'
efficiencyCoucheney=function(O,P){
  sum.error<-0
  sum.mean.dev<-0
  O.mean=mean(O)
  for(k in 1:length(P)){
    dif<-(O[k]-P[k])
    error<-dif^2
    sum.error<-sum.error+error
    mean.dev<-(O[k]-O.mean)^2
    sum.mean.dev<-sum.mean.dev+mean.dev
  }
  EF=1-sum.error/sum.mean.dev
  EF[!is.finite(EF)]<-NA
  return(EF)
}
