
##############
## fonctions criteres stat eval model
## as in Coucheney et al. 2015 (https://doi.org/10.1016/j.envsoft.2014.11.024)
##############



#' Calculate the Root Mean Square Error (RMSE) of predicted model outputs
#'
#' @param O An observed vector of values
#' @param P A vector of values predicted by a model
#' @return The Root Mean Square Error between \code{O} and \code{P}
#' @examples
#' obs <- c(1, 2, 3, 5)
#' sim <- c(0.9, 2.2, 3.1, 5)
#' rmse(obs, sim)
rmse <- function(O,P){
  sqrt(mean((O-P)**2,na.rm=T))
}



#' Calculate the systematic component of Root Mean Square Error (RMSEs) of predicted model outputs
#'
#' @param O An observed vector of values
#' @param P A vector of values predicted by a model
#' @return The systematic Root Mean Square Error between \code{O} and \code{P}
#' @examples
#' obs <- c(1, 2, 3, 5)
#' sim <- c(0.9, 2.2, 3.1, 5)
#' rmsesCoucheney(obs, sim)
#' @family rmse
rmsesCoucheney=function(O,P){
  Nb.O=length(O)
  sum.error.egs<-0
  resreg<-lm(P~O)
  Preg<-fitted.values(resreg)
  for(k in 1:length(P)){
    error.egs<-(Preg[k]-O[k])^2 # difference entre predit par la regression lineaire et observe
    sum.error.egs<-sum.error.egs+error.egs
  }
  sum.error.egs=unname(sum.error.egs) #permet de rendre le vecteur sans-nom (sans quoi il est avec-nom "1" pour une raison qui m'echappe...)
  return(((1/Nb.O)*sum.error.egs)^0.5)
}






#' Calculate the unsystematic component of Root Mean Square Error (RMSEu) of predicted model outputs
#'
#' @param O An observed vector of values
#' @param P A vector of values predicted by a model
#' @return The unsystematic Root Mean Square Error between \code{O} and \code{P}
#' @examples
#' obs <- c(1, 2, 3, 5)
#' sim <- c(0.9, 2.2, 3.1, 5)
#' rmseuCoucheney(obs, sim)
#' @family rmse
rmseuCoucheney=function(O,P){
  Nb.O=length(O)
  sum.error.egu<-0
  resreg<-lm(P~O)
  Preg<-fitted.values(resreg)
  for(k in 1:length(P)){
    error.egu<-(P[k]-Preg[k])^2 # difference entre predit par la regression lineaire et simule
    sum.error.egu<-sum.error.egu+error.egu
  }
  sum.error.egu=unname(sum.error.egu) #permet de rendre le vecteur sans-nom (sans quoi il est avec-nom "1" pour une raison qui m'echappe...)
  return(((1/Nb.O)*sum.error.egu)^0.5)
}





#' Calculate the relative Root Mean Square Error (rRMSE) of predicted model outputs
#'
#' @param O An observed vector of values
#' @param P A vector of values predicted by a model
#' @return The relative Root Mean Square Error between \code{O} and \code{P} to the mean of observed values
#' @examples
#' obs <- c(1, 2, 3, 5)
#' sim <- c(0.9, 2.2, 3.1, 5)
#' rrmseCoucheney(obs, sim)
rrmseCoucheney=function(O,P){
  O.mean=mean(O)
  return(100*rmse(O,P)/O.mean)
}




#' Calculate the partial systematic component of Root Mean Square Error (pRMSEs) of predicted model outputs
#'
#' @param O An observed vector of values
#' @param P A vector of values predicted by a model
#' @return The partial systematic Root Mean Square Error between \code{O} and \code{P}
#' @examples
#' obs <- c(1, 2, 3, 5)
#' sim <- c(0.9, 2.2, 3.1, 5)
#' rmse_val <- rmse(obs, sim)
#' rmses_val <- rmsesCoucheney(obs, sim)
#' rmseu_val <- rmseuCoucheney(obs, sim)
#' pRMSEs(rmse_val, rmses_val)
#' pRMSEu(rmse_val, rmseu_val)
#' pRMSEs(rmse_val, rmses_val) + pRMSEu(rmse_val, rmseu_val)
pRMSEs=function(rmse_val, rmses_val)
{
  rmses_val**2 / rmse_val**2
}



#' Calculate the partial unsystematic component of Root Mean Square Error (pRMSEu) of predicted model outputs
#'
#' @param O An observed vector of values
#' @param P A vector of values predicted by a model
#' @return The partial unsystematic Root Mean Square Error between \code{O} and \code{P}
#' @examples
#' obs <- c(1, 2, 3, 5)
#' sim <- c(0.9, 2.2, 3.1, 5)
#' rmse_val <- rmse(obs, sim)
#' rmses_val <- rmsesCoucheney(obs, sim)
#' rmseu_val <- rmseuCoucheney(obs, sim)
#' pRMSEs(rmse_val, rmses_val)
#' pRMSEu(rmse_val, rmseu_val)
#' pRMSEs + pRMSEu
pRMSEu=function(rmse_val, rmseu_val)
{
  rmseu_val**2 / rmse_val**2
}




#' Calculate the modelling Efficicency (EF) of predicted model outputs
#'
#' @param O An observed vector of values
#' @param P A vector of values predicted by a model
#' @return The modelling Efficicency of \code{P} with respect to \code{O}
#' @examples
#' obs <- c(1, 2, 3, 5)
#' sim <- c(0.9, 2.2, 3.1, 5)
#' efficiencyCoucheney(obs, sim)
efficiencyCoucheney=function(O,P){
  sum.error<-0
  sum.mean.dev<-0
  O.mean=mean(O)
  for(k in 1:length(P)){
    dif<-(O[k]-P[k]) # difference entre valeur simulee et observee
    error<-dif^2
    sum.error<-sum.error+error
    mean.dev<-(O[k]-O.mean)^2 # deviation des observes / moyenne des observes
    sum.mean.dev<-sum.mean.dev+mean.dev
  }
  EF=1-sum.error/sum.mean.dev # efficience du modele ; max 1; the closer to 1 the better
  EF[!is.finite(EF)]<-NA
  return(EF)
}

