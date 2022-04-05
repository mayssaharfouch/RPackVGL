#' This function determines the average of a variable 'var' for all the plants in a canopy
#'
#' @param ltoto A large list of three data frames (seeds 0, 1 and 2)
#' @param lsusm names(ltoto)
#' @param var Output variable of L-egume
#' @param esp Species
#' @param optSD Optional parameter.Returns standard deviation of the sum of the individuals
#'
#' @return The mean value of 'var' for each DOY
moysimval <- function(ltoto, lsusm, var,esp=NA, optSD=F)
{
  res <- vector("list",length(lsusm))
  names(res) <- lsusm
  for (usm in lsusm)
  {

    if (is.na(esp))
    {dat <- ltoto[[usm]]
    } else
    {
      nomcol <- names(ltoto[[usm]])
      idcols <- grepl(esp, nomcol)
      dat <- cbind(ltoto[[usm]][,c(1:2)], ltoto[[usm]][,idcols])
    }

    nbplt <- length(dat)-2
    xplt <- as.matrix(dat[dat$V1==var,3:(3+nbplt-1)], ncol=nbplt)
    xsum <- rowSums(xplt)
    res[[usm]] <- xsum
  }
  if (optSD==F)
  {
    xav <- rowSums(as.data.frame(res))/length(lsusm)
  }else
  {
    xav <- apply(as.data.frame(res),MARGIN=1,stats::sd)
  }

  xav
}


#' This function groups all the mean of the output variables in one data frame
#'
#' @param ltoto A large list of three data frames (seeds 0, 1 and 2)
#' @param lsusm names(ltoto)
#' @param esp Species
#' @param optSD Optional parameter.Returns standard deviation of the sum of the individuals
#'
#' @return A data frame grouping all the outputs variables and their mean values for each DOY
#' @export
#'
build_simmoy <- function(ltoto, lsusm, esp=NA, optSD=F)
{
  if (is.na(esp))
  {dat <- ltoto[[lsusm[1]]]
  } else
  {
    nomcol <- names(ltoto[[lsusm[1]]])
    idcols <- grepl(esp, nomcol)
    dat <- cbind(ltoto[[lsusm[1]]][,c(1:2)], ltoto[[lsusm[1]]][,idcols])
  }

  TT <- dat[dat$V1=='TT',3]
  STEPS <- dat[dat$V1=='TT',2]
  nbplt <- length(dat)-2
  surfsolref <- dat[dat$V1=='pattern',3]

  LAI <- moysimval(ltoto, lsusm, var='SurfPlante', esp, optSD)/ surfsolref
  MSA_esp_canopy <- moysimval(ltoto,lsusm, var='MSaerien', esp, optSD)/ surfsolref
  MSpiv <- moysimval(ltoto,lsusm, var='MS_pivot', esp, optSD)/ surfsolref
  MSracfine <- moysimval(ltoto,lsusm, var='MS_rac_fine', esp, optSD)/ surfsolref
  MSrac <- MSpiv + MSracfine
  NBI_plt <- moysimval(ltoto,lsusm, var='NBI', esp, optSD)/ nbplt
  NBI_plt <- pmax(0, NBI_plt - 0.75)
  NBphyto <- moysimval(ltoto, lsusm, var='NBphyto', esp, optSD)/ surfsolref
  Nbapex <- moysimval(ltoto, lsusm, var='NBapexAct', esp, optSD)/ surfsolref
  NBphyto <- pmax(0,NBphyto - 0.5*Nbapex)
  NBsh_esp_canopy <- moysimval(ltoto, lsusm, var='NBsh', esp, optSD)/ surfsolref

  RDepth <- moysimval(ltoto,lsusm, var='RDepth', esp, optSD)/ nbplt
  Hmax_canopy <- moysimval(ltoto,lsusm, var='Hplante', esp, optSD)/ nbplt
  FTSW <- moysimval(ltoto,lsusm, var='FTSW', esp, optSD)/ nbplt
  NNI <- moysimval(ltoto,lsusm, var='NNI', esp, optSD)/ nbplt
  R_DemandC_Root <- moysimval(ltoto,lsusm, var='R_DemandC_Root', esp, optSD)/ nbplt
  cutNB <- moysimval(ltoto,lsusm, var='cutNB', esp, optSD)/ nbplt
  Npc_aer <- moysimval(ltoto,lsusm, var='Npc_aer', esp, optSD)/ nbplt
  Ndfa <- moysimval(ltoto,lsusm, var='Ndfa', esp, optSD)/ surfsolref
  Epsi <- moysimval(ltoto,lsusm, var='epsi', esp, optSD)
  time <- moysimval(ltoto,lsusm, var='time', esp, optSD)
  aliveB <- moysimval(ltoto,lsusm, var='aliveB', esp, optSD)/ nbplt

  simmoy <- data.frame(nbplt,surfsolref, time,STEPS, TT, NBI_plt, NBphyto, LAI, MSA_esp_canopy, MSpiv, MSracfine, MSrac, RDepth, Hmax_canopy, FTSW, NNI, R_DemandC_Root, cutNB, Npc_aer,Ndfa,Epsi,NBsh_esp_canopy,aliveB)
  simmoy
}
