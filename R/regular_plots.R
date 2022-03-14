
##############
## Regular plot functions
##
##############



#' Plot 2D spatio-temporal change of soil variables from the outHR file
#'
#' @param dat A dataframe with outHR file data
#' @param var_ Name of the soil variable to plot
#' @param epcouche Soil voxel (cm)
#' @param Boundarycols A vector of two color names defining the extremes of the gradient of values
#' @return 2D spatio-temporal plot of soil variables from the outHR file
#' @examples
#' dat <- outHR_0
#' Plot2D_soilVar(dat, var_='FTSW')
#' Plot2D_soilVar(dat, var_='HRp')
#' Plot2D_soilVar(dat, var_='m_NO3', Boundarycols =c("brown", "green"))
#' Plot2D_soilVar(dat, var_='m_NN4', Boundarycols =c("brown", "green"))
Plot2D_soilVar <- function(dat, var_='FTSW', epcouche=5., Boundarycols = c("blue", "red"))
{

  #plot spatio-temporel ftsw

  #var_ <- 'FTSW'
  #epcouche = 5.
  lscol <- colorRampPalette(Boundarycols)( 101 ) #palette de couleur
  sdat <- split(dat, dat$var)


  DOYs <- sdat[[var_]]$DOY
  vals <- sdat[[var_]][,c(-1,-2)]
  nbcouches <- dim(vals)[2]
  #normalisation
  if (var_ == 'm_NO3' | var_ == 'm_NN4')
  {
    vals <- vals/max(vals) #normalisation / max
  }

  #plot
  plot(-10,-10,ylim=c(-1*epcouche*nbcouches,-0), xlim=c(min(DOYs),max(DOYs)), main=var_, xlab='DOY', ylab='soil depth')

  #draw a sequence of recatngle
  #jour 1
  j <- 1
  for (j in 1:dim(vals)[1])
  {
    ftswj <- as.numeric(vals[j,])
    cols <- rev(lscol[round((1-ftswj)*100,0)+1])

    xleft <- rep(DOYs[j], nbcouches)
    ybottom <- seq(-1*epcouche*nbcouches, -1*epcouche, epcouche)
    xright <- rep(DOYs[j+1], nbcouches)
    ytop <- seq(-1*epcouche*nbcouches, -1*epcouche, epcouche)+epcouche
    rect(xleft, ybottom, xright, ytop, col=cols, border=cols)
  }

}







col100 <- function(valrel100, lscols)
{
  #pour gestion des couleur: vecteur 100
  # fonction pour definir un vecteur de couleur a partir de valeur relative et d'une liste de 101 couleur
  #lscols = vecteur de 100 couleurs
  # valrel100 = position dans ce vecteur (% du max)

  #lscols[rdtrel]#pas bon!
  cols_ <- NULL
  for(i in valrel100)
  {
    cols_ <- rbind(cols_, lscols[i+1])
  }
  cols_ <- as.vector(cols_)
  cols_
}






