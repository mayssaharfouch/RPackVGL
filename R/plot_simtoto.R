#' This functions plots the dynamics of the output variables of L-egume
#'
#' @param varsim Simulated variable
#' @param simmoy output of build_simmoy function (mean values)
#' @param simsd output of build_simmoy function (standard deviation)
#' @param name name of the plot
#' @param col color of the plot
#' @param colref color of the second plot of the same graph
#'
#' @return A 2D graph
#' @export
#' @importFrom ggplot2 aes element_text labs theme ylim
#' @importFrom rlang .data
#'
gg_plotsim <- function(varsim, simmoy, simsd, name = " ", col="blue", colref="red")
{
  var_ <- varsim

  min <- 0
  max <- 1.5*max(simmoy[,var_])

  plot_var <- ggplot2::ggplot(NULL, aes(x = simmoy$STEPS)) +
    ggplot2::geom_line(aes(y = simmoy[,var_]), color=col)+
    ggplot2::geom_ribbon(aes(ymin=simmoy[,var_]-simsd[,var_],ymax=simmoy[,var_]+simsd[,var_]),fill=col,alpha=0.2)+
    ggplot2::geom_hline(yintercept=0)+
    ylim(min,max)+
    ggplot2::geom_text(x=1.20*min(simmoy$STEPS), y=0.98*max, size=4, label=name)+
    theme(axis.text.x = element_text(size=6),axis.text.y = element_text(size=6))+
    labs(title = "obs",subtitle = "sim",x = "DOY", y = var_)+
    theme(plot.title=element_text(size=10,color = colref),plot.subtitle = element_text(size=10,color = col))

  plot_var
}


#' This function adds the observed points to plot_var
#'
#' @param plot_var plot of the dynamics of the output variables of L-egume
#' @param var_ Simulated variable
#' @param obsOK observation values of the same DOY range of the simulation values
#' @param obsMerge data frame containing the simulation and observation values
#'
#' @return A 2D graph
#' @export
#'
gg_addplotobs <- function(plot_var, var_, obsOK, obsMerge)
{
  if(var_ %in% names(obsOK)) {plot_var2 <- plot_var + ggplot2::geom_point(aes(obsMerge$DOY, obsMerge[,var_]), fill="red",color="red" , size=2)}
  else {plot_var2 <- plot_var}

  plot_var2
}
