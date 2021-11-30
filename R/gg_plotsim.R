#' Plotting simulation outputs
#'
#' @param varsim name of the output variable
#' @param simmoy data frame containing the simulated data of all the output variables
#' @param onglet excel sheet name, used as a title for the plot
#' @param opt_cut optional parameter, displays the dates when cuts where done when activated as well as the height of the cuts
#' @param opt_Irrig optional parameter, displays dates and quantity of irrigation when activated
#' @param opt_Fert optional parameter, displays dates and quantity of fertilizers used when activated
#'
#' @return A graph showing the dynamics of the selected output variable
#' @export
#'
#' @importFrom ggplot2 aes scale_x_continuous scale_y_continuous geom_line ggtitle geom_text geom_point sec_axis geom_col labs element_line element_blank element_text theme theme_bw
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom data.table last
#'

# Update this function call
utils::globalVariables(c("STEPS", "FertNH4_DOY", "hcut", "FertNO3_DOY", "FertNH4", "Irrig_DOY",
                         "cut_DOY", "Irrig", "FertNO3" ))

gg_plotsim <- function(varsim, simmoy,onglet, opt_cut = NULL, opt_Irrig = NULL, opt_Fert = NULL) {
  var_ <- varsim
  min <- 0
  max <- 1.5*max(simmoy[,var_])
  digits <- floor(log10(max)) + 1 #determine the nb of digits of the max value of 1st y-axis
  # loop to determine the coeff used for the second axis scale
  if (digits == 2 && is.nan(digits) == FALSE){
    coeff = 1
  } else if (digits == 3 && is.nan(digits) == FALSE){
    coeff = 0.1
  } else {
    coeff = 0.01
  }
  ## plot without cuts
  plot_var <- ggplot2::ggplot(NULL, aes(x,y)) +
    geom_line(data = simmoy,aes(x = STEPS, y = simmoy[,var_]), color="blue", size = 0.7)+
    scale_x_continuous(breaks = seq(0, last(simmoy$STEPS), by = 100)) +
    theme_ipsum() +
    geom_text(x=3*min(simmoy$STEPS), y=0.98*max, size=4, label=name)+
    labs(x = "DOY", y = var_)+
    ggtitle(onglet)+
    theme_bw()+
    theme(legend.text=element_text(family = "serif", size=10),
          axis.text.x = element_text(angle = 90, size = 7),
          axis.text.y = element_text(size = 7),
          axis.title.x =element_text(family = "serif", size=10),
          axis.title.y=element_text(family = "serif", size=10),
          axis.line = element_line(colour = "black"), panel.background = element_blank(),
          panel.grid.major = element_blank())
  # Cuts
  if (!is.null(opt_cut)) {
    # eliminate overlayed x-axis labels
    break_init <- sort(c(seq(0, last(simmoy$STEPS), by = 100), cut_DOY$x), decreasing = FALSE)
    bln10 <- diff(break_init) < 10
    breaks_x <- replace(break_init, (c(FALSE, bln10) | c(bln10, FALSE)) & break_init %% 100 == 0, NA)
    plot_var +
      geom_point(data = cut_DOY, aes(x = x, y = hcut/coeff), shape = 4, color = "black", size = 2)+
      scale_x_continuous(breaks = breaks_x) +
      scale_y_continuous(
        # Add a second axis and specify its features
        sec.axis = sec_axis(~.*coeff, name="hcut (m)"))
  }
  # Irrigation
  else if  (!is.null(opt_Irrig)) {
    plot_var +
      geom_col(data = Irrig_DOY, aes(x = x, y = Irrig/coeff), fill = "#8C00E1", alpha = 0.6)+
      scale_y_continuous(
        # Add a second axis and specify its features
        sec.axis = sec_axis(~.*coeff, name="Irrigation"))
  }
  # Fertilization
  else if  (!is.null(opt_Fert)) {
    plot_var +
      geom_col(data = FertNO3_DOY, aes(x = x, y = FertNO3/coeff), fill = "#B969B5", alpha = 0.6)+
      geom_col(data = FertNH4_DOY, aes(x = x, y = FertNH4/coeff), fill = "#E69F00", alpha = 0.6)+
      scale_y_continuous(
        # Add a second axis and specify its features
        sec.axis = sec_axis(~.*coeff, name="Fertilization"))
  }
  else {
    plot_var
  }
}
