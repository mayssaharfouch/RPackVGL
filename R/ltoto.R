#' This function builds dtoto
#'
#' @param ls_toto list of all toto files of the same experiment
#'
#' @return A dataframe. Each row of this list correspond to an IDusm
#' @export
#'
build_dtoto <- function(ls_toto)
{
  cols_ <- strsplit(ls_toto, '_')
  test_long <- as.numeric(lapply(cols_, length))
  dtoto <- as.data.frame(t(as.data.frame(cols_[test_long==11])))
  row.names(dtoto) <- 1: length(dtoto[,1])
  dtoto <- dtoto[,c(2,3,4,5,6,7,8,10)]
  names(dtoto) <- c('usm','lsystem','mix','damier','scenario','Mng', 'seed','sd')
  dtoto$name <- ls_toto[test_long==11]
  dtoto$seed <- as.numeric(as.character(dtoto$seed))
  dtoto$scenario <- substr(as.character(dtoto$scenario), 9, nchar(as.character(dtoto$scenario)))
  dtoto$keysc <- paste(dtoto$scenario, dtoto$mix, dtoto$Mng, dtoto$sd)
  dtoto
}

#' This function takes the group of USMs of a specific Key, reads the content of the simulation file "toto" of each USM and stocks the output in a list
#'
#' @param ls_toto list of toto files belonging to the same key
#'
#' @return A list. Each element of this list correspond to a single USM
#' @export
#'
#' @importFrom utils read.table
#'
read_ltoto <- function(ls_toto)
{
  ltoto <- vector('list', length(ls_toto))
  names(ltoto) <- ls_toto

  for (i in 1:length(ls_toto))
  {
    name <- ls_toto[i]
    ltoto[[name]] <- read.table(name, header=T, sep=';')
  }
  ltoto
}
