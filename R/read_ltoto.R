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
