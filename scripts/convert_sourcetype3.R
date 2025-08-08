#' convert source types of outgroup ST93
#'
#' @param x sequence type
#' @param y source type
#'
#' @return a column with source types of outgr being NA
#' @export
#'
#' @examples
sourceconvert3 <- function(x, y){  
  ifelse(x == 10, y, NA)
}