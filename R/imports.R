#' @importFrom Rcpp evalCpp
NULL

#' @useDynLib sarabc
NULL

# Importing from the Matrix pkg ------------------------------------------------

#' @import methods
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Matrix Matrix rowSums colSums
#' @importMethodsFrom Matrix t
NULL

#' @importFrom mvtnorm rmvnorm
#' @importFrom stats runif
NULL

#' @importFrom parallel stopCluster clusterMap clusterSplit
NULL