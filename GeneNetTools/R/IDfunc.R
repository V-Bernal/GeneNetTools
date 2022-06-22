#' id.func
#'
#' assign names.
#'
#' @docType package
#'
#' @author Victor Bernal \email{v.a.bernal.arzola@rug.nl; victor.arturo.bernal@gmail.com}
#'
#' @name id.func
#'
#' @import GeneNet
#' @import corpcor
#' @import stats4
#'
#' @references \url{https://doi.org/10.2202/1544-6115.1175}
#'
#' @param x correlation matrix
#'
#' @return names (identifiers) for the correlations
#'
#' @examples
#'
#' @export
id.func = function(x){

  if(is.null(colnames(x))){

    id <- lower.tri(diag(ncol(x)))
    id <- which(id == T, arr.ind = T)
    ctemp <- 1:ncol(x)
    id.names <- paste(ctemp[id[,2]],
                      ctemp[id[,1]], sep = '-')

  } else{

    id <- lower.tri(diag(ncol(x)))
    id <- which(id == T, arr.ind = T)
    id.names <- paste(colnames(x)[id[,2]],
                     colnames(x)[id[,1]], sep = '-')
  }

  return(id.names)
}

NULL
