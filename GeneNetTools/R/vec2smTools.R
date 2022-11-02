#' vec2sm.tools
#'
#' Symm matrix to a vector
#'
#' @docType package
#'
#' @author Victor Bernal \email{v.a.bernal.arzola@rug.nl; victor.arturo.bernal@gmail.com}
#'
#' @name vec2sm.tools
#'
#' @import GeneNet
#' @import corpcor
#' @import stats4
#'
#' @references \url{https://doi.org/10.2202/1544-6115.1175}
#'
#' @param vec vector
#'
#' @return symmetric matrix
#'
# #' @examples
#'
#' @export
vec2sm.tools = function(vec){

  p <-  (sqrt(1 + 8 * length(vec)) + 1)/2

  id <- which(x = lower.tri(x = diag(p), diag = F) , arr.ind = T)

  m <- matrix(NA, p,p)

  m[id] <- vec

  m <- (m + t(m))

  return(m)
}

NULL
