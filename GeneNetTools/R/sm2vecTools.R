#' sm2vec.tools
#'
#' Symm matrix to a vector
#'
#' @docType package
#'
#' @author Victor Bernal \email{v.a.bernal.arzola@rug.nl; victor.arturo.bernal@gmail.com}
#'
#' @name sm2vec.tools
#'
#' @import GeneNet
#' @import corpcor
#' @import stats4
#'
#' @references \url{https://doi.org/10.2202/1544-6115.1175}
#'
#' @param m symmetric matrix
#' @param diag include the main diagonal
#'
#' @return vector
#'
# #' @examples
#'
sm2vec.tools = function (m, diag = FALSE) {

  id <- lower.tri(x = m, diag = diag)

  vec <- m[id]

  return(vec)
}

NULL
