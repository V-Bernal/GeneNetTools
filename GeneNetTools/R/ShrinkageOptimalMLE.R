#' optimal.shrinkage
#'
#' Computes the optimal shrinkage value
#'
#' @docType package
#'
#' @author Victor Bernal \email{v.a.bernal.arzola@rug.nl; victor.arturo.bernal@gmail.com}
#'
#' @name optimal.shrinkage
#'
#' @import GeneNet
#' @import corpcor
#' @import stats4
#'
#' @references \url{https://doi.org/10.2202/1544-6115.1175}
#'
#' @param x data matrix
#'
#' @return optimal shrinkage value \code{lambda}.
#'
#' @examples
#' #---------------------------------
#' # 1- Simulate a 40 x 40 matrix of partial correlations
#' true.pcor <- GeneNet::ggm.simulate.pcor(40)
#'
#' # 2- Simulate random data from true.pcor
#' dat.sim <- GeneNet::ggm.simulate.data(40, true.pcor)
#'
#' # 3- Estimate the optimal shrinkage value
#' lambda <- optimal.shrinkage(x = dat.sim)
#' lambda
#' #---------------------------------
#'
#' @export
#'
optimal.shrinkage <- function ( x = NULL) {

  demean.x <- scale(x)

  p <- ncol(demean.x)
  n <- nrow(demean.x)

  ids <- which(x = lower.tri(x = diag(ncol(demean.x)), diag = F), arr.ind = T)


  func = function(ids){
    w_kij <- demean.x[,ids[1]]*demean.x[,ids[2]]
    w_ij <- t(t(rep(1, n))) %*% mean(w_kij) #t(t(rep(1, n))) %*% colMeans(w_kij)
    var_w <- sum( (w_kij - w_ij)^2 ) * (n/(n-1)^3)
    return(var_w)

  }


  var_w <- apply(X = ids, MARGIN = 1, FUN = func)

  shrinkage <- (sum(var_w)/sum(cor(demean.x)[lower.tri(diag(p))]^2))

  return( max(0,min(1,shrinkage)) )

  #--------------------------------

}

NULL



