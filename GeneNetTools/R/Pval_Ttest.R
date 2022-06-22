#' ttest.shrunk
#'
#' P-values for the partial correlations with shrinkage.
#' Test whether the partial correlation is different from zero.
#'
#' @docType package
#'
#' @author Victor Bernal \email{v.a.bernal.arzola@rug.nl; victor.arturo.bernal@gmail.com}
#'
#' @name ttest.shrunk
#'
#' @import GeneNet
#' @import corpcor
#' @import stats4
#'
#' @references \url{to be submitted}
#'
#' @param x matrix of partial correlations
#' @param lambda (optional) shrinkage value (see \code{\link[optimal.shrinkage]{par}})
#'
#' @return vector of p-values for the partial correlations \code{P} with shrinkage \code{lambda}.
#'
#'
#' @examples
#' #---------------------------------
#' # 1- Simulate a 40 x 40 matrix of partial correlations
#' true.pcor <- GeneNet::ggm.simulate.pcor(40)
#'
#' # 2- Simulate random data from true.pcor
#' dat.sim <- GeneNet::ggm.simulate.data(40, true.pcor)
#'
#' # 3- Estimate the partial correlations with shrinkage
#' estimated.pcor <- GGM.shrunk(dat.sim)
#' pval  <- ttest.shrunk (x = estimated.pcor)
#'
#'
#' @export

ttest.shrunk <- function( x = NULL, lambda = NULL){


  #-----------------------
  # Meta data
  p <- attr(x, 'number.variables')

  n <- attr(x, 'sample.size')

  lambda <- attr(x, 'lambda')

  k <- attr( x, 'degrees.freedom')

  #--------------------------------------------
  # Estimate p values using the shrunk t-test

  # Rescale pcor
  rr <-  unlist(x) / (1- lambda)

  # T-test
  tt <- rr * sqrt( (k-1)/(1 - rr^2 ) )

  # p-vals student
  pval.shrunk <- 2*pt(q = -abs(tt),
                  df = (k - 1),
                  lower.tail = TRUE, log.p = FALSE)

  pval.shrunk <- as.matrix( pval.shrunk)

  #------------------------------------
  # Add meta data to the pval object
  attr( pval.shrunk, 'lambda') <- lambda

  attr( pval.shrunk, 'sample.size') <- n

  attr( pval.shrunk, 'degrees. freedom') <- k

  attr( pval.shrunk, 'number.variables') <- p

  rownames(pval.shrunk) <- rownames(x)

  colnames(pval.shrunk) <- 'p-values'

  return( pval.shrunk )

}

NULL

