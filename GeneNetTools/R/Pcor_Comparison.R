#' compare.GGM
#'
#' Test whether two (partial) correlations are equal.
#'
#' @docType package
#'
#' @author Victor Bernal \email{v.a.bernal.arzola@rug.nl; victor.arturo.bernal@gmail.com}
#'
#' @name compare.GGM
#'
#' @import GeneNet
#' @import corpcor
#' @import stats4
#'
#' @references \url{to be submitted}
#'
#' @param x1 data matrix
#' @param x2 shrinkage value \code{lambda}
#'
#' @return z-score that compares two (partial) correlations accounting for the shrinkage \code{lambda}.
#'
#'
#' @examples
#' #---------------------------------
#' # 1- Simulate a 40 x 40 matrix of partial correlations
#' true.pcor <- GeneNet::ggm.simulate.pcor(40)
#'
#' # 2- Simulate random data from true.pcor
#' dat.sim1 <- GeneNet::ggm.simulate.data(40, true.pcor)
#' dat.sim2 <- GeneNet::ggm.simulate.data(40, true.pcor)
#'
#' # 3- Estimate the partial correlations with shrinkage
#' estimated.pcor1 <- GGM.shrunk(dat.sim1)
#' estimated.pcor2 <- GGM.shrunk(dat.sim2)
#'
#' # 4- Compare
#' z <- compare.GGM(estimated.pcor1, estimated.pcor2)
#'
#' # 5- Compare the estimate with the true value
#' sum( abs(z) >=1.96)
#' plot( abs(z), pch=20, ylim=c(0,4) , cex=0.5)
#' abline(h=1.96)
#' #---------------------------------
#'
#' @export
compare.GGM <- function ( x1 = NULL , x2 = NULL) {

  #------------------------------
  p1 <- attr( x1, 'number.variables')
  p2 <- attr( x2, 'number.variables')

  n1 <- attr(x1, 'sample.size')
  n2 <- attr(x2, 'sample.size')

  lambda1 <- attr(x1, 'lambda')
  lambda2 <- attr(x2, 'lambda')

  k1 <- attr( x1, 'degrees.freedom')
  k2 <- attr( x2, 'degrees.freedom')

  #------------------------------
  # Rescale
  new.pcor1 <- unlist(x1)/(1- lambda1)
  new.pcor2 <- unlist(x2) /(1-lambda2)

  # Fisher transform is normal
  FisherT11 <- atanh( new.pcor1 )
  FisherT22 <- atanh( new.pcor2 )

  # z score for the difference
  diff.effect <- ( FisherT11 - FisherT22 )

  diff.SE <- sqrt( (1 /( k1 - 2 )) + (1 / ( k2 - 2 ) ))

  Zobserved <- diff.effect / diff.SE

  # P-values
  pval <- 2*pnorm(-abs( Zobserved ),
                mean = 0, sd = 1,
                lower.tail = TRUE, log.p = FALSE)

  # Return
  Zobserved <- c(unlist(Zobserved))

  Zobserved <- as.matrix( cbind(Zobserved, pval))

  colnames(Zobserved) <-  c('z-score', 'p-value')



  #
  return( Zobserved )

}

NULL
