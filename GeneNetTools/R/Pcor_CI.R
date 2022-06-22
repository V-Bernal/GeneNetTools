#' confint.GGM
#'
#' Confidence intervals for the (partial) correlation with shrinkage.
#'
#' @docType package
#'
#' @author Victor Bernal \email{v.a.bernal.arzola@rug.nl; victor.arturo.bernal@gmail.com}
#'
#' @name confint.GGM
#'
#' @import GeneNet
#' @import corpcor
#' @import stats4
#'
#' @references \url{to be submitted}
#'
#' @param x vector of (partial) correlations
#' @param alpha significance level (two-sided test)
#'
#' @return Confidence intervals of the (partial) correlations with shrinkage \code{lambda}.
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
#'
#' # 4- Compare the estimate with the true value
#' confint.GGM( x = estimated.pcor, alpha = 0.05 )
#'
#' #---------------------------------
#' @export
#'
confint.GGM <- function( x = NULL, alpha = 0.05 ){


  #--------------------------------
  # 0.- Error handling
  if ( missing( x ) ) {

    stop()

  } else {

  #--------------------------------
  # Meta data
  p <- attr(x, 'number.variables')

  n <- attr(x, 'sample.size')

  lambda <- attr(x, 'lambda')

  k <- attr( x, 'degrees.freedom')

  }

  #--------------------------------
  if ( lambda != 0 ){

    #--------------------------------
    # CI Confidence interval for the shrunk
    #---------------------
    zr <- atanh( unlist(x) / (1 - lambda) )
    z1 <- qnorm( ( 1 - alpha*0.5 ) , mean = 0 , sd = 1)

    L <- zr  - ( z1 / sqrt( k - 2 ) ) #n-3
    U <- zr  + ( z1 / sqrt( k - 2 ) ) #n-3

    L.ci <- (1 - lambda) * tanh( L )
    U.ci <- (1 - lambda) * tanh( U )

    CI <- cbind( 'shrunk_pcor'= unlist(x),
                  'low.conf.int' = L.ci,
                  'upper.conf.int' = U.ci )

    CI <- as.matrix(CI)

    colnames(CI) <- c( 'shrunk_pcor',
                 'low.conf.int',
                 'upper.conf.int')

    rownames(CI) <- rownames(x)

    return( CI )
    #--------------------------------

  } else {

    print( 'Error: object not reconigzed' )

    stop()

    }

  }

NULL
