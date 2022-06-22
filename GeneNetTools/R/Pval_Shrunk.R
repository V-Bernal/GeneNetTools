#' pval.shrunk
#'
#' P-values for the partial correlations with shrinkage.
#' Test whether the partial correlation is different from zero.
#'
#'
#' @docType package
#'
#' @author Victor Bernal \email{v.a.bernal.arzola@rug.nl; victor.arturo.bernal@gmail.com}
#'
#' @name pval.shrunk
#'
#' @import GeneNet
#' @import corpcor
#' @import stats4
#'
#' @references \url{https://doi.org/10.1093/bioinformatics/btz357}
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
#'
#' # 4.- P-values
#' pval <- pval.shrunk( x = estimated.pcor)
#'
#' #---------------------------------
#'
#' @export

pval.shrunk <- function( x = NULL, lambda = NULL){

  #-----------------------
  # Meta data
  p <- attr(x, 'number.variables')

  n <- attr(x, 'sample.size')

  lambda <- attr(x, 'lambda')

  k <- attr( x, 'degrees.freedom')

  # -----------------------------------
  if ( lambda >= 1 ){

    #------------------------------------
    cat( paste( "Singular case: Shrinkage lambda =", round(lambda, 3), "\n"))

    cat( paste( "=> The GGM is a zero matrix, and p-values are 1"))

    pval.shrunk <- rep(1, 0.5*p*(p-1))

    #------------------------------------
    # Add meta data to the GGM object
    attr( pval.shrunk, "lambda") <- lambda

    attr( pval.shrunk, "sample.size") <- n

    attr( pval.shrunk, "degrees.freedom") <- NA

    attr( pval.shrunk, 'number.variables') <- p

    #-------------------------------------------
    return( pval.shrunk  )

    stop()

    } else {

    cat( "NON-Singular case")

  }


  #--------------------------------------------
  # Estimate p values using the shrunk prob  density
  r <- unlist( x )

  if( lambda != 0 ){

   density.shrunk <- function(r) {(  ((1)^2-(r/(1-lambda))^2) ^((k-3)*0.5)  )/( beta(0.5, 0.5*(k-1))*(1-lambda)) }

   pval.shrunk <- matrix(Inf, length(r), 1)

   pval.shrunk <- sapply(X = r, FUN = function(X){

                         int <- integrate( density.shrunk ,
                                           lower = -( 1 - lambda ),
                                           upper = - abs( X )   )

                         return( 2*( int$value ) )

                         }
   )
  }

  pval.shrunk <- as.matrix( pval.shrunk )

  #------------------------------------
  # Add meta data to the pval object
  attr( pval.shrunk, 'lambda') <- lambda

  attr( pval.shrunk, 'sample.size') <- n

  attr( pval.shrunk, 'number.variables') <- p

  rownames(pval.shrunk) <- rownames(x)

  colnames(pval.shrunk) <- 'p-values'

return( pval.shrunk )

}


NULL

