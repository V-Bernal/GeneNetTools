#' corr.shrunk
#'
#' Computes the (Pearson's) correlation matrix with shrinkage.
#'
#' @docType package
#'
#' @author Victor Bernal \email{v.a.bernal.arzola@rug.nl; victor.arturo.bernal@gmail.com}
#'
#' @name corr.shrunk
#'
#' @import GeneNet
#' @import corpcor
#' @import stats4
#'
#' @references \url{https://doi.org/10.2202/1544-6115.1175}
#'
#' @param x data matrix
#' @param lambda shrinkage value
#' @param verbose print details
#'
#' @return The matrix of (Pearson's) correlations \code{R} with shrinkage.
#'
#' @examples
#'
#' @export
corr.shrunk <- function( x = NULL, lambda, verbose = TRUE ) {

  #================================
  # 0.- Error handling
  if ( missing( x ) ) {

    stop() }

  else {

    p <- ncol( x )
    n <- nrow( x )

    cat( paste( 'Number of samples = ', n , '\n') )
    cat( paste( 'Number of variables = ', p , '\n') )

   }

  #================================
  # 1.- Estimate optimal shrinkage lambda

  if ( missing( lambda ) ) {

    lambda <- optimal.shrinkage(x)
    cat( paste( 'Shrinkage (lambda) optimal = ', round( lambda, 3 ), '\n'))


  } else { # user gives lambda

    if ( lambda < 0 )

      lambda <- 0

    if ( lambda > 1 )

      lambda <- 1

    if ( verbose )

      cat( paste( 'Shrinkage (lambda) from the user = ', round( lambda, 3 ), '\n'))

    }

  if ( lambda >= 1 ){

    cat( paste( 'Singular case: Shrinkage lambda =', round(lambda, 3), '\n'))
    cat( paste( '=> The correlations are zero, and their p-values are 1'))

    COR <- matrix(data = 0, nrow = p, ncol = p)

    attr( COR, 'lambda') <- lambda
    attr( COR, 'sample.size') <- n

    colnames(COR) <- rownames(COR)

    return( COR )

    stop()

  } else {

    cat( 'NON-Singular case')

  }

  #--------------------------------
  # 2.- Estimate the partial correlations with shrinkage lambda

  #----------------------
  # Estimate the Correlation matrix
  R = cor(x)

  #--------------------------------
  # Use the Ledoit-Wolf shrinkage on the correlation matrix
  R_lambda = ( 1 - lambda )*( R ) + ( lambda ) * diag( diag( R ) ) # shrinkage

  R_lambda = R_lambda[lower.tri(R_lambda, diag = F)]

  COR <- as.matrix( R_lambda )

  #--------------------------------
  # 4.1 - Add meta-data
  attr( COR, 'lambda' ) <- lambda

  attr( COR, 'sample.size' ) <- n

  attr( COR, 'number.variables' ) <- p

  attr( COR, 'degrees.freedom' ) <- (n - 1)

  rownames(COR) <- id.func(x)

  colnames(COR) <- 'shrunk_cors'

  return( COR )

}

NULL
