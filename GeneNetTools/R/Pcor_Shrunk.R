#' GGM.shrunk
#'
#' Computes the partial correlations with shrinkage.
#'
#' @docType package
#'
#' @author Victor Bernal \email{v.a.bernal.arzola@rug.nl; victor.arturo.bernal@gmail.com}
#'
#' @name GGM.shrunk
#'
#' @import GeneNet
#' @import corpcor
#' @import stats4
#'
#' @references \url{https://doi.org/10.2202/1544-6115.1175; to be submitted}
#'
#' @param x data matrix
#' @param lambda shrinkage value (see \code{\link[optimal.shrinkage]{par}})
#' @param verbose print details
#'
#' @return The matrix of partial correlations \code{P} with shrinkage \code{lambda}.
#'
#' @examples
#' #---------------------------------
#' @export

GGM.shrunk <- function( x = NULL, lambda, verbose = TRUE ) {

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


      } else { # user gives lambda

        lambda <- lambda
      }

  #     if ( lambda < 0 )
  #
  #       lambda <- 0
  #
  #     if ( lambda > 1 )
  #
  #       lambda <- 1
  #
  #     if ( verbose )
  #
  #       cat( paste( 'Shrinkage lambda from user = ', round( lambda, 3 ), '\n'))
  #       #lambda.estimated = FALSE
  #   }
  #
  # if ( lambda >= 1 ){
  #
  #   cat( paste( 'Singular case: Shrinkage lambda =', round(lambda, 3), '\n'))
  #   cat( paste( '=> The partial correlations are zero, and their p-values are 1'))
  #
  #   GGM <- matrix(data = 0, nrow = p, ncol = p)
  #
  #   attr( GGM, 'lambda') <- lambda
  #   attr( GGM, 'sample.size') <- n
  #
  #   colnames(GGM) <- rownames(GGM)
  #
  #   return( GGM )
  #
  #   stop()
  #
  #   } else {
  #
  #     cat( 'NON-Singular case')
  #
  #   }

  #--------------------------------
  # 2.- Estimate the partial correlations with shrinkage lambda

  # Estimate the Correlation matrix
  R = cor(x)

  #--------------------------------
  # Use the Ledoit-Wolf shrinkage on the correlation matrix
  R_lambda = ( 1 - lambda )*( R ) + ( lambda ) * diag( diag( R ) ) # shrinkage

  #--------------------------------
  # 3.- Compute the Precision matrix (Omega)
  D = chol( x = R_lambda ) # Cholesky decomposition

  omega.chol = solve( D, diag(ncol(x)) ) %*% solve( t(D), diag(ncol(x)) )# precision matrix

  #--------------------------------
  # 4.- Compute the Partial correlations
  chol.PCOR = -cov2cor( omega.chol )

  diag( chol.PCOR ) = 1

  GGM <- chol.PCOR

  GGM <- GGM[ lower.tri( GGM , diag = F ) ]# symm matrix to long vector

  GGM <- as.matrix(GGM)

  #--------------------------------
  # degrees of freedom k using MLE
  k <- sapply(X = c(1:100), FUN = function(X){ k.shrunk(p = p,
                                                        n = n,
                                                        lambda = lambda)
  }
  )

  k <- mean(k)

  cat(paste('degrees of freedom k =' , k))

  #--------------------------------
  # 4.1 - Add meta-data

  attr( GGM, 'lambda' ) <- lambda

  attr( GGM, 'sample.size' ) <- n

  attr( GGM, 'number.variables' ) <- p

  attr( GGM, 'degrees.freedom' ) <- k

  rownames(GGM) <- id.func(x)

  colnames(GGM) <- 'shrunk_pcors'

  return( GGM )

  }

NULL
