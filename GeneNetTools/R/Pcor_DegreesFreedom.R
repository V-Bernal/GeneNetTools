#' k.shrunk
#'
#' Degrees of freedom k via Maximum Likelihood Estimation
#'
#' @docType package
#'
#' @author Victor Bernal \email{v.a.bernal.arzola@rug.nl; victor.arturo.bernal@gmail.com}
#'
#' @name k.shrunk
#'
#' @import GeneNet
#' @import corpcor
#' @import stats4
#'
#' @references \url{https://doi.org/10.2202/1544-6115.1175}
#'
#' @param p number of variables
#' @param n number of samples
#' @param lambda (optional) shrinkage value (see \code{\link[lambda]{par}})
#'
#' @return degrees of freedom with shrinkage \code{k}
#'
#' @examples
#' #---------------------------------
#' # 1- Simulate a 40 x 40 matrix of partial correlations
#' true.pcor <- GeneNet::ggm.simulate.pcor(40)
#'
#' # 2- Simulate random data from true.pcor
#' dat.sim <- GeneNet::ggm.simulate.data(40, true.pcor)
#'
#' # 3- Estimate the degrees of freedom
#' k.shrunk(p = ncol(dat.sim) , n = nrow(dat.sim) , lambda = 0.1)
#' #---------------------------------
#'
#' @export
k.shrunk <- function( p , n , lambda){

  #----------------------------
  # 1.- Simulate data under the null hypothesis: zero partial correlation.
  sim.pcor.S <- diag(p)

  sim.data.S <- GeneNet::ggm.simulate.data(sample.size =  n, pcor = sim.pcor.S) # data from identity GGM (sim.pcor.S)

  #----------------------
  # Estimate the Correlation matrix
  R <- cor(sim.data.S)

  R_lambda <- ( 1 - lambda )*( R ) + ( lambda ) * diag( diag( R ) ) # shrinkage

  D <- chol( x = R_lambda ) # Cholesky decomposition

  omega.chol <- solve( D, diag(ncol(sim.data.S)) ) %*% solve( t(D), diag(ncol(sim.data.S)) )# precision matrix

  chol.PCOR <- -cov2cor( omega.chol )

  #----------------------
  GGM.S <- chol.PCOR #GGM.shrunk( x= sim.data.S, lambda = lambda, verbose = F ) #corpcor::pcor.shrink(x = sim.data.S, lambda = lambda)

  r.S <- sm2vec(GGM.S) # vectorize

  #----------------------------
  # 2.- Compute the negative log-likelihood [Bernal et al.]
  nlogL.shrunk <- function(k) {

    log.density.shrunk <- function(x) {  ((k-3)*0.5)*log(1-(x/(1-lambda))^2) - log( beta(0.5, 0.5*(k-1)))-log(1-lambda) }

    log.f <- log.density.shrunk(r.S)

    return( -sum( log.f ) )
  }

  #--------------------------
  # 3.- Find the degrees of freedom k by minimizing the negative log-likelihood
  k.fit.shrunk <-  mle( nlogL.shrunk,
                        start = list( k = 100 ),
                        method = "L-BFGS-B", lower = c( 5 ),
                        upper = c( Inf ) )

  #--------------------------
  return( k.fit.shrunk@coef[1] )

}

NULL
