#' pval.montecarlo
#'
#' P-values for the partial correlations with shrinkage using MonteCarlo.
#' It test whether the partial correlation is different from zero.
#'
#' @docType package
#'
#' @author Victor Bernal \email{v.a.bernal.arzola@rug.nl; victor.arturo.bernal@gmail.com}
#'
#' @name pval.montecarlo
#'
#' @import GeneNet
#' @import corpcor
#' @import stats4
#'
#' @references \url{https://doi.org/10.1093/bioinformatics/btz357}
#'
#' @param r matrix of partial correlations
#' @param number number of simulations
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
#' pval <- pval.montecarlo(r = estimated.pcor)
#'
#' #---------------------------------
#'
#' @export

pval.montecarlo <- function( r = NULL, number = 10 ){

  #-----------------------
  # Meta data
  p <- attr(r, 'number.variables')

  n <- attr(r, 'sample.size')

  lambda <- attr(r, 'lambda')

  #-----------------------
  # Initialize variables
  r.monte <- matrix(NA, nrow(r), 1)

  p.values <- matrix(NA, nrow(r), ncol(r))

  cum.pv <- matrix(0, nrow(r) ,1)

  #-----------------------
  # Simulate null hypothetic GGM coefficients for "number" times
  for (i in 1 : number){

    #-----------------------
    sim.data.S <- GeneNet::ggm.simulate.data(sample.size = n, pcor = diag(p))

    #----------------------
    # Estimate the Correlation matrix
    R = cor(sim.data.S)

    R_lambda = ( 1 - lambda )*( R ) + ( lambda ) * diag( diag( R ) ) # shrinkage

    D = chol( x = R_lambda ) # Cholesky decomposition

    omega.chol = solve( D, diag( p) ) %*% solve( t(D), diag(p) )# precision matrix

    chol.PCOR = -cov2cor( omega.chol )

    #----------------------
    GGM.S <- chol.PCOR #GGM.shrunk( x= sim.data.S, lambda = lambda, verbose = F ) #corpcor::pcor.shrink(x = sim.data.S, lambda = lambda)

    r.monte <- sm2vec(GGM.S) # vectorize

    #-----------------------
    # compare the real coefficients against r.monte
    pv <- apply(X = r, MARGIN = 1,  FUN = function(x) sum( abs( r.monte ) >= abs(x) )/length( r.monte ) )

    cum.pv <- cum.pv+pv

  }

  return('pval.MC'= cum.pv/number)

}

NULL
