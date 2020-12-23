################################################################################################
####                                                                                        ####
####             Selection of Estimators for the LKC of Gaussian fields                     ####
####                                                                                        ####
################################################################################################
#
#  Implemented Estimators:
#     - direct LKC estimator by estimating the Riemannian metrix and computing the LKCs directly (1D and 2D)
#       as in Telschow Schwartzmann (2019)
#     - hermite projection estimator as in Schwartzman et al ()
#     - warping estimator as in Taylor (2007)
#
# Currently the code only works for rectangular domains. Extensions are possible and might follow
#
#################################### direct LKC estimator ####################################
#' Estimates the Lipschitz killing curvatures of a Gaussian field over a rectangular domain using a finite difference estimator for the differential process. It is numerical (tested up to d=2) identical to the warping estimator proposed in  TaylorWorsley (2007).
#'
#' IMPORTANT NOTES: the process must have mean 0 and variance 1 and current implementation only works for 1d and 2d domains of the random field.
#'
#' @param R Array (dimension K_1 x ... x K_d x N containing N-realizations of a standardaized Gaussian random field over a d-dimensional domain.)
#' @param subdiv Vector (containing the breakpoints to split the array in independent fields, e.g. for 2 sample tests)
#' @return LKC Vector containing the Lipschitz killing curvatures of dimension greater than 1. Note that the 0-th LKC must be calculated seperately as the Euler characteristic of the domain.
#' @export
LKC_estim_direct = function( R, subdiv=dim(R)[length(dim(R))], heavy = FALSE ){
  #### Get dimension of the domain of the data
  dimR = dim(R); N = diff(c(0,subdiv)); D = length(dimR)-1;

  #### Check user input
  if( !all(diff(c(1,subdiv))>1) || subdiv[length(subdiv)]!=dimR[length(dimR)] ){
    stop("The subdiff vector must be strictly increasing and the last entry must match the last dimension of the array") ;
  }
  #### Make the subdiff vector usable for the for-loop
  Nsubdiv = length(subdiv)
  subdiv  = c( 1, subdiv );

  #### Initiate the LKC vector
  LKC = rep( NA, D ) ;

  #### Switch cases depending on dimension
  if( D == 1 ){
    if( heavy ){
      dVar = apply( diff( R ), 1,
                    function(x) fit_mvt( x )$cov )
      LKC = sum( sqrt( dVar ) )
    }else{
      if( !all( subdiv == dim( R )[ length( dim( R ) ) ] ) ){
      VardR = rep(0, dimR[1]-1)
      # Estimate dSigma^2 assuming independence across the subdivisions of the fields
      for( u in 1:Nsubdiv ){
            dtmp    = ( R[2:dimR[1], subdiv[u]:subdiv[u+1]] - R[1:(dimR[1]-1), subdiv[u]:subdiv[u+1]] )^2 ;
            VardR   = VardR + apply( dtmp, 1, function(row) sum(row) ) / (N[u]-1) ;
      }
      # Integrate dSigma over the domain. Note that dt cancels, since we did not divide the finite differences.
      LKC    = sum( sqrt(VardR) ) ;
    }else{
      x_scl    <- t( apply( R, 1, scale ) ) # pointwise standardization (mean=0, sd=1)
      x_scl    <- x_scl
      x_p      <- apply( x_scl, 2,
                         FUN = function( yy ){     # differentiation
                           xx <- seq( 0, 1, len = length( yy ) )
                           fn <- stats::splinefun( x = xx,
                                                   y = yy,
                                                   method = "natural")
                           pracma::fderiv( f = fn,
                                           x = xx,
                                           n = 1,
                                           h = diff(xx)[1],
                                           method = "central")
                           } )
      tau_t    <- apply( x_p, 1, stats::sd )           # pointwise sd
    }
    }
  }else if( D==2 ){
    # Compute length of horizontal/vertical/diagonal edges, the
    # triangulation of the rectangle is given by the following pattern for
    # 3 x 3 square
    #
    # o---o---o---o---o
    # | \ | \ | \ | \ |
    # o---o---o---o---o
    # | \ | \ | \ | \ |
    # o---o---o---o---o

    ####### Compute estimator L1 = integral over the boundary of the standard deviation of the differentiated process
    edges_vert = matrix(0, dimR[1]-1, dimR[2]) ;
    edges_horz = matrix(0, dimR[1], dimR[2]-1) ;
    # Estimate dSigma^2 assuming independence across the subdivisions of the fields
    for( u in 1:Nsubdiv ){
      tmp_edges_vert = ( R[ 2:dimR[1], , subdiv[u]:subdiv[u+1]] - R[1:(dimR[1]-1), , subdiv[u]:subdiv[u+1]] )^2 ;
      tmp_edges_horz = ( R[ , 2:dimR[2], subdiv[u]:subdiv[u+1]] - R[, 1:(dimR[2]-1), subdiv[u]:subdiv[u+1]] )^2 ;

      edges_vert = edges_vert + rowSums( tmp_edges_vert, dims=2 ) / (N[u]-1) ;
      edges_horz = edges_horz + rowSums( tmp_edges_horz, dims=2 ) / (N[u]-1) ;
    }

    edges_vert_sd = sqrt( edges_vert ) ;
    edges_horz_sd = sqrt( edges_horz ) ;

    # integrate the standard deviation. Note that dt cancels, since we did not divide the finite differences
    LKC[1] = 0.5*(    sum( edges_horz_sd[1,] ) + sum( edges_horz_sd[dimR[1],] ) + sum( edges_vert_sd[,1] ) + sum( edges_vert_sd[,dimR[2]] )     ) ;

    ####### Compute estimator L2 = integral over the domain of the determinant of the standard deviation of the differentiated process
    edges_vert = R[ 2:dimR[1],,] - R[1:(dimR[1]-1),, ] ;
    edges_horz = R[,2:dimR[2], ] - R[ ,1:(dimR[2]-1),] ;
    ## Compute sample variance matrix for forward derivatives
    sampCov11 <- sampCov22 <- sampCov12 <- 0;
    for( u in 1:Nsubdiv ){
        gradRx1      <- edges_vert[,1:(dimR[2]-1), subdiv[u]:subdiv[u+1]]
        gradRy1      <- edges_horz[1:(dimR[1]-1),, subdiv[u]:subdiv[u+1]]
        sampCov11    <- sampCov11 + rowSums( (gradRx1)^2, dims=2 ) / (N[u]-1)
        sampCov22    <- sampCov22 + rowSums( (gradRy1)^2, dims=2 ) / (N[u]-1)
        sampCov12    <- sampCov12 + rowSums( gradRx1*gradRy1, dims=2 ) / (N[u]-1)
    }
    sqrtdetCov1 = sqrt( sampCov11*sampCov22 - sampCov12^2 ) ;

    ## Compute sample variance matrix for backward derivatives
    sampCov11 <- sampCov22 <- sampCov12 <- 0;
    for( u in 1:Nsubdiv ){
      gradRx1      <- edges_vert[,2:(dimR[2]), subdiv[u]:subdiv[u+1]]
      gradRy1      <- edges_horz[2:(dimR[1]),, subdiv[u]:subdiv[u+1]]
      sampCov11    <- sampCov11 + rowSums( (gradRx1)^2, dims=2 ) / (N[u]-1)
      sampCov22    <- sampCov22 + rowSums( (gradRy1)^2, dims=2 ) / (N[u]-1)
      sampCov12    <- sampCov12 + rowSums( gradRx1*gradRy1, dims=2 ) / (N[u]-1)
    }
    sqrtdetCov2 = sqrt( sampCov11*sampCov22 - sampCov12^2 ) ;

    # integrate the square root of the determinant of the covariance estimate. Note that dt^2 cancels, since we did not divide the finite differences!
    LKC[2]      = ( sum(sqrtdetCov1) + sum(sqrtdetCov2) ) / 2

  }else{
    stop("Estimation of LKC for Domains of Dimension greater than 2 is currently not supported.") ;
  }

  return( LKC ) ;
}

#################################### Hermite polynomial LKC estimator ####################################
#################### Estimation Euler characterisitic curves from the data
#' Estimates the Euler characteristic of the exceedance set A(u)={ s in S | f(s) >u } of a function f:S->R.
#' 4 or 8 connectivity  is currently used and it is implemented up to dimension 2.
#'
#' @param f Vector/Matrix observation of a random function on a grid
#' @param u Vector exceedance values to be computed
#' @param connectivity Integer amount of points considered as neighbors
#' @return Vector of Euler characteristics of the exceedance sets for the values u
#' @export
EulerChar <- function( f, u, connectivity=8 ){
  # Make vector input a matrix
  if( is.vector(f) ){
    f <- matrix( f, 1, length(f) )
  }

  sz = dim(f)
  N  = ifelse( sz[1]==1, 1, length(sz) )

  EC <- rep( 0, length(u) )

  for( j in length(u):1 ){
    A = ( f >= u[j] )

    if( N == 1 ){
      vertices <- sum( A )
      edges    <- sum( A[1:sz[2]-1] & A[2:sz[2]] )
      EC[j]    <- vertices - edges

    }else if( N == 2 ){
      if( connectivity == 4 ){
        vertices <- sum( A )
        edges    <- sum( A[1:(sz[1]-1),] & A[2:sz[1],] ) + sum( A[,1:(sz[2]-1)] & A[,2:sz[2]] )
        faces    <- sum( A[1:(sz[1]-1),1:(sz[2]-1)] & A[2:(sz[1]),1:(sz[2]-1)] & A[1:(sz[1]-1),2:sz[2]] & A[2:sz[1],2:sz[2]] )
        EC[j]    <- vertices - edges + faces ;

      }else if( connectivity == 8 ){
        vertices <- sum( A )
        edges    <- sum( A[1:(sz[1]-1), ] & A[2:sz[1], ] ) + sum( A[ ,2:sz[2]] & A[ ,1:(sz[2]-1)] ) +
          sum( A[1:(sz[1]-1),1:(sz[2]-1)] & A[2:sz[1],2:sz[2]] ) + sum( A[2:sz[1],1:(sz[2]-1)] & A[1:(sz[1]-1),2:sz[2]] )
        # Summands left to right correspond to the following faces
        #   o--o    o  x   o--o    x  o
        #   \ /     \ \      \\     / \
        #   o  x    o--o   x  o    o--o
        p.faces  <- sum( A[1:(sz[1]-1),1:(sz[2]-1)] & A[2:sz[1],1:(sz[2]-1)] &  A[1:(sz[1]-1),2:sz[2]] ) +
          sum( A[1:(sz[1]-1),1:(sz[2]-1)] & A[2:sz[1],1:(sz[2]-1)] &  A[2:sz[1],2:sz[2]] ) +
          sum( A[1:(sz[1]-1),2:sz[2]]     & A[2:sz[1],2:sz[2]]     &  A[1:(sz[1]-1),1:(sz[2]-1)] ) +
          sum( A[1:(sz[1]-1),2:sz[2]]     & A[2:sz[1],2:sz[2]]     &  A[2:sz[1],1:(sz[2]-1)] )
        # Note that we only count 4 instead of 5 verices and 6 instead of 8 edges. Thus, we need to subtract a face if the following scheme appears
        # o---o
        # \ x \
        # o---o
        n.faces  <- sum( A[1:(sz[1]-1),1:(sz[2]-1)] & A[2:(sz[1]),1:(sz[2]-1)] & A[1:(sz[1]-1),2:sz[2]] & A[2:sz[1],2:sz[2]] )
        EC[j]    <- vertices - edges + p.faces - n.faces ;

      }else{
        stop("For 2D fields 'connectivity' has to have the value 4 or 8.")
      }

    }else{
      stop("N >= 2 currently not implemented")
    }
  }
  return( EC )
}

#' Computes the LKCs of the boundary of a domain.
#'
#' @param x x-Coordinates of the grid on which the data is observed.
#' @param y y-Coordinates of the grid on which the data is observed.
#' @param cont The contour of f at value level
#' @param R An array of dimension c(length(x),length(y),n) containing the
#'          realizations of the field.
#' @importFrom stats na.omit pnorm
#' @return A function g that computes for u>0 the probility that the supremum of
#'         the field exceeds u.
LKCcontDomain = function( x, y, bdry, R, h = 0.01 ){

  # Number of contour elements
  nCont <- length(bdry$bdry)
  cont  <- bdry$bdry
  if( nCont == 0) return( function(u) return(0) )

  # number of samples
  nn = dim(R)[3]

  #---- Compute the contour line integrals
  SC = vector( "list", nCont )

  # interpolate the residuals to the contour line
  for( i in 1:nCont ){
    for( j in 1:nn ){
      if( j == 1 ){
        SC[[i]] = fields::interp.surface(
                      list( x = x, y = y, z = R[,,j] ),
                      cbind( cont[[i]]$x, cont[[i]]$y )
                      )
      }else{
        SC[[i]] = cbind( SC[[i]], fields::interp.surface(
                          list( x = x, y = y, z = R[,,j] ),
                                cbind( cont[[i]]$x, cont[[i]]$y ) ) )
      }
    }
  }
  # clean up NAs due to extrapolation
#  SC = sapply( SC, na.omit )

  #Gives the length of one component.
  L1 = function( X ){
    if( !is.matrix( X ) ) return( 0 ) # Pathological cases.
    if( nrow(X) == 1 ) return( 0 )
    X = X / sqrt( rowSums( X^2 ) )
    sum( sqrt( rowSums( diff(X)^2 ) ) )
  }

  # sums the length of all components wrt the Riemannian metric
  L1 = sum( sapply( SC, L1 ) )


  ##---- Compute L2
  volform <- function( ss ){
    N = dim( R )[3]
    normR = R / sqrt( rowSums( R^2 ) / ( N - 1 )  )
    dRx <- dRy <- matrix( NaN, dim(ss)[1], N )
    for( j in 1:N ){
      dRx[,j] <- ( fields::interp.surface(
                        list( x = x, y = y, z = normR[,,j] ),
                        cbind( ss[,1] + h, ss[,2] )
                      ) -
                   fields::interp.surface(
                      list( x = x, y = y, z = normR[,,j] ),
                      cbind( ss[,1] - h, ss[,2] )
                   ) ) / 2 / h

      dRy[,j] <- ( fields::interp.surface(
                      list( x = x, y = y, z = R[,,j] ),
                      cbind( ss[,1], ss[,2] + h )
                    ) -
                    fields::interp.surface(
                      list( x = x, y = y, z = R[,,j] ),
                      cbind( ss[,1], ss[,2] - h )
                    ) ) / 2 / h
    }

    sqrt( rowSums( dRx^2 ) * rowSums( dRy^2 ) - rowSums( dRx*dRy )^2 ) / ( N-1 )

  }

  volumes <- NaN * ( 1:nCont )
  for ( i in 1:nCont ) {

    volumes[i] <- polyCub.midpoint( bdry, f = volform, plot = TRUE )
  }

 L2 = sum( volumes )

  c( L1 = L1, L2 = L2 )
}
