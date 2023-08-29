mvnorm <- function(n,vs,rhos){
  Sigma <- ( sqrt(vs) %o% sqrt(vs) ) * ( rhos )
  matrix( rnorm(n*length(vs)), n, length(vs) ) %*% phenix:::mat.sqrt(Sigma)
}

simfxn <- function(N, S, h2s, rhog, rhoe, r2fgrs, prevs ){
  G     <- scale(matrix( rnorm(N*S), N, S ))
  betas <- mvnorm(n=S,vs=h2s,rhos=rhog) / sqrt(S) 
  liabs <- G %*% betas

  rhofgrs   <- ( sqrt(h2s) %o% sqrt(h2s) ) * ( rhog ) + ( sqrt(1-h2s) %o% sqrt(1-h2s) ) * ( rhoe )
  sig2fgrs  <- h2s^2 / r2fgrs - h2s
  fgrs      <- liabs + mvnorm(n=N,vs=sig2fgrs,rho=rhofgrs) 
  fgrs      <- scale( fgrs )

  y     <- liabs + mvnorm(n=N,vs=1-h2s,rhos=rhoe)
  y     <- sapply( 1:ncol(y), function(j){ y[,j]>quantile(y[,j],1-prevs[j]) })

  list( fgrs=fgrs, liabs=liabs, y=y )
}

run_one <- function( N = 3e3, S = 30, r2fgrs = c(.03,.03), prevs=c(.1,.1), h2s=c(.5,.5), niter=1e2, rhog=NA, rhoe=NA, savefile ){
  deltas  <- array( NA, dim=c( length(phens), length(phens), niter ), dimnames=list( phens, phens, 1:niter ) )
  for( it in 1:niter ) try({
    print(it)
    simdat <- simfxn(N, S, h2s, rhog=rhog, rhoe=rhoe, r2fgrs, prevs ) 
    for( phen1 in phens )
      for( phen2 in phens )
    {
      score <- simdat$fgrs [,which(phens==phen2)]
      case  <- simdat$y    [,which(phens==phen2)]==1
      if( phen1 == phen2 ){
        #base <- rep( TRUE, nrow(simdat$y) )
        base <- rowSums( simdat$y[,-which(phens==phen1)] ) == 0
      } else {
        #base <- simdat$y[,which(phens==phen1)]==1
        base <- ( simdat$y[,which(phens==phen1)]==1 ) & ( rowSums( simdat$y[,-which(phens%in%c(phen1,phen2))] )==0 )
      }
      deltas[phen1,phen2,it]<- mean( score[  case & base] ) -
                               mean( score[(!case)& base] )
    }
    save( deltas, file=savefile )
  })
}
