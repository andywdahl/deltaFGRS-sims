rm( list=ls() ) 
load( 'data/deltas.Rdata' )
source( 'sim_fxn.R')  

N     <- 1e5
S     <- 100
niter <- 1e2

rhoetype  <- 'rhoe=rhog'
rhoe  <- rhogs

rhoetype  <- 'rhoe=0'
rhoe  <- rhogs * diag(nrow(rhogs))

savefile <- paste0( 'data/', rhoetype, '_N=', N, '_S=', S, '.Rdata' )
if( !file.exists(savefile) )
run_one( N=N, S=S, niter=niter, savefile=savefile,
  r2fgrs = r2fgrs, prevs = prevs, h2s = h2s, rhog = rhogs, rhoe = rhoe )
load( savefile ) 

deltameans  <- apply( deltas, 1:2, mean, na.rm=T )
deltases    <- apply( deltas, 1:2, function(x) sd(x,na.rm=T)/sqrt(sum(!is.na(x))) )

pdf( paste0( 'deltaplot_', rhoetype, '_N=', N, '_S=', S, '.pdf' ), w=5, h=5 )
par( mar=c(4.5,4.5,.5,.5) )

x   <- c(deltameans)
y   <- c(deltahats[,,1])
xse <- c(deltases)
yse <- c(deltahats[,,2])

phencols  <- cbind(phens,phens,phens,phens)
phenlabs  <- rbind(phens,phens,phens,phens)
diag( phencols ) <- 'None'
phencols  <- c( phencols )
phenlabs  <- c( phenlabs )

cols  <- 1:4
names(cols)  <- phens

lims   <- c(0.12,0.92)
plot(   x, y , xlim=lims, ylim=lims, pch=16, ylab='Observed Delta-FGRS', xlab='Simulated Delta-FGRS', type='n' )
abline(a=0,b=1,col='grey',lty=2, lwd=2)

arrows( x, y - 2*yse, x, y + 2*yse, cex=.5, code=3, angle=90, length=.04, lwd=0.3, col=cols[phenlabs] )
text(   x, y, cex=1.0, lab=phencols, col=cols[phenlabs] )

legend( 'bottomright' , cex=1.0, bty='n', title='FGRS', fill=c(cols), leg=c(names(cols) ) )
legend( 'right'       , cex=1.0, bty='n', leg=paste0( 'r=', round( cor( x, y  ), 2 ) ) )

dev.off()
