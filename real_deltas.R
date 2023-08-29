rm( list=ls() ) 
phens   <- c( 'MD','AD','AUD','DUD' )

prevs       <- as.numeric(as.character(read.table( 'prev_table.txt', fill=T )[9:12,3]))/100

r2fgrs<- c( 0.0393, 0.0422, 0.0464, 0.0505 )
#MD: 0.0393 (0.0381-0.0405)
#AD: 0.0422 (0.0410-0.0435)
#AUD: 0.0464 (0.0446-0.0481)
#DUD: 0.0505 (0.0490-0.0521)

deltas    <- read.csv( 'deltas.csv', head=F )
deltahats <- array( NA, dim=c( 4, 4, 2 ), dimnames=list( phens, paste0( phens, '-FGRS' ), c( 'est', 'sd' ) ) )
for( i in 1:4 ){
  fgrs1  <- as.numeric(as.character(deltas[(i-1)*15 + c(2,4,6,8)  ,3]))
  fgrs2  <- as.numeric(as.character(deltas[(i-1)*15 + c(2,4,6,8)  ,4]))
  fgrs1se<- as.numeric(as.character(deltas[(i-1)*15 + c(2,4,6,8)+1,3]))
  fgrs2se<- as.numeric(as.character(deltas[(i-1)*15 + c(2,4,6,8)+1,4]))

  fgrs1  [1]  <- 0
  fgrs1se[1]  <- 0
  deltahats[c(i,setdiff(1:4,i)),i,1]  <- fgrs2-fgrs1
  deltahats[c(i,setdiff(1:4,i)),i,2]  <- sqrt( fgrs2se^2 + fgrs1se^2 )
}
deltahats[,,1]
deltahats[,,2]

rhogs <- as.matrix( read.table( 'rhogs.txt', head=T ) )
h2s   <- diag(rhogs)
diag(rhogs) <- 1

names(h2s)      <- phens
rownames(rhogs) <- phens
colnames(rhogs) <- phens
names(prevs)    <- phens
names(r2fgrs)   <- phens

save.image( file='deltas.Rdata' )
