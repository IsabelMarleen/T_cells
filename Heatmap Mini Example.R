#Matrix CD3
m <- data$DLBCL2[ , c("smooth_CD3E", "smooth_CD4", "smooth_CD8B" ) ] ^ .1

geneplotter::multidensity( m, ylim=c(0,10) )

sigmoid <- function( x, m=0, h=0.025 )
  1 / ( 1 + exp( -(x-m)/h ) )

m[,1] <- sigmoid( m[,1], 0.4 )
m[,2] <- sigmoid( m[,2], 0.25 )
m[,3] <- sigmoid( m[,3], 0.2 )

pheatmap( m[ rowSums(m) > .1, ], cluster_cols=FALSE )
