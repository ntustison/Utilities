#setwd( '/Users/nick/Desktop/Data/ConnieDogData' );
files <- list.files( "/Users/nick/Desktop/Data/ConnieDogData/",
  full.names = TRUE, pattern = "*.csv" );

LMvalues <- mat.or.vec( length( files ), 3 );
LUvalues <- mat.or.vec( length( files ), 3 );
RUvalues <- mat.or.vec( length( files ), 3 );

for( i in 1:length( files ) )
  {
  print( files[i] );

  d <- read.csv( files[i] );
  mydata <- data.matrix( d[,-1] );

  LMvalues[i,1] <- mydata[1,3];
  LMvalues[i,2] <- mydata[10,3];
  LMvalues[i,3] <- mydata[19,3];
  LUvalues[i,1] <- mydata[1,5];
  LUvalues[i,2] <- mydata[10,5];
  LUvalues[i,3] <- mydata[19,5];
  RUvalues[i,1] <- mydata[1,13];
  RUvalues[i,2] <- mydata[10,13];
  RUvalues[i,3] <- mydata[19,13];
  }


colors <- c( "red", "green", "blue", "brown", "black", "darkgreen" );
pchNumbers <- c( 19, 21, 22, 23, 24, 25 );


timePoints <- c( 0.0, 3.0, 15.0 );
plot( c( 0, 15 ), c( min( LMvalues ), max( LMvalues ) ), type = "n",
  xlab = expression( paste( bold( "Time (Months)" ) ) ),
  ylab = expression( paste( bold( "Strain" ) ) ),
  frame.plot = FALSE );
title( main = "Lung Principal Strain", col.main = "black", font.main = 2 );
for( i in 1:length( files ) )
  {
  lines( timePoints, LMvalues[i,], col = colors[i], lty = i, lwd = 3 );
  points( timePoints, LMvalues[i,], col = colors[i], pch = pchNumbers[i], lwd = 3 );
  }
legend( "topright", c( "Achilles", "Agamemnon", "Hector", "Paris" ),
  pch = pchNumbers[1:length( files )], lty = 1:4,
  col = colors[1:length(files)],
  title = expression( paste( bold( "Canine Subject" ) ) ),
  bty = "n",
  lwd = 1.5 );
axis( 1, tick = TRUE, lwd = 2, font = 2 );
axis( 2, tick = TRUE, lwd = 2, font = 2 );
