setwd( '/Users/ntustison/Desktop/' );

volumes <- read.csv( "vdv.csv", header = FALSE );
volumes <- volumes / 1e6;
rowmeans <- rowMeans( volumes );

readerNames <- c( "Atropos", "Carlos", "JAM", "JG" );
colors <- c( "red", "green", "blue", "orange" );

for( i in 1:length( readerNames ) )
  {
  difference <- ( volumes[,i] - rowmeans );

  mu <- mean( difference );
  std <- sd( difference );

  insubjects <- c();
  outsubjects <- c();

  for( k in 1:length( difference ) )
    {
    if( difference[[k]] <= mu + 2*std && difference[[k]] >= mu - 2*std )
      {
      insubjects <- c( insubjects, k );
      }
    else
      {
      outsubjects <- c( outsubjects, k );
      }
    }


  pdf( paste( "~/Desktop/BA_", readerNames[i], ".pdf", sep = "" ) );
  average <- 0.5 * ( volumes[insubjects,i] + rowmeans[insubjects] );
  plot( average, difference[insubjects],
        xlim = c(0,3.5), ylim = c(-1., 1.),
        xlab = expression( paste( "Mean" ) ),
        ylab = expression( paste( "Difference" ) ),
        lwd = 2.0, col = "green" );

  average <- 0.5 * ( volumes[outsubjects,i] + rowmeans[outsubjects] );
  points( average, difference[outsubjects], lwd = 2.0, col = "red" );
  abline( a = 0, b = 0, col = "black", lwd = 0.5 );
  abline( a = mu, b = 0, col = "blue", lty = "dotted", lwd = 1.5 );
  text( 3, mu, expression(mu), col = "blue", adj = c(0, -.1) )
  abline( a = mu + 2.0*std, b = 0, col = "red", lty = "dashed", lwd = 1.5 );
  text( 3, mu + 2*std, expression(mu + 2 * sigma), col = "red", adj = c(0, -.1) );
  abline( a = mu - 2.0*std, b = 0, col = "red", lty = "dashed", lwd = 1.5 );
  text( 3, mu-2*std, expression(mu + 2 * sigma), col = "red", adj = c(0, -.1) );
  title( main = paste( "Bland-Altman (", readerNames[i], ")", sep = "" ), col.main = "black", font.main = 2 );
  axis( 1, lwd = 2, font = 2, labels = FALSE, tick = TRUE );
  axis( 2, tick = TRUE, lwd = 2, font = 2 );
  dev.off();
  }
