setwd( '/Users/nick/Desktop/' );
r1 <- read.csv( "reader1.csv" );
r2 <- read.csv( "reader2.csv" );
at <- read.csv( "defectsMNDS.csv" );
r1_mnds <- r1[,3]/r1[,2];
r2_mnds <- r2[,3]/r2[,2];
mnds <- cbind( r1_mnds, r2_mnds, at[,2] );

colnames( mnds ) <- c( "reader1", "reader2", "atropos" );

ICC( mnds[,1:2] );
ICC( mnds[,2:3] );

assay1 = 1;
assay2 = 3;

difference <- mnds[,assay1] - mnds[,assay2];

mu <- mean( difference );
std <- sd( difference );

insubjects = c();
outsubjects = c();

for( i in 1:length( difference ) )
  {
  if( difference[[i]] <= mu + 2*std && difference[[i]] >= mu - 2*std )
    {
    insubjects <- c( insubjects, i );
    }
  else
    {
    outsubjects <- c( outsubjects, i );
    }
  }

#pdf( "Reader1_Atropos.pdf" );
average <- 0.5 * ( mnds[insubjects,assay1] + mnds[insubjects,assay2] );
plot( average, difference[insubjects],
      xlim = c(0,4), ylim = c(-4, 4),
      xlab = expression( paste( "Mean" ) ),
      ylab = expression( paste( "Difference" ) ),
      lwd = 2.0, col = "green" );
average <- 0.5 * ( mnds[outsubjects,assay1] + mnds[outsubjects,assay2] );
points( average, difference[outsubjects],
      lwd = 2.0, col = "red" );
abline( a = 0, b = 0, col = "black", lwd = 0.5 );
abline( a = mu, b = 0, col = "blue", lty = "dotted", lwd = 1.5 );
text( 3.5, mu, expression(mu), col = "blue", adj = c(0, -.1) )
abline( a = mu + 2.0*std, b = 0, col = "red", lty = "dashed", lwd = 1.5 );
text( 3.5, mu + 2*std, expression(mu + 2 * sigma), col = "red", adj = c(0, -.1) )
abline( a = mu - 2.0*std, b = 0, col = "red", lty = "dashed", lwd = 1.5 );
text( 3.5, mu-2*std, expression(mu + 2 * sigma), col = "red", adj = c(0, -.1) )
#title( main = "Bland-Altman Plot (ICC = 0.85)", col.main = "black", font.main = 2 );
axis( 1, lwd = 2, font = 2, labels = FALSE, tick = TRUE );
axis( 2, tick = TRUE, lwd = 2, font = 2 );
#dev.off();

#BlandAltman( mnds[,1], mnds[,2], gui = FALSE, bandsOn = TRUE, biasOn = FALSE, regionOn = FALSE );
#
#BlandAltman( mnds[,2], mnds[,3], gui = FALSE, bandsOn = TRUE, biasOn = FALSE, regionOn = FALSE );
