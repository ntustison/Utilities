regions = c(
  expression( paste( bold( "LM Lobe" ) ) ),
  expression( paste( bold( "LU Lobe" ) ) ),
  expression( paste( bold( "RU Lobe" ) ) ),
  expression( paste( bold( "LM Vessels" ) ) ),
  expression( paste( bold( "LU Vessels" ) ) ),
  expression( paste( bold( "RU Vessels" ) ) ),
  expression( paste( bold( "Airways" ) ) )
  );


pdf( "/Users/nick/Desktop/PSR4_15mo.pdf" )
dice <- read.csv(
  file = "/Users/nick/Desktop/Data/ConnieDogData/Results/psrLabelOverlapMeasures15moOut.txt",
  sep = ",", head = TRUE );

par( mar = c( 7, 4, 4, 2 ) + 0.1 );
boxplot( dice, varwidth = TRUE, notch = FALSE,
  col = rainbow(7), pch = 'o', lwd = 2,
  outlwd = 2, names = FALSE,
  ylim = c(0.0, 1.0), ylab = expression( paste( bold( "Dice Coefficient" ) ) ) );
#  names = regions,
#  xlab = expression( paste( bold( "Lung Regions" ) ) ) );
title( main = "Point Set Registration Results (3 mo)", col.main = "black", font.main = 2 );
axis( 1, lwd = 2, font = 2, labels = FALSE, tick = FALSE );
axis( 2, tick = TRUE, lwd = 2, font = 2 );

text( 1:7, par( "usr" )[3] - 0.075, srt = 45, adj = 1,
  labels = regions, xpd = TRUE, col = rainbow( 7 ) );
mtext( 1, text = expression( paste( bold( "Lung Regions" ) ) ), line = 5.5 );
dev.off()
