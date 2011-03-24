regions = c(
  expression( paste( bold( "LL Lobe" ) ) ),
  expression( paste( bold( "LM Lobe" ) ) ),
  expression( paste( bold( "LU Lobe" ) ) ),
  expression( paste( bold( "RC Lobe" ) ) ),
  expression( paste( bold( "RL Lobe" ) ) ),
  expression( paste( bold( "RM Lobe" ) ) ),
  expression( paste( bold( "RU Lobe" ) ) ),
  expression( paste( bold( "LL Vessels" ) ) ),
  expression( paste( bold( "LM Vessels" ) ) ),
  expression( paste( bold( "LU Vessels" ) ) ),
  expression( paste( bold( "RC Vessels" ) ) ),
  expression( paste( bold( "RL Vessels" ) ) ),
  expression( paste( bold( "RM Vessels" ) ) ),
  expression( paste( bold( "RU Vessels" ) ) ),
  expression( paste( bold( "Airways" ) ) )
  );


pdf( "/Users/nick/Desktop/PSR4_pre.pdf" )
dice <- read.csv(
  file = "/Users/nick/Desktop/Data/ConnieDogData/Results/psrLabelOverlapMeasurespreOUT.txt",
  sep = ",", head = TRUE );

par( mar = c( 7, 4, 4, 2 ) + 0.1 );
boxplot( dice, varwidth = TRUE, notch = FALSE,
  col = rainbow(15), pch = 'o', lwd = 2,
  outlwd = 2, names = FALSE,
  ylim = c(0.0, 1.0), ylab = expression( paste( bold( "Dice Coefficient" ) ) ) );
#  names = regions,
#  xlab = expression( paste( bold( "Lung Regions" ) ) ) );
title( main = "Point Set Registration Results (Pre)", col.main = "black", font.main = 2 );
#title( main = "DMFFD Registration Results (Pre)", col.main = "black", font.main = 2 );
axis( 1, lwd = 2, font = 2, labels = FALSE, tick = FALSE );
axis( 2, tick = TRUE, lwd = 2, font = 2 );

text( 1:15, par( "usr" )[3] - 0.075, srt = 45, adj = 1,
  labels = regions, xpd = TRUE, col = rainbow( 15 ) );
mtext( 1, text = expression( paste( bold( "Lung Regions" ) ) ), line = 5.5 );
dev.off()
