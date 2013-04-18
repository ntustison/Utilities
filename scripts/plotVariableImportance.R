library( ggplot2 )
library( randomForest )

stopQuietly <- function(...)
  {
  blankMsg <- sprintf( "\r%s\r", paste( rep(" ", getOption( "width" ) - 1L ), collapse = " ") );
  stop( simpleError( blankMsg ) );
  } # stopQuietly()

args <- commandArgs( trailingOnly = TRUE )

###############################################
#
# Selected parameters
#
###############################################

if( length( args ) < 2 )
  {
  cat( "Usage: Rscript plotVariableImportance.R inputModel outputFile <title>", sep = "" )
  stopQuietly()
  }

inputModel <- args[1]
outputFile <- args[2]
plotTitle <- ''
if( length( args ) > 2 )
  {
  plotTitle <- args[3]
  }

load( inputModel )

forestImp <- importance( modelForest )
forestImp.df <- data.frame( Statistic = names( forestImp[,1] ), Importance = as.numeric( forestImp[,1] )  )
forestImp.df <- forestImp.df[order( forestImp.df$Importance ),]

forestImp.df$Statistic <- factor( x = forestImp.df$Statistic, levels = forestImp.df$Statistic )

vPlot <- ggplot( data = forestImp.df, aes( x = Importance, y = Statistic ) ) +
         geom_point( aes( color = Importance ) ) +
         labs( title = plotTitle ) +
         ylab( "Variable" ) +
         xlab( "Importance" ) +
         scale_color_continuous( low = "navyblue", high = "darkred" ) +
         theme( axis.text.y = element_text( size = 5 ) ) +
         theme( legend.position = "none" )

ggsave( file = outputFile, plot = vPlot, width = 8, height = 8 )



