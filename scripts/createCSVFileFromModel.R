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
  cat( "Usage: Rscript createCSVFileFromModel.R inputModel maskImage imagePrefix outputCSVFile", sep = "" )
  stopQuietly()
  }

inputModelName <- args[1]
maskImage <- args[2]
imagePrefix <- args[3]
outputFile <- args[4]

###############################################
#
# Load model:  contained in the variable "modelForest"
#
###############################################

load( inputModelName )

featureNames <- attr( modelForest$terms, "term.labels" )
fileNames <- c();
for( i in 1:length( featureNames ) )
  {
  fileNames[i] <- paste0( imagePrefix, featureNames[i], ".nii.gz" )
  }

featureNames <- append( featureNames, "MASK", after = 0 )
fileNames <- append( fileNames, maskImage, after = 0 )

write.table( rbind( featureNames, fileNames ), file = outputFile,
  append = FALSE, col.names = FALSE, row.names = FALSE, sep = ",")
