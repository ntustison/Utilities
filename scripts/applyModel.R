library( ANTsR )
library( randomForest )
library( snowfall )

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

if( length( args ) < 3 )
  {
  cat( "Usage: Rscript applyModel.R dimension inputModel inputCSVFile ",
       "outputProbabilityImagePrefix <numberOfThreads=4>", sep = "" )
  stopQuietly()
  }

dimension <- as.numeric( args[1] )
inputModelName <- args[2]
fileList <- read.csv( args[3] )
probImagePrefix <- args[4]

numberOfThreads <- 1
if( length( args ) >= 5 )
  {
  numberOfThreads <- as.numeric( args[5] )
  }

###############################################
#
# Load model:  contained in the variable "modelForest"
#
###############################################

load( inputModelName )

###############################################
#
# Put the image data into a data frame (modelData)
#
###############################################

maskName <- fileList[1,1]
featureImages <- fileList[1,2:ncol( fileList )]
featureNames <- colnames( featureImages )

## Create the model data frame

maskImage <- antsImageRead( as.character( maskName ), dimension = dimension, pixeltype = 'unsigned int' )
mask <- as.array( maskImage )

maskIndices <- which( mask != 0 )

subjectData <- matrix( NA, nrow = length( maskIndices ), ncol = length( featureNames ) )
for( j in 1:length( featureNames ) )
  {
  cat( "  Reading feature image ", featureNames[j], ".\n", sep = "" )
  featureImage <- as.array( antsImageRead( as.character( featureImages[1,j] ), dimension = dimension, pixeltype = 'float' ) )

  values <- featureImage[maskIndices]
  subjectData[, j] <- values
  }

colnames( subjectData ) <- c( featureNames )
subjectData <- as.data.frame( subjectData )

# If the subject data has NA's, we need to get rid of them
# since predict.randomForest will return NA's otherwise.
# Setting NA's to 0 is a complete hack.
subjectData[is.na( subjectData )] <- 0

###############################################
#
# Predict using the model (in parallel)
#
###############################################

# Start the clock!
ptm <- proc.time()

# the function each thread calls
parallelPredict <- function( i ) {
  numberOfSamplesPerThread <- as.integer( nrow( subjectData ) / numberOfThreads )
  threadIndexRange <- ( ( i - 1 ) * numberOfSamplesPerThread + 1 ):( i * numberOfSamplesPerThread )
  if( i == numberOfThreads )
    {
    threadIndexRange <- ( ( i - 1 ) * numberOfSamplesPerThread + 1 ):( nrow( subjectData ) )
    }
  return( predict( modelForest, subjectData[threadIndexRange,], type = "prob" ) )
}

if( numberOfThreads == 1 )
  {

  subjectProbabilities <- predict( modelForest, subjectData, type = "prob" )

  # Stop the clock
  elapsedTime <- proc.time() - ptm
  cat( "Prediction took ", as.numeric( elapsedTime[3] ), " seconds.\n", sep = "" )

  ###############################################
  #
  # Write the probability images to disk
  #
  ###############################################

  for( i in 1:ncol( subjectProbabilities ) )
    {
    probImage <- antsImageClone( maskImage, "float" )
    probImage[maskImage != 0] <- subjectProbabilities[,i];
    probFileName <- paste( probImagePrefix, i, ".nii.gz", sep = "" )
    cat( "Writing ", probFileName, ".\n" )
    antsImageWrite( probImage, probFileName )
    }
  } else {
  # Initialize the cluster
  sfInit( parallel = TRUE, cpus = numberOfThreads, type = 'SOCK' )

  # Make data available to each R instance / node
  sfExport( list = c( "modelForest", "subjectData", "numberOfThreads" ) )

  # Load library on each R instance / node
  sfClusterEval( library( randomForest ) )

  # Use a parallel random number generator to avoid correlated random numbers
  # this requires rlecuyer (which is default)
  sfClusterSetupRNG()

  # build the random forests
  parallelProbabilities <- sfClusterApply( 1:numberOfThreads, parallelPredict )

  sfStop()

  # everything finished so merge all forests into one
  subjectProbabilities <- parallelProbabilities[[1]]
  for( i in 2:numberOfThreads )
    {
    subjectProbabilities <- rbind( subjectProbabilities, parallelProbabilities[[i]] )
    }

  # Stop the clock
  elapsedTime <- proc.time() - ptm
  cat( "Prediction took ", as.numeric( elapsedTime[3] ), " seconds.\n", sep = "" )

  ###############################################
  #
  # Write the probability images to disk
  #
  ###############################################

  for( i in 1:ncol( subjectProbabilities ) )
    {
    probImage <- antsImageClone( maskImage, "float" )
    probImage[maskImage != 0] <- subjectProbabilities[,i];
    probFileName <- paste( probImagePrefix, i, ".nii.gz", sep = "" )
    cat( "Writing ", probFileName, ".\n" )
    antsImageWrite( probImage, probFileName )
    }
  }

