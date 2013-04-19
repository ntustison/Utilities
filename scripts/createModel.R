library( ANTsR )
library( randomForest )
library( snowfall )
library( rlecuyer )

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
  cat( "Usage: Rscript createModel.R inputFileList outputModelPrefix ",
       "<numberOfThreads=4> <trainingPortion=1.0> <numberOfTreesPerThread=1000> ",
       "<numberOfSamplesPerLabel=1000>", sep = "" )
  stopQuietly()
  }

fileList <- read.csv( args[1] )
outputModelName <- paste0( args[2], ".RData" )

numberOfThreads <- 4
if( length( args ) >= 3 )
  {
  numberOfThreads <- as.numeric( args[3] )
  }
trainingPortion <- 1.0
if( length( args ) >= 4 )
  {
  trainingPortion <- as.numeric( args[4] )
  }
numberOfSamplesPerLabel <- 1000
if( length( args ) >= 5 )
  {
  numberOfSamplesPerLabel <- as.numeric( args[5] )
  }
numberOfTreesPerThread <- 1000
if( length( args ) >= 6 )
  {
  numberOfTreesPerThread <- as.numeric( args[6] )
  }

###############################################
#
# Put the image data into a data frame (modelData)
#
###############################################

truthLabels <- fileList[,1]
masks <- fileList[,2]
featureImages <- fileList[,3:ncol( fileList )]
featureNames <- colnames( featureImages )

totalNumberOfSubjects <- length( masks )
modelNumberOfSubjects <- floor( trainingPortion * totalNumberOfSubjects )

## Create the model data frame

modelData <- matrix()
indices <- sort( sample.int( totalNumberOfSubjects, modelNumberOfSubjects, replace = FALSE ) )
for( i in indices )
  {
  cat( as.character( truthLabels[i] ), "\n" )

  mask <- as.array( antsImageRead( as.character( masks[i] ), dimension = 3, pixeltype = 'unsigned int' ) )
  truth <- as.array( antsImageRead( as.character( truthLabels[i] ), dimension = 3, pixeltype = 'unsigned int' ) )

  uniqueTruthLabels <- sort( unique( truth[which( mask == 1 )] ) )
  uniqueTruthLabels <- uniqueTruthLabels[which( uniqueTruthLabels != 0 )]
  cat( "Unique truth labels: ", uniqueTruthLabels, "\n", sep = " " )

  truthLabelIndices <- list()
  numberOfSamplesPerLabelInSubjectData <- rep( 0, length( uniqueTruthLabels ) )
  for( n in 1:length( uniqueTruthLabels ) )
    {
    labelIndices <- which( truth == uniqueTruthLabels[n] )
    numberOfSamplesPerLabelInSubjectData[n] <- min( numberOfSamplesPerLabel, length( labelIndices ) )
    truthLabelIndices[[n]] <- labelIndices[sample.int( length( labelIndices ), numberOfSamplesPerLabelInSubjectData[n], replace = FALSE )]
    }

  subjectData <- matrix( NA, nrow = sum( numberOfSamplesPerLabelInSubjectData ), ncol = length( featureNames ) + 1 )
  for( j in 1:length( featureNames ) )
    {
    cat( "  Reading feature image ", featureNames[j], ".\n", sep = "" )
    featureImage <- as.array( antsImageRead( as.character( featureImages[i,j] ), dimension = 3, pixeltype = 'float' ) )
    for( n in 1:length( uniqueTruthLabels ) )
      {
      values <- featureImage[truthLabelIndices[[n]]]

      startIndex <- 1
      if( n > 1 )
        {
        startIndex <- sum( numberOfSamplesPerLabelInSubjectData[1:(n-1)] ) + 1
        }
      endIndex <- startIndex + length( values ) - 1

      subjectData[startIndex:endIndex, j] <- values
      if( j == 1 )
        {
        subjectData[startIndex:endIndex, length( featureNames ) + 1] <- rep.int( uniqueTruthLabels[n], length( truthLabelIndices[[n]] ) )
        }
      }
    }
  if( i == indices[1] )
    {
    modelData <- subjectData;
    }
  else
    {
    modelData <- rbind( modelData, subjectData )
    }
  }
colnames( modelData ) <- c( featureNames, "Labels" )
modelData <- as.data.frame( modelData )
modelData$Labels <- as.factor( modelData$Labels )

###############################################
#
# Create the random forest model in parallel
#
###############################################

# Start the clock!
ptm <- proc.time()

modelFormula <- as.formula( "Labels ~ . " )

#the function each thread calls
parallelRF <- function( i ) {
  return( randomForest( modelFormula, modelData, ntree = numberOfTreesPerThread, type = classification ) )
}

if( numberOfThreads == 1 )
  {
  modelForest <- randomForest( modelFormula, modelData, ntree = numberOfTreesPerThread, type = classification )

  # Stop the clock
  elapsedTime <- proc.time() - ptm
  cat( "Model creation took ", as.numeric( elapsedTime[3] ), " seconds.\n", sep = "" )

  ###############################################
  #
  # Save the model
  #
  ###############################################

  save( modelForest, file = outputModelName )
  }
else
  {
  # Initialize the cluster
  sfInit( parallel = TRUE, cpus = numberOfThreads, type = 'SOCK' )

  # Make data available to each R instance / node
  sfExport( list = c( "modelData", "modelFormula", "numberOfTreesPerThread" ) )

  # Load library on each R instance / node
  sfClusterEval( library( randomForest ) )

  # Use a parallel random number generator to avoid correlated random numbers
  # this requires rlecuyer (which is default)
  sfClusterSetupRNG()

  # build the random forests
  parallelForests <- sfClusterApply( 1:numberOfThreads, parallelRF )

  sfStop()

  # everything finished so merge all forests into one
  modelForest <- parallelForests[[1]]
  if( numberOfThreads > 1 )
    {
    for( i in 2:numberOfThreads )
      {
      modelForest <- combine( modelForest, parallelForests[[i]] )
      }
    }

  # Stop the clock
  elapsedTime <- proc.time() - ptm
  cat( "Model creation took ", as.numeric( elapsedTime[3] ), " seconds.\n", sep = "" )

  ###############################################
  #
  # Save the model
  #
  ###############################################

  save( modelForest, file = outputModelName )
  }
