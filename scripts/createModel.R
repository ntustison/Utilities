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
  cat( "Usage: Rscript createModel.R dimension inputFileList outputModelPrefix ",
       "<numberOfThreads=4> <trainingPortion=1.0> <numberOfTreesPerThread=1000> ",
       "<numberOfSamplesPerLabel=1000> <numberOfUniqueLabels=NA>", sep = "" )
  stopQuietly()
  }

dimension <- as.numeric( args[1] )
fileList <- read.csv( args[2] )
outputModelName <- paste0( args[3], ".RData" )

numberOfThreads <- 4
if( length( args ) >= 4 )
  {
  numberOfThreads <- as.numeric( args[4] )
  }
trainingPortion <- 1.0
if( length( args ) >= 5 )
  {
  trainingPortion <- as.numeric( args[5] )
  }
numberOfSamplesPerLabel <- 1000
if( length( args ) >= 6 )
  {
  numberOfSamplesPerLabel <- as.numeric( args[6] )
  }
numberOfTreesPerThread <- 1000
if( length( args ) >= 7 )
  {
  numberOfTreesPerThread <- as.numeric( args[7] )
  }
numberOfUniqueLabels <- NA
if( length( args ) >= 8 )
  {
  numberOfUniqueLabels <- as.numeric( args[8] )
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

  mask <- as.array( antsImageRead( as.character( masks[i] ), dimension = dimension, pixeltype = 'unsigned int' ) )
  truth <- as.array( antsImageRead( as.character( truthLabels[i] ), dimension = dimension, pixeltype = 'unsigned int' ) )
  if( is.na( numberOfUniqueLabels ) )
    {
    uniqueTruthLabels <- sort( unique( truth[which( mask == 1 )] ) )
    uniqueTruthLabels <- uniqueTruthLabels[which( uniqueTruthLabels != 0 )]
    }
  else
    {
    uniqueTruthLabels <- 1:numberOfUniqueLabels
    }
  cat( "Unique truth labels: ", uniqueTruthLabels, "\n", sep = " " )

  truthLabelIndices <- list()
  numberOfSamplesPerLabelInSubjectData <- rep( 0, length( uniqueTruthLabels ) )
  for( n in 1:length( uniqueTruthLabels ) )
    {
    labelIndices <- which( truth == uniqueTruthLabels[n] )
    numberOfSamplesPerLabelInSubjectData[n] <- min( numberOfSamplesPerLabel, length( labelIndices ) )
    if( length( labelIndices ) > 0 )
      {
      truthLabelIndices[[n]] <- labelIndices[sample.int( length( labelIndices ), numberOfSamplesPerLabelInSubjectData[n], replace = FALSE )]
      }
    }

  subjectData <- matrix( NA, nrow = sum( numberOfSamplesPerLabelInSubjectData ), ncol = length( featureNames ) + 1 )
  for( j in 1:length( featureNames ) )
    {
    cat( "  Reading feature image ", featureNames[j], ".\n", sep = "" )
    featureImage <- as.array( antsImageRead( as.character( featureImages[i,j] ), dimension = dimension, pixeltype = 'float' ) )
    for( n in 1:length( uniqueTruthLabels ) )
      {
      if( numberOfSamplesPerLabelInSubjectData[n] == 0 )
        {
        next
        }

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

# outputModelDataFileName <- paste0( args[2], "ModelData.RData" )
# save( modelData, file = outputModelDataFileName )

###############################################
#
# Create the random forest model in parallel
#
###############################################

cat( "\nCreating the RF model.  ", sep = "" )

# Start the clock!
ptm <- proc.time()

modelFormula <- as.formula( "Labels ~ . " )

#the function each thread calls
parallelRF <- function( i ) {
#   modelData.imputed <- rfImpute( modelFormula, modelData )
#   return( randomForest( modelFormula, modelData.imputed, ntree = numberOfTreesPerThread, type = classification ) )
  return( randomForest( modelFormula, modelData, ntree = numberOfTreesPerThread, type = classification, na.action = na.omit ) )
}

if( numberOfThreads == 1 )
  {
#   modelData.imputed <- rfImpute( modelFormula, modelData )
#   modelForest <- randomForest( modelFormula, modelData.imputed, ntree = numberOfTreesPerThread, type = classification )
  modelForest <- randomForest( modelFormula, modelData, ntree = numberOfTreesPerThread, type = classification, na.action = na.omit )

  # Stop the clock
  elapsedTime <- proc.time() - ptm
  cat( "Model creation took ", as.numeric( elapsedTime[3] ), " seconds.\n", sep = "" )

  ###############################################
  #
  # Save the model
  #
  ###############################################

  save( modelForest, file = outputModelName )

  } else {
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
