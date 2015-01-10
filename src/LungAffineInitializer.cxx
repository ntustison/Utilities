#include "itkAffineTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkGradientDescentOptimizerv4.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkJointHistogramMutualInformationImageToImageMetricv4.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkRegistrationParameterScalesFromPhysicalShift.h"
#include "itkResampleImageFilter.h"
#include "itkTransformFileWriter.h"

#include <stdio.h>
#include <string>
#include <vector>

#include "Common.h"

int main( int argc, char *argv[] )
{

  if( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0] << " fixedImage movingImage transformName <lowerScalexUpperScalexNumberOfScaleSamples=0.5x2.0x5> <whichMetric[GCorMIorMattes]=GC> <numberOfIterations=10>" << std::endl;
    return EXIT_FAILURE;
    }

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

  const unsigned int ImageDimension = 3;

  typedef double PixelType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::Vector<PixelType, ImageDimension> VectorType;

  // read in fixed and moving images

  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer fixedReader = ReaderType::New();
  fixedReader->SetFileName( argv[1] );

  ImageType::Pointer fixedImage = ImageType::New();
  fixedImage = fixedReader->GetOutput();
  fixedImage->Update();
  fixedImage->DisconnectPipeline();

  ReaderType::Pointer movingReader = ReaderType::New();
  movingReader->SetFileName( argv[2] );

  ImageType::Pointer movingImage = ImageType::New();
  movingImage = movingReader->GetOutput();
  movingImage->Update();
  movingImage->DisconnectPipeline();

  // Get scale parameters

  std::string scaleParametersString = std::string( "0.5x2.0x5" );
  if( argc > 4 )
    {
    scaleParametersString = std::string( argv[4] );
    }
  std::vector<double> scaleParameters = ConvertVector<double>( scaleParametersString );
  if( scaleParameters.size() != 3 )
    {
    std::cerr << "The scale parameters were improperly specified.  See usage." << std::endl;
    return EXIT_FAILURE;
    }
  double scaleLowerBoundLog = vcl_log( scaleParameters[0] );
  double scaleUpperBoundLog = vcl_log( scaleParameters[1] );
  unsigned int scaleNumberOfSamples = static_cast<unsigned int>( scaleParameters[2] );
  double scaleDelta = ( scaleUpperBoundLog - scaleLowerBoundLog ) / static_cast<double>( scaleNumberOfSamples - 1 );

  // Get number of iterations

  unsigned int numberOfIterations = 10;
  if( argc > 6 )
    {
    numberOfIterations = atoi( argv[6] );
    }

  // Set up metric

  std::string metricString = "GC";
  if( argc > 5 )
    {
    metricString = std::string( argv[5] );
    }

  typedef itk::ImageToImageMetricv4<ImageType, ImageType, ImageType> ImageMetricType;
  typedef ImageMetricType::FixedSampledPointSetType PointSetType;

  ImageMetricType::Pointer imageMetric = ITK_NULLPTR;

  if( std::strcmp( metricString.c_str(), "Mattes" ) == 0 )
    {
    typedef itk::MattesMutualInformationImageToImageMetricv4<ImageType, ImageType, ImageType> MattesMetricType;
    MattesMetricType::Pointer mattesMetric = MattesMetricType::New();
    mattesMetric->SetNumberOfHistogramBins( 20 );

    imageMetric = mattesMetric;
    }
  else if( std::strcmp( metricString.c_str(), "GC" ) == 0 )
    {
    typedef itk::CorrelationImageToImageMetricv4<ImageType, ImageType, ImageType> GCMetricType;
    GCMetricType::Pointer gcMetric = GCMetricType::New();

    imageMetric = gcMetric;
    }
  else if( std::strcmp( metricString.c_str(), "MI" ) == 0 )
    {
    typedef itk::JointHistogramMutualInformationImageToImageMetricv4<ImageType, ImageType> MIMetricType;
    MIMetricType::Pointer miMetric = MIMetricType::New();
    miMetric->SetNumberOfHistogramBins( 20 );

    imageMetric = miMetric;
    }
  else
    {
    std::cerr << "Unrecognized metric option." << std::endl;
    return EXIT_FAILURE;
    }

  imageMetric->SetFixedImage( fixedImage );
  imageMetric->SetMovingImage( movingImage );
  imageMetric->SetVirtualDomainFromImage( fixedImage );
  imageMetric->SetUseMovingImageGradientFilter( false );
  imageMetric->SetUseFixedImageGradientFilter( false );

  // identity transform for fixed image

  typedef itk::IdentityTransform<double, ImageDimension> IdentityTransformType;
  IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();
  identityTransform->SetIdentity();

  imageMetric->SetFixedTransform( identityTransform );

  // Do a random sampling

  unsigned int index = 0;
  unsigned int count = 0;

  PointSetType::Pointer pointSet = PointSetType::New();
  pointSet->Initialize();

  itk::ImageRegionIteratorWithIndex<ImageType> It( fixedImage, fixedImage->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    // take every N^th point
    if ( count % 20 == 0  )
      {
      PointSetType::PointType point;
      fixedImage->TransformIndexToPhysicalPoint( It.GetIndex(), point );
      pointSet->SetPoint( index++, point );
      }
    count++;
    }
  imageMetric->SetFixedSampledPointSet( pointSet );
  imageMetric->SetUseFixedSampledPointSet( true );

  // Now go through the rotations + scalings to find the optimal pose.

  typedef itk::AffineTransform<double, ImageDimension> AffineTransformType;

  double optimalMetricValue = itk::NumericTraits<double>::max();
  AffineTransformType::Pointer optimalTransform = ITK_NULLPTR;

  // Initialize centered transform (based on the center of the image)

  AffineTransformType::Pointer initialTransform = AffineTransformType::New();

  typedef itk::CenteredTransformInitializer<AffineTransformType, ImageType, ImageType> TransformInitializerType;
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
  initializer->SetTransform( initialTransform );
  initializer->SetFixedImage( fixedImage );
  initializer->SetMovingImage( movingImage );
  initializer->GeometryOn();

  initializer->InitializeTransform();

  for( double scaleLog = scaleLowerBoundLog; scaleLog <= scaleUpperBoundLog; scaleLog += scaleDelta )
    {
    double scale = vcl_exp( scaleLog );

    AffineTransformType::MatrixType matrix = initialTransform->GetMatrix();

    AffineTransformType::Pointer affineTransform = AffineTransformType::New();
    affineTransform->SetCenter( initialTransform->GetCenter() );
    affineTransform->SetTranslation( initialTransform->GetTranslation() );
    affineTransform->SetMatrix( matrix * scale );

    imageMetric->SetMovingTransform( affineTransform );
    imageMetric->Initialize();

    if( numberOfIterations > 0 )
      {
      typedef itk::RegistrationParameterScalesFromPhysicalShift<ImageMetricType> ScalesEstimatorType;
      ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
      scalesEstimator->SetMetric( imageMetric );
      scalesEstimator->SetTransformForward( true );

      typedef itk::ConjugateGradientLineSearchOptimizerv4Template<double> ConjugateGradientDescentOptimizerType;
      ConjugateGradientDescentOptimizerType::Pointer optimizer = ConjugateGradientDescentOptimizerType::New();
      optimizer->SetLowerLimit( 0 );
      optimizer->SetUpperLimit( 2 );
      optimizer->SetEpsilon( 0.2 );
      optimizer->SetLearningRate( 0.15 );
      optimizer->SetMaximumStepSizeInPhysicalUnits( 0.15 );
      optimizer->SetNumberOfIterations( numberOfIterations );
      optimizer->SetScalesEstimator( scalesEstimator );
      optimizer->SetMinimumConvergenceValue( 1e-10 );
      optimizer->SetConvergenceWindowSize( 10 );
      optimizer->SetDoEstimateLearningRateAtEachIteration( false );
      optimizer->SetDoEstimateLearningRateOnce( true );
      optimizer->SetMetric( imageMetric );

//         typedef itk::GradientDescentOptimizerv4Template<double> GradientDescentOptimizerType;
//         GradientDescentOptimizerType::Pointer optimizer2 = GradientDescentOptimizerType::New();
//         optimizer2->SetLearningRate( 0.15 );
//         optimizer2->SetMaximumStepSizeInPhysicalUnits( 0.15 );
//         optimizer2->SetNumberOfIterations( numberOfIterations );
//         optimizer2->SetScalesEstimator( scalesEstimator );
//         optimizer2->SetMinimumConvergenceValue( 1e-6 );
//         optimizer2->SetConvergenceWindowSize( 10 );
//         optimizer2->SetDoEstimateLearningRateAtEachIteration( true );
//         optimizer2->SetDoEstimateLearningRateOnce( false );
//         optimizer2->SetMetric( imageMetric );

      try
        {
        optimizer->StartOptimization();
        }
      catch( ... )
        {
        continue;
        }
      }

    double metricValue = imageMetric->GetValue();

    if( metricValue < optimalMetricValue )
      {
      optimalMetricValue = metricValue;
      optimalTransform = affineTransform;
      }

    }

  typedef itk::TransformFileWriter TransformWriterType;
  TransformWriterType::Pointer transformWriter = TransformWriterType::New();
  transformWriter->SetInput( optimalTransform );
  transformWriter->SetFileName( argv[3] );
  transformWriter->Update();

  return EXIT_SUCCESS;
}
