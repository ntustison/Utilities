#include "itkImage.h"
#include "itkIdentityTransform.h"
#include "itkGeodesicActiveContourShapePriorLevelSetImageFilter.h"
#include "itkPCAShapeSignedDistanceFunction.h"
#include "itkShapePriorMAPCostFunction.h"
#include "itkOnePlusOneEvolutionaryOptimizer.h"
#include "itkNormalVariateGenerator.h"
#include "itkNumericSeriesFileNames.h"

#include "itkFastMarchingImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSpatialFunctionImageEvaluatorFilter.h"


#include "vnl/vnl_sample.h"

#include "itkCommand.h"

template<class TFilter>
class CommandIterationUpdate : public itk::Command
{
public:
  typedef CommandIterationUpdate   Self;
  typedef itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() {};
public:

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *) caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
    const TFilter * filter =
      dynamic_cast< const TFilter * >( object );
    if( typeid( event ) != typeid( itk::IterationEvent ) )
      { return; }

    std::cout << filter->GetElapsedIterations() << ": ";
    std::cout << filter->GetRMSChange() << " ";
    std::cout << filter->GetCurrentParameters() << std::endl;
    }

};

template<unsigned int ImageDimension>
int ShapePriorLevelSet( int argc, char *argv[] )
{
  typedef float RealType;
  typedef RealType PixelType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef unsigned char LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;

  typedef  itk::ImageFileReader<RealImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[4] );
  reader->Update();

  /**
   * Set up the PCA shape function
   */

  // Read in the pca images.  startIndex is meanImage and startIndex+1:endIndex
  //   are the principal components.
  unsigned int numberOfComponents = atoi( argv[3] );

  typedef itk::PCAShapeSignedDistanceFunction
    <double, ImageDimension, RealImageType> ShapeFunctionType;
  typename ShapeFunctionType::Pointer shape = ShapeFunctionType::New();
  shape->SetNumberOfPrincipalComponents( numberOfComponents );

  itk::NumericSeriesFileNames::Pointer fileNamesCreator =
          itk::NumericSeriesFileNames::New();
  fileNamesCreator->SetStartIndex( 0 );
  fileNamesCreator->SetEndIndex( numberOfComponents );
  fileNamesCreator->SetSeriesFormat( argv[2] );
  const std::vector<std::string> & pcaNames = fileNamesCreator->GetFileNames();

  std::vector<typename RealImageType::Pointer>
    shapeModeImages( pcaNames.size() - 1 );
  for ( unsigned int k = 0; k < pcaNames.size(); k++ )
    {
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( pcaNames[k].c_str() );
    reader->Update();
    if( k == 0 )
      {
      shape->SetMeanImage( reader->GetOutput() );
      }
    else
      {
      shapeModeImages[k-1] = reader->GetOutput();
      }
    }
  shape->SetPrincipalComponentImages( shapeModeImages );

  typename ShapeFunctionType::ParametersType
    pcaStandardDeviations( numberOfComponents );
  pcaStandardDeviations.Fill( 1.0 );
  shape->SetPrincipalComponentStandardDeviations( pcaStandardDeviations );

  typedef itk::IdentityTransform<double, ImageDimension> TransformType;
  typename TransformType::Pointer transform = TransformType::New();
  shape->SetTransform( transform );

  typedef itk::ShapePriorMAPCostFunction<RealImageType, PixelType> CostFunctionType;
  typename CostFunctionType::Pointer costFunction = CostFunctionType::New();
  typename CostFunctionType::WeightsType weights;

  weights[0] = ( argc > 10 ) ? atof( argv[10] ) : 1.0; // weight for contour fit term
  weights[1] = ( argc > 11 ) ? atof( argv[11] ) : 20.0; // weight for image fit term
  weights[2] = ( argc > 12 ) ? atof( argv[12] ) : 1.0;  // weight for shape prior term
  weights[3] =  0.0;  // weight for pose prior term

  costFunction->SetWeights( weights );

  typename CostFunctionType::ArrayType mean(
    shape->GetNumberOfShapeParameters() );
  typename CostFunctionType::ArrayType stddev(
    shape->GetNumberOfShapeParameters() );
  mean.Fill( 0.0 );
  stddev.Fill( 1.0 );
  costFunction->SetShapeParameterMeans( mean );
  costFunction->SetShapeParameterStandardDeviations( stddev );

  /**
   * Set up the optimizer
   */

  typedef itk::OnePlusOneEvolutionaryOptimizer OptimizerType;
  typename OptimizerType::Pointer optimizer = OptimizerType::New();

  typedef itk::Statistics::NormalVariateGenerator GeneratorType;
  GeneratorType::Pointer generator = GeneratorType::New();
  generator->Initialize( 20020702 );

  optimizer->SetNormalVariateGenerator( generator );
  typename OptimizerType::ScalesType scales( shape->GetNumberOfParameters() );
  scales.Fill( 1.0 );
  for( unsigned int k = 0; k < numberOfComponents; k++ )
    {
    scales[k] = 20.0;  // scales for the pca mode multiplier
    }
  optimizer->SetScales( scales );

  double initRadius = 1.05;
  double grow = 1.1;
  double shrink = pow(grow, -0.25);

  optimizer->Initialize( initRadius, grow, shrink );
  optimizer->SetEpsilon( 1.0e-6 );
  optimizer->SetMaximumIteration( 15 );

  typename ShapeFunctionType::ParametersType
    parameters( shape->GetNumberOfParameters() );
  parameters.Fill( 0.0 );

  /**
   * now put everything in the prior shape level set filter.
   */
  typename ReaderType::Pointer levelSetReader = ReaderType::New();
  levelSetReader->SetFileName( argv[5] );
  levelSetReader->Update();

  unsigned int numberOfIterations = ( argc > 7 ) ? atoi( argv[7] ): 100;
  const double shapePriorScaling = ( argc > 8 ) ? atof( argv[8] ) : 5.0;
  const double propagationScaling = ( argc > 9 ) ? atof( argv[9] ) : 1.0;

  typedef itk::GeodesicActiveContourShapePriorLevelSetImageFilter
    <RealImageType, RealImageType>   GeodesicActiveContourFilterType;
  typename GeodesicActiveContourFilterType::Pointer priorFilter
    = GeodesicActiveContourFilterType::New();
  priorFilter->SetPropagationScaling( propagationScaling );
  priorFilter->SetShapePriorScaling( shapePriorScaling );
  priorFilter->SetCurvatureScaling( 1.0 );
  priorFilter->SetAdvectionScaling( 1.0 );
  priorFilter->SetMaximumRMSError( 0.005 );
  priorFilter->SetNumberOfIterations( numberOfIterations );
  priorFilter->SetNumberOfLayers( 4 );
  priorFilter->SetInput( levelSetReader->GetOutput() );
  priorFilter->SetFeatureImage( reader->GetOutput() );
  priorFilter->SetShapeFunction( shape );
  priorFilter->SetCostFunction( costFunction );
  priorFilter->SetOptimizer( optimizer );
  priorFilter->SetInitialParameters( parameters );

  typedef CommandIterationUpdate<GeodesicActiveContourFilterType> CommandType;
  typename CommandType::Pointer observer = CommandType::New();
  priorFilter->AddObserver( itk::IterationEvent(), observer );

  try
    {
    priorFilter->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

  std::cout << std::endl;
  std::cout << "Max. no. iterations: " << priorFilter->GetNumberOfIterations() << std::endl;
  std::cout << "Max. RMS error: " << priorFilter->GetMaximumRMSError() << std::endl;
  std::cout << std::endl;
  std::cout << "No. elpased iterations: " << priorFilter->GetElapsedIterations() << std::endl;
  std::cout << "RMS change: " << priorFilter->GetRMSChange() << std::endl;
  std::cout << "Parameters: " << priorFilter->GetCurrentParameters() << std::endl;

  shape->SetParameters( priorFilter->GetCurrentParameters() );

  typedef itk::SpatialFunctionImageEvaluatorFilter
    <ShapeFunctionType, RealImageType, RealImageType> EvaluatorFilterType;
  typename EvaluatorFilterType::Pointer evaluator = EvaluatorFilterType::New();
  evaluator->SetInput( priorFilter->GetOutput() );
  evaluator->SetFunction( shape );
  evaluator->Update();

  typedef itk::BinaryThresholdImageFilter<RealImageType, LabelImageType> ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( evaluator->GetOutput() );
  thresholder->SetLowerThreshold( -1000.0 );
  thresholder->SetUpperThreshold( 0.0 );
  thresholder->SetOutsideValue( 0 );
  thresholder->SetInsideValue( 1 );
  thresholder->Update();

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( thresholder->GetOutput() );
  writer->SetFileName( argv[6] );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 7 )
    {
    std::cerr << argv[0] << " imageDimension pcaFilesFormat "
      << "numberOfPCAComponents featureImage initialLevelSetImage "
      << "outputImage [numberOfIterations] [shapePriorScaling] "
      << "[propogationScaling] [contourWeighting] "
      << "[imageWeighting] [shapePriorWeighting] "
      << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     ShapePriorLevelSet<2>( argc, argv );
     break;
   case 3:
     ShapePriorLevelSet<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
