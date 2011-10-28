#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkSimpleImageRegistrationMethod.h"

#include "itkANTSNeighborhoodCorrelationImageToImageObjectMetric.h"
#include "itkCompositeTransform.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransform.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor.h"
#include "itkGradientDescentObjectOptimizer.h"
#include "itkIdentityTransform.h"
#include "itkRegistrationParameterScalesFromShift.h"
#include "itkResampleImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkTransform.h"
#include "itkVector.h"

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

    unsigned int currentLevel = filter->GetCurrentLevel();
    typename TFilter::NumberOfIterationsArrayType numIterations = filter->GetNumberOfIterationsPerLevel();
    typename TFilter::LearningRatesArrayType learningRates = filter->GetLearningRatesPerLevel();
    typename TFilter::ShrinkFactorsArrayType shrinkFactors = filter->GetShrinkFactorsPerLevel();
    typename TFilter::SmoothingSigmasArrayType smoothingSigmas = filter->GetSmoothingSigmasPerLevel();
    typename TFilter::TransformParametersAdaptorsContainerType adaptors = filter->GetTransformParametersAdaptorsPerLevel();

    std::cout << "  Current level = " << currentLevel << std::endl;
    std::cout << "    number of iterations = " << numIterations[currentLevel] << std::endl;
    std::cout << "    learning rate = " << learningRates[currentLevel] << std::endl;
    std::cout << "    shrink factor = " << shrinkFactors[currentLevel] << std::endl;
    std::cout << "    smoothing sigma = " << smoothingSigmas[currentLevel] << std::endl;
    std::cout << "    required fixed parameters = " << adaptors[currentLevel]->GetRequiredFixedParameters() << std::endl;
    }
};

template <unsigned int ImageDimension>
int SimpleRegistration( unsigned int argc, char *argv[] )
{

  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> FixedImageType;
  typedef itk::Image<PixelType, ImageDimension> MovingImageType;

  typedef itk::ImageFileReader<FixedImageType> ImageReaderType;

  typename ImageReaderType::Pointer fixedImageReader = ImageReaderType::New();
  fixedImageReader->SetFileName( argv[2] );
  fixedImageReader->Update();
  typename FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
  fixedImage->Update();
  fixedImage->DisconnectPipeline();

  typename ImageReaderType::Pointer movingImageReader = ImageReaderType::New();
  movingImageReader->SetFileName( argv[3] );
  movingImageReader->Update();
  typename MovingImageType::Pointer movingImage = movingImageReader->GetOutput();
  movingImage->Update();
  movingImage->DisconnectPipeline();

  typedef itk::SimpleImageRegistrationMethod<FixedImageType, MovingImageType> AffineRegistrationType;
  typename AffineRegistrationType::Pointer affineSimple = AffineRegistrationType::New();
  affineSimple->SetFixedImage( fixedImage );
  affineSimple->SetMovingImage( movingImage );

  typedef CommandIterationUpdate<AffineRegistrationType> AffineCommandType;
  typename AffineCommandType::Pointer affineObserver = AffineCommandType::New();
  affineSimple->AddObserver( itk::IterationEvent(), affineObserver );

  try
    {
    std::cout << "Affine transform" << std::endl;
    affineSimple->StartRegistration();
    }
  catch( itk::ExceptionObject &e )
    {
    std::cerr << "Exception caught: " << e << std::endl;
    return EXIT_FAILURE;
    }

  //
  // Now do the displacement field transform with gaussian smoothing using
  // the composite transform.
  //

  typedef typename AffineRegistrationType::RealType RealType;

  typedef itk::CompositeTransform<RealType, ImageDimension> CompositeTransformType;
  typename CompositeTransformType::Pointer compositeTransform = CompositeTransformType::New();
  compositeTransform->AddTransform( const_cast<typename AffineRegistrationType::TransformType *>( affineSimple->GetOutput()->Get() ) );

  typedef itk::Vector<RealType, ImageDimension> VectorType;
  VectorType zeroVector( 0.0 );
  typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;
  typename DisplacementFieldType::Pointer displacementField = DisplacementFieldType::New();
  displacementField->CopyInformation( fixedImageReader->GetOutput() );
  displacementField->SetRegions( fixedImage->GetBufferedRegion() );
  displacementField->Allocate();
  displacementField->FillBuffer( zeroVector );

  typedef itk::GaussianSmoothingOnUpdateDisplacementFieldTransform<RealType, ImageDimension> DisplacementFieldTransformType;
  typename DisplacementFieldTransformType::Pointer fieldTransform = DisplacementFieldTransformType::New();
  fieldTransform->SetGaussianSmoothingVarianceForTheUpdateField( 0 );
  fieldTransform->SetGaussianSmoothingVarianceForTheTotalField( 1.5 );
  fieldTransform->SetDisplacementField( displacementField );

  typedef itk::ANTSNeighborhoodCorrelationImageToImageObjectMetric<FixedImageType, MovingImageType> CorrelationMetricType;
  typename CorrelationMetricType::Pointer correlationMetric = CorrelationMetricType::New();
  typename CorrelationMetricType::RadiusType radius;
  radius.Fill( 2 );
  correlationMetric->SetRadius( radius );
  correlationMetric->SetDoFixedImagePreWarp( true );
  correlationMetric->SetDoMovingImagePreWarp( true );
  correlationMetric->SetUseMovingImageGradientFilter( false );
  correlationMetric->SetUseFixedImageGradientFilter( false );

  typedef itk::SimpleImageRegistrationMethod<FixedImageType, MovingImageType, DisplacementFieldTransformType> DisplacementFieldRegistrationType;
  typename DisplacementFieldRegistrationType::Pointer displacementFieldSimple = DisplacementFieldRegistrationType::New();
  displacementFieldSimple->SetFixedImage( fixedImage );
  displacementFieldSimple->SetMovingImage( movingImage );
  displacementFieldSimple->SetNumberOfLevels( 3 );
  displacementFieldSimple->SetCompositeTransform( compositeTransform );
  displacementFieldSimple->SetTransform( fieldTransform );
  displacementFieldSimple->SetScalesEstimator( NULL );
  displacementFieldSimple->SetMetric( correlationMetric );

  typename DisplacementFieldRegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
  shrinkFactorsPerLevel.SetSize( 3 );
  shrinkFactorsPerLevel[0] = 2;
  shrinkFactorsPerLevel[1] = 2;
  shrinkFactorsPerLevel[2] = 1;
  displacementFieldSimple->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );

  typename DisplacementFieldRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
  smoothingSigmasPerLevel.SetSize( 3 );
  smoothingSigmasPerLevel[0] = 2;
  smoothingSigmasPerLevel[1] = 1;
  smoothingSigmasPerLevel[2] = 0;
  displacementFieldSimple->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );

  typename DisplacementFieldRegistrationType::NumberOfIterationsArrayType numberOfIterationsPerLevel;
  numberOfIterationsPerLevel.SetSize( 3 );
  numberOfIterationsPerLevel[0] = 10;
  numberOfIterationsPerLevel[1] = 100;
  numberOfIterationsPerLevel[2] = 200;
  displacementFieldSimple->SetNumberOfIterationsPerLevel( numberOfIterationsPerLevel );

  typename DisplacementFieldRegistrationType::LearningRatesArrayType learningRatesPerLevel;
  learningRatesPerLevel.SetSize( 3 );
  learningRatesPerLevel[0] = 50;
  learningRatesPerLevel[1] = 50;
  learningRatesPerLevel[2] = 50;
  displacementFieldSimple->SetLearningRatesPerLevel( learningRatesPerLevel );

  typedef itk::GaussianSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor<DisplacementFieldTransformType> DisplacementFieldTransformAdaptorType;

  typename FixedImageType::SpacingType spacingAtFullResolution = fixedImage->GetSpacing();
  typename FixedImageType::SizeType sizeAtFullResolution = fixedImage->GetBufferedRegion().GetSize();

  typename DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;

  for( unsigned int level = 0; level < shrinkFactorsPerLevel.Size(); level++ )
    {
    typename DisplacementFieldType::SpacingType spacing;
    typename DisplacementFieldType::SizeType size;

    typedef itk::ShrinkImageFilter<DisplacementFieldType, DisplacementFieldType> ShrinkFilterType;
    typename ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
    shrinkFilter->SetShrinkFactors( shrinkFactorsPerLevel[level] );
    shrinkFilter->SetInput( displacementField );
    shrinkFilter->Update();

    typename DisplacementFieldTransformAdaptorType::Pointer fieldTransformAdaptor = DisplacementFieldTransformAdaptorType::New();
    fieldTransformAdaptor->SetRequiredSpacing( shrinkFilter->GetOutput()->GetSpacing() );
    fieldTransformAdaptor->SetRequiredSize( shrinkFilter->GetOutput()->GetBufferedRegion().GetSize() );
    fieldTransformAdaptor->SetRequiredDirection( shrinkFilter->GetOutput()->GetDirection() );
    fieldTransformAdaptor->SetRequiredOrigin( shrinkFilter->GetOutput()->GetOrigin() );

    adaptors.push_back( fieldTransformAdaptor.GetPointer() );
    }
  displacementFieldSimple->SetTransformParametersAdaptorsPerLevel( adaptors );

  typedef CommandIterationUpdate<DisplacementFieldRegistrationType> DisplacementFieldRegistrationCommandType;
  typename DisplacementFieldRegistrationCommandType::Pointer displacementFieldObserver = DisplacementFieldRegistrationCommandType::New();
  displacementFieldSimple->AddObserver( itk::IterationEvent(), displacementFieldObserver );

  try
    {
    std::cout << "Displacement field transform (gaussian update)" << std::endl;
    displacementFieldSimple->StartRegistration();
    }
  catch( itk::ExceptionObject &e )
    {
    std::cerr << "Exception caught: " << e << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::ResampleImageFilter<MovingImageType, FixedImageType> ResampleFilterType;
  typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetTransform( compositeTransform );
  resampler->SetInput( movingImage );
  resampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
  resampler->SetOutputSpacing( fixedImage->GetSpacing() );
  resampler->SetOutputDirection( fixedImage->GetDirection() );
  resampler->SetDefaultPixelValue( 0 );
  resampler->Update();

  typedef itk::ImageFileWriter<FixedImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( resampler->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << argv[0] << " imageDimension fixedImage movingImage outputImage" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     SimpleRegistration<2>( argc, argv );
     break;
   case 3:
     SimpleRegistration<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

