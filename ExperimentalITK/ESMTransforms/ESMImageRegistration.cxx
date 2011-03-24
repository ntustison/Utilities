#include "itkESMRigid2DTransform.h"
#include "itkESMMeanSquaresImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkESMDogLegOptimizer.h"
#include "itkImage.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSubtractImageFilter.h"



class CommandIterationUpdate : public itk::Command 
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate() {};

public:

  static const unsigned int  Dimension = 2;
  typedef  float             PixelType;

  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  typedef itk::Image< PixelType, Dimension >  MovingImageType;
  
  typedef itk::ESMDogLegOptimizer<
    FixedImageType,MovingImageType>           OptimizerType;
  typedef const OptimizerType                 *OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    OptimizerPointer optimizer = 
                         dynamic_cast< OptimizerPointer >( object );

    if( ! itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }

    std::cout << optimizer->GetCurrentIteration() << " = ";
    std::cout << optimizer->GetValue() << " : ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
  }
   
};


int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile  movingImageFile ";
    std::cerr << "outputImagefile [differenceImageAfter]";
    std::cerr << "[differenceImageBefore]" << std::endl;
    return EXIT_FAILURE;
    }

  const    unsigned int    Dimension = 2;
  typedef  float           PixelType;

  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  typedef itk::Image< PixelType, Dimension >  MovingImageType;

  typedef itk::ESMRigid2DTransform< double >  TransformType;
  
  typedef itk::ESMDogLegOptimizer<
    FixedImageType,MovingImageType>           OptimizerType;
  
  typedef itk::ESMMeanSquaresImageToImageMetric< 
    FixedImageType, MovingImageType >         MetricType;
  
  typedef itk:: LinearInterpolateImageFunction< 
    MovingImageType, double >                 InterpolatorType;
  
  MetricType::Pointer         metric        = MetricType::New();
  TransformType::Pointer      transform     = TransformType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();

  optimizer->SetESMCostFunction( metric );
  
  metric->SetESMTransform( transform );
  metric->SetInterpolator(  interpolator  );


  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName(  argv[1] );
  movingImageReader->SetFileName( argv[2] );
  
  metric->SetFixedImage(    fixedImageReader->GetOutput()    );
  metric->SetMovingImage(   movingImageReader->GetOutput()   );
  
  fixedImageReader->Update();
  metric->SetFixedImageRegion( 
                    fixedImageReader->GetOutput()->GetBufferedRegion() );
  
  typedef OptimizerType::ParametersType ParametersType;
  ParametersType initialParameters( transform->GetNumberOfParameters() );

  for ( unsigned int i=0; i<transform->GetNumberOfParameters(); ++i )
    {
    initialParameters[i] = 0.0;
    }

  optimizer->SetInitialPosition( initialParameters );
  
  optimizer->SetMaximumNumberOfIterations( 10 );

  // Connect an observer
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );
  
  try 
    { 
    optimizer->StartOptimization();
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    }

  ParametersType finalParameters = optimizer->GetOptimalParameters();
  
  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  
  const double bestValue = optimizer->GetOptimalValue();
  
  // Print out results
  //
  std::cout << "Result = " << std::endl;
  std::cout << " Parameters    = " << finalParameters    << std::endl;
  std::cout << " Stop code     = " << optimizer->GetStopConditionDescription() << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;
  
  typedef itk::ResampleImageFilter< 
                            MovingImageType, 
                            FixedImageType >    ResampleFilterType;
  
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetInput( movingImageReader->GetOutput() );
  
  resampler->SetTransform( transform );
  
  FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
  resampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
  resampler->SetOutputSpacing( fixedImage->GetSpacing() );
  resampler->SetOutputDirection( fixedImage->GetDirection() );
  resampler->SetDefaultPixelValue( 100 );
  
  typedef unsigned char OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::CastImageFilter< 
                        FixedImageType,
                        OutputImageType > CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  
  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();

  writer->SetFileName( argv[3] );
  
  caster->SetInput( resampler->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();
  
  typedef itk::SubtractImageFilter< 
                                  FixedImageType, 
                                  FixedImageType, 
                                  FixedImageType > DifferenceFilterType;

  DifferenceFilterType::Pointer difference = DifferenceFilterType::New();

  difference->SetInput1( fixedImageReader->GetOutput() );
  difference->SetInput2( resampler->GetOutput() );
  
  typedef itk::RescaleIntensityImageFilter< 
                                  FixedImageType, 
                                  OutputImageType >   RescalerType;

  RescalerType::Pointer intensityRescaler = RescalerType::New();
  
  intensityRescaler->SetInput( difference->GetOutput() );
  intensityRescaler->SetOutputMinimum(   0 );
  intensityRescaler->SetOutputMaximum( 255 );

  resampler->SetDefaultPixelValue( 1 );
  
  WriterType::Pointer writer2 = WriterType::New();
  writer2->SetInput( intensityRescaler->GetOutput() );

  if( argc > 4 )
    {
    writer2->SetFileName( argv[4] );
    writer2->Update();
    }
  
  TransformType::Pointer identityTransform = TransformType::New();
  identityTransform->SetIdentity();
  resampler->SetTransform( identityTransform );


  if( argc > 5 )
    {
    writer2->SetFileName( argv[5] );
    writer2->Update();
    }
  
  return EXIT_SUCCESS;
}

