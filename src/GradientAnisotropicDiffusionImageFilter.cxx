#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

template <unsigned int ImageDimension>
int GradientAnisotropicDiffusionImageFilter( int argc, char * argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType; 
  
  typename itk::ImageFileReader <ImageType>::Pointer reader;
  reader= itk::ImageFileReader <ImageType>::New();
  reader->SetFileName ( argv[2] );
  reader->Modified();
  reader->Update();

  const unsigned int numberOfIterations = atoi( argv[4] );
  const double       timeStep = atof( argv[5] );
  const double       conductance = atof( argv[6] );

  typedef itk::GradientAnisotropicDiffusionImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetNumberOfIterations( numberOfIterations );
  filter->SetTimeStep( timeStep );
  filter->SetConductanceParameter( conductance );
  filter->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();
  
  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 7 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " imageDimension inputImage  outputImage ";
    std::cerr << "numberOfIterations  timeStep conductance" << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     GradientAnisotropicDiffusionImageFilter<2>( argc, argv );
     break;
   case 3:
     GradientAnisotropicDiffusionImageFilter<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

