#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "global.h"

int main ( const int argc, const char * argv[] )
{
  if ( argc < 6 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  inputImageFile  outputImageFile ";
    std::cerr << "numberOfIterations  timeStep  conductance" << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType; 
  
  itk::ImageFileReader <ImageType>::Pointer reader;
  reader= itk::ImageFileReader <ImageType>::New();
  reader->SetFileName ( argv[1] );
  reader->Modified();
  reader->Update();

  const unsigned int numberOfIterations = atoi( argv[3] );
  const double       timeStep = atof( argv[4] );
  const double       conductance = atof( argv[5] );

  typedef itk::GradientAnisotropicDiffusionImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetNumberOfIterations( numberOfIterations );
  filter->SetTimeStep( timeStep );
  filter->SetConductanceParameter( conductance );
  filter->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[2] );
  writer->Update();
  
  return 0;
}
