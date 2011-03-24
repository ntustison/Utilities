#include "itkAdaptiveHistogramEqualizationImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

template <unsigned int ImageDimension>
int AdaptiveHistogramEqualizeImage( int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType; 
  
  typename itk::ImageFileReader <ImageType>::Pointer reader;
  reader= itk::ImageFileReader <ImageType>::New();
  reader->SetFileName ( argv[2] );
  reader->Modified();
  reader->Update();

  typename ImageType::SizeType radius;
  radius.Fill( atoi( argv[6] ) );
  
  typename itk::AdaptiveHistogramEqualizationImageFilter<ImageType>::Pointer filter;
  filter = itk::AdaptiveHistogramEqualizationImageFilter<ImageType>::New();
  filter->SetInput ( reader->GetOutput() );
  filter->SetAlpha( atof( argv[4] ) );
  filter->SetBeta( atof( argv[5] ) );
  filter->SetRadius( radius );
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
      std::cout << "Usage " << argv[0] << " imageDimension" 
                << " inputImage outputImage"
                << " alpha(0 = classical histogram equalization -> 1 = unsharp mask)"
                << " beta(0 = unsharp mask -> 1 = pass through)"
                << " radius"
                << std::endl;
      return 1;
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     AdaptiveHistogramEqualizeImage<2>( argc, argv );
     break;
   case 3:
     AdaptiveHistogramEqualizeImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}



