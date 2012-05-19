#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkHistogramMatchingImageFilter.h"
#include "itkDynamicHistogramWarpingImageFilter.h"

template <unsigned int ImageDimension>
int HistogramMatchImages( int argc, char * argv[] )
{
  typedef float PixelType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[2] );
  reader1->Update();

  typename ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[3] );
  reader2->Update();

  typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetSourceImage( reader1->GetOutput() );
  filter->SetReferenceImage( reader2->GetOutput() );
  filter->ThresholdAtMeanIntensityOn();
  filter->SetNumberOfHistogramLevels( ( argc > 5 ) ? atoi( argv[5] ) : 255 );
  filter->SetNumberOfMatchPoints( ( argc > 6 ) ? atoi( argv[6] ) : 12 );
  filter->Update();

//   typedef itk::DynamicHistogramWarpingImageFilter<ImageType, ImageType> FilterType;
//   typename FilterType::Pointer filter = FilterType::New();
//   filter->SetSourceImage( reader1->GetOutput() );
//   filter->SetReferenceImage( reader2->GetOutput() );
//   filter->SetNumberOfHistogramLevels( ( argc > 5 ) ? atoi( argv[5] ) : 255 );
//   filter->SetMaximumSourceBinCompressionSize(
//     ( argc > 6 ) ? atoi( argv[6] ) : 10 );
//   filter->SetMaximumReferenceBinCompressionSize(
//     ( argc > 7 ) ? atoi( argv[7] ) : 10 );
//   filter->Update();

  typedef itk::ImageFileWriter<ImageType>  WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[4] );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
//     std::cout << "Usage: " << argv[0] << " imageDimension "
//       << "sourceImage referenceImage outputImage [histLevels] "
//       << "[sourceCompressionSize] [referenceCompressionSize]" << std::endl;
    std::cout << "Usage: " << argv[0] << " imageDimension "
      << "sourceImage referenceImage outputImage [histLevels] [numberOfMatchPoints] "
      << std::endl;
    exit( 1 );

    }


  switch( atoi( argv[1] ) )
   {
   case 2:
     HistogramMatchImages<2>( argc, argv );
     break;
   case 3:
     HistogramMatchImages<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
