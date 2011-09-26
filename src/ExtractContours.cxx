#include "itkBinaryContourImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLabelContourImageFilter.h"

template <unsigned int ImageDimension>
int ExtractContours( int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  if( argc < 7 )
    {
    typedef itk::LabelContourImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( reader->GetOutput() );
    if( argc > 4 )
      {
      filter->SetFullyConnected( static_cast<PixelType>( atof( argv[4] ) ) );
      }
    if( argc > 5 )
      {
      filter->SetBackgroundValue( static_cast<PixelType>( atof( argv[5] ) ) );
      }

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[3] );
    writer->SetInput( filter->GetOutput() );
    writer->Update();
    }
  else
    {
    typedef itk::BinaryContourImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( reader->GetOutput() );
    if( argc > 4 )
      {
      filter->SetFullyConnected( static_cast<PixelType>( atof( argv[4] ) ) );
      }
    if( argc > 5 )
      {
      filter->SetBackgroundValue( static_cast<PixelType>( atof( argv[5] ) ) );
      }
    if( argc > 6 )
      {
      filter->SetForegroundValue( static_cast<PixelType>( atof( argv[6] ) ) );
      }

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[3] );
    writer->SetInput( filter->GetOutput() );
    writer->Update();
    }


  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension inputImage outputImage"
      << " [fullyConnected] "
      << "[backgroundValue] [foregroundValue] " << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     ExtractContours<2>( argc, argv );
     break;
   case 3:
     ExtractContours<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
