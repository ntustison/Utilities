#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkPermuteAxesImageFilter.h"

#include "Common.h"

template<unsigned int ImageDimension>
int PermuteAxesImage( int argc, char *argv[] )
{

  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::PermuteAxesImageFilter<ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  typename FilterType::PermuteOrderArrayType array;

  std::vector<unsigned int> which
    = ConvertVector<unsigned int>( std::string( argv[4] ) );
  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    array[d] = which[d];
    }

  filter->SetInput( reader->GetOutput() );
  filter->SetOrder( array );

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << argv[0] << " imageDimension inputImage outputImage axesOrder[e.g.1x2x0]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     PermuteAxesImage<2>( argc, argv );
     break;
   case 3:
     PermuteAxesImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

