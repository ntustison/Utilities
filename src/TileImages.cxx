#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkTileImageFilter.h"

#include "Common.h"

template <unsigned int ImageDimension>
int TileImages( unsigned int argc, char *argv[] )
{

  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::TileImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  typename FilterType::LayoutArrayType array;

  std::vector<unsigned int> layout = ConvertVector<unsigned int>( std::string( argv[3] ) );
  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    array[d] = layout[d];
    }
  filter->SetLayout( array );

  for ( unsigned int n = 4; n < argc; n++ )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[n] );
    reader->Update();

    filter->SetInput( n-4, reader->GetOutput() );
    }
  filter->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension outputImage layout inputImage1 ... inputImageN" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     TileImages<2>( argc, argv );
     break;
   case 3:
     TileImages<3>( argc, argv );
     break;
   case 4:
     TileImages<4>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

