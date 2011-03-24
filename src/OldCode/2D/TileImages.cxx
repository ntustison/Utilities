#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkTileImageFilter.h"
#include "global.h"

int main( unsigned int argc, char *argv[] )
{
  if ( argc == 1 )
    {
    std::cout << argv[0] << " outputImage layout[0] ... layout[d-1] inputImage1 ... inputImageN" << std::endl;
    exit( 1 );
    }
  
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::TileImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  FilterType::LayoutArrayType array;
  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    array[d] = atoi( argv[2+d] );
    }
  filter->SetLayout( array ); 

  for ( unsigned int n = 2+ImageDimension; n < argc; n++ )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[n] );
    reader->Update();
 
    filter->SetInput( n-(2+ImageDimension), reader->GetOutput() );
    }
  filter->Update();
  
  
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[1] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return 0;
}
