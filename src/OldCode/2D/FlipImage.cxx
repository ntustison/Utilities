#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkFlipImageFilter.h"
#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc != 6 )
    {
    std::cout << "Usage: FlipImage input_filename output_filename which_axes(e.g. 0 1 1)" << std::endl;
    exit( 1 );
    }
  
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  
  typedef itk::FlipImageFilter<ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  FilterType::FlipAxesArrayType array;
  array[0] = atoi( argv[3] );
  array[1] = atoi( argv[4] );
  array[2] = atoi( argv[5] );
  filter->SetInput( reader->GetOutput() );
  filter->SetFlipAxes( array ); 
  
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return 0;
}
