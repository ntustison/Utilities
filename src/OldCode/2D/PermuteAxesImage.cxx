#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkPermuteAxesImageFilter.h"
#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc != ImageDimension+3 )
    {
    std::cout << "Usage: PermuteAxesImage input_filename output_filename axes_order( e.g. 1 2 0 for 3-D )" << std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  
  typedef itk::PermuteAxesImageFilter<ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  FilterType::PermuteOrderArrayType array;
  array[0] = atoi( argv[3] );
  array[1] = atoi( argv[4] );
  if ( ImageDimension == 3 )
    {
    array[2] = atoi( argv[5] );
    }
  filter->SetInput( reader->GetOutput() );
  filter->SetOrder( array ); 
  
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return 0;
}
