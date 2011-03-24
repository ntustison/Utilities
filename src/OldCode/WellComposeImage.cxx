#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkWellComposedImageFilter.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 || argc > 4 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage outputImage numberOfLabels(optional) " << std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  int numberOfLabels = ( argc == 4 ) ? atoi( argv[3] ) : 2;

  typedef itk::WellComposedImageFilter<ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();

  filter->SetInput( reader->GetOutput() );
  filter->SetTotalNumberOfLabels( numberOfLabels );
  
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[2] );                                          
  writer->Update();
    
  return 0;
}
