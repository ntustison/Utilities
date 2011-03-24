#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkZeroCrossingBasedEdgeDetectionImageFilter.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage outputImage [variance dim1] .. [variance dimN]" << std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::ZeroCrossingBasedEdgeDetectionImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();

  filter->SetInput( reader->GetOutput() );
  filter->SetBackgroundValue( 0 );
  filter->SetForegroundValue( 1 );

  FilterType::ArrayType array;

  if ( argc > 3 )
    {
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      array[d] = atof( argv[3+d] );
      }
    filter->SetVariance( array );
    }
  filter->Update();
  
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[2] );                                          
  writer->Update();
    
  return 0;
}
