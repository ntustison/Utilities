#include <stdio.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3  )
    {
    std::cout << "Usage: " << argv[0] << " intputImage outputImage [constant]" << std::endl;
    exit( 1 );
    }
  typedef itk::Image<PixelType, ImageDimension> ImageType;   

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  ImageType::Pointer image = ImageType::New();
  image->SetOrigin( reader->GetOutput()->GetOrigin() );
  image->SetSpacing( reader->GetOutput()->GetSpacing() );
  image->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  image->Allocate();
  if ( argc == 3 )
    {
    image->FillBuffer( 0 );
    } 
  else
    {
    image->FillBuffer( atof( argv[3] ) );
    } 
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( image );
  writer->Update();
  
  return 0;
}
