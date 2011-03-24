#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

int main( int argc, char *argv[] )
{
  if ( argc < 2 )
    {
    std::cout << "Usage: " << argv[0] << " image [label]" << std::endl;
    exit( 1 );
    }

  typedef int PixelType;
  const unsigned int ImageDimension = 3;
  typedef float RealType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  PixelType label = itk::NumericTraits<PixelType>::One;
  if ( argc > 2 )
    {
    label = static_cast<PixelType>( atoi( argv[2] ) );
    } 

  itk::ImageRegionIterator<ImageType> It( reader->GetOutput(),
    reader->GetOutput()->GetLargestPossibleRegion() );

  unsigned long numberOfPixelsInsideMask = 0;

  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( It.Get() == label )
      {
      numberOfPixelsInsideMask++;
      }
    }   

  double volume = static_cast<double>( numberOfPixelsInsideMask );
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    volume *= reader->GetOutput()->GetSpacing()[i];
    }
 
  std::cout << "[" << argv[0] << "]" << std::endl;
  std::cout << volume << std::endl;

  return 0;
}
