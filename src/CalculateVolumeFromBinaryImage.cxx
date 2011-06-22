#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

template <unsigned int ImageDimension>
int CalculateVolumeFromBinaryImage( int argc, char *argv[] )
{
  typedef int PixelType;
  typedef float RealType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  PixelType label = itk::NumericTraits<PixelType>::One;
  if ( argc > 3 )
    {
    label = static_cast<PixelType>( atoi( argv[3] ) );
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
  if( argc > 4 )
    {
    PixelType label = argc > 5 ? static_cast<PixelType>( atof( argv[5] ) ) : 1; 
     
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[4] );
    reader->Update();
     
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
     
    double volume2 = static_cast<double>( numberOfPixelsInsideMask );
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      volume2 *= reader->GetOutput()->GetSpacing()[i];
      }
      
    volume /= volume2;  
    }
   
  std::cout << "[" << argv[0] << "]" << std::endl;
  std::cout << volume << std::endl;

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension inputImage [label] " 
      "[image2] [label2]" << std::endl;
    std::cerr << "If image2 is specified, we calculate volume(inputImage)/volume(image2)." << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     CalculateVolumeFromBinaryImage<2>( argc, argv );
     break;
   case 3:
     CalculateVolumeFromBinaryImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

