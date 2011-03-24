#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

template <unsigned int ImageDimension>
int ExtractMaskImage( int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[2] );
  reader1->Update();
  typename ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[3] );
  reader2->Update();

  PixelType insideValue = 1;
  if ( argc > 5 )
    {
    insideValue = static_cast<PixelType>( atof( argv[5] ) );
    }
  PixelType outsideValue = 0;
  if ( argc > 6 )
    {
    outsideValue = static_cast<PixelType>( atof( argv[6] ) );
    }

  itk::ImageRegionIterator<ImageType> It1( reader1->GetOutput(),
    reader1->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> It2( reader2->GetOutput(),
    reader2->GetOutput()->GetLargestPossibleRegion() );

  for ( It1.GoToBegin(), It2.GoToBegin(); !It1.IsAtEnd();
        ++It1, ++It2 )
    {
    if ( It2.Get() != insideValue )
      {
      It1.Set( outsideValue );
      }
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( reader1->GetOutput() );
  writer->Update();
  
  return 0;
}


int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension inputImage "
      << "inputMaskImage outputImage [insideValue] [outsideValue]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     ExtractMaskImage<2>( argc, argv );
     break;
   case 3:
     ExtractMaskImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
