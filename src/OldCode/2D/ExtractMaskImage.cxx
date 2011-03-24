#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage inputMaskImage outputImage [insideValue] [outsideValue]" << std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[1] );
  reader1->Update();
  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  reader2->Update();

  PixelType insideValue = 1;
  if ( argc > 4 )
    {
    insideValue = static_cast<PixelType>( atof( argv[4] ) );
    }
  PixelType outsideValue = 0;
  if ( argc > 5 )
    {
    outsideValue = static_cast<PixelType>( atof( argv[5] ) );
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
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( reader1->GetOutput() );
  writer->Update();


  return 0;
}
