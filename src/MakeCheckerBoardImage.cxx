#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCheckerBoardImageFilter.h"

int main( int argc, char *argv[] )
{
  if ( argc != 5 )
    {
    std::cout << "Usage: MakeCheckerBoardImage image1_filename image2_filename number_of_boxes_per_dimension checkerimage_filename" << std::endl;
    exit( 1 );
    }

  const unsigned int ImageDimension = 2;
  typedef float PixelType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[1] );
  reader1->Update();

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  reader2->Update();

  typedef itk::CheckerBoardImageFilter<ImageType> CheckerType;
  CheckerType::Pointer checker = CheckerType::New();
  CheckerType::PatternArrayType pattern;
  pattern.Fill( atoi( argv[3] ) );
  checker->SetInput1( reader1->GetOutput() );
  checker->SetInput2( reader2->GetOutput() );
  checker->SetCheckerPattern( pattern );
  checker->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetInput( checker->GetOutput() );
  writer->SetFileName( argv[4] );
  writer->Update();

  return 0;
}
