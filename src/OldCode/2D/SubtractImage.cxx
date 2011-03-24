#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkSubtractImageFilter.h"

#include "global.h"

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  minuend subtrahend remainder" << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[1] );
  reader1->Update();

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  reader2->Update();

  typedef itk::SubtractImageFilter<
                            ImageType, 
                            ImageType,
                            ImageType >  SubtracterType;
  SubtracterType::Pointer subtracter = SubtracterType::New();
  subtracter->SetInput1( reader1->GetOutput() );
  subtracter->SetInput2( reader2->GetOutput() );
  subtracter->Update();

  typedef itk::ImageFileWriter<ImageType>  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( subtracter->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();

  return EXIT_SUCCESS;
}

