#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkAbsImageFilter.h"
#include "itkSubtractImageFilter.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc != 4 )
    {
    std::cout << "Usage: " << argv[0] << " image1 image2 outputImage" << std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[1] );  
  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );  

  typedef itk::SubtractImageFilter<ImageType, ImageType, ImageType> SubtracterType;
  SubtracterType::Pointer subtracter = SubtracterType::New();
  subtracter->SetInput1( reader1->GetOutput() );
  subtracter->SetInput2( reader2->GetOutput() );
  subtracter->Update(); 

  typedef itk::AbsImageFilter<ImageType, ImageType> AbsFilterType;
  AbsFilterType::Pointer abs = AbsFilterType::New();
  abs->SetInput( subtracter->GetOutput() );
  abs->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( abs->GetOutput() );
  writer->Update();

  return 0;
}
