#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc != 4 )
    {
    std::cout << "Usage: " << argv[0] << "inputImage outputImage ballRadius" << std::endl;
    exit( 1 );
    }
  
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension> KernelType;
  KernelType ball;
  KernelType::SizeType size;
  size.Fill( atoi( argv[3] ) );
  ball.SetRadius( size );
  ball.CreateStructuringElement();
  
  typedef itk::BinaryMorphologicalClosingImageFilter<ImageType, ImageType, KernelType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetKernel( ball );
  filter->SetSafeBorder( true );
  filter->SetForegroundValue( 1 );

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return 0;
}
