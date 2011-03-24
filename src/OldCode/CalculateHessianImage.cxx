#include <stdio.h>

#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkImage.h"

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << "inputImage outputImage [sigma]"<< std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::HessianRecursiveGaussianImageFilter<ImageType> HessianFilterType;    
  HessianFilterType::Pointer hessian = HessianFilterType::New();
  hessian->SetInput( reader->GetOutput() );
  hessian->SetNormalizeAcrossScale( false );
  if ( argv[3] )
    {
    hessian->SetSigma( atof( argv[3] ) );   
    }  
  hessian->SetNumberOfThreads( 1 );
  hessian->Update();

  typedef itk::ImageFileWriter<HessianFilterType::OutputImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( hessian->GetOutput() );
  writer->SetFileName( argv[2] );
  writer->Update();

  return 0;
}
