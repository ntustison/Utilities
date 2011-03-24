#include <stdio.h>

#include "itkGaussianImageSource.h"
#include "itkImageFileWriter.h"
#include "itkTimeProbe.h"
#include "itkBSplineKernelFunction.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc != 3 + ImageDimension*5 )
    {
    std::cout << "Usage: " << argv[0] << " outputImageFilename scale ";
    std::cout << "origin[0] ... origin[n-1] ";  
    std::cout << "size[0] ... size[n-1] ";  
    std::cout << "spacing[0] ... spacing[n-1] ";  
    std::cout << "sigma[0] ... sigma[n-1] ";
    std::cout << "mean[0] .. mean[n-1] " << std::endl;
    exit( 1 );
    }
  typedef itk::Image<PixelType, ImageDimension> ImageType;   

  itk::TimeProbe timer;
  timer.Start();

  typedef itk::GaussianImageSource<ImageType> GaussianSourceType;
  GaussianSourceType::Pointer gaussian = GaussianSourceType::New();

  double scale = atof( argv[2] );
  ImageType::SizeType size;
  ImageType::PointType origin;
  ImageType::SpacingType spacing;
  GaussianSourceType::ArrayType mean; 
  GaussianSourceType::ArrayType sigma; 
  
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin[i] = atof( argv[3+i] );
    size[i] = atoi( argv[3+2*ImageDimension+i] );
    spacing[i] = atof( argv[3+1*ImageDimension+i] );
    sigma[i] = atof( argv[3+3*ImageDimension+i] );
    mean[i] = atoi( argv[3+4*ImageDimension+i] );
    }

  gaussian->SetSpacing( spacing );
  gaussian->SetOrigin( origin );
  gaussian->SetSize( size );
  gaussian->SetSigma( sigma );
  gaussian->SetScale( scale );
  gaussian->SetMean( mean );
  gaussian->SetNormalized( false );
  gaussian->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[1] );
  writer->SetInput( gaussian->GetOutput() );
  writer->Update();
  
  return 0;
}
