#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc != 3 + ImageDimension )
    {
    std::cout << "Usage: " << argv[0] << " image_filename output_filename spacing_dim_0 spacing_dim_1 (optional)spacing_dim_2" << std::endl;
    exit( 1 );
    }
  
  typedef double RealType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  
  typedef itk::IdentityTransform<RealType, ImageDimension> TransformType;
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  
  typedef itk::LinearInterpolateImageFunction<ImageType, RealType> LinearInterpolatorType;
  LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
  interpolator->SetInputImage( reader->GetOutput() );

  typedef itk::ResampleImageFilter<ImageType, ImageType, RealType> ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();
  ResamplerType::SpacingType spacing; 
  spacing[0] = atof( argv[3] );
  spacing[1] = atof( argv[4] );
  if ( ImageDimension == 3 )
    {
    spacing[2] = atof( argv[5] );
    }
    
  ResamplerType::SizeType size;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    RealType spacing_old = reader->GetOutput()->GetSpacing()[i];
    RealType size_old = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[i];
    size[i] = static_cast<int>( ( spacing_old * size_old ) / spacing[i] + 0.5 ); 
    }
      
  resampler->SetTransform( transform );
  resampler->SetInterpolator( interpolator );
  resampler->SetInput( reader->GetOutput() );
  resampler->SetOutputSpacing( spacing );
  resampler->SetOutputOrigin( reader->GetOutput()->GetOrigin() );
  resampler->SetSize( size );
  resampler->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( resampler->GetOutput() );
  writer->Update();

 return 0;
}
