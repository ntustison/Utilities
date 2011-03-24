#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkEuler3DTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc != 12 )
    {
    std::cout << "Usage: " << argv[0] << " image_filename output_filename center_x center_y center_z rotation_x rotation_y rotation_z translation_x translation_y translation_z" << std::endl;
    exit( 1 );
    }
    
  const unsigned int ImageDimension = 3;  
  
  typedef double RealType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  
  typedef itk::Euler3DTransform<RealType> TransformType;
  TransformType::Pointer transform = TransformType::New();

  transform->SetRotation( static_cast<RealType>( atof( argv[6] ) ),
                          static_cast<RealType>( atof( argv[7] ) ),
                          static_cast<RealType>( atof( argv[8] ) ) );
  TransformType::OutputVectorType translation;
  translation[0] = static_cast<RealType>( atof( argv[9] ) );
  translation[1] = static_cast<RealType>( atof( argv[10] ) );
  translation[2] = static_cast<RealType>( atof( argv[11] ) );
  transform->SetTranslation( translation );
  TransformType::InputPointType center;
  center[0] = static_cast<RealType>( atof( argv[3] ) );
  center[1] = static_cast<RealType>( atof( argv[4] ) );
  center[2] = static_cast<RealType>( atof( argv[5] ) );
  transform->SetCenter( center );
  
  typedef itk::LinearInterpolateImageFunction<ImageType, RealType> LinearInterpolatorType;
  LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
  interpolator->SetInputImage( reader->GetOutput() );

  typedef itk::ResampleImageFilter<ImageType, ImageType, RealType> ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();
  resampler->SetTransform( transform );
  resampler->SetInterpolator( interpolator );
  resampler->SetInput( reader->GetOutput() );
  resampler->SetOutputSpacing( reader->GetOutput()->GetSpacing() );
  resampler->SetOutputOrigin( reader->GetOutput()->GetOrigin() );
  resampler->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  resampler->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( resampler->GetOutput() );
  writer->Update();

 return 0;
}
