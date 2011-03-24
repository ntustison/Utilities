#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkExtractImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkAdaptiveHistogramEqualizationImageFilter.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " image_filename output_filename x_1 y_1 z_1 x_2 y_2 z_2 size" << std::endl;
    exit( 1 );
    }
  
  typedef double RealType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  
  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescalerType;
  RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetInput( reader->GetOutput() );
  rescaler->SetOutputMinimum( 0.0 );
  rescaler->SetOutputMaximum( 255.0 );
  rescaler->Update();

/*
  typedef itk::AdaptiveHistogramEqualizationImageFilter<ImageType> EqualizerType;
  EqualizerType::Pointer equalizer = EqualizerType::New();
  equalizer->SetInput( rescaler->GetOutput() );
  equalizer->SetAlpha( 1 );
  equalizer->SetBeta( 1 );
  equalizer->Update();  
*/
  typedef itk::ExtractImageFilter<ImageType, ImageType> CropperType;
  CropperType::Pointer cropper = CropperType::New();

  ImageType::RegionType region;
  ImageType::RegionType::SizeType size;
  ImageType::RegionType::IndexType index;
  if ( argc > 3 )
    {
    size[0] = atoi( argv[6] ) - atoi( argv[3] ) + 1;
    size[1] = atoi( argv[7] ) - atoi( argv[4] ) + 1;
    size[2] = atoi( argv[8] ) - atoi( argv[5] ) + 1;
    index[0] = atoi( argv[3] );
    index[1] = atoi( argv[4] );
    index[2] = atoi( argv[5] );
    }
  else
   {
   for ( unsigned int d = 0; d < ImageDimension; d++ )
     {
     size[d] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[d]; 
     index[d] = reader->GetOutput()->GetLargestPossibleRegion().GetIndex()[d];
     }
   }
  region.SetSize( size );
  region.SetIndex( index );
  
  cropper->SetInput( rescaler->GetOutput() );
  cropper->SetExtractionRegion( region );
  cropper->Update();  
  
  index.Fill( 0 );
  region.SetIndex( index );
  cropper->GetOutput()->SetRegions( region ); 
  
  typedef itk::IdentityTransform<RealType, ImageDimension> TransformType;
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  
  typedef itk::LinearInterpolateImageFunction<ImageType, RealType> InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  typedef itk::ResampleImageFilter<ImageType, ImageType, RealType> ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();
  ResamplerType::SpacingType spacing;
  spacing[0] = vnl_math_min( cropper->GetOutput()->GetSpacing()[0], 
    vnl_math_min( cropper->GetOutput()->GetSpacing()[1], cropper->GetOutput()->GetSpacing()[2] ) );
  spacing[1] = spacing[2] = spacing[0];

  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    RealType spacing_old = cropper->GetOutput()->GetSpacing()[i];
    RealType size_old = cropper->GetOutput()->GetLargestPossibleRegion().GetSize()[i];
    size[i] = static_cast<int>( ( spacing_old * size_old ) / spacing[i] + 0.5 ); 
    }

  resampler->SetTransform( transform );
  resampler->SetInterpolator( interpolator );
  resampler->SetInput( cropper->GetOutput() );
  resampler->SetOutputSpacing( spacing );
  resampler->SetOutputOrigin( cropper->GetOutput()->GetOrigin() );
  resampler->SetSize( size );
  resampler->SetDefaultPixelValue( 1.0 );
  resampler->Update();

  int max_size = vnl_math_max( size[0], vnl_math_max( size[1], size[2] ) );

  unsigned long upperfactors[3];
  unsigned long lowerfactors[3];

  for ( unsigned int i = 0; i < 3; i++ )
    {
    upperfactors[i] = static_cast<int>( 0.5*( max_size-size[i] ) + 0.5 );
    lowerfactors[i] = max_size - size[i] - upperfactors[i];
    }

  typedef itk::ConstantPadImageFilter<ImageType, ImageType> PadderType;
  PadderType::Pointer padder = PadderType::New();
  padder->SetInput( resampler->GetOutput() );
  padder->SetConstant( 0 );
  padder->SetPadUpperBound( upperfactors );
  padder->SetPadLowerBound( lowerfactors );
  padder->Update();

  index.Fill( 0 );
  region.SetSize( padder->GetOutput()->GetLargestPossibleRegion().GetSize() );
  region.SetIndex( index );
  padder->GetOutput()->SetRegions( region );
  
  size.Fill( atoi(argv[9]) );
  spacing.Fill( padder->GetOutput()->GetSpacing()[0] *
    padder->GetOutput()->GetLargestPossibleRegion().GetSize()[0]/size[0] );
    
  ResamplerType::Pointer resampler2 = ResamplerType::New();
  resampler2->SetInput( padder->GetOutput() );
  resampler2->SetSize( size );
  resampler2->SetOutputOrigin( padder->GetOutput()->GetOrigin() );
  resampler2->SetOutputSpacing( spacing );
  resampler2->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( resampler2->GetOutput() );
  writer->Update();

 return 0;
}
