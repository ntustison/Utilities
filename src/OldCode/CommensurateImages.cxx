#include <stdio.h>

#include "itkConstantPadImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"

#include "global.h"
#include <string>

int main( int argc, char *argv[] )
{
  /**
   * Take two images -> resample the images to the finest resolution
   *    -> pad the images so that they are the same size and have the same
   *    index.
   */

  if ( argc < 5 )
    {
    std::cout << "Usage: " << argv[0]
              << " inputImage1 inputImage2 outputImage1 outputImage2 [padValue1] [padValue2]" << std::endl;
    exit( 1 );
    }

  typedef float RealType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[1] );
  reader1->Update();

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  reader2->Update();

  // Resample the images to the finest resolution

  typedef itk::IdentityTransform<RealType, ImageDimension> TransformType;
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  typedef itk::LinearInterpolateImageFunction<ImageType, RealType> LinearInterpolatorType;

  LinearInterpolatorType::Pointer interpolator1 = LinearInterpolatorType::New();
  interpolator1->SetInputImage( reader1->GetOutput() );

  typedef itk::ResampleImageFilter<ImageType, ImageType, RealType> ResamplerType;

  ResamplerType::SpacingType spacing;
  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    spacing[d] = vnl_math_min( reader1->GetOutput()->GetSpacing()[d],
      reader2->GetOutput()->GetSpacing()[d] );
    }

  ResamplerType::SizeType size1;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    RealType spacing_old = reader1->GetOutput()->GetSpacing()[i];
    RealType size_old = reader1->GetOutput()->GetLargestPossibleRegion().GetSize()[i];
    size1[i] = static_cast<int>( ( spacing_old * size_old ) / spacing[i] + 0.5 );
    }
  ResamplerType::Pointer resampler1 = ResamplerType::New();
  resampler1->SetTransform( transform );
  resampler1->SetInterpolator( interpolator1 );
  resampler1->SetInput( reader1->GetOutput() );
  resampler1->SetOutputSpacing( spacing );
  resampler1->SetOutputOrigin( reader1->GetOutput()->GetOrigin() );
  resampler1->SetSize( size1 );
  resampler1->Update();

  LinearInterpolatorType::Pointer interpolator2 = LinearInterpolatorType::New();
  interpolator2->SetInputImage( reader2->GetOutput() );

  ResamplerType::SizeType size2;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    RealType spacing_old = reader2->GetOutput()->GetSpacing()[i];
    RealType size_old = reader2->GetOutput()->GetLargestPossibleRegion().GetSize()[i];
    size2[i] = static_cast<int>( ( spacing_old * size_old ) / spacing[i] + 0.5 );
    }
  ResamplerType::Pointer resampler2 = ResamplerType::New();
  resampler2->SetTransform( transform );
  resampler2->SetInterpolator( interpolator2 );
  resampler2->SetInput( reader2->GetOutput() );
  resampler2->SetOutputSpacing( spacing );
  resampler2->SetOutputOrigin( reader2->GetOutput()->GetOrigin() );
  resampler2->SetSize( size2 );
  resampler2->Update();

  // Pad the images so that they're the same size

  unsigned long lowerBound1[ImageDimension];
  unsigned long upperBound1[ImageDimension];
  unsigned long lowerBound2[ImageDimension];
  unsigned long upperBound2[ImageDimension];
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    int sizeDifference = resampler1->GetOutput()->GetLargestPossibleRegion().GetSize()[i]
      - resampler2->GetOutput()->GetLargestPossibleRegion().GetSize()[i];
    if ( sizeDifference > 0 )
      {
      lowerBound2[i] = static_cast<unsigned int>( vcl_floor( 0.5*static_cast<RealType>( sizeDifference ) ) );
      upperBound2[i] = static_cast<unsigned int>( vcl_ceil( 0.5*static_cast<RealType>( sizeDifference ) ) );
      lowerBound1[i] = 0;
      upperBound1[i] = 0;
      }
    else
      {
      sizeDifference *= -1;
      lowerBound2[i] = 0;
      upperBound2[i] = 0;
      lowerBound1[i] = static_cast<unsigned int>( vcl_floor( 0.5*static_cast<RealType>( sizeDifference ) ) );
      upperBound1[i] = static_cast<unsigned int>( vcl_ceil( 0.5*static_cast<RealType>( sizeDifference ) ) );
      }
    }

  typedef itk::ConstantPadImageFilter<ImageType, ImageType> PadderType;

  PadderType::Pointer padder1 = PadderType::New();
  padder1->SetInput( resampler1->GetOutput() );
  padder1->SetPadLowerBound( lowerBound1 );
  padder1->SetPadUpperBound( upperBound1 );
  padder1->SetConstant( 0 );
  if ( argc > 5 )
    {
    padder1->SetConstant( atof( argv[5] ) );
    }
  padder1->Update();

  PadderType::Pointer padder2 = PadderType::New();
  padder2->SetInput( resampler2->GetOutput() );
  padder2->SetPadLowerBound( lowerBound2 );
  padder2->SetPadUpperBound( upperBound2 );
  padder2->SetConstant( 0 );
  if ( argc > 6 )
    {
    padder2->SetConstant( atof( argv[6] ) );
    }
  padder2->Update();

  padder1->GetOutput()->SetOrigin( padder2->GetOutput()->GetOrigin() );
  padder1->GetOutput()->SetDirection( padder2->GetOutput()->GetDirection() );

  typedef itk::ImageFileWriter<ImageType> WriterType;

  WriterType::Pointer writer1 = WriterType::New();
  writer1->SetInput( padder1->GetOutput() );
  writer1->SetFileName( argv[3] );
  writer1->Update();

  WriterType::Pointer writer2 = WriterType::New();
  writer2->SetInput( padder2->GetOutput() );
  writer2->SetFileName( argv[4] );
  writer2->Update();

  return 0;
}