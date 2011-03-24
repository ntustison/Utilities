#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericTraits.h"
#include "itkAdaptiveHistogramEqualizationImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkFlipImageFilter.h"

#include "global.h"
#include <string>

int main( int argc, char *argv[] )
{
  if ( argc != 3 )
    {
    std::cout << "Usage: " << argv[0] << " input_image output_image" << std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  // If the requested output image is a .png or .jpg, optimize for visualization
  if ( strstr( argv[2], ".png" ) != NULL || strstr( argv[2], ".jpg" ) != NULL )
    {
    typedef itk::Image<unsigned char, ImageDimension> ShortImageType;

/*
    typedef itk::AdaptiveHistogramEqualizationImageFilter<ImageType> EqualizerType;
    EqualizerType::Pointer equalizer = EqualizerType::New();
    equalizer->SetInput( reader->GetOutput() );
    equalizer->SetAlpha( 0.0 );
    equalizer->SetBeta( 0.0 );
    EqualizerType::ImageSizeType size;
    size.Fill( 5 );
    equalizer->SetRadius( size );
*/

    typedef itk::RescaleIntensityImageFilter<ImageType, ShortImageType> FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput( reader->GetOutput() );
    filter->SetOutputMinimum( 0 );
    filter->SetOutputMaximum( 255 );

//    typedef itk::FlipImageFilter<ShortImageType> FlipperType;
//    FlipperType::Pointer flipper = FlipperType::New();
//    flipper->SetInput( filter->GetOutput() );
//    FlipperType::FlipAxesArrayType array;
//    array.Fill( false );
//    array[1] = true;
//    flipper->SetFlipAxes( array );
    
    typedef itk::ImageFileWriter<ShortImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( filter->GetOutput() );
    writer->SetFileName( argv[2] );                                          
    writer->Update();
    }
  else
    {  

/*
    ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();
    spacing[0] = 0.390625;
    spacing[2] = 1.0;
    spacing[1] = 0.390625;
    reader->GetOutput()->SetSpacing( spacing );    
    ImageType::PointType origin = reader->GetOutput()->GetOrigin();
    origin[0] = 13;
    origin[1] = -78;
    origin[2] = 0;
    reader->GetOutput()->SetOrigin( origin );    
*/ 

    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( reader->GetOutput() );
    writer->SetFileName( argv[2] );                                          
    writer->Update();
    }
    
  return 0;
}
