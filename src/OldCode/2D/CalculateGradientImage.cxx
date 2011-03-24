#include <stdio.h>

#include "itkGradientImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage outputImage [magnitude]" << std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::GradientImageFilter<ImageType, RealType, RealType> GradientFilterType;
  GradientFilterType::Pointer gradient = GradientFilterType::New();
  gradient->SetInput( reader->GetOutput() );
  gradient->SetUseImageSpacing( true );
  gradient->Update();

  typedef GradientFilterType::OutputImageType      GradientImageType;

  ImageType::Pointer gradientMagnitudeImage = ImageType::New();
  gradientMagnitudeImage->SetOrigin( gradient->GetOutput()->GetOrigin() );
  gradientMagnitudeImage->SetSpacing( gradient->GetOutput()->GetSpacing() );
  gradientMagnitudeImage->SetRegions( gradient->GetOutput()->GetRequestedRegion() );
  gradientMagnitudeImage->Allocate();

  itk::ImageRegionIterator<GradientImageType> 
    ItG( gradient->GetOutput(), gradient->GetOutput()->GetRequestedRegion() );
  itk::ImageRegionIterator<RealImageType> 
    ItM( gradientMagnitudeImage, gradientMagnitudeImage->GetRequestedRegion() );
  for ( ItG.GoToBegin(), ItM.GoToBegin(); !ItG.IsAtEnd(); ++ItG, ++ItM )
    {
    ItM.Set( ( ItG.Get() ).GetNorm() );
    } 
  
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( gradientMagnitudeImage );
  writer->SetFileName( argv[2] );                                          
  writer->Update();
    
  return 0;
}
