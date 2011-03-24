#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << " image outputImage min max [outputLabel] [labelImage] [label]" << std::endl;
    exit( 1 );
    }

  typedef int PixelType;
  const unsigned int ImageDimension = 3;
  typedef float RealType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[1] );
  imageReader->Update();

  PixelType outputLabel = itk::NumericTraits<PixelType>::One;
  if ( argc > 5 )
    {
    outputLabel = static_cast<PixelType>( atoi( argv[5] ) );
    }
  
  ImageType::Pointer maskImage = ImageType::New();
  if ( argc > 6 )
    {
    ReaderType::Pointer labelImageReader = ReaderType::New();
    labelImageReader->SetFileName( argv[6] );
    labelImageReader->Update();
    maskImage = labelImageReader->GetOutput();
    } 
  else
    {
    maskImage->SetOrigin( imageReader->GetOutput()->GetOrigin() );
    maskImage->SetSpacing( imageReader->GetOutput()->GetSpacing() );
    maskImage->SetRegions( imageReader->GetOutput()->GetLargestPossibleRegion() );
    maskImage->Allocate();
    maskImage->FillBuffer( itk::NumericTraits<PixelType>::One );
    } 
  PixelType label = itk::NumericTraits<PixelType>::One;
  if ( argc > 7 )
    {
    label = static_cast<PixelType>( atoi( argv[7] ) );
    } 

  itk::ImageRegionIterator<ImageType> It( imageReader->GetOutput(),
    imageReader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> ItM( maskImage,
    maskImage->GetLargestPossibleRegion() );
  
  PixelType minVoxel = static_cast<PixelType>( atoi( argv[3] ) );
  PixelType maxVoxel = static_cast<PixelType>( atoi( argv[4] ) );

  for ( It.GoToBegin(), ItM.GoToBegin(); !It.IsAtEnd(); ++It, ++ItM )
    {
    if ( ItM.Get() == label && It.Get() >= minVoxel && It.Get() <= maxVoxel )
      {
      ItM.Set( outputLabel );
      }
    else
      { 
      ItM.Set(  0 );
      }
    }   

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( maskImage );
  writer->SetFileName( argv[2] );
  writer->Update();

  return 0;
}
