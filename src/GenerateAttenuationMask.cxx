#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

template <unsigned int ImageDimension>
int GenerateAttenuationMask( int argc, char *argv[] )
{

  typedef int PixelType;
  typedef float RealType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[2] );
  imageReader->Update();

  PixelType outputLabel = itk::NumericTraits<PixelType>::One;
  if ( argc > 6 )
    {
    outputLabel = static_cast<PixelType>( atoi( argv[6] ) );
    }
  
  typename ImageType::Pointer maskImage = ImageType::New();
  if ( argc > 7 )
    {
    typename ReaderType::Pointer labelImageReader = ReaderType::New();
    labelImageReader->SetFileName( argv[7] );
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
  if ( argc > 8 )
    {
    label = static_cast<PixelType>( atoi( argv[8] ) );
    } 

  itk::ImageRegionIterator<ImageType> It( imageReader->GetOutput(),
    imageReader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> ItM( maskImage,
    maskImage->GetLargestPossibleRegion() );
  
  PixelType minVoxel = static_cast<PixelType>( atoi( argv[4] ) );
  PixelType maxVoxel = static_cast<PixelType>( atoi( argv[5] ) );

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
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( maskImage );
  writer->SetFileName( argv[3] );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage outputImage min max [outputLabel] [labelImage] [label]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     GenerateAttenuationMask<2>( argc, argv );
     break;
   case 3:
     GenerateAttenuationMask<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

