#include <stdio.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkCastImageFilter.h"
#include "itkScalarToFractalImageFilter.h"

int main( int argc, char *argv[] )
{
  if ( argc < 2 )
    {
    std::cout << "Usage: " << argv[0] << "inputImage outputImage [labelImage] [label]" << std::endl;
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
  
  ImageType::Pointer mask = ImageType::New();
  mask->SetRegions( imageReader->GetOutput()->GetLargestPossibleRegion() );
  mask->SetSpacing( imageReader->GetOutput()->GetSpacing() );
  mask->SetOrigin( imageReader->GetOutput()->GetOrigin() );
  mask->Allocate();

  if ( argc < 4 )
    {
    mask->FillBuffer( itk::NumericTraits<PixelType>::One );
    }
  else
    { 
    ReaderType::Pointer labelImageReader = ReaderType::New();
    labelImageReader->SetFileName( argv[3] );
    labelImageReader->Update();
  
    PixelType label = static_cast<PixelType>( atoi( argv[4] ) );

    mask->FillBuffer( itk::NumericTraits<PixelType>::Zero );

    itk::ImageRegionIterator<ImageType> ItI( imageReader->GetOutput(),
      imageReader->GetOutput()->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<ImageType> ItL( labelImageReader->GetOutput(),
      labelImageReader->GetOutput()->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<ImageType> ItM( mask,
      mask->GetLargestPossibleRegion() );
  
    for ( ItM.GoToBegin(), ItL.GoToBegin(); !ItM.IsAtEnd(); ++ItM, ++ItL )
      {
      if ( ItL.Get() == label )
        {
        ItM.Set( itk::NumericTraits<PixelType>::One );
        }
      }   
    }       

  typedef itk::ScalarToFractalImageFilter<ImageType, RealImageType> FractalFilterType;

  typedef itk::CastImageFilter<ImageType, FractalFilterType::MaskImageType> CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput( mask );
  caster->Update();

  FractalFilterType::Pointer fractal = FractalFilterType::New();
  fractal->SetInput( imageReader->GetOutput() );
  fractal->SetMaskImage( caster->GetOutput() );
  fractal->Update();
  
  typedef itk::ImageFileWriter<RealImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( fractal->GetOutput() );
  writer->SetFileName( argv[2] );
  writer->Update();



  return 0;
}
