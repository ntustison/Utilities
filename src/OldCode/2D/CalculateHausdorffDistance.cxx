#include <stdio.h>

#include "itkHausdorffDistanceImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage1 inputImage2 label(optional)"<< std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[1] );
  reader1->Update();

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  reader2->Update();

  long unsigned int differences = 0;

  typedef itk::HausdorffDistanceImageFilter<ImageType, ImageType> HausdorfferType;
  HausdorfferType::Pointer hausdorffer = HausdorfferType::New();

  if ( argc == 3 )
    {
    hausdorffer->SetInput1( reader1->GetOutput() );
    hausdorffer->SetInput2( reader2->GetOutput() );

    itk::ImageRegionIterator<ImageType> It1( reader1->GetOutput(), reader1->GetOutput()->GetLargestPossibleRegion() );  
    itk::ImageRegionIterator<ImageType> It2( reader2->GetOutput(), reader2->GetOutput()->GetLargestPossibleRegion() );  
    for ( It2.GoToBegin(), It1.GoToBegin(); !It2.IsAtEnd(); ++It2, ++It1 )
      {
      if ( It1.Get() != It2.Get() )
        {
        differences++;
        }
      }
    }
  else
    {

    ImageType::Pointer image1 = ImageType::New();
    image1->SetRegions( reader1->GetOutput()->GetLargestPossibleRegion() );
    image1->SetOrigin( reader1->GetOutput()->GetOrigin() );
    image1->SetSpacing( reader1->GetOutput()->GetSpacing() );
    image1->Allocate();
    image1->FillBuffer( 0 );

    itk::ImageRegionIterator<ImageType> ItR1( reader1->GetOutput(), reader1->GetOutput()->GetLargestPossibleRegion() );  
    itk::ImageRegionIterator<ImageType> It1( image1, image1->GetLargestPossibleRegion() );  
    for ( It1.GoToBegin(), ItR1.GoToBegin(); !It1.IsAtEnd(); ++It1, ++ItR1 )
      {
      if ( ItR1.Get() == static_cast<PixelType>( atoi( argv[3] ) ) )
        {
        It1.Set( 1 );
        }
      }

    ImageType::Pointer image2 = ImageType::New();
    image2->SetRegions( reader2->GetOutput()->GetLargestPossibleRegion() );
    image2->SetOrigin( reader2->GetOutput()->GetOrigin() );
    image2->SetSpacing( reader2->GetOutput()->GetSpacing() );
    image2->Allocate();
    image2->FillBuffer( 0 );

    itk::ImageRegionIterator<ImageType> ItR2( reader2->GetOutput(), reader2->GetOutput()->GetLargestPossibleRegion() );  
    itk::ImageRegionIterator<ImageType> It2( image2, image2->GetLargestPossibleRegion() );  
    for ( It2.GoToBegin(), ItR2.GoToBegin(); !It2.IsAtEnd(); ++It2, ++ItR2 )
      {
      if ( ItR2.Get() == static_cast<PixelType>( atoi( argv[3] ) ) )
        {
        It2.Set( 1 );
        }
      }
    for ( It2.GoToBegin(), It1.GoToBegin(); !It2.IsAtEnd(); ++It2, ++It1 )
      {
      if ( It1.Get() != It2.Get() )
        {
        differences++;
        }
      }

    hausdorffer->SetInput1( image1 );
    hausdorffer->SetInput2( image2 );
    }        

  hausdorffer->Update();
  std::cout << "Pixel-wise difference = " << differences << std::endl;
  std::cout << "Hausdorff distance = " << hausdorffer->GetHausdorffDistance() << std::endl;  

  return 0;
}
