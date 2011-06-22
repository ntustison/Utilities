#include <stdio.h>

#include "itkHausdorffDistanceImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"


template <unsigned int ImageDimension>
int CalculateHausdorffDistance( int argc, char *argv[] )
{

  typedef int PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[2] );
  reader1->Update();

  typename ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[3] );
  reader2->Update();

  long unsigned int differences = 0;

  typedef itk::HausdorffDistanceImageFilter<ImageType, ImageType> HausdorfferType;
  typename HausdorfferType::Pointer hausdorffer = HausdorfferType::New();

  if ( argc <= 4 )
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

    typename ImageType::Pointer image1 = ImageType::New();
    image1->SetRegions( reader1->GetOutput()->GetLargestPossibleRegion() );
    image1->SetOrigin( reader1->GetOutput()->GetOrigin() );
    image1->SetSpacing( reader1->GetOutput()->GetSpacing() );
    image1->Allocate();
    image1->FillBuffer( 0 );

    itk::ImageRegionIterator<ImageType> ItR1( reader1->GetOutput(), reader1->GetOutput()->GetLargestPossibleRegion() );  
    itk::ImageRegionIterator<ImageType> It1( image1, image1->GetLargestPossibleRegion() );  
    for ( It1.GoToBegin(), ItR1.GoToBegin(); !It1.IsAtEnd(); ++It1, ++ItR1 )
      {
      if ( ItR1.Get() == static_cast<PixelType>( atoi( argv[4] ) ) )
        {
        It1.Set( 1 );
        }
      }

    typename ImageType::Pointer image2 = ImageType::New();
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

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension inputImage1 inputImage2 [label]"<< std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     CalculateHausdorffDistance<2>( argc, argv );
     break;
   case 3:
     CalculateHausdorffDistance<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

