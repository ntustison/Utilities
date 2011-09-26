#include <stdio.h>

#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

template <unsigned int ImageDimension>
int FindClosestPointsBetweenTwoObjects( int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef int LabelType;
  typedef itk::Image<LabelType, ImageDimension> MaskImageType;

  LabelType label1 = 1;
  if( argc > 5 )
    {
    label1 = static_cast<LabelType>( atof( argv[5] ) );
    }

  LabelType label2 = 1;
  if( argc > 6 )
    {
    label2 = static_cast<LabelType>( atof( argv[6] ) );
    }

  typedef itk::ImageFileReader<MaskImageType> ReaderType;

  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[2] );
  reader1->Update();

  typename ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[3] );
  reader2->Update();

  typedef itk::BinaryThresholdImageFilter<MaskImageType, MaskImageType> ThresholderType;

  typename ThresholderType::Pointer thresholder1 = ThresholderType::New();
  thresholder1->SetInput( reader1->GetOutput() );
  thresholder1->SetLowerThreshold( label1 );
  thresholder1->SetUpperThreshold( label1 );
  thresholder1->SetInsideValue( 1 );
  thresholder1->SetOutsideValue( 0 );

  typename ThresholderType::Pointer thresholder2 = ThresholderType::New();
  thresholder2->SetInput( reader2->GetOutput() );
  thresholder2->SetLowerThreshold( label2 );
  thresholder2->SetUpperThreshold( label2 );
  thresholder2->SetInsideValue( 1 );
  thresholder2->SetOutsideValue( 0 );

  typedef itk::SignedMaurerDistanceMapImageFilter<MaskImageType, ImageType> FilterType;

  typename FilterType::Pointer filter1 = FilterType::New();
  filter1->SetInput( thresholder1->GetOutput() );
  filter1->SetSquaredDistance( false );
  filter1->SetUseImageSpacing( true );
  filter1->SetInsideIsPositive( false );
  filter1->Update();

  typename FilterType::Pointer filter2 = FilterType::New();
  filter2->SetInput( thresholder2->GetOutput() );
  filter2->SetSquaredDistance( false );
  filter2->SetUseImageSpacing( true );
  filter2->SetInsideIsPositive( false );
  filter2->Update();

  typename ImageType::PixelType minDistance = itk::NumericTraits<typename ImageType::PixelType>::max();

  itk::ImageRegionIteratorWithIndex<ImageType> ItF(
    filter1->GetOutput(), filter1->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<MaskImageType> ItM(
    thresholder2->GetOutput(), thresholder2->GetOutput()->GetLargestPossibleRegion() );
  for( ItF.GoToBegin(), ItM.GoToBegin(); !ItF.IsAtEnd(); ++ItF, ++ItM )
    {
    if( ItF.Get() < minDistance && ItM.Get() == 1 )
      {
      minDistance = ItF.Get();
      }
    }

  typename MaskImageType::Pointer output = MaskImageType::New();
  output->CopyInformation( reader2->GetOutput() );
  output->SetRegions( reader2->GetOutput()->GetLargestPossibleRegion() );
  output->Allocate();
  output->FillBuffer( 0 );

  itk::ImageRegionIteratorWithIndex<MaskImageType> ItO(
    output, output->GetLargestPossibleRegion() );
  for( ItF.GoToBegin(), ItM.GoToBegin(), ItO.GoToBegin(); !ItF.IsAtEnd(); ++ItF, ++ItM, ++ItO )
    {
    if( ItF.Get() <= minDistance && ItM.Get() == 1 )
      {
      ItO.Set( 1 );
      }
    }
  std::cout << "Min distance = " << minDistance << std::endl;

  typedef itk::ImageFileWriter<MaskImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( output );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension inputMaskImage1 "
      << "inputMaskImage2 outputImage [label1=1] [label2=1]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     FindClosestPointsBetweenTwoObjects<2>( argc, argv );
     break;
   case 3:
     FindClosestPointsBetweenTwoObjects<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
