#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkAddImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"


template <unsigned int ImageDimension>
int FindOverlappingLabels( int argc, char * argv[] )
{
  typedef int LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;
  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;

  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::ImageFileReader<RealImageType> ReaderType;

  typename LabelReaderType::Pointer sourceLabelReader = LabelReaderType::New();
  sourceLabelReader->SetFileName( argv[2] );
  sourceLabelReader->Update();

  typename LabelReaderType::Pointer targetLabelReader = LabelReaderType::New();
  targetLabelReader->SetFileName( argv[3] );
  targetLabelReader->Update();

  float thresholdPercentage = 0.90;
  if( argc > 5 )
    {
    thresholdPercentage = atof( argv[5] );
    }

  typedef itk::LabelGeometryImageFilter<LabelImageType, RealImageType> FilterType;
  typename FilterType::Pointer sourceFilter = FilterType::New();
  sourceFilter->SetInput( sourceLabelReader->GetOutput() );
  sourceFilter->CalculatePixelIndicesOff();
  sourceFilter->CalculateOrientedBoundingBoxOff();
  sourceFilter->CalculateOrientedLabelRegionsOff();
  sourceFilter->Update();

  typename FilterType::Pointer targetFilter = FilterType::New();
  targetFilter->SetInput( targetLabelReader->GetOutput() );
  targetFilter->CalculatePixelIndicesOff();
  targetFilter->CalculateOrientedBoundingBoxOff();
  targetFilter->CalculateOrientedLabelRegionsOff();
  targetFilter->Update();

  typename RealImageType::Pointer outputImage = RealImageType::New();
  outputImage->SetRegions( targetLabelReader->GetOutput()->GetRequestedRegion() );
  outputImage->Allocate();
  outputImage->FillBuffer( 0 );

  typename FilterType::LabelsType sourceLabels = sourceFilter->GetLabels();
  typename FilterType::LabelsType::iterator sourceLabelsIt;
  for( sourceLabelsIt = sourceLabels.begin(); sourceLabelsIt != sourceLabels.end(); sourceLabelsIt++ )
    {
    if( *sourceLabelsIt == 0 )
      {
      continue;
      }

    typedef itk::BinaryThresholdImageFilter<LabelImageType, RealImageType> ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput( sourceLabelReader->GetOutput() );
    thresholder->SetLowerThreshold( *sourceLabelsIt );
    thresholder->SetUpperThreshold( *sourceLabelsIt );
    thresholder->SetInsideValue( 1 );
    thresholder->SetOutsideValue( 1 );

    typedef itk::MultiplyImageFilter<LabelImageType, RealImageType> MultiplierType;
    typename MultiplierType::Pointer multiplier = MultiplierType::New();
    multiplier->SetInput1( targetLabelReader->GetOutput() );
    multiplier->SetInput2( thresholder->GetOutput() );

    typename FilterType::Pointer localFilter = FilterType::New();
    localFilter->SetInput( multiplier->GetOutput() );
    localFilter->CalculatePixelIndicesOff();
    localFilter->CalculateOrientedBoundingBoxOff();
    localFilter->CalculateOrientedLabelRegionsOff();
    localFilter->Update();

    typename FilterType::LabelsType localLabels = localFilter->GetLabels();
    typename FilterType::LabelsType::iterator localLabelsIt;
    for( localLabelsIt = localLabels.begin(); localLabelsIt != localLabels.end(); localLabelsIt++ )
      {
      if( *localLabelsIt == 0 )
        {
        continue;
        }

      float volumeRatio = static_cast<float>( localFilter->GetVolume( *localLabelsIt ) ) /
        static_cast<float>( localFilter->GetVolume( *localLabelsIt ) );
      if( volumeRatio >= thresholdPercentage )
        {
        typename ThresholderType::Pointer thresholder2 = ThresholderType::New();
        thresholder2->SetInput( targetLabelReader->GetOutput() );
        thresholder2->SetLowerThreshold( *localLabelsIt );
        thresholder2->SetUpperThreshold( *localLabelsIt );
        thresholder2->SetInsideValue( *sourceLabelsIt );
        thresholder2->SetOutsideValue( 0 );

        typedef itk::AddImageFilter<RealImageType, RealImageType> AdderType;
        typename AdderType::Pointer adder = AdderType::New();
        adder->SetInput1( outputImage );
        adder->SetInput2( thresholder2->GetOutput() );
        adder->Update();
        adder->GetOutput()->DisconnectPipeline();

        outputImage = adder->GetOutput();
        }
      }
    }

  typedef itk::ImageFileWriter<RealImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( outputImage );
  writer->SetFileName( argv[4] );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension sourceLabelImage targetLabelImage outputImage [thresholdPercentOverlap=0.9]"
      << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     FindOverlappingLabels<2>( argc, argv );
     break;
   case 3:
     FindOverlappingLabels<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

