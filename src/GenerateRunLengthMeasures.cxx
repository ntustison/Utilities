#include <stdio.h>

#include "itkBoundingBox.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPointSet.h"

#include "itkScalarImageToRunLengthFeaturesFilter.h"
#include "itkDenseFrequencyContainer2.h"

template <unsigned int ImageDimension>
int GenerateRunLengthMeasures( int argc, char *argv[] )
{

  typedef float PixelType;
  typedef float RealType;


  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[2] );
  imageReader->Update();

  typedef itk::Statistics::DenseFrequencyContainer2 HistogramFrequencyContainerType;

  typedef itk::Statistics::ScalarImageToRunLengthFeaturesFilter
    <RealImageType, HistogramFrequencyContainerType> RunLengthFilterType;
  typename RunLengthFilterType::Pointer runLengthFilter = RunLengthFilterType::New();
  runLengthFilter->SetInput( imageReader->GetOutput() );

  typename ImageType::Pointer mask = NULL;
  PixelType label = itk::NumericTraits<PixelType>::One;
  if ( argc > 4 )
    {
    typename ReaderType::Pointer labelImageReader = ReaderType::New();
    labelImageReader->SetFileName( argv[4] );
    labelImageReader->Update();
    mask = labelImageReader->GetOutput();
    runLengthFilter->SetInput( mask );

    if ( argc > 5 )
      {
      label = static_cast<PixelType>( atoi( argv[5] ) );
      }
    runLengthFilter->SetInsidePixelValue( label );
    }


  unsigned int numberOfBins = 256;
  if ( argc > 3 )
    {
    numberOfBins = static_cast<PixelType>( atoi( argv[3] ) );
    }
  runLengthFilter->SetNumberOfBinsPerAxis( numberOfBins );


  itk::ImageRegionIteratorWithIndex<ImageType> ItI( imageReader->GetOutput(),
    imageReader->GetOutput()->GetLargestPossibleRegion() );

  PixelType maxValue = itk::NumericTraits<PixelType>::NonpositiveMin();
  PixelType minValue = itk::NumericTraits<PixelType>::max();

  typedef itk::BoundingBox<unsigned long,
       ImageDimension, RealType> BoundingBoxType;
  typename BoundingBoxType::Pointer bbox = BoundingBoxType::New();
  typename BoundingBoxType::PointsContainerPointer points
       = BoundingBoxType::PointsContainer::New();
  itk::Point<RealType, ImageDimension> point;

  unsigned int idx = 0;

  for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
    {
    if ( !mask || ( mask->GetPixel( ItI.GetIndex() ) == label ) )
      {
      if ( ItI.Get() < minValue )
        {
        minValue = ItI.Get();
        }
      else if ( ItI.Get() > maxValue )
        {
        maxValue = ItI.Get();
        }
      imageReader->GetOutput()->TransformIndexToPhysicalPoint( ItI.GetIndex(), point );
      points->InsertElement( idx++, point );
      }
    }
  bbox->SetPoints( points );
  bbox->ComputeBoundingBox();
  typename BoundingBoxType::PointType pointMin = bbox->GetMinimum();
  typename BoundingBoxType::PointType pointMax = bbox->GetMaximum();

  runLengthFilter->SetPixelValueMinMax( minValue, maxValue );
  runLengthFilter->SetDistanceValueMinMax( 0, pointMin.EuclideanDistanceTo( pointMax ) );
  runLengthFilter->SetNumberOfBinsPerAxis( numberOfBins );
  runLengthFilter->FastCalculationsOff();

  runLengthFilter->Update();

  typename RunLengthFilterType::FeatureValueVectorPointer means =
    runLengthFilter->GetFeatureMeans();
  const typename RunLengthFilterType::FeatureNameVector* names =
    runLengthFilter->GetRequestedFeatures();

  typename RunLengthFilterType::FeatureValueVector::ConstIterator mIt =
    means->Begin();
  typename RunLengthFilterType::FeatureNameVector::ConstIterator nIt =
    names->Begin();

  std::cout << "ShortRunEmphasis,LongRunEmphasis,GreyLevelNonuniformity,RunLengthNonuniformity,LowGreyLevelRunEmphasis,HighGreyLevelRunEmphasis,ShortRunLowGreyLevelEmphasis,ShortRunHighGreyLevelEmphasis,LongRunLowGreyLevelEmphasis,LongRunHighGreyLevelEmphasis" << std::endl;
  while( mIt != means->End() )
    {
//    std::cout << nIt.Value() << ": " << mIt.Value() << std::endl;
    std::cout << mIt.Value() << " ";
    ++mIt;
    ++nIt;
    }
  std::cout << std::endl;



  return 0;
}


int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage "
     << "[numberOfBinsPerAxis=256] [maskImage] [maskLabel=1]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GenerateRunLengthMeasures<2>( argc, argv );
     break;
   case 3:
     GenerateRunLengthMeasures<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


