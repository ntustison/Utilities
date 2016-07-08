#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryThresholdImageFunction.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkExtractImageFilter.h"
#include "itkFloodFilledImageFunctionConditionalIterator.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include "Common.h"

template <unsigned int ImageDimension>
int ConnectedSegmentImage(int argc, char *argv[])
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<int, ImageDimension> LabelImageType;
  typedef itk::Image<PixelType, ImageDimension-1> SliceType;
  typedef itk::Image<int, ImageDimension-1> LabelSliceType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  typename ImageType::Pointer inputImage = reader->GetOutput();
  inputImage->Update();
  inputImage->DisconnectPipeline();

  typename LabelImageType::Pointer outputImage = LabelImageType::New();
  outputImage->CopyInformation( inputImage );
  outputImage->SetRegions( inputImage->GetRequestedRegion() );
  outputImage->Allocate();
  outputImage->FillBuffer( 0 );

  typedef int LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;
  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;

  typename LabelReaderType::Pointer seedReader = LabelReaderType::New();
  seedReader->SetFileName( argv[2] );

  typename LabelImageType::Pointer seedImage = seedReader->GetOutput();
  seedImage->Update();
  seedImage->DisconnectPipeline();

  float lowerMultiplier = 2.5;
  if ( argc > 4 )
    {
    lowerMultiplier = atof( argv[4] );
    }

  float upperMultiplier = 2.5;
  if ( argc > 5 )
    {
    upperMultiplier = atof( argv[5] );
    }

  unsigned int radius = 1;
  if ( argc > 6 )
    {
    radius = atoi( argv[6] );
    }

  typedef itk::LabelGeometryImageFilter<LabelImageType, ImageType> LabelFilterType;
  typename LabelFilterType::Pointer geometryFilter = LabelFilterType::New();
  geometryFilter->SetInput( seedImage );
  geometryFilter->CalculatePixelIndicesOff();
  geometryFilter->CalculateOrientedBoundingBoxOff();
  geometryFilter->CalculateOrientedLabelRegionsOff();
  geometryFilter->Update();

  typename LabelImageType::PointType centroid = geometryFilter->GetCentroid( 1 );
  typename ImageType::IndexType centroidIndex;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    centroidIndex[d] = centroid[d];
    }
  seedImage->TransformIndexToPhysicalPoint( centroidIndex, centroid );

  typename ImageType::IndexType inputIndex;
  inputImage->TransformPhysicalPointToIndex( centroid, inputIndex );

  // Temporarily set the seed pixel to perform binary dilation to get statistical info.

  outputImage->SetPixel( inputIndex, 1 );

  typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension> StructuringElementType;
  StructuringElementType  element;
  element.SetRadius( radius );
  element.CreateStructuringElement();

  typedef itk::BinaryDilateImageFilter<LabelImageType, LabelImageType, StructuringElementType>  DilateFilterType;
  typename DilateFilterType::Pointer  dilateFilter = DilateFilterType::New();
  dilateFilter->SetKernel( element );
  dilateFilter->SetInput( outputImage );
  dilateFilter->SetBackgroundValue( 0 );
  dilateFilter->SetForegroundValue( 1 );
  dilateFilter->Update();

  typedef itk::LabelStatisticsImageFilter<ImageType, LabelImageType> StatsFilterType;
  typename StatsFilterType::Pointer stats = StatsFilterType::New();
  stats->SetLabelInput( dilateFilter->GetOutput() );
  stats->SetInput( inputImage );
  stats->Update();

  float meanValue = stats->GetMean( 1 );
  float sigmaValue = stats->GetSigma( 1 );

  outputImage->FillBuffer( 0 );

  typedef itk::BinaryThresholdImageFunction<SliceType, float> FunctionType;
  typedef itk::FloodFilledImageFunctionConditionalIterator<SliceType, FunctionType> IteratorType;

  // Set up the image function used for connectivity
  typename FunctionType::Pointer function = FunctionType::New();

  typename SliceType::IndexType seedIndex;
  seedIndex[0] = inputIndex[0];
  seedIndex[1] = inputIndex[1];

  typedef std::vector<typename SliceType::IndexType> SeedsContainerType;
  SeedsContainerType seedsContainer;
  seedsContainer.push_back( seedIndex );

  // Do segmentations on all the slices

  typename ImageType::RegionType::SizeType size = inputImage->GetRequestedRegion().GetSize();
  typename ImageType::RegionType::IndexType index = inputImage->GetLargestPossibleRegion().GetIndex();
  unsigned int numberOfSlices = size[2];
  unsigned int beginIndex = index[2];

  size[2] = 0;

  for( unsigned int s = 0; s < numberOfSlices; s++ )
    {
    typename ImageType::RegionType region;

    index[2] = beginIndex + s;
    region.SetIndex( index );
    region.SetSize( size );

    typedef itk::ExtractImageFilter<ImageType, SliceType> ExtracterType;
    typename ExtracterType::Pointer extracter = ExtracterType::New();
    extracter->SetInput( inputImage );
    extracter->SetExtractionRegion( region );
    extracter->SetDirectionCollapseToIdentity();
    extracter->Update();

    function->SetInputImage( extracter->GetOutput() );
    function->ThresholdBetween( static_cast<PixelType>( meanValue - sigmaValue * lowerMultiplier ),
                                static_cast<PixelType>( meanValue + sigmaValue * upperMultiplier ) );

    IteratorType It = IteratorType( extracter->GetOutput(), function, seedsContainer );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      typename SliceType::IndexType sliceIndex = It.GetIndex();

      typename LabelImageType::IndexType outputIndex;
      outputIndex[0] = sliceIndex[0];
      outputIndex[1] = sliceIndex[1];
      outputIndex[2] = index[2];

      outputImage->SetPixel( outputIndex, 1 );
      }
    }

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( outputImage );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
  {
  if ( argc < 5 )
    {
    std::cerr << "Usage: "<< argv[0] << " inputImage seedImage outputImage "
      "<lowerMultiplier=2.5> <upperMultiplier=2.5> <radius=1>" << std::endl;
    return EXIT_FAILURE;
    }

  ConnectedSegmentImage<3>( argc, argv );
  }
