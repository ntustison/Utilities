#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkConfidenceConnectedImageFilter.h"

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

  float multiplier = 2.5;
  if ( argc > 4 )
    {
    multiplier = atof( argv[4] );
    }

  unsigned int numberOfIterations = 1;
  if ( argc > 5 )
    {
    numberOfIterations = atoi( argv[5] );
    }

  typedef itk::LabelGeometryImageFilter<LabelImageType, ImageType> LabelFilterType;
  typename LabelFilterType::Pointer filter = LabelFilterType::New();
  filter->SetInput( seedImage );
  filter->CalculatePixelIndicesOff();
  filter->CalculateOrientedBoundingBoxOff();
  filter->CalculateOrientedLabelRegionsOff();
  filter->Update();

  typename LabelImageType::PointType centroid = filter->GetCentroid( 1 );
  typename ImageType::IndexType centroidIndex;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    centroidIndex[d] = centroid[d];
    }
  seedImage->TransformIndexToPhysicalPoint( centroidIndex, centroid );

  typename ImageType::IndexType inputIndex;
  inputImage->TransformPhysicalPointToIndex( centroid, inputIndex );

  typename SliceType::IndexType seedIndex;
  seedIndex[0] = inputIndex[0];
  seedIndex[1] = inputIndex[1];

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

    typedef itk::ConfidenceConnectedImageFilter<SliceType, LabelSliceType> SegmenterType;
    typename SegmenterType::Pointer segmenter = SegmenterType::New();
    segmenter->SetInput( extracter->GetOutput() );
    segmenter->SetReplaceValue( 1 );
    segmenter->SetInitialNeighborhoodRadius( 1 );
    segmenter->SetSeed( seedIndex );
    segmenter->SetNumberOfIterations( numberOfIterations );
    segmenter->SetMultiplier( multiplier );
    segmenter->Update();

    itk::ImageRegionIteratorWithIndex<LabelSliceType> It( segmenter->GetOutput(),
      segmenter->GetOutput()->GetRequestedRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      typename SliceType::IndexType sliceIndex = It.GetIndex();

      typename LabelImageType::IndexType outputIndex;
      outputIndex[0] = sliceIndex[0];
      outputIndex[1] = sliceIndex[1];
      outputIndex[2] = index[2];

      outputImage->SetPixel( outputIndex, It.Get() );
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
      "<multiplier=2.5> <numberOfIterations=1>" << std::endl;
    return EXIT_FAILURE;
    }

  ConnectedSegmentImage<3>( argc, argv );
  }
