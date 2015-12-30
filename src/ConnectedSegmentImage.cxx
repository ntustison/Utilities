#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkConfidenceConnectedImageFilter.h"

template <unsigned int ImageDimension>
int ConnectedSegmentImage(int argc, char *argv[])
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<int, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );

  typename ImageType::Pointer inputImage = reader->GetOutput();
  inputImage->Update();
  inputImage->DisconnectPipeline();

  typedef int LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;
  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;

  typename LabelReaderType::Pointer seedReader = LabelReaderType::New();
  seedReader->SetFileName( argv[3] );

  typename LabelImageType::Pointer seedImage = seedReader->GetOutput();
  seedImage->Update();
  seedImage->DisconnectPipeline();

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

  typedef itk::ConfidenceConnectedImageFilter<ImageType, LabelImageType> SegmenterType;
  typename SegmenterType::Pointer segmenter = SegmenterType::New();
  segmenter->SetInput( inputImage );
  segmenter->SetReplaceValue( 1 );
  segmenter->SetInitialNeighborhoodRadius( 1 );

  segmenter->SetSeed( inputIndex );

  float multiplier = 2.5;
  if ( argc > 5 )
    {
    multiplier = atof( argv[5] );
    }
  segmenter->SetMultiplier( multiplier );

  unsigned int numberOfIterations = 1;
  if ( argc > 6 )
    {
    numberOfIterations = atoi( argv[6] );
    }
  segmenter->SetNumberOfIterations( numberOfIterations );
  segmenter->Update();

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( segmenter->GetOutput() );
  writer->Update();

 return 0;
}

int main( int argc, char *argv[] )
  {
  if ( argc < 5 )
    {
    std::cerr << "Usage: "<< argv[0] << " imageDimension inputImage seedImage outputImage "
      "<multiplier=2.5> <numberOfIterations=1>" << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      ConnectedSegmentImage<2>( argc, argv );
      break;
    case 3:
      ConnectedSegmentImage<3>( argc, argv );
      break;
    case 4:
      ConnectedSegmentImage<4>( argc, argv );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
  }
