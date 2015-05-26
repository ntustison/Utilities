#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkExtractImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include <string>
#include <vector>

typedef float RealType;

std::vector<RealType> FitRegressionLine(
  std::vector<RealType> X, std::vector<RealType> Y )
{
  if( X.size() != Y.size() )
    {
    std::cerr << "Vectors are not the same size" << std::endl;
    exit( 1 );
    }

  RealType N = X.size();

  std::vector<RealType>::const_iterator itX;
  std::vector<RealType>::const_iterator itY;

  RealType sumX = 0.0;
  RealType sumY = 0.0;
  RealType sumX2 = 0.0;
  RealType sumXY = 0.0;

  for( itX = X.begin(), itY = Y.begin(); itX != X.end(); ++itX, ++itY )
    {
    sumX  += (*itX);
    sumY  += (*itY);
    sumXY += (*itX) * (*itY);
    sumX2 += (*itX) * (*itX);
    }

  std::vector<RealType> line( 2 );
  line[0] = ( N * sumXY - sumX*sumY ) / ( N *sumX2 - sumX*sumX );
  if( sumX2 == 0 )
    {
    line[1] = sumY / N;
    }
  else
    {
    line[1] = ( sumY - line[0] * sumX ) / N;
    }

  return line;
}


template <unsigned int ImageDimension>
int DirectionalBiasCorrection( int argc, char *argv[] )
{
  typedef float                                      PixelType;
  typedef int                                        MaskPixelType;

  typedef itk::Image<PixelType, ImageDimension>      ImageType;
  typedef itk::Image<MaskPixelType, ImageDimension>  MaskImageType;

  typedef itk::Image<PixelType, ImageDimension-1>      SliceType;
  typedef itk::Image<MaskPixelType, ImageDimension-1>  MaskSliceType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );

  typename ImageType::Pointer inputImage = reader->GetOutput();
  inputImage->Update();
  inputImage->DisconnectPipeline();

  typename MaskImageType::Pointer maskImage = MaskImageType::New();
  if( argc > 3 )
    {
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    typename MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName( argv[3] );
    maskImage = maskReader->GetOutput();

    maskImage->Update();
    maskImage->DisconnectPipeline();
    }
  else
    {
    maskImage->CopyInformation( inputImage );
    maskImage->SetRegions( inputImage->GetRequestedRegion() );
    maskImage->Allocate();
    maskImage->FillBuffer( 1 );
    }

  unsigned int direction = 0;
  if( argc > 4 )
    {
    atoi( argv[4] );
    }
  else
    {
    typename ImageType::SpacingType spacing = inputImage->GetSpacing();
    direction = 0;
    RealType maxSpacing = spacing[0];
    for( unsigned int d = 1; d < ImageDimension; d++ )
      {
      if( spacing[d] > maxSpacing )
        {
        direction = d;
        maxSpacing = spacing[d];
        }
      }
    }

  unsigned int model = 0;
  if( argc > 5 )
    {
    model = atoi( argv[5] );
    }

  unsigned int numberOfSlices = inputImage->GetRequestedRegion().GetSize()[direction];

  std::vector<RealType> X;
  std::vector<RealType> Y;

  typename ImageType::RegionType region;
  typename ImageType::RegionType::SizeType size = inputImage->GetRequestedRegion().GetSize();
  size[direction] = 0;

  typename ImageType::IndexType index;
  index.Fill( 0 );

  RealType totalAverageIntensityValue = 0.0;
  unsigned int nonZeroSliceCount = 0;

  for( unsigned int n = 0; n < numberOfSlices; n++ )
    {
    index[direction] = n;
    region.SetIndex( index );
    region.SetSize( size );

    typedef itk::ExtractImageFilter<ImageType, SliceType> ExtracterType;
    typename ExtracterType::Pointer extracter = ExtracterType::New();
    extracter->SetInput( inputImage );
    extracter->SetExtractionRegion( region );
    extracter->SetDirectionCollapseToIdentity();
    extracter->Update();

    typedef itk::ExtractImageFilter<MaskImageType, MaskSliceType> MaskExtracterType;
    typename MaskExtracterType::Pointer maskExtracter = MaskExtracterType::New();
    maskExtracter->SetInput( maskImage );
    maskExtracter->SetExtractionRegion( region );
    maskExtracter->SetDirectionCollapseToIdentity();
    maskExtracter->Update();

    typedef itk::LabelStatisticsImageFilter<SliceType, MaskSliceType> StatsFilterType;
    typename StatsFilterType::Pointer stats = StatsFilterType::New();
    stats->SetInput( extracter->GetOutput() );
    stats->SetLabelInput( maskExtracter->GetOutput() );
    stats->Update();

    if( stats->GetCount( 1 ) > 0 )
      {
      X.push_back( static_cast<RealType>( n ) );
      if( model == 0 )
        {
        Y.push_back( std::log( stats->GetMean( 1 ) ) );
        }
      totalAverageIntensityValue = totalAverageIntensityValue +
        ( Y[nonZeroSliceCount] - totalAverageIntensityValue ) /
        static_cast<RealType>( nonZeroSliceCount + 1 );
      nonZeroSliceCount++;
      }
    }

  /////////////////////////////
  //
  // Perform fitting
  //
  /////////////////////////////

  std::vector<RealType> line = FitRegressionLine( X, Y );

  typename ImageType::Pointer outputImage = ImageType::New();
  outputImage->CopyInformation( inputImage );
  outputImage->SetRegions( inputImage->GetRequestedRegion() );
  outputImage->Allocate();
  outputImage->FillBuffer( 0 );

  itk::ImageRegionIteratorWithIndex<ImageType> ItO( outputImage,
    outputImage->GetRequestedRegion() );
  itk::ImageRegionConstIterator<ImageType> ItI( inputImage,
    inputImage->GetRequestedRegion() );
  for( ItI.GoToBegin(), ItO.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItO )
    {
    typename ImageType::IndexType index = ItO.GetIndex();

    RealType predictedValue = line[0] + line[1] * static_cast<RealType>( index[direction] );
    if( model == 0 )
      {
      predictedValue = std::exp( line[0] ) * std::exp( line[1] * static_cast<RealType>( index[direction] ) );
      }
    RealType residual = ItI.Get() - predictedValue;
    ItO.Set( residual );
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( outputImage );
  writer->SetFileName( argv[2] );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0]
      << "inputImage outputImage <maskImage> <direction=largest spacing> <model=0>" << std::endl;

    std::cerr << "  model types: " << std::endl;
    std::cerr << "      0: exponential " << std::endl;
    std::cerr << "      1: linear " << std::endl;

    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
   {
   case 3:
     DirectionalBiasCorrection<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
   }
}
