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
  reader->SetFileName( argv[1] );

  typename ImageType::Pointer inputImage = reader->GetOutput();
  inputImage->Update();
  inputImage->DisconnectPipeline();

  typename MaskImageType::Pointer maskImageRegression = MaskImageType::New();
  if( argc > 3 )
    {
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    typename MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName( argv[3] );
    maskImageRegression = maskReader->GetOutput();

    maskImageRegression->Update();
    maskImageRegression->DisconnectPipeline();
    }
  else
    {
    maskImageRegression->CopyInformation( inputImage );
    maskImageRegression->SetRegions( inputImage->GetRequestedRegion() );
    maskImageRegression->Allocate();
    maskImageRegression->FillBuffer( 1 );
    }

  typename MaskImageType::Pointer maskImageUpdate = MaskImageType::New();
  if( argc > 4 )
    {
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    typename MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName( argv[4] );
    maskImageUpdate = maskReader->GetOutput();

    maskImageUpdate->Update();
    maskImageUpdate->DisconnectPipeline();
    }
  else
    {
    maskImageUpdate->CopyInformation( inputImage );
    maskImageUpdate->SetRegions( inputImage->GetRequestedRegion() );
    maskImageUpdate->Allocate();
    maskImageUpdate->FillBuffer( 1 );
    }


  unsigned int physicalCoordinateComponent = 1;
  if( argc > 5 )
    {
    if( std::strcmp( argv[5], "x" ) == 0 )
      {
      physicalCoordinateComponent = 0;
      }
    else if( std::strcmp( argv[5], "y" ) == 0 )
      {
      physicalCoordinateComponent = 1;
      }
    else if( std::strcmp( argv[5], "z" ) == 0 )
      {
      physicalCoordinateComponent = 2;
      }
    else
      {
      std::cerr << "Unrecognized direction option.  Must choose 'x', 'y', or 'z'." << std::endl;
      return EXIT_FAILURE;
      }
    }

  unsigned int direction = 0;
  float maxComponentValue = 0.0;

  typename ImageType::IndexType index;
  index.Fill( 0 );
  typename ImageType::PointType pointOrigin;
  inputImage->TransformIndexToPhysicalPoint( index, pointOrigin );

  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    typename ImageType::PointType point;
    index.Fill( 0 );
    index[d] = 1;
    inputImage->TransformIndexToPhysicalPoint( index, point );

    typename ImageType::PointType::VectorType directionalVector = point - pointOrigin;

    if( vnl_math_abs( directionalVector[physicalCoordinateComponent] ) > maxComponentValue )
      {
      direction = d;
      }
    }
  index.Fill( 0 );

  unsigned int model = 1;
  if( argc > 6 )
    {
    model = atoi( argv[6] );
    }

  unsigned int numberOfSlices = inputImage->GetRequestedRegion().GetSize()[direction];

  std::vector<RealType> X;
  std::vector<RealType> Y;

  typename ImageType::RegionType region;
  typename ImageType::RegionType::SizeType size = inputImage->GetRequestedRegion().GetSize();
  size[direction] = 0;


  RealType totalAverageIntensityValue = 0.0;
  RealType maxIntensityValue = itk::NumericTraits<RealType>::NonpositiveMin();
  RealType minIntensityValue = itk::NumericTraits<RealType>::max();
  unsigned int nonZeroSliceCount = 0;

  std::cout << "direction = " << direction << std::endl;
  std::cout << "slice,intensity" << std::endl;
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
    maskExtracter->SetInput( maskImageRegression );
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
      if( stats->GetMaximum( 1 ) > maxIntensityValue )
        {
        maxIntensityValue = stats->GetMaximum( 1 );
        }
      if( stats->GetMinimum( 1 ) < minIntensityValue )
        {
        minIntensityValue = stats->GetMinimum( 1 );
        }

      X.push_back( static_cast<RealType>( n ) );
      if( model == 0 )
        {
        Y.push_back( std::log( stats->GetMean( 1 ) ) );
        }
      else
        {
        Y.push_back( stats->GetMean( 1 ) );
        }
      totalAverageIntensityValue = totalAverageIntensityValue +
        ( stats->GetMean( 1 ) - totalAverageIntensityValue ) /
        static_cast<RealType>( nonZeroSliceCount + 1 );
      std::cout << X[nonZeroSliceCount] << "," << stats->GetMean( 1 ) << std::endl;
      nonZeroSliceCount++;
      }
    }

  /////////////////////////////
  //
  // Perform fitting
  //
  /////////////////////////////

  std::vector<RealType> line = FitRegressionLine( X, Y );

  std::cout << std::endl << std::endl;
  if( model == 0 )
    {
    std::cout << "Model:  y = " << std::exp( line[1] ) << " * exp(" << line[0] << "x)" << std::endl;;
    }
  else
    {
    std::cout << "Model:  y = " << line[1] << " + " << line[0] << "x" << std::endl;;
    }

  typename ImageType::Pointer outputImage = ImageType::New();
  outputImage->CopyInformation( inputImage );
  outputImage->SetRegions( inputImage->GetLargestPossibleRegion() );
  outputImage->Allocate();
  outputImage->FillBuffer( 0 );

  RealType maxResidualValue = itk::NumericTraits<RealType>::NonpositiveMin();
  RealType minResidualValue = itk::NumericTraits<RealType>::max();

  itk::ImageRegionIteratorWithIndex<ImageType> ItO( outputImage,
    outputImage->GetLargestPossibleRegion() );
  itk::ImageRegionConstIterator<ImageType> ItI( inputImage,
    inputImage->GetLargestPossibleRegion() );
  for( ItI.GoToBegin(), ItO.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItO )
    {
    typename ImageType::IndexType index = ItO.GetIndex();
    if( maskImageUpdate->GetPixel( index ) == 1 )
      {
      RealType predictedValue = line[1] + line[0] * static_cast<RealType>( index[direction] );
      if( model == 0 )
        {
        predictedValue = std::exp( line[1] ) * std::exp( line[0] * static_cast<RealType>( index[direction] ) );
        }
      RealType residual = ItI.Get() - predictedValue;
      ItO.Set( residual );

      if( residual > maxResidualValue )
        {
        maxResidualValue = residual;
        }
      if( residual < minResidualValue )
        {
        minResidualValue = residual;
        }
      }
    else
      {
      ItO.Set( ItI.Get() );
      }
    }

  ///////////
  //
  // now do a global rescale
  //

  RealType slope =  ( maxIntensityValue - minIntensityValue ) / ( maxResidualValue - minResidualValue );

  for( ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO )
    {
    typename ImageType::IndexType index = ItO.GetIndex();
    if( maskImageUpdate->GetPixel( index ) == 1 )
      {
      RealType rescaledValue = maxIntensityValue + slope * ( ItO.Get() - maxResidualValue );
      ItO.Set( rescaledValue );
      }
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
      << " inputImage outputImage <maskImageToDoRegression> <maskImageForUpdate> <direction=x,(y),z> <model=1>" << std::endl;

    std::cerr << "  model types: " << std::endl;
    std::cerr << "      0: exponential " << std::endl;
    std::cerr << "      1: linear " << std::endl;

    return EXIT_FAILURE;
    }

//   switch( atoi( argv[1] ) )
//    {
//    case 3:
     DirectionalBiasCorrection<3>( argc, argv );
//      break;
//    default:
//       std::cerr << "Unsupported dimension" << std::endl;
//       return EXIT_FAILURE;
//    }
}
