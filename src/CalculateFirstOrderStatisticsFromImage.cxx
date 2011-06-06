#include <stdio.h>

#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkLabelStatisticsImageFilter.h"

#include <fstream.h>

template <unsigned int ImageDimension>
int CalculateFirstOrderStatistics( int argc, char *argv[] )
{

  itk::MultiThreader::SetGlobalDefaultNumberOfThreads( 1 );

  typedef int PixelType;
  typedef float RealType;

  unsigned int numberOfBins = 200;
  if ( argc > 5 )
    {
    numberOfBins = atoi( argv[5] );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  typename ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[2] );
  imageReader->Update();

  typename ImageType::Pointer mask = ImageType::New();
  if ( argc > 3 )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer labelImageReader = ReaderType::New();
    labelImageReader->SetFileName( argv[3] );
    labelImageReader->Update();
    mask = labelImageReader->GetOutput();
    }
  else
    {
    mask->SetOrigin( imageReader->GetOutput()->GetOrigin() );
    mask->SetSpacing( imageReader->GetOutput()->GetSpacing() );
    mask->SetRegions( imageReader->GetOutput()->GetLargestPossibleRegion() );
    mask->SetDirection( imageReader->GetOutput()->GetDirection() );
    mask->Allocate();
    mask->FillBuffer( itk::NumericTraits<PixelType>::One );
    }
  PixelType label = itk::NumericTraits<PixelType>::One;
  if ( argc > 4 )
    {
    label = static_cast<PixelType>( atoi( argv[4] ) );
    }

  itk::ImageRegionIterator<RealImageType> ItI( imageReader->GetOutput(),
    imageReader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> ItM( mask,
    mask->GetLargestPossibleRegion() );

  RealType maxValue = itk::NumericTraits<RealType>::NonpositiveMin();
  RealType minValue = itk::NumericTraits<RealType>::max();

  for ( ItM.GoToBegin(), ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItM, ++ItI )
    {
    if ( ItM.Get() == label )
      {
      if ( ItI.Get() < minValue )
        {
        minValue = ItI.Get();
        }
      else if ( ItI.Get() > maxValue )
        {
        maxValue = ItI.Get();
        }
      ItM.Set( itk::NumericTraits<PixelType>::One );
      }
    else
      {
      ItM.Set( itk::NumericTraits<PixelType>::Zero );
      }
    if( std::isnan( ItI.Get() ) || std::isinf( ItI.Get() ) )
      {
      ItM.Set( itk::NumericTraits<PixelType>::Zero );
      }
    }

  /**
   * Zeroth and first order measurements
   * These include:
   *   1. mean
   *   2. variance
   *   3. kurtosis
   *   4. skewness
   *   5. entropy
   *   6. fifth percentile value
   *   7. ninety-fifth percentile value
   *   8. mean of lower fifth percentile
   *   9. mean of upper fifth percentile
   */

  RealType mean;
  RealType sigma;
  RealType sum;
  RealType variance;
  RealType skewness;
  RealType kurtosis;
  RealType median;
  RealType entropy = 0;

  typedef itk::LabelStatisticsImageFilter<RealImageType, ImageType> HistogramGeneratorType;
  typename HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
  stats->SetInput( imageReader->GetOutput() );
  stats->SetLabelInput( mask );
  stats->SetUseHistograms( argc > 5 );
  if( stats->GetUseHistograms() )
    {
    stats->SetHistogramParameters( numberOfBins, minValue, maxValue );
    }
  stats->Update();

  mean = stats->GetMean( 1 );
  sum = stats->GetSum( 1 );
  sigma = stats->GetSigma( 1 );
  variance = sigma * sigma;
  median = stats->GetMedian( 1 );

  kurtosis = 0.0;
  skewness = 0.0;

  RealType N = 0.0;
  for ( ItI.GoToBegin(), ItM.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItM )
    {
    if ( ItM.Get() == itk::NumericTraits<PixelType>::One )
      {
      RealType value = ItI.Get();

      RealType diff = value - mean;
      skewness += ( diff * diff * diff );
      kurtosis += ( diff * diff * diff * diff );

      N += 1.0;
      }
    }
  skewness /= ( ( N - 1 ) * variance * sigma );
  kurtosis /= ( ( N - 1 ) * variance * variance );

		double fifthPercentileValue = 0.0;
		double ninetyFifthPercentileValue = 0.0;

		double fifthPercentileMean = 0.0;
		double fifthN = 0.0;
		double ninetyFifthPercentileMean = 0.0;
		double ninetyFifthN = 0.0;

  if ( argc > 6 )
    {
				typedef typename HistogramGeneratorType::HistogramType  HistogramType;
				const HistogramType *histogram = stats->GetHistogram( 1 );

				if( !histogram )
						{
						std::cerr << "ERROR:  No histogram created." << std::endl;
						exit( 0 );
						}

//     ofstream str( argv[6] );
//     for( unsigned int i = 0; i < histogram->Size(); i++ )
//       {
//       str << histogram->GetMeasurement( i, 0 ) << " "
//           << histogram->GetFrequency( i, 0 ) << std::endl;
//       }
//     str.close();

				ofstream str( argv[6] );
				for ( ItM.GoToBegin(), ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItM, ++ItI )
						{
						if ( ItM.Get() == 1 )
								{
								str << ItI.Get() << std::endl;
								}
						}
				str.close();

				fifthPercentileValue = histogram->Quantile( 0, 0.05 );
				ninetyFifthPercentileValue = histogram->Quantile( 0, 0.95 );

				entropy = 0.0;
				for( unsigned int i = 0; i < histogram->Size(); i++ )
						{
						RealType p = static_cast<RealType>( histogram->GetFrequency( i, 0 )  )
								/ static_cast<RealType>( histogram->GetTotalFrequency() );
						if ( p > 0 )
								{
								entropy += ( -p * vcl_log( p ) / vcl_log( 2.0 ) );
								}
						}

				for ( ItM.GoToBegin(), ItI.GoToBegin(); !ItM.IsAtEnd(); ++ItM, ++ItI )
						{
						RealType value = ItI.Get();
						if ( value <= fifthPercentileValue )
								{
								fifthPercentileMean += value;
								fifthN++;
								}
						else if ( value >= ninetyFifthPercentileValue )
								{
								ninetyFifthPercentileMean += value;
								ninetyFifthN++;
								}
						}

				fifthPercentileMean /= fifthN;
				ninetyFifthPercentileMean /= ninetyFifthN;
    }

  std::cout << "[" << argv[0] << "]" << std::endl;
  std::cout << mean << " "
            << sigma << " "
            << sum << " "
            << skewness << " "
            << kurtosis << " "
            << entropy << " "
            << fifthPercentileValue << " "
            << ninetyFifthPercentileValue << " "
            << fifthPercentileMean << " "
            << ninetyFifthPercentileMean << " "
            << minValue << " "
            << maxValue << " "
            << median << std::endl;

/*
  std::cout << "mean:        " << mean << std::endl;
  std::cout << "sigma:       " << sigma << std::endl;
  std::cout << "sum:         " << sum << std::endl;
  std::cout << "skewness:    " << skewness << std::endl;
  std::cout << "kurtosis:    " << kurtosis << std::endl;
  std::cout << "entropy:     " << entropy << std::endl;
  std::cout << "5th %:       " << fifthPercentileValue << std::endl;
  std::cout << "95th %:      " << ninetyFifthPercentileValue << std::endl;
  std::cout << "5th % mean:  " << fifthPercentileMean << std::endl;
  std::cout << "95th % mean: " << ninetyFifthPercentileMean << std::endl;
  std::cout << "Min value:   " << minValue << std::endl;
  std::cout << "Max value:   " << maxValue << std::endl;
  std::cout << "Median:   " << median << std::endl;
*/


  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension inputImage "
              << "[labelImage] [label] [numberOfBins=100] [histogramFile]" << std::endl;

    std::cerr << "  Output:  mean sigma sum skewness kurtosis entropy 5th% "
              << "95th% 5th%mean 95th%mean min max median"<< std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     CalculateFirstOrderStatistics<2>( argc, argv );
     break;
   case 3:
     CalculateFirstOrderStatistics<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
