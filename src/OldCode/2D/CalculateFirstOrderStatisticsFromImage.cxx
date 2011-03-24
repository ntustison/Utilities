#include <stdio.h>

#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkLabelStatisticsImageFilter.h"

#include <fstream.h>

int main( int argc, char *argv[] )
{
  if ( argc < 2 )
    {
    std::cout << "Usage: " << argv[0] << " image [labelImage] [label] [numberOfBins=100] [histogramFile]" << std::endl;
    exit( 1 );
    }

  typedef int PixelType;
  const unsigned int ImageDimension = 3;
  typedef float RealType;

  unsigned int numberOfBins = 100;
  if ( argc > 4 )
    {
    numberOfBins = atoi( argv[4] );
    }


  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[1] );
  imageReader->Update();

  ImageType::Pointer mask = ImageType::New();
  if ( argc > 2 )
    {
    ReaderType::Pointer labelImageReader = ReaderType::New();
    labelImageReader->SetFileName( argv[2] );
    labelImageReader->Update();
    mask = labelImageReader->GetOutput();
    }
  else
    {
    mask->SetOrigin( imageReader->GetOutput()->GetOrigin() );
    mask->SetSpacing( imageReader->GetOutput()->GetSpacing() );
    mask->SetRegions( imageReader->GetOutput()->GetLargestPossibleRegion() );
    mask->Allocate();
    mask->FillBuffer( itk::NumericTraits<PixelType>::Zero );
    }
  PixelType label = itk::NumericTraits<PixelType>::One;
  if ( argc > 3 )
    {
    label = static_cast<PixelType>( atoi( argv[3] ) );
    }

  itk::ImageRegionIterator<ImageType> ItI( imageReader->GetOutput(),
    imageReader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> ItM( mask,
    mask->GetLargestPossibleRegion() );

  PixelType maxValue = itk::NumericTraits<PixelType>::min();
  PixelType minValue = itk::NumericTraits<PixelType>::max();

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
  RealType entropy;

  typedef itk::LabelStatisticsImageFilter<ImageType, ImageType> HistogramGeneratorType;
  HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
  stats->SetInput( imageReader->GetOutput() );
  stats->SetLabelInput( mask );
  stats->UseHistogramsOn();
  stats->SetHistogramParameters( numberOfBins, minValue, maxValue );
  stats->Update();

  mean = stats->GetMean( 1 );
  sum = stats->GetSum( 1 );
  sigma = stats->GetSigma( 1 );
  variance = sigma * sigma;

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

      N++;
      }
    }
  skewness /= ( ( N - 1 ) * variance * sigma );
  kurtosis /= ( ( N - 1 ) * variance * variance );

  typedef HistogramGeneratorType::HistogramType  HistogramType;
  const HistogramType *histogram = stats->GetHistogram( 1 );

  if ( argc > 5 )
    {
    ofstream str( argv[5] );
    for( unsigned int i = 0; i < histogram->Size(); i++ )
      {
      str << histogram->GetMeasurement( i, 0 ) << " "
          << histogram->GetFrequency( i, 0 ) << std::endl;
      }
    str.close();
    }

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

  double fifthPercentileValue = histogram->Quantile( 0, 0.05 );
  double ninetyFifthPercentileValue = histogram->Quantile( 0, 0.95 );

  double fifthPercentileMean = 0.0;
  double fifthN = 0.0;
  double ninetyFifthPercentileMean = 0.0;
  double ninetyFifthN = 0.0;

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
            << ninetyFifthPercentileMean << std::endl;

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
*/


  return 0;
}
