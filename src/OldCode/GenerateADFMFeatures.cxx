#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkCastImageFilter.h"
#include "itkGreyLevelCooccurrenceMatrixTextureCoefficientsCalculator.h"
#include "itkGreyLevelRunLengthMatrixTextureCoefficientsCalculator.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMaskedScalarImageToGreyLevelCooccurrenceMatrixGenerator.h"
#include "itkMaskedScalarImageToGreyLevelRunLengthMatrixGenerator.h"
#include "itkScalarToFractalImageFilter.h"

int main( int argc, char *argv[] )
{
  if ( argc != 5 )
    {
    std::cout << "Usage: " << argv[0] << "image labelImage label" << std::endl;
    exit( 1 );
    }

  typedef int PixelType;
  const unsigned int ImageDimension = 3;
  typedef float RealType;

  unsigned int numberOfBins = 100;
  PixelType pixelMax = -1000 ;
  PixelType pixelMin = 1000;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[1] );
  imageReader->Update();

  ReaderType::Pointer labelImageReader = ReaderType::New();
  labelImageReader->SetFileName( argv[2] );
  labelImageReader->Update();

  PixelType label = static_cast<PixelType>( atoi( argv[3] ) );
  
  ImageType::Pointer mask = ImageType::New();
  mask->SetRegions( imageReader->GetOutput()->GetLargestPossibleRegion() );
  mask->SetSpacing( imageReader->GetOutput()->GetSpacing() );
  mask->SetOrigin( imageReader->GetOutput()->GetOrigin() );
  mask->Allocate();
  mask->FillBuffer( itk::NumericTraits<PixelType>::Zero );

  itk::ImageRegionIterator<ImageType> ItI( imageReader->GetOutput(),
    imageReader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> ItL( labelImageReader->GetOutput(),
    labelImageReader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> ItM( mask,
    mask->GetLargestPossibleRegion() );

  unsigned long numberOfPixelsInsideMask = 0;

  for ( ItM.GoToBegin(), ItL.GoToBegin(); !ItM.IsAtEnd(); ++ItM, ++ItL )
    {
    if ( ItL.Get() == label )
      {
      ItM.Set( itk::NumericTraits<PixelType>::One );
      numberOfPixelsInsideMask++;
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
  stats->SetHistogramParameters( numberOfBins, pixelMin, pixelMax );
  stats->Update();

  mean = stats->GetMean( 1 );
  sum = stats->GetSum( 1 );
  sigma = stats->GetSigma( 1 );
  variance = sigma * sigma;
  
  kurtosis = 0.0;
  skewness = 0.0; 
 
  RealType N = 0.0;
  for ( ItI.GoToBegin(), ItL.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItL )
    {
    if ( ItL.Get() == label )
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

  entropy = 0.0;
  for( unsigned int i = 0; i < histogram->Size(); i++ )
    {
    RealType p = static_cast<RealType>( histogram->GetFrequency( i, 0 )  )
      / static_cast<RealType>( histogram->GetTotalFrequency() );
    entropy += ( -p * vcl_log( p ) / vcl_log( 2.0 ) );  
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
  





  /**
   * Second order measurements
   * These include:
   *   1. energy (f1) *cth, *amfm
   *   2. entropy (f2) *cth, *amfm
   *   3. correlation (f3) *amfm
   *   4. inverse difference moment (f4) *cth, *amfm
   *   5. inertia (f5) *cth, *amfm
   *   6. cluster shade (f6) *cth
   *   7. cluster prominence (f7) *cth
   *   8. haralick's correlation (f8)
   */

  typedef itk::Statistics::MaskedScalarImageToGreyLevelCooccurrenceMatrixGenerator<ImageType>
    CooccurrenceMatrixGeneratorType;
  CooccurrenceMatrixGeneratorType::Pointer generator = CooccurrenceMatrixGeneratorType::New();
  generator->SetInput( imageReader->GetOutput() );
  generator->SetImageMask( mask );
  generator->SetNumberOfBinsPerAxis( numberOfBins );
  generator->SetPixelValueMinMax( pixelMin, pixelMax );
//  generator->SetOffset();
  generator->Compute();   

  typedef itk::Statistics::GreyLevelCooccurrenceMatrixTextureCoefficientsCalculator
    <CooccurrenceMatrixGeneratorType::HistogramType> CalculatorType;
  CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetHistogram( generator->GetOutput() );
  calculator->Compute();

  RealType energy = calculator->GetEnergy();
  RealType entropy2 = calculator->GetEntropy();
  RealType correlation = calculator->GetCorrelation();
  RealType inverseDifferenceMoment = calculator->GetInverseDifferenceMoment();
  RealType inertia = calculator->GetInertia();
  RealType clusterShade = calculator->GetClusterShade();
  RealType clusterProminence = calculator->GetClusterProminence();
  RealType haralickCorrelation = calculator->GetHaralickCorrelation();


  /**
   * Calculate measurements from fractal image
   * These include:
   *   1. fractal mean
   *   2. fractal variance
   *   3. fractal skewness
   *   4. fractal kurtosis
   *   5. fractal entropy
   */

  RealType meanFractal;
  RealType sigmaFractal;
  RealType sumFractal;
  RealType varianceFractal;
  RealType skewnessFractal;
  RealType kurtosisFractal;
  RealType entropyFractal;

  typedef itk::ScalarToFractalImageFilter<ImageType, RealImageType> FractalFilterType;

  typedef itk::CastImageFilter<ImageType, FractalFilterType::MaskImageType> CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput( mask );
  caster->Update();

  FractalFilterType::Pointer fractal = FractalFilterType::New();
  fractal->SetInput( imageReader->GetOutput() );
  fractal->SetMaskImage( caster->GetOutput() );
  fractal->Update();
  
  typedef itk::LabelStatisticsImageFilter<RealImageType, ImageType> HistogramGeneratorType2;
  HistogramGeneratorType2::Pointer statsFractal = HistogramGeneratorType2::New();
  statsFractal->SetInput( fractal->GetOutput() );
  statsFractal->SetLabelInput( mask );
  statsFractal->UseHistogramsOn();
  statsFractal->SetHistogramParameters( numberOfBins, pixelMin, pixelMax );
  statsFractal->Update();

  meanFractal = statsFractal->GetMean( 1 );
  sumFractal = statsFractal->GetSum( 1 );
  sigmaFractal = statsFractal->GetSigma( 1 );
  varianceFractal = sigmaFractal * sigmaFractal;
  
  kurtosisFractal = 0.0;
  skewnessFractal = 0.0; 
 
  N = 0.0;
  for ( ItI.GoToBegin(), ItL.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItL )
    {
    if ( ItL.Get() == label )
      {
      RealType value = ItI.Get();

      RealType diff = value - mean; 
      skewnessFractal += ( diff * diff * diff );        
      kurtosisFractal += ( diff * diff * diff * diff );
 
      N++;
      }
    }   
  skewnessFractal /= ( ( N - 1 ) * varianceFractal * sigmaFractal );
  kurtosisFractal /= ( ( N - 1 ) * varianceFractal * varianceFractal );

  typedef HistogramGeneratorType2::HistogramType  HistogramType2;
  const HistogramType2 *histogramFractal = statsFractal->GetHistogram( 1 );

  entropy = 0.0;
  for( unsigned int i = 0; i < histogram->Size(); i++ )
    {
    RealType p = static_cast<RealType>( histogramFractal->GetFrequency( i, 0 )  )
      / static_cast<RealType>( histogramFractal->GetTotalFrequency() );
    entropyFractal += ( -p * vcl_log( p ) / vcl_log( 2.0 ) );   
    }




  /**
   * Second order measurements
   * These include:
   *   1. energy (f1) *cth, *amfm
   *   2. entropy (f2) *cth, *amfm
   *   3. correlation (f3) *amfm
   *   4. inverse difference moment (f4) *cth, *amfm
   *   5. inertia (f5) *cth, *amfm
   *   6. cluster shade (f6) *cth
   *   7. cluster prominence (f7) *cth
   *   8. haralick's correlation (f8)
   */

  {
  typedef itk::Statistics::MaskedScalarImageToGreyLevelRunLengthMatrixGenerator<ImageType>
    RunLengthMatrixGeneratorType;
  RunLengthMatrixGeneratorType::Pointer generator = RunLengthMatrixGeneratorType::New();
  generator->SetInput( imageReader->GetOutput() );
  generator->SetImageMask( mask );
  generator->SetNumberOfBinsPerAxis( numberOfBins );
  generator->SetPixelValueMinMax( pixelMin, pixelMax );
//  generator->SetOffset();
  generator->Compute();   

  typedef itk::Statistics::GreyLevelRunLengthMatrixTextureCoefficientsCalculator
    <RunLengthMatrixGeneratorType::HistogramType> CalculatorType;
  CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetHistogram( generator->GetOutput() );
  calculator->Compute();

  RealType sre = calculator->GetShortRunEmphasis();
  RealType lre = calculator->GetLongRunEmphasis();
  RealType gln = calculator->GetGreyLevelNonuniformity();
  RealType rln = calculator->GetRunLengthNonuniformity();
  RealType rp = static_cast<RealType>( calculator->GetTotalNumberOfRuns() ) 
    / static_cast<RealType>( numberOfPixelsInsideMask );
  RealType lgre = calculator->GetLowGreyLevelRunEmphasis();
  RealType hgre = calculator->GetHighGreyLevelRunEmphasis();
  RealType srlge = calculator->GetShortRunLowGreyLevelEmphasis();
  RealType srhge = calculator->GetShortRunHighGreyLevelEmphasis();
  RealType lrlge = calculator->GetLongRunLowGreyLevelEmphasis();
  RealType lrhge = calculator->GetLongRunHighGreyLevelEmphasis();
  }


/*

  std::string histogram_file = std::string( "histogram_" )
    + std::string( argv[4] ) + std::string( ".dat" );
  ofstream str( histogram_file.c_str() );

  unsigned int i;
  for( i = 0; i < histogram->Size(); i++ )
    {
    str << histogram->GetMeasurement( i, 0 ) << " "
        << histogram->GetFrequency( i, 0 ) << " "
        << histogram->GetTotalFrequency() << std::endl;
    }
    
  std::string stats_file = std::string( "stats_" )
    + std::string( argv[4] ) + std::string( ".dat" );
  ofstream str1( stats_file.c_str() );
*/


  return 0;
}
