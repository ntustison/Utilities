#include <stdio.h>

#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkGreyLevelCooccurrenceMatrixTextureCoefficientsCalculator.h"
#include "itkHistogram.h"
#include "itkLabelStatisticsImageFilter.h"

template <unsigned int ImageDimension>
int CalculateStatistics( int argc, char *argv[] )
{

  typedef int MaskPixelType;
  typedef float RealType;

  unsigned int numberOfBins = 200;
  if ( argc > 6 )
    {
    numberOfBins = atoi( argv[6] );
    }

  typedef itk::Image<MaskPixelType, ImageDimension> MaskImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  typename ReaderType::Pointer imageReader1 = ReaderType::New();
  imageReader1->SetFileName( argv[2] );
  imageReader1->Update();

  typename ReaderType::Pointer imageReader2 = ReaderType::New();
  imageReader2->SetFileName( argv[3] );
  imageReader2->Update();

  if( imageReader1->GetOutput()->GetLargestPossibleRegion().GetSize()
    != imageReader2->GetOutput()->GetLargestPossibleRegion().GetSize() )
    {
    std::cerr << "Error: Images are not the same size." << std::endl;
    return EXIT_FAILURE;
    }

  typename MaskImageType::Pointer mask = MaskImageType::New();
  if ( argc > 4 )
    {
    typedef itk::ImageFileReader<MaskImageType> ReaderType;
    typename ReaderType::Pointer labelImageReader = ReaderType::New();
    labelImageReader->SetFileName( argv[4] );
    labelImageReader->Update();
    mask = labelImageReader->GetOutput();
    }
  else
    {
    mask->SetOrigin( imageReader1->GetOutput()->GetOrigin() );
    mask->SetSpacing( imageReader1->GetOutput()->GetSpacing() );
    mask->SetRegions( imageReader1->GetOutput()->GetLargestPossibleRegion() );
    mask->SetDirection( imageReader1->GetOutput()->GetDirection() );
    mask->Allocate();
    mask->FillBuffer( itk::NumericTraits<MaskPixelType>::One );
    }

  MaskPixelType label = itk::NumericTraits<MaskPixelType>::One;
  if ( argc > 5 )
    {
    label = static_cast<MaskPixelType>( atoi( argv[5] ) );
    }

  itk::ImageRegionIterator<MaskImageType> ItM( mask,
    mask->GetLargestPossibleRegion() );

  /**
   * Calculate Pearson's correlation coefficient
   */

  /**
   * Calculate mean, sigma for image 1
   */
  itk::ImageRegionIterator<RealImageType> It1( imageReader1->GetOutput(),
    imageReader1->GetOutput()->GetLargestPossibleRegion() );

  RealType maxValue1 = itk::NumericTraits<RealType>::NonpositiveMin();
  RealType minValue1 = itk::NumericTraits<RealType>::max();

  for ( It1.GoToBegin(), ItM.GoToBegin(); !It1.IsAtEnd(); ++It1, ++ItM )
    {
    if( ItM.Get() == label )
      {
      if ( It1.Get() < minValue1 )
        {
        minValue1 = It1.Get();
        }
      else if ( It1.Get() > maxValue1 )
        {
        maxValue1 = It1.Get();
        }
      }
    }

  typedef itk::LabelStatisticsImageFilter<RealImageType, MaskImageType> HistogramGeneratorType;
  typename HistogramGeneratorType::Pointer stats1 = HistogramGeneratorType::New();
  stats1->SetInput( imageReader1->GetOutput() );
  stats1->SetLabelInput( mask );
  stats1->UseHistogramsOn();
  stats1->SetHistogramParameters( numberOfBins, minValue1, maxValue1 );
  stats1->Update();

  RealType mean1 = stats1->GetMean( label );
  RealType sigma1 = stats1->GetSigma( label );

  /**
   * Calculate mean, sigma for image 2
   */
  itk::ImageRegionIterator<RealImageType> It2( imageReader2->GetOutput(),
    imageReader2->GetOutput()->GetLargestPossibleRegion() );

  RealType maxValue2 = itk::NumericTraits<RealType>::NonpositiveMin();
  RealType minValue2 = itk::NumericTraits<RealType>::max();

  for ( It2.GoToBegin(), ItM.GoToBegin(); !It2.IsAtEnd(); ++It2, ++ItM )
    {
    if( ItM.Get() == label )
      {
      if ( It2.Get() < minValue2 )
        {
        minValue2 = It2.Get();
        }
      else if ( It2.Get() > maxValue2 )
        {
        maxValue2 = It2.Get();
        }
      }
    }

  typename HistogramGeneratorType::Pointer stats2 = HistogramGeneratorType::New();
  stats2->SetInput( imageReader2->GetOutput() );
  stats2->SetLabelInput( mask );
  stats2->UseHistogramsOn();
  stats2->SetHistogramParameters( numberOfBins, minValue2, maxValue2 );
  stats2->Update();

  RealType mean2 = stats2->GetMean( label );
  RealType sigma2 = stats2->GetSigma( label );

  RealType pearson = 0.0;
  RealType numberOfVoxels = 0.0;

  It1.GoToBegin();
  It2.GoToBegin();
  ItM.GoToBegin();
  while( !ItM.IsAtEnd() )
    {
    if( ItM.Get() == label )
      {
      pearson += ( ( It1.Get() - mean1 ) / sigma1 )
        * ( ( It2.Get() - mean2 ) / sigma2 );

      numberOfVoxels++;
      }
    ++ItM;
    ++It1;
    ++It2;
    }
  pearson /= numberOfVoxels;

  /**
   * Calculate cooccurrence matrix between the two images
   */

  typedef itk::Statistics::Histogram<RealType, 2> HistogramType;
  typename HistogramType::Pointer histogram = HistogramType::New();

  typename HistogramType::SizeType histogramSize;
  typename HistogramType::MeasurementVectorType histogramLowerBounds;
  typename HistogramType::MeasurementVectorType histogramUpperBounds;

  histogramSize[0] = numberOfBins;
  histogramSize[1] = numberOfBins;
  histogramLowerBounds[0] = minValue1;
  histogramLowerBounds[1] = minValue2;
  histogramUpperBounds[0] = maxValue1;
  histogramUpperBounds[1] = maxValue2;

  histogram->Initialize( histogramSize, histogramLowerBounds, histogramUpperBounds );
  histogram->SetClipBinsAtEnds( false );

  It1.GoToBegin();
  It2.GoToBegin();
  ItM.GoToBegin();
  while( !ItM.IsAtEnd() )
    {
    if( ItM.Get() == label )
      {
      typename HistogramType::MeasurementVectorType measurement;
      measurement[0] = It1.Get();
      measurement[1] = It2.Get();
      histogram->IncreaseFrequency( measurement,
        static_cast<typename HistogramType::FrequencyType>( 1 ) );
      }
    ++ItM;
    ++It1;
    ++It2;
    }

  typedef itk::Statistics::GreyLevelCooccurrenceMatrixTextureCoefficientsCalculator
    <HistogramType> CalculatorType;
  typename CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetHistogram( histogram );
  calculator->Compute();

  RealType energy = calculator->GetEnergy();
  RealType entropy = calculator->GetEntropy();
  RealType correlation = calculator->GetCorrelation();
  RealType inverseDifferenceMoment = calculator->GetInverseDifferenceMoment();
  RealType inertia = calculator->GetInertia();
  RealType clusterShade = calculator->GetClusterShade();
  RealType clusterProminence = calculator->GetClusterProminence();
  RealType haralickCorrelation = calculator->GetHaralickCorrelation();

  std::cout << "[" << argv[0] << "]" << std::endl;
  std::cout << pearson << " "
            << energy << " "
            << entropy << " "
            << correlation << " "
            << inverseDifferenceMoment << " "
            << inertia << " "
            << clusterShade << " "
            << clusterProminence << " "
            << haralickCorrelation << std::endl;

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension inputImage1 inputImage2 "
              << "[labelImage] [label]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     CalculateStatistics<2>( argc, argv );
     break;
   case 3:
     CalculateStatistics<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
