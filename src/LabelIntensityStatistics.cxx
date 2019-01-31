#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkLabelStatisticsImageFilter.h"

#include <vector>
#include <algorithm>
#include <iomanip>

template <unsigned int ImageDimension>
int LabelIntensityStatistics( int argc, char *argv[] )
{
  typedef int LabelType;
  typedef float RealType;

  typedef itk::Image<LabelType, ImageDimension> LabelImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  typename ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[2] );
  imageReader->Update();

		typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
		typename LabelReaderType::Pointer labelImageReader = LabelReaderType::New();
		labelImageReader->SetFileName( argv[3] );
		labelImageReader->Update();

  std::vector<LabelType> labels;
  itk::ImageRegionIterator<LabelImageType> It( labelImageReader->GetOutput(),
    labelImageReader->GetOutput()->GetLargestPossibleRegion() );

  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( It.Get() != 0 &&
      std::find( labels.begin(), labels.end(), It.Get() ) == labels.end() )
      {
      labels.push_back( It.Get() );
      }
    }
  std::sort( labels.begin(), labels.end() );

  RealType maxValue = itk::NumericTraits<RealType>::NonpositiveMin();
  RealType minValue = itk::NumericTraits<RealType>::max();

  itk::ImageRegionConstIterator<RealImageType> ItI( imageReader->GetOutput(),
    imageReader->GetOutput()->GetRequestedRegion() );
  for( ItI.GoToBegin(), It.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++It )
    {
    if( It.Get() == 0 )
      {
      continue;
      }
    RealType value = ItI.Get();

    if( value < minValue )
      {
      minValue = value;
      }
    else if( value > maxValue )
      {
      maxValue = value;
      }
    }

  typedef itk::LabelStatisticsImageFilter<RealImageType, LabelImageType> HistogramGeneratorType;
  typename HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
  stats->SetInput( imageReader->GetOutput() );
  stats->SetLabelInput( labelImageReader->GetOutput() );
  stats->UseHistogramsOn();
  stats->SetHistogramParameters( 200, minValue, maxValue );
  stats->Update();

  // Calculate moments for skewness and kurtosis calculations

  std::vector<RealType> m3( labels.size(), 0.0 );
  std::vector<RealType> m4( labels.size(), 0.0 );
  std::vector<RealType> N( labels.size(), 0.0 );
  for( ItI.GoToBegin(), It.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++It )
    {
    LabelType label = It.Get();
    if( label == 0 )
      {
      continue;
      }
    typename std::vector<LabelType>::iterator it =
      std::find( labels.begin(), labels.end(), label );
    if( it == labels.end() )
      {
      std::cerr << "Label not found.  Shouldn't get here." << std::endl;
      return EXIT_FAILURE;
      }
    RealType difference = ItI.Get() - stats->GetMean( label );

    unsigned long index = it - labels.begin();

    m3[index] += ( difference * difference * difference );
    m4[index] += ( m3[index] * difference );
    N[index] += 1.0;
    }
  for( unsigned int n = 0; n < N.size(); n++ )
    {
    m3[n] /= N[n];
    m4[n] /= N[n];
    }

//   std::cout << "                                       "
//             << "************ Individual Labels *************" << std::endl;
  std::cout << std::setw( 8 )  << "Label"
            << std::setw( 14 ) << "Mean"
            << std::setw( 14 ) << "Sigma"
            << std::setw( 14 ) << "Skewness"
            << std::setw( 14 ) << "Kurtosis"
            << std::setw( 14 ) << "Entropy"
            << std::setw( 14 ) << "Sum"
            << std::setw( 14 ) << "5th%"
            << std::setw( 14 ) << "95th%"
            << std::setw( 14 ) << "Min"
            << std::setw( 14 ) << "Max" << std::endl;

  std::vector<LabelType>::iterator it;
  for( it = labels.begin(); it != labels.end(); ++it )
    {
    typedef typename HistogramGeneratorType::HistogramType HistogramType;
    const HistogramType *histogram = stats->GetHistogram( *it );

    RealType fifthPercentileValue = histogram->Quantile( 0, 0.05 );
    RealType ninetyFifthPercentileValue = histogram->Quantile( 0, 0.95 );
    RealType entropy = 0.0;
    for( unsigned int i = 0; i < histogram->Size(); i++ )
      {
      RealType p = static_cast<RealType>( histogram->GetFrequency( i, 0 )  )
        / static_cast<RealType>( histogram->GetTotalFrequency() );
      if ( p > 0 )
        {
        entropy += ( -p * std::log( p ) / std::log( 2.0 ) );
        }
      }

    typename std::vector<LabelType>::iterator it2 =
      std::find( labels.begin(), labels.end(), *it );
    unsigned long index = it2 - labels.begin();

    RealType m2 = std::pow( stats->GetSigma( *it ), 2 );
    RealType k2 = ( N[index] ) / ( N[index] - 1.0 ) * m2;

    RealType prefactor3 = std::pow( N[index], 2 ) / ( ( N[index] - 1.0 ) * ( N[index] - 2.0 ) );
    RealType k3 = prefactor3 * m3[index];

    RealType prefactor4 = std::pow( N[index], 2 ) / ( ( N[index] - 1.0 ) * ( N[index] - 2.0 ) * ( N[index] - 3.0 ) );
    RealType k4 = prefactor4 * ( ( N[index] + 1 ) * m4[index] - 3 * ( N[index] - 1 ) ) * std::pow( m2, 2 );

    RealType skewness = k3 / std::sqrt( k2 * k2 * k2 );
    RealType kurtosis = k4 / std::pow( k2, 2 );

    std::cout << std::setw( 8  ) << *it;
    std::cout << std::setw( 14 ) << stats->GetMean( *it );
    std::cout << std::setw( 14 ) << stats->GetSigma( *it );
    std::cout << std::setw( 14 ) << skewness;
    std::cout << std::setw( 14 ) << kurtosis;
    std::cout << std::setw( 14 ) << entropy;
    std::cout << std::setw( 14 ) << stats->GetSum( *it );
    std::cout << std::setw( 14 ) << fifthPercentileValue;
    std::cout << std::setw( 14 ) << ninetyFifthPercentileValue;
    std::cout << std::setw( 14 ) << stats->GetMinimum( *it );
    std::cout << std::setw( 14 ) << stats->GetMaximum( *it );
    std::cout << std::endl;
    }

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension inputImage labelImage" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     LabelIntensityStatistics<2>( argc, argv );
     break;
   case 3:
     LabelIntensityStatistics<3>( argc, argv );
     break;
   case 4:
     LabelIntensityStatistics<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
