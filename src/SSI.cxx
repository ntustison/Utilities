#include "itkArray.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkListSample.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodIterator.h"
#include "itkWeightedCovarianceSampleFilter.h"

#include <string>
#include <vector>

template<class TValue>
TValue Convert( std::string optionString )
			{
			TValue value;
			std::istringstream iss( optionString );
			iss >> value;
			return value;
			}

template<class TValue>
std::vector<TValue> ConvertVector( std::string optionString )
			{
			std::vector<TValue> values;
			std::string::size_type crosspos = optionString.find( 'x', 0 );

			if ( crosspos == std::string::npos )
					{
					values.push_back( Convert<TValue>( optionString ) );
					}
			else
					{
					std::string element = optionString.substr( 0, crosspos ) ;
					TValue value;
					std::istringstream iss( element );
					iss >> value;
					values.push_back( value );
					while ( crosspos != std::string::npos )
							{
							std::string::size_type crossposfrom = crosspos;
							crosspos = optionString.find( 'x', crossposfrom + 1 );
							if ( crosspos == std::string::npos )
									{
									element = optionString.substr( crossposfrom + 1, optionString.length() );
									}
							else
									{
									element = optionString.substr( crossposfrom + 1, crosspos ) ;
									}
							std::istringstream iss( element );
							iss >> value;
							values.push_back( value );
							}
					}
			return values;
			}

template <unsigned int ImageDimension>
int SSI( unsigned int argc, char *argv[] )
{
  typedef double                                     RealType;
  typedef itk::Image<RealType, ImageDimension>        ImageType;
  typedef itk::Image<int, ImageDimension>             MaskImageType;
  typedef itk::Array<RealType>                        MeasurementVectorType;
  typedef typename itk::Statistics::ListSample
    <MeasurementVectorType>                           SampleType;


  RealType K1 = 0.01;
  RealType K2 = 0.03;

  RealType C1 = vnl_math_sqr( K1 * 1.0 );
  RealType C2 = vnl_math_sqr( K2 * 1.0 );
  RealType C3 = 0.5 * C2;

  RealType alpha = 1.0;
  RealType beta = 1.0;
  RealType gamma = 1.0;


  typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 4 );

  typename MaskImageType::Pointer mask = NULL;
  if( argc > 5 )
    {
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    typename MaskReaderType::Pointer reader = MaskReaderType::New();
    reader->SetFileName( argv[5] );
    reader->Update();
    mask = reader->GetOutput();
    }

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[2] );
  reader1->Update();
  typename ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[3] );
  reader2->Update();

  typename ImageType::Pointer ssiImage = ImageType::New();
  ssiImage->SetOrigin( reader1->GetOutput()->GetOrigin() );
  ssiImage->SetSpacing( reader1->GetOutput()->GetSpacing() );
  ssiImage->SetRegions( reader1->GetOutput()->GetLargestPossibleRegion() );
  ssiImage->SetDirection( reader1->GetOutput()->GetDirection() );
  ssiImage->Allocate();
  ssiImage->FillBuffer( 0 );

  typedef itk::Array<double> WeightArrayType;

  NeighborhoodIteratorType It1( radius, reader1->GetOutput(),
    reader1->GetOutput()->GetLargestPossibleRegion() );
  NeighborhoodIteratorType It2( radius, reader2->GetOutput(),
    reader2->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> ItS( ssiImage,
    ssiImage->GetLargestPossibleRegion() );

  It1.GoToBegin();
  It2.GoToBegin();
  ItS.GoToBegin();

  while( !ItS.IsAtEnd() )
    {
    if( mask && mask->GetPixel( It1.GetIndex() ) == 0 )
      {
      ++ItS;
      ++It1;
      ++It2;

      continue;
      }

    typename SampleType::Pointer list = SampleType::New();
    list->SetMeasurementVectorSize( 2 );

    std::vector<RealType> w;

    for( unsigned int n = 0; n < ( It1.GetNeighborhood() ).Size(); n++ )
      {
      bool isInside = false;
      RealType val = It1.GetPixel( n, isInside );

      if( isInside )
        {
        MeasurementVectorType measurement;
        measurement.SetSize( 2 );
        measurement[0] = val;
        measurement[1] = It2.GetPixel( n );
        list->PushBack( measurement );

        typename NeighborhoodIteratorType::OffsetType offset =
          It1.GetOffset( n );

        RealType distance = 0.0;
        RealType variance = 0.0;
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          variance += ssiImage->GetSpacing()[d];
          distance += vnl_math_sqr( offset[d] * ssiImage->GetSpacing()[d] );
          }
        w.push_back( vcl_exp( -0.5 * distance / ( 1.5 * variance ) ) );
        }
      }

    WeightArrayType weights( w.size() );
    for( unsigned int n = 0; n < w.size(); n++ )
      {
      weights[n] = w[n];
      }

    // Calculate the statistics in the neighborhood of both images
    typedef typename itk::Statistics::
      WeightedCovarianceSampleFilter<SampleType> CovarianceCalculatorType;
    typename CovarianceCalculatorType::Pointer covarianceCalculator =
      CovarianceCalculatorType::New();
    covarianceCalculator->SetWeights( weights );
    covarianceCalculator->SetInput( list );
    covarianceCalculator->Update();

    MeasurementVectorType mean = covarianceCalculator->GetMean();
    typename CovarianceCalculatorType::MatrixType covariance =
      covarianceCalculator->GetCovarianceMatrix();

    RealType luminance = ( 2.0 * mean[0] * mean[1] + C1 ) /
      ( vnl_math_sqr( mean[0] ) + vnl_math_sqr( mean[1] ) + C1 );
    RealType contrast = ( 2.0 * vcl_sqrt( covariance( 0, 0 ) ) * vcl_sqrt( covariance( 1, 1 ) ) + C2 ) /
      ( covariance( 0, 0 ) + covariance( 1, 1 ) + C2 );
    RealType correlation = ( covariance( 1, 0 ) + C3 ) /
      ( vcl_sqrt( covariance( 0, 0 ) ) * vcl_sqrt( covariance( 1, 1 ) ) + C3 );

    ItS.Set(
      vcl_pow( static_cast<double>( luminance ), static_cast<double>( alpha ) ) *
      vcl_pow( static_cast<double>( contrast ), static_cast<double>( beta ) ) *
      vcl_pow( static_cast<double>( correlation ), static_cast<double>( gamma ) )
      );

    ++ItS;
    ++It1;
    ++It2;
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( ssiImage );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cerr << argv[0] << " imageDimension inputImage1 inputImage2 outputImage "
      << "[mask] [radius=4] [alpha=1.0] [beta=1.0] [gamma=1.0]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      SSI<2>( argc, argv );
      break;
    case 3:
      SSI<3>( argc, argv );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}

