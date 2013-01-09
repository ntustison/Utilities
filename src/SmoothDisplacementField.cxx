#include "itkDisplacementFieldToBSplineImageFilter.h"

#include "itkGaussianOperator.h"
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTimeProbe.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkVectorNeighborhoodOperatorImageFilter.h"

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
int SmoothDisplacementField( int argc, char *argv[] )
{

  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;

  typedef itk::ImageFileReader<DisplacementFieldType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename DisplacementFieldType::Pointer field = reader->GetOutput();

  typename DisplacementFieldType::Pointer smoothField = DisplacementFieldType::New();

  float elapsedTime = 0.0;

  std::vector<float> var = ConvertVector<float>( std::string( argv[4] ) );
  if( var.size() == 1 )
    {
    float variance = var[0];

    typedef itk::GaussianOperator<float, ImageDimension>                                             GaussianSmoothingOperatorType;
    typedef itk::VectorNeighborhoodOperatorImageFilter<DisplacementFieldType, DisplacementFieldType> GaussianSmoothingSmootherType;

    GaussianSmoothingOperatorType gaussianSmoothingOperator;

    typedef itk::ImageDuplicator<DisplacementFieldType> DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage( field );
    duplicator->Update();

    smoothField = duplicator->GetOutput();

    itk::TimeProbe timer;
    timer.Start();
    typename GaussianSmoothingSmootherType::Pointer smoother = GaussianSmoothingSmootherType::New();

    for( unsigned int dimension = 0; dimension < ImageDimension; ++dimension )
      {
      // smooth along this dimension
      gaussianSmoothingOperator.SetDirection( dimension );
      gaussianSmoothingOperator.SetVariance( variance );
      gaussianSmoothingOperator.SetMaximumError( 0.001 );
      gaussianSmoothingOperator.SetMaximumKernelWidth( smoothField->GetRequestedRegion().GetSize()[dimension] );
      gaussianSmoothingOperator.CreateDirectional();

      // todo: make sure we only smooth within the buffered region
      smoother->SetOperator( gaussianSmoothingOperator );
      smoother->SetInput( smoothField );
      smoother->Update();

      smoothField = smoother->GetOutput();
      smoothField->Update();
      smoothField->DisconnectPipeline();
      }

    const VectorType zeroVector( 0.0 );

    //make sure boundary does not move
    float weight1 = 1.0;
    if (variance < 0.5)
      {
      weight1 = 1.0 - 1.0 * ( variance / 0.5);
      }
    float weight2 = 1.0 - weight1;

    const typename DisplacementFieldType::RegionType region = field->GetLargestPossibleRegion();
    const typename DisplacementFieldType::SizeType size = region.GetSize();
    const typename DisplacementFieldType::IndexType startIndex = region.GetIndex();

    itk::ImageRegionConstIteratorWithIndex<DisplacementFieldType> fieldIt( field, field->GetLargestPossibleRegion() );
    itk::ImageRegionIteratorWithIndex<DisplacementFieldType> smoothedFieldIt( smoothField, smoothField->GetLargestPossibleRegion() );
    for( fieldIt.GoToBegin(), smoothedFieldIt.GoToBegin(); !fieldIt.IsAtEnd(); ++fieldIt, ++smoothedFieldIt )
      {
      typename DisplacementFieldType::IndexType index = fieldIt.GetIndex();
      bool isOnBoundary = false;
      for ( unsigned int dimension = 0; dimension < ImageDimension; ++dimension )
        {
        if( index[dimension] == startIndex[dimension] || index[dimension] == static_cast<int>( size[dimension] ) - startIndex[dimension] - 1 )
          {
          isOnBoundary = true;
          break;
          }
        }
      if( isOnBoundary )
        {
        smoothedFieldIt.Set( zeroVector );
        }
      else
        {
        smoothedFieldIt.Set( smoothedFieldIt.Get() * weight1 + fieldIt.Get() * weight2 );
        }
      }
    timer.Stop();
    elapsedTime = timer.GetMean();

    }
  else if( var.size() == ImageDimension )
    {
    typedef itk::DisplacementFieldToBSplineImageFilter<DisplacementFieldType, DisplacementFieldType> BSplineFilterType;

    unsigned int numberOfLevels = 1;
    if( argc > 5 )
      {
      numberOfLevels = atoi( argv[5] );
      }

    unsigned int splineOrder = 3;
    if( argc > 6 )
      {
      splineOrder = atoi( argv[6] );
      }

    typename BSplineFilterType::ArrayType ncps;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      ncps[d] = var[d] + splineOrder;
      }

    typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
    bspliner->SetDisplacementField( field );
    bspliner->SetNumberOfControlPoints( ncps );
    bspliner->SetSplineOrder( splineOrder );
    bspliner->SetNumberOfFittingLevels( numberOfLevels );
    bspliner->SetEnforceStationaryBoundary( true );
    if( argc > 7 )
      {
      typedef itk::ImageFileReader<typename BSplineFilterType::RealImageType> ConfidenceImageReaderType;
      typename ConfidenceImageReaderType::Pointer cReader = ConfidenceImageReaderType::New();
      cReader->SetFileName( argv[7] );
      cReader->Update();

      bspliner->SetConfidenceImage( cReader->GetOutput() );
      }
    bspliner->SetEstimateInverse( false );

    itk::TimeProbe timer;
    timer.Start();
    bspliner->Update();
    timer.Stop();
    elapsedTime = timer.GetMean();

    smoothField = bspliner->GetOutput();
    smoothField->DisconnectPipeline();
    }
  else
    {
    std::cerr << "Error: unexpected variance format." << std::endl;
    return EXIT_FAILURE;
    }

  typename RealImageType::Pointer rmseImage = RealImageType::New();
  rmseImage->SetOrigin( field->GetOrigin() );
  rmseImage->SetDirection( field->GetDirection() );
  rmseImage->SetSpacing( field->GetSpacing() );
  rmseImage->SetRegions( field->GetLargestPossibleRegion() );
  rmseImage->Allocate();
  rmseImage->FillBuffer( 0.0 );

  itk::ImageRegionConstIterator<DisplacementFieldType> fieldIt( field, field->GetLargestPossibleRegion() );
  itk::ImageRegionConstIterator<DisplacementFieldType> smoothedFieldIt( smoothField, smoothField->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<RealImageType> ItR( rmseImage, rmseImage->GetLargestPossibleRegion() );

  float rmse = 0.0;
  vnl_vector<float> rmse_comp( 2 );
  rmse_comp.fill( 0.0 );
  float N = 0.0;
  for( fieldIt.GoToBegin(), smoothedFieldIt.GoToBegin(), ItR.GoToBegin(); !fieldIt.IsAtEnd(); ++fieldIt, ++smoothedFieldIt, ++ItR )
    {
    ItR.Set( ( fieldIt.Get() - smoothedFieldIt.Get() ).GetSquaredNorm() );
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      rmse_comp[d] += vnl_math_sqr( fieldIt.Get()[d] - smoothedFieldIt.Get()[d] );
      }
    rmse += ( fieldIt.Get() - smoothedFieldIt.Get() ).GetSquaredNorm();
    N += 1.0;
    }
  rmse = vcl_sqrt( rmse / N );

  std::cout << "Elapsed time: " << elapsedTime << std::endl;
  std::cout << "RMSE = " << rmse << std::endl;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    std::cout << "  rmse[" << d << "] = " << vcl_sqrt( rmse_comp[d] / N ) << std::endl;
    }

  {
  typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( smoothField );
  writer->Update();
  }

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << argv[0] << " imageDimension inputField outputField variance_or_mesh_size_base_level "
              << "[numberOfevels=1] [splineOrder=3] [confidenceImage]" << std::endl;
    exit( 0 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     SmoothDisplacementField<2>( argc, argv );
     break;
   case 3:
     SmoothDisplacementField<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

