#include "itkJensenHavrdaCharvatTsallisLabeledPointSetMetric.h"
#include "itkEuclideanDistancePointMetric.h"

#include "itkLabeledPointSetFileReader.h"
#include "itkIdentityTransform.h"
#include "itkPointSet.h"
#include "itkVector.h"

template <unsigned int ImageDimension>
int JHCT( int argc, char *argv[] )
{
  const unsigned int Dimension = ImageDimension;

  typedef float RealType;

  typedef long PointDataType;
  typedef itk::PointSet<PointDataType, Dimension> PointSetType;

  typedef itk::LabeledPointSetFileReader<PointSetType> ReaderType;

  typename ReaderType::Pointer fixedPointSetReader = ReaderType::New();
  fixedPointSetReader->SetFileName( argv[2] );
  fixedPointSetReader->Update();

  typename ReaderType::Pointer movingPointSetReader = ReaderType::New();
  movingPointSetReader->SetFileName( argv[3] );
  movingPointSetReader->Update();

  typedef itk::JensenHavrdaCharvatTsallisLabeledPointSetMetric<PointSetType>
    PointSetMetricType;
  typename PointSetMetricType::Pointer pointSetMetric = PointSetMetricType::New();

  pointSetMetric->SetFixedPointSet( fixedPointSetReader->GetOutput() );
  pointSetMetric->SetFixedPointSetSigma( atof( argv[5] ) );
  pointSetMetric->SetFixedEvaluationKNeighborhood( atoi( argv[7] ) );

  pointSetMetric->SetMovingPointSet( movingPointSetReader->GetOutput() );
  pointSetMetric->SetMovingPointSetSigma( atof( argv[5] ) );
  pointSetMetric->SetMovingEvaluationKNeighborhood( atoi( argv[7] ) );

  pointSetMetric->SetUseRegularizationTerm( atoi( argv[12] ) );
  pointSetMetric->SetUseInputAsSamples( atoi( argv[11] ) );
  pointSetMetric->SetAlpha( atof( argv[10] ) );
  pointSetMetric->SetUseAnisotropicCovariances( atoi( argv[9] ) );

  if( pointSetMetric->GetUseAnisotropicCovariances() )
    {
    pointSetMetric->SetFixedCovarianceKNeighborhood( atoi( argv[8] ) );
    pointSetMetric->SetFixedKernelSigma( atof( argv[6] ) );
    pointSetMetric->SetMovingCovarianceKNeighborhood( atoi( argv[8] ) );
    pointSetMetric->SetMovingKernelSigma( atof( argv[6] ) );
    }

  if( !pointSetMetric->GetUseInputAsSamples() )
    {
    pointSetMetric->SetNumberOfFixedSamples( atoi( argv[13] ) );
    pointSetMetric->SetNumberOfMovingSamples( atoi( argv[14] ) );
    }

  pointSetMetric->Initialize();
  typename PointSetMetricType::DefaultTransformType::ParametersType parameters;
  parameters.Fill( 0 );

  typename PointSetMetricType::MeasureType measureFixed
    = pointSetMetric->GetValue( parameters );

  std::cout << measureFixed[0] << std::endl << std::endl;

  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension>
int ICP( int argc, char *argv[] )
{
  const unsigned int Dimension = ImageDimension;

  typedef float RealType;

  typedef long PointDataType;
  typedef itk::PointSet<PointDataType, Dimension> PointSetType;
  typedef itk::Image<RealType, Dimension> RealImageType;

  typedef itk::LabeledPointSetFileReader<PointSetType> ReaderType;

  typename ReaderType::Pointer fixedPointSetReader = ReaderType::New();
  fixedPointSetReader->SetFileName( argv[2] );
  fixedPointSetReader->Update();

  typename ReaderType::Pointer movingPointSetReader = ReaderType::New();
  movingPointSetReader->SetFileName( argv[3] );
  movingPointSetReader->Update();

  typedef itk::EuclideanDistancePointMetric
    <PointSetType,PointSetType,RealImageType> PointSetMetricType;
  typename PointSetMetricType::Pointer pointSetMetric = PointSetMetricType::New();

  typedef itk::IdentityTransform<double, Dimension> DefaultTransformType;
  typename DefaultTransformType::Pointer transform = DefaultTransformType::New();
  transform->SetIdentity();

  pointSetMetric->SetFixedPointSet( fixedPointSetReader->GetOutput() );
  pointSetMetric->SetMovingPointSet( movingPointSetReader->GetOutput() );
  pointSetMetric->SetComputeSquaredDistance( false );
  pointSetMetric->SetTransform( transform );

  pointSetMetric->Initialize();
  typename PointSetMetricType::TransformParametersType parameters;
  parameters.Fill( 0 );

  typename PointSetMetricType::MeasureType measureFixed
    = pointSetMetric->GetValue( parameters );

  RealType measure = 0.0;
  for( unsigned int d = 0; d < measureFixed.GetSize(); d++ )
    {
    measure += measureFixed[d];
    }
  std::cout << measure / static_cast<float>( measureFixed.GetSize() ) << std::endl << std::endl;

  return EXIT_SUCCESS;
}


int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Usage: " << argv[0] << " Dimension fixedPointSet movingPointSet metric"
              << "   Metric options: " << std::endl
              << "     0: ICP" << std::endl
              << "     1: JHCT -> pointSetSigma kernelSigma evaluationKNeighborhood covarianceKNeighborhood "
              << "useAnisotropicCovariances alpha useInputAsSamples useRegularization "
              << "[numberofFixedSamples] [numberOfMovingSamples]"
              << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     {
     switch( atoi( argv[4] ) )
       {
       case 0:
         ICP<2>( argc, argv );
         break;
       case 1:
         JHCT<2>( argc, argv );
         break;
       default:
         std::cerr << "Unrecognized option." << std::endl;
         break;
       }
     break;
     }
   case 3:
     {
     switch( atoi( argv[4] ) )
       {
       case 0:
         ICP<3>( argc, argv );
         break;
       case 1:
         JHCT<3>( argc, argv );
         break;
       default:
         std::cerr << "Unrecognized option." << std::endl;
         break;
       }
     break;
     }
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }

}

