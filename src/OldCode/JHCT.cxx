#include "itkJensenHavrdaCharvatTsallisPointSetMetric.h"

#include "itkLabeledPointSetFileReader.h"
#include "itkPointSet.h"

#include <fstream.h>

#include "global.h"

int itkJensenHavrdaCharvatTsallisPointSetMetricTest(
  int argc, char *argv[] )
{
  const unsigned int Dimension = ImageDimension;

  typedef float RealType;
  typedef itk::PointSet<RealType, Dimension> PointSetType;

  typedef itk::LabeledPointSetFileReader<PointSetType> ReaderType;

  ReaderType::Pointer fixedPointSetReader = ReaderType::New();
  fixedPointSetReader->SetFileName( argv[2] );
  fixedPointSetReader->Update();

  ReaderType::Pointer movingPointSetReader = ReaderType::New();
  movingPointSetReader->SetFileName( argv[3] );
  movingPointSetReader->Update();

  typedef itk::JensenHavrdaCharvatTsallisPointSetMetric<PointSetType>
    PointSetMetricType;
  PointSetMetricType::Pointer pointSetMetric = PointSetMetricType::New();

  pointSetMetric->SetFixedPointSet( fixedPointSetReader->GetOutput() );
  pointSetMetric->SetFixedPointSetSigma( atof( argv[4] ) );
  pointSetMetric->SetFixedEvaluationKNeighborhood( atoi( argv[7] ) );

  pointSetMetric->SetMovingPointSet( movingPointSetReader->GetOutput() );
  pointSetMetric->SetMovingPointSetSigma( atof( argv[9] ) );
  pointSetMetric->SetMovingEvaluationKNeighborhood( atoi( argv[12] ) );

  pointSetMetric->SetUseRegularizationTerm( atoi( argv[17] ) );
  pointSetMetric->SetUseInputAsSamples( atoi( argv[16] ) );
  pointSetMetric->SetUseAnisotropicCovariances( atoi( argv[15] ) );
  pointSetMetric->SetAlpha( atof( argv[14] ) );

  if( pointSetMetric->GetUseAnisotropicCovariances() )
    {
    pointSetMetric->SetFixedCovarianceKNeighborhood( atoi( argv[8] ) );
    pointSetMetric->SetFixedKernelSigma( atof( argv[5] ) );
    pointSetMetric->SetMovingCovarianceKNeighborhood( atoi( argv[13] ) );
    pointSetMetric->SetMovingKernelSigma( atof( argv[10] ) );
    }

  if( !pointSetMetric->GetUseInputAsSamples() )
    {
    pointSetMetric->SetNumberOfFixedSamples( atoi( argv[6] ) );
    pointSetMetric->SetNumberOfMovingSamples( atoi( argv[11] ) );
    }

  try
    {
    pointSetMetric->Initialize();

    PointSetMetricType::DefaultTransformType::ParametersType parameters;
    parameters.Fill( 0 );

    pointSetMetric->SetUseWithRespectToTheMovingPointSet( false );
    PointSetMetricType::DerivativeType gradientFixed;
    pointSetMetric->GetDerivative( parameters, gradientFixed );
    
    {
    std::string source = std::string( argv[1] ) 
      + std::string( "FixedSource.txt" );
    std::string target = std::string( argv[1] ) 
      + std::string( "FixedTarget.txt" );

    ofstream sourceStr( source.c_str() );
    ofstream targetStr( target.c_str() );

    sourceStr << "0 0 0 0" << std::endl;
    targetStr << "0 0 0 0" << std::endl;

    for( unsigned int n = 0; n < 
         fixedPointSetReader->GetOutput()->GetNumberOfPoints(); n++ )
      {
      PointSetType::PointType point;
      fixedPointSetReader->GetOutput()->GetPoint( n, &point );

      sourceStr << point[0] << " " << point[1] << " ";
      targetStr << point[0] + gradientFixed(n, 0)/1e-10 << " "
                << point[1] + gradientFixed(n, 1)/1e-10 << " ";
      if ( Dimension == 2 )
        {
        sourceStr << "0 " << n+1 << std::endl; 
        targetStr << "0 " << n+1 << std::endl; 
        }           
      else
        {
        sourceStr << point[2] << " " << n+1 << std::endl; 
        targetStr << point[2] + gradientFixed(n, 2)/1e-10 << " " << n+1 << std::endl; 
        }  
      } 
    sourceStr << "0 0 0 0" << std::endl;
    targetStr << "0 0 0 0" << std::endl;
    
    sourceStr.close();
    targetStr.close();
    }
    PointSetMetricType::MeasureType measureFixed
      = pointSetMetric->GetValue( parameters );

    std::cout << "Fixed value: " << std::endl;
    std::cout << measureFixed << std::endl << std::endl;

    PointSetMetricType::MeasureType measureFixedTest;
    PointSetMetricType::DerivativeType gradientFixedTest;
    pointSetMetric->GetValueAndDerivative( parameters, 
      measureFixedTest, gradientFixedTest );
    if ( measureFixedTest != measureFixed )
      {
      std::cout << "Warning: fixed values from GetValue() and GetValueAndDerivative() "
                << "differ. " << measureFixed << " != " 
                << measureFixedTest << std::endl;
      }
    if ( gradientFixedTest != gradientFixed )
      {
      std::cout << "Warning: fixed values from GetDerivative() and GetValueAndDerivative() "
                << "differ. " << gradientFixed << " != " 
                << gradientFixedTest << std::endl;
      }
  
  
    /**
     * With respect to moving point set
     */

    pointSetMetric->SetUseWithRespectToTheMovingPointSet( true );
    PointSetMetricType::DerivativeType gradientMoving;
    pointSetMetric->GetDerivative( parameters, gradientMoving );

    {
    std::string source = std::string( argv[1] ) 
      + std::string( "MovingSource.txt" );
    std::string target = std::string( argv[1] ) 
      + std::string( "MovingTarget.txt" );

    ofstream sourceStr( source.c_str() );
    ofstream targetStr( target.c_str() );

    sourceStr << "0 0 0 0" << std::endl;
    targetStr << "0 0 0 0" << std::endl;

    for( unsigned int n = 0; n < 
         movingPointSetReader->GetOutput()->GetNumberOfPoints(); n++ )
      {
      PointSetType::PointType point;
      movingPointSetReader->GetOutput()->GetPoint( n, &point );

      sourceStr << point[0] << " " << point[1] << " ";
      targetStr << point[0] + gradientMoving(n, 0)/1e-10 << " "
                << point[1] + gradientMoving(n, 1)/1e-10 << " ";
      if ( Dimension == 2 )
        {
        sourceStr << "0 " << n+1 << std::endl; 
        targetStr << "0 " << n+1 << std::endl; 
        }           
      else
        {
        sourceStr << point[2] << " " << n+1 << std::endl; 
        targetStr << point[2] + gradientMoving(n, 2)/1e-10 << " " << n+1 << std::endl; 
        }  
      } 
    sourceStr << "0 0 0 0" << std::endl;
    targetStr << "0 0 0 0" << std::endl;
    
    sourceStr.close();
    targetStr.close();
    }

    PointSetMetricType::MeasureType measureMoving
      = pointSetMetric->GetValue( parameters );
    std::cout << "Moving value: " << std::endl;
    std::cout << measureMoving << std::endl;
    
    pointSetMetric->GetValueAndDerivative( parameters, 
      measureMoving, gradientMoving );

    PointSetMetricType::MeasureType measureMovingTest;
    PointSetMetricType::DerivativeType gradientMovingTest;
    pointSetMetric->GetValueAndDerivative( parameters, 
      measureMovingTest, gradientMovingTest );
    if ( measureMovingTest != measureMoving )
      {
      std::cout << "Warning: moving values from GetValue() and GetValueAndDerivative() "
                << "differ. " << measureMoving << " != " 
                << measureMovingTest << std::endl;
      }
    if ( gradientMovingTest != gradientMoving )
      {
      std::cout << "Warning: moving values from GetDerivative() and GetValueAndDerivative() "
                << "differ. " << gradientMoving << " != " 
                << gradientMovingTest << std::endl;
      }


    }
  catch(...)
    {
    std::cerr << "JensenHavrdaCharvatTsallisPointSetMetricTest: "
              << " Exception thrown." << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 18 )
    {
    std::cerr << "Usage: " << argv[0] << " outputPrefix fixedPointSet movingPointSet "
               << "fixedPointSetSigma fixedKernelSigma numberOfFixedSamples "
               << "fixedEvaluationKNeighborhood fixedCovarianceKNeighborhood "
               << "movingPointSetSigma movingKernelSigma numberOfMovingSamples "
               << "movingEvaluationKNeighborhood movingCovarianceKNeighborhood "
               << "alpha useAnisotropicCovariances "
               << "useInputAsSamples useRegularization" << std::endl;
    return EXIT_FAILURE;
    }
  return itkJensenHavrdaCharvatTsallisPointSetMetricTest(
    argc, argv );

}

