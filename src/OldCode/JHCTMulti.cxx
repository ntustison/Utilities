#include "itkJensenHavrdaCharvatTsallisMultiplePointSetMetric.h"

#include "itkLabeledPointSetFileReader.h"
#include "itkPointSet.h"

#include <fstream.h>

#include "global.h"

int itkJensenHavrdaCharvatTsallisMultiplePointSetMetricTest(
  unsigned int argc, char *argv[] )
{
  if( argc % 6 != 0 )
    {
    std::cerr << "Incorrect number of parameters." << std::endl;
    exit( 0 ); 
    }

  const unsigned int Dimension = ImageDimension;

  typedef float RealType;
  typedef itk::PointSet<RealType, Dimension> PointSetType;

  typedef itk::JensenHavrdaCharvatTsallisMultiplePointSetMetric<PointSetType>
    PointSetMetricType;
  PointSetMetricType::Pointer pointSetMetric = PointSetMetricType::New();

  pointSetMetric->SetAlpha( atof( argv[2] ) );
  pointSetMetric->SetUseAnisotropicCovariances( atoi( argv[3] ) );
  pointSetMetric->SetUseInputAsSamples( atoi( argv[4] ) );
  pointSetMetric->SetUseRegularizationTerm( atoi( argv[5] ) );

  for ( unsigned int n = 0; n < static_cast<RealType>( argc ) / 6.0 - 1.0; n++ )
    {
    typedef itk::LabeledPointSetFileReader<PointSetType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[6+n*6+0] );
    reader->Update();
    
    pointSetMetric->SetInput( n, reader->GetOutput() );
    pointSetMetric->SetPointSetSigma( n, atof( argv[6+n*6+1] ) );
    pointSetMetric->SetKernelSigma( n, atof( argv[6+n*6+2] ) );
    pointSetMetric->SetNumberOfSamples( n, atoi( argv[6+n*6+3] ) );
    pointSetMetric->SetEvaluationKNeighborhood( n, atoi( argv[6+n*6+4] ) );
    pointSetMetric->SetCovarianceKNeighborhood( n, atoi( argv[6+n*6+5] ) );
    }  

  std::cout << pointSetMetric << std::endl;

  pointSetMetric->Initialize();
  
  PointSetMetricType::MeasureType measure;
  PointSetMetricType::DerivativeType gradient;
  pointSetMetric->GetValueAndDerivative( measure, gradient );

  RealType avgNorm = 0;
  for ( unsigned int n = 0; n < pointSetMetric->GetNumberOfValues(); n++ )
    {
    RealType norm = 0.0;
    for ( unsigned int d = 0; d < Dimension; d++ )
      {
      norm += ( gradient(n, d)*gradient(n, d) ); 
      }
    avgNorm += vcl_sqrt( norm );
    }
  avgNorm /= static_cast<RealType>( pointSetMetric->GetNumberOfValues() );  

  for ( unsigned int n = 0; n < static_cast<RealType>( argc ) / 6.0 - 1.0; n++ )
    {
    itk::OStringStream buf;
    buf << n;

    std::string source = std::string( argv[1] ) + std::string( "_" ) + 
      buf.str() + std::string( "_source.txt" );
    std::string target = std::string( argv[1] ) + std::string( "_" ) + 
      buf.str() + std::string( "_target.txt" );

    ofstream sourceStr( source.c_str() );
    ofstream targetStr( target.c_str() );

    sourceStr << "0 0 0 0" << std::endl;
    targetStr << "0 0 0 0" << std::endl;

    unsigned int index = 0;
    for ( unsigned int m = 0; m < n; m++ )
      {
      index += pointSetMetric->GetInput( m )->GetNumberOfPoints(); 
      }

    for( unsigned int i = 0; i < 
         pointSetMetric->GetInput( n )->GetNumberOfPoints(); i++ )
      {
      PointSetType::PointType point;
      pointSetMetric->GetInput( n )->GetPoint( i, &point );

      sourceStr << point[0] << " " << point[1] << " ";
      targetStr << point[0] + gradient(index+i, 0)/avgNorm << " "
                << point[1] + gradient(index+i, 1)/avgNorm << " ";
      if ( Dimension == 2 )
        {
        sourceStr << "0 " << n+1 << std::endl; 
        targetStr << "0 " << n+1 << std::endl; 
        }           
      else
        {
        sourceStr << point[2] << " " << n+1 << std::endl; 
        targetStr << point[2] + gradient(index+i, 2)/avgNorm << " " << n+1 << std::endl; 
        }  
      } 
    sourceStr << "0 0 0 0" << std::endl;
    targetStr << "0 0 0 0" << std::endl;
    
    sourceStr.close();
    targetStr.close();
    }
  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc == 1 )
    {
    std::cerr << "Usage: " << argv[0] << " outputPrefix alpha "
              << "useAnisotropicCovariances useInputAsSamples "
              << "useRegularizationTerm "
              << "pointSet[0] pointSetSigma[0] kernelSigma[0] numberOfSamples[0] "
              << "evaluationKNeighborhood[0] covarianceKNeighborhood[0] ..."
              << "pointSet[n] pointSetSigma[n] kernelSigma[n] numberOfSamples[n] "
              << "evaluationKNeighborhood[n] covarianceKNeighborhood[n] "
              << std::endl;
    return EXIT_FAILURE;
    }
  return itkJensenHavrdaCharvatTsallisMultiplePointSetMetricTest(
    argc, argv );

}

