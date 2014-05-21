#include "itkJensenHavrdaCharvatTsallisLabeledPointSetMetric.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImageFileReader.h"
#include "itkLabeledPointSetFileReader.h"
#include "itkPointSet.h"
#include "itkVector.h"

#include <fstream>

#include "Common.h"

template <unsigned int ImageDimension>
int JHCT( int argc, char *argv[] )
{

  const unsigned int Dimension = ImageDimension;

  typedef float RealType;

  typedef long PointDataType;
  typedef itk::PointSet<PointDataType, Dimension> PointSetType;

  typedef itk::LabeledPointSetFileReader<PointSetType> ReaderType;

  typename ReaderType::Pointer fixedPointSetReader = ReaderType::New();
  fixedPointSetReader->SetFileName( argv[3] );
  fixedPointSetReader->Update();

  typename ReaderType::Pointer movingPointSetReader = ReaderType::New();
  movingPointSetReader->SetFileName( argv[4] );
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
    std::cout << "HERE" << std::endl;
    pointSetMetric->SetNumberOfFixedSamples( atoi( argv[13] ) );
    pointSetMetric->SetNumberOfMovingSamples( atoi( argv[14] ) );
    }


  pointSetMetric->Initialize();

  typename PointSetMetricType::DefaultTransformType::ParametersType parameters;
  parameters.Fill( 0 );

  {
  pointSetMetric->SetUseWithRespectToTheMovingPointSet( false );
  typename PointSetMetricType::DerivativeType gradientFixed;
  pointSetMetric->GetDerivative( parameters, gradientFixed );

  typedef itk::Image<RealType, Dimension> ImageType;
  typedef itk::Vector<RealType, Dimension> VectorType;
  typedef itk::Image<VectorType, Dimension> DeformationFieldType;
  typedef itk::PointSet<VectorType, Dimension> BSplinePointSetType;
  typedef itk::BSplineScatteredDataPointSetToImageFilter
    <BSplinePointSetType, DeformationFieldType> BSplineFilterType;

  typename BSplinePointSetType::Pointer fieldPoints = BSplinePointSetType::New();
  fieldPoints->Initialize();

  for( unsigned int n = 0; n < pointSetMetric->GetNumberOfValues(); n++ )
    {
    typename BSplinePointSetType::PointType fieldPoint;
    VectorType gradient;

    typename PointSetType::PointType point;
    fixedPointSetReader->GetOutput()->GetPoint( n, &point );

    PointDataType data = 1;
    fixedPointSetReader->GetOutput()->GetPointData( n, &data );

    for( unsigned d = 0; d < Dimension; d++ )
      {
      gradient[d] = gradientFixed(n, d);
      fieldPoint[d] = point[d];
      }
    fieldPoints->SetPoint( n, fieldPoint );
    fieldPoints->SetPointData( n, gradient );
    }

  RealType sumSquaredNorm = 0.0;
  for( unsigned int i = 0; i < fieldPoints->GetNumberOfPoints(); i++ )
    {
    VectorType gradient;
    fieldPoints->GetPointData( i, &gradient );
    sumSquaredNorm += gradient.GetSquaredNorm();
    }
  VectorType V;
  RealType sigma;
  for( unsigned int i = 0; i < Dimension; i++ )
    {
    V[i] = 1.0;
    }
  if( Dimension == 2 )
    {
    sigma = V.GetNorm();
    }
  else if( Dimension == 3 )
    {
    sigma = V.GetNorm()/vcl_sqrt( 2.0 );
    }
  RealType gradientScalingFactor = sigma*vcl_sqrt
    ( static_cast<RealType>( Dimension *
    fieldPoints->GetNumberOfPoints() ) / sumSquaredNorm );

//   std::cout << "Fixed gradient scaling factor: " << gradientScalingFactor
//     << std::endl;

  for( unsigned int i = 0; i < fieldPoints->GetNumberOfPoints(); i++ )
    {
    VectorType gradient;
    fieldPoints->GetPointData( i, &gradient );
    fieldPoints->SetPointData( i, gradient*gradientScalingFactor );
    }

  std::string source = std::string( argv[2] )
    + std::string( "FixedSource.txt" );
  std::string target = std::string( argv[2] )
    + std::string( "FixedTarget.txt" );

  std::ofstream sourceStr( source.c_str() );
  std::ofstream targetStr( target.c_str() );

  sourceStr << "0 0 0 0" << std::endl;
  targetStr << "0 0 0 0" << std::endl;

  if( argc < 16 )
    {
    for( unsigned int n = 0; n < fieldPoints->GetNumberOfPoints(); n++ )
      {
      typename PointSetType::PointType point;
      point.Fill( 0.0 );
      fieldPoints->GetPoint( n, &point );

      PointDataType data = 1;
      fixedPointSetReader->GetOutput()->GetPointData( n, &data );
      VectorType gradient;
      gradient.Fill( 0.0 );
      fieldPoints->GetPointData( n, &gradient );

      sourceStr << point[0] << " " << point[1] << " ";
      targetStr << point[0] + gradient[0] << " "
                << point[1] + gradient[1] << " ";
      if ( Dimension == 2 )
        {
        sourceStr << "0 " << data << std::endl;
        targetStr << "0 " << data << std::endl;
        }
      else
        {
        sourceStr << point[2] << " " << data << std::endl;
        targetStr << point[2] + gradient[2] << " " << data << std::endl;
        }
      }
    }
//   else
//     {
//     typedef itk::ImageFileReader<ImageType> ReaderType;
//     typename ReaderType::Pointer reader = ReaderType::New();
//     reader->SetFileName( argv[17] );
//     reader->Update();
//
//     int expansionFactor = 0;
//     if( argc > 17 )
//       {
//       expansionFactor = atoi( argv[18] );
//       }
//     typename ImageType::PointType origin = reader->GetOutput()->GetOrigin();
//     typename ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();
//     typename ImageType::SizeType size
//       = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
//     for( unsigned int d = 0; d < Dimension; d++ )
//       {
//       origin[d] -= spacing[d]*static_cast<RealType>( expansionFactor );
//       size[d] += 2*size[d]*expansionFactor;
//       }
//
//     std::vector<unsigned int> ncps
//       = ConvertVector<unsigned int>( std::string( argv[15] ) );
//     typename BSplineFilterType::ArrayType numberOfControlPoints;
//     for( unsigned int d = 0; d < Dimension; d++ )
//       {
//       numberOfControlPoints[d] = ncps[d];
//       }
//
//     typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
//     bspliner->SetInput( fieldPoints );
//     bspliner->SetOrigin( origin );
//     bspliner->SetSpacing( spacing );
//     bspliner->SetSize( size );
//     bspliner->SetNumberOfLevels( 1 );
//     bspliner->SetSplineOrder( atoi( argv[16] ) );
//     bspliner->SetNumberOfControlPoints( numberOfControlPoints );
//     bspliner->SetGenerateOutputImage( false );
//     bspliner->Update();
//
//     for( unsigned int n = 0; n < fieldPoints->GetNumberOfPoints(); n++ )
//       {
//       typename PointSetType::PointType point;
//       point.Fill( 0.0 );
//       fieldPoints->GetPoint( n, &point );
//
//       PointDataType data = 1;
//       fixedPointSetReader->GetOutput()->GetPointData( n, &data );
//       VectorType gradient;
//       gradient.Fill( 0.0 );
//       bspliner->EvaluateAtPoint( point, gradient );
//
//       sourceStr << point[0] << " " << point[1] << " ";
//       targetStr << point[0] + gradient[0] << " "
//                 << point[1] + gradient[1] << " ";
//       if ( Dimension == 2 )
//         {
//         sourceStr << "0 " << data << std::endl;
//         targetStr << "0 " << data << std::endl;
//         }
//       else
//         {
//         sourceStr << point[2] << " " << data << std::endl;
//         targetStr << point[2] + gradient[2] << " " << data << std::endl;
//         }
//       }
//     }
  sourceStr << "0 0 0 0" << std::endl;
  targetStr << "0 0 0 0" << std::endl;
  sourceStr.close();
  targetStr.close();

  typename PointSetMetricType::MeasureType measureFixed
    = pointSetMetric->GetValue( parameters );

  std::cout << "Fixed value: " << std::endl;
  std::cout << measureFixed << std::endl << std::endl;

  typename PointSetMetricType::MeasureType measureFixedTest;
  typename PointSetMetricType::DerivativeType gradientFixedTest;
  pointSetMetric->GetValueAndDerivative( parameters,
    measureFixedTest, gradientFixedTest );
  if ( measureFixedTest != measureFixed )
    {
    std::cout << "Warning: fixed values from GetValue() and GetValueAndDerivative() "
              << "differ. "
              << measureFixed << " != " << measureFixedTest << std::endl;
    }
  if ( gradientFixedTest != gradientFixed )
    {
    std::cout << "Warning: fixed values from GetDerivative() and GetValueAndDerivative() "
              << "differ. "
//                << gradientFixed << " != " << gradientFixedTest
              << std::endl;
    }
  }

  /**
    * With respect to moving point set
    */
  {
  pointSetMetric->SetUseWithRespectToTheMovingPointSet( true );
  typename PointSetMetricType::DerivativeType gradientMoving;
  pointSetMetric->GetDerivative( parameters, gradientMoving );

  typedef itk::Image<RealType, Dimension> ImageType;
  typedef itk::Vector<RealType, Dimension> VectorType;
  typedef itk::Image<VectorType, Dimension> DeformationFieldType;
  typedef itk::PointSet<VectorType, Dimension> BSplinePointSetType;
  typedef itk::BSplineScatteredDataPointSetToImageFilter
    <BSplinePointSetType, DeformationFieldType> BSplineFilterType;

  typename BSplinePointSetType::Pointer fieldPoints = BSplinePointSetType::New();
  fieldPoints->Initialize();

  for( unsigned int n = 0; n < pointSetMetric->GetNumberOfValues(); n++ )
    {
    typename BSplinePointSetType::PointType fieldPoint;
    VectorType gradient;

    typename PointSetType::PointType point;
    movingPointSetReader->GetOutput()->GetPoint( n, &point );

    PointDataType data = 1;
    movingPointSetReader->GetOutput()->GetPointData( n, &data );

    for( unsigned d = 0; d < Dimension; d++ )
      {
      gradient[d] = gradientMoving(n, d);
      fieldPoint[d] = point[d];
      }
    fieldPoints->SetPoint( n, fieldPoint );
    fieldPoints->SetPointData( n, gradient );
    }

  RealType sumSquaredNorm = 0.0;
  for( unsigned int i = 0; i < fieldPoints->GetNumberOfPoints(); i++ )
    {
    VectorType gradient;
    fieldPoints->GetPointData( i, &gradient );
    sumSquaredNorm += gradient.GetSquaredNorm();
    }
  VectorType V;
  RealType sigma;
  for( unsigned int i = 0; i < Dimension; i++ )
    {
    V[i] = 1.0;
    }
  if( Dimension == 2 )
    {
    sigma = V.GetNorm();
    }
  else if( Dimension == 3 )
    {
    sigma = V.GetNorm()/vcl_sqrt( 2.0 );
    }
  RealType gradientScalingFactor = sigma*vcl_sqrt
    ( static_cast<RealType>( Dimension *
    fieldPoints->GetNumberOfPoints() ) / sumSquaredNorm );

//   std::cout << "Moving gradient scaling factor: " << gradientScalingFactor
//     << std::endl;

  for( unsigned int i = 0; i < fieldPoints->GetNumberOfPoints(); i++ )
    {
    VectorType gradient;
    fieldPoints->GetPointData( i, &gradient );
    fieldPoints->SetPointData( i, gradient*gradientScalingFactor );
    }

  std::string source = std::string( argv[2] )
    + std::string( "MovingSource.txt" );
  std::string target = std::string( argv[2] )
    + std::string( "MovingTarget.txt" );

  std::ofstream sourceStr( source.c_str() );
  std::ofstream targetStr( target.c_str() );

  sourceStr << "0 0 0 0" << std::endl;
  targetStr << "0 0 0 0" << std::endl;

  if( argc < 16 )
    {
    for( unsigned int n = 0; n < fieldPoints->GetNumberOfPoints(); n++ )
      {
      typename PointSetType::PointType point;
      point.Fill( 0.0 );
      fieldPoints->GetPoint( n, &point );

      PointDataType data = 1;
      movingPointSetReader->GetOutput()->GetPointData( n, &data );
      VectorType gradient;
      gradient.Fill( 0.0 );
      fieldPoints->GetPointData( n, &gradient );

      sourceStr << point[0] << " " << point[1] << " ";
      targetStr << point[0] + gradient[0] << " "
                << point[1] + gradient[1] << " ";
      if ( Dimension == 2 )
        {
        sourceStr << "0 " << data << std::endl;
        targetStr << "0 " << data << std::endl;
        }
      else
        {
        sourceStr << point[2] << " " << data << std::endl;
        targetStr << point[2] + gradient[2] << " " << data << std::endl;
        }
      }
    }
//   else
//     {
//     typedef itk::ImageFileReader<ImageType> ReaderType;
//     typename ReaderType::Pointer reader = ReaderType::New();
//     reader->SetFileName( argv[17] );
//     reader->Update();
//
//     int expansionFactor = 0;
//     if( argc > 17 )
//       {
//       expansionFactor = atoi( argv[18] );
//       }
//     typename ImageType::PointType origin = reader->GetOutput()->GetOrigin();
//     typename ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();
//     typename ImageType::SizeType size
//       = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
//     for( unsigned int d = 0; d < Dimension; d++ )
//       {
//       origin[d] -= spacing[d]*static_cast<RealType>( expansionFactor );
//       size[d] += 2*expansionFactor;
//       }
//
//     std::vector<unsigned int> ncps
//       = ConvertVector<unsigned int>( std::string( argv[15] ) );
//     typename BSplineFilterType::ArrayType numberOfControlPoints;
//     for( unsigned int d = 0; d < Dimension; d++ )
//       {
//       numberOfControlPoints[d] = ncps[d];
//       }
//
//     typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
//     bspliner->SetInput( fieldPoints );
//     bspliner->SetOrigin( origin );
//     bspliner->SetSpacing( spacing );
//     bspliner->SetSize( size );
//     bspliner->SetNumberOfLevels( 1 );
//     bspliner->SetSplineOrder( atoi( argv[16] ) );
//     bspliner->SetNumberOfControlPoints( numberOfControlPoints );
//     bspliner->SetGenerateOutputImage( false );
//     bspliner->Update();
//
//     for( unsigned int n = 0; n < fieldPoints->GetNumberOfPoints(); n++ )
//       {
//       typename PointSetType::PointType point;
//       point.Fill( 0.0 );
//
//       fieldPoints->GetPoint( n, &point );
//
//       PointDataType data = 1;
//       movingPointSetReader->GetOutput()->GetPointData( n, &data );
//       VectorType gradient;
//       gradient.Fill( 0.0 );
//       bspliner->EvaluateAtPoint( point, gradient );
//
//       sourceStr << point[0] << " " << point[1] << " ";
//       targetStr << point[0] + gradient[0] << " "
//                 << point[1] + gradient[1] << " ";
//       if ( Dimension == 2 )
//         {
//         sourceStr << "0 " << data << std::endl;
//         targetStr << "0 " << data << std::endl;
//         }
//       else
//         {
//         sourceStr << point[2] << " " << data << std::endl;
//         targetStr << point[2] + gradient[2] << " " << data << std::endl;
//         }
//       }
//     }
  sourceStr << "0 0 0 0" << std::endl;
  targetStr << "0 0 0 0" << std::endl;

  sourceStr.close();
  targetStr.close();

  typename PointSetMetricType::MeasureType measureMoving
    = pointSetMetric->GetValue( parameters );
  std::cout << "Moving value: " << std::endl;
  std::cout << measureMoving << std::endl;

  pointSetMetric->GetValueAndDerivative( parameters,
    measureMoving, gradientMoving );

  typename PointSetMetricType::MeasureType measureMovingTest;
  typename PointSetMetricType::DerivativeType gradientMovingTest;
  pointSetMetric->GetValueAndDerivative( parameters,
    measureMovingTest, gradientMovingTest );
//   if ( measureMovingTest != measureMoving )
//     {
//     std::cout << "Warning: moving values from GetValue() and GetValueAndDerivative() "
//               << "differ. " << measureMoving << " != "
//               << measureMovingTest << std::endl;
//     }
//   if ( gradientMovingTest != gradientMoving )
//     {
//     std::cout << "Warning: moving values from GetDerivative() and GetValueAndDerivative() "
//               << "differ. "
// //                << gradientMoving << " != " << gradientMovingTest
//               << std::endl;
//     }
  }


  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 13 )
    {
    std::cerr << "Usage: " << argv[0] << " Dimension outputPrefix fixedPointSet movingPointSet "
               << "pointSetSigma kernelSigma evaluationKNeighborhood covarianceKNeighborhood "
               << "useAnisotropicCovariances alpha useInputAsSamples useRegularization "
               << "[numberofFixedSamples] [numberOfMovingSamples] "
//                << " [ncps] [order] [domainImage] [expansionFactor]"
               << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     JHCT<2>( argc, argv );
     break;
   case 3:
     JHCT<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }

}

