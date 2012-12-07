#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryThinning3DImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkPointSet.h"

#include <fstream>

#include "vnl/vnl_cross.h"
#include "vcl_cmath.h"

int main( int argc, char *argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << argv[0] << " inputAvantsSampledCurveFile.txt" << std::endl;

    return 1;
    }

  const unsigned int ImageDimension = 3;

  typedef unsigned int LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;
  typedef LabelImageType::PointType PointType;
  typedef PointType::VectorType     VectorType;

  typedef float RealType;

  std::vector<PointType> points;
  std::vector<LabelType> labels;
  std::vector<LabelType> distinctLabels;

  RealType x[3];
  LabelType l;

  std::fstream str( argv[1] );

  while ( str >> x[1] >> x[0] >> x[2] >> l )
    {
    if ( x[0] == 0 && x[1] == 0 && x[2] == 0 && l == 0 )
      {
      continue;
      }

    PointType point;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      point[d] = x[d];
      }
    points.push_back( point );

    if( std::find( distinctLabels.begin(), distinctLabels.end(), l ) == distinctLabels.end() )
      {
      distinctLabels.push_back( l );
      }
    labels.push_back( l );
    }
  str.close();

  // Calculate the DM metric
  float dm = 0.0;

  for( unsigned int n = 0; n < points.size()-1; n++ )
    {
    dm += ( points[n] - points[n+1] ).GetNorm();
    }


  // Calculate the ICM metric

  std::vector<VectorType> Ns;

  for( unsigned int n = 1; n < points.size()-1; n++ )
    {
    VectorType T1 = points[n] - points[n-1];
    VectorType T2 = points[n+1] - points[n];
    VectorType V = points[n+1] - points[n-1];
    VectorType A = T2 - T1;

    VectorType T = V / V.GetNorm();

    VectorType VxA;
    VxA.SetVnlVector( vnl_cross_3d( V.GetVnlVector(), A.GetVnlVector() ) );

    VectorType N;
    N.SetVnlVector( vnl_cross_3d( VxA.GetVnlVector(), V.GetVnlVector() ) );
    N /= N.GetNorm();

    VectorType B;
    B.SetVnlVector( vnl_cross_3d( T.GetVnlVector(), N.GetVnlVector() ) );

    Ns.push_back( N );
    }

  std::vector<RealType> deltaNMag;
  for( unsigned int i = 1; i < Ns.size(); ++i )
    {
    VectorType deltaN = Ns[i] - Ns[i-1];
    deltaNMag.push_back( deltaN.GetSquaredNorm() );
    }

  unsigned int icmCount = 1;
  for( unsigned int i = 1; i < deltaNMag.size()-1; ++i )
    {
    if( deltaNMag[i] > 1.0 && deltaNMag[i] > deltaNMag[i-1] && deltaNMag[i] > deltaNMag[i-1] )
      {
      icmCount++;
      }
    }

  // Calculate the SOAM metric
  float CP = 0.0;

  for( unsigned int n = 1; n < points.size()-2; n++ )
    {
    VectorType T1 = points[n] - points[n-1];
    VectorType T2 = points[n+1] - points[n];
    VectorType T3 = points[n+2] - points[n+1];

    VectorType T1xT2;
    T1xT2.SetVnlVector( vnl_cross_3d( T1.GetVnlVector(), T2.GetVnlVector() ) );
    VectorType T2xT3;
    T2xT3.SetVnlVector( vnl_cross_3d( T2.GetVnlVector(), T3.GetVnlVector() ) );

    VectorType T1norm = T1 / T1.GetNorm();
    VectorType T2norm = T2 / T2.GetNorm();
    VectorType T3norm = T3 / T3.GetNorm();
    VectorType T1xT2norm = T1xT2 / T1xT2.GetNorm();
    VectorType T2xT3norm = T2xT3 / T2xT3.GetNorm();

    float IP = vcl_acos( dot_product( T1norm.GetVnlVector(), T2norm.GetVnlVector() ) );
    float TP = vcl_acos( dot_product( T1xT2norm.GetVnlVector(), T2xT3norm.GetVnlVector() ) );

    CP += vcl_sqrt( vnl_math_sqr( IP ) + vnl_math_sqr( TP ) );
    }

  RealType soam = CP / dm;
  dm /= ( points[0] - points[points.size()-1] ).GetNorm();

  std::cout << "DistanceMetric,InflectionCountMetric,SumOfAnglesMetric" << std::endl;
  std::cout << dm << "," << icmCount << "," << soam << std::endl;

  return 0;
}

