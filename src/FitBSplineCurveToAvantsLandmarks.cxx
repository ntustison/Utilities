#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkPointSet.h"

#include <stdio.h>
#include <vector>
#include <fstream.h>
#include <string>

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " inputLandmarks outputLandmarksPrefix [order] [nlevels] "
      << "[ncps] [sampleSpacing] [closed]" << std::endl;
    std::cout << "  Note2:  1. Points are assumed to be parametrically ordered. " << std::endl
              << "          2. The fourth column is used for weights." << std::endl;
    exit( 1 );
    }

  typedef float RealType;

  typedef itk::Vector<RealType, 3> VectorType;
  typedef itk::Image<VectorType, 1> CurveImageType;

  typedef itk::PointSet<VectorType, 1> PointSetType;
  PointSetType::Pointer pointSet = PointSetType::New();
  pointSet->Initialize();

  typedef itk::BSplineScatteredDataPointSetToImageFilter
     <PointSetType, CurveImageType>  FilterType;
  FilterType::Pointer filter = FilterType::New();

  FilterType::WeightsContainerType::Pointer weights = FilterType::WeightsContainerType::New();

  RealType x, y, z, l;
  RealType totalDistance = 0.0;

  ifstream str( argv[1] );

  unsigned int count = 0;
  while ( str >> x >> y >> z >> l )
    {
    if ( x == 0 && y == 0 && z == 0 && l == 0 )
      {
      continue;
      }
    VectorType vector;
    vector[0] = x;
    vector[1] = y;
    vector[2] = z;
    pointSet->SetPointData( count, vector );

    if ( count > 0 )
      {
      VectorType previous;
      pointSet->GetPointData( count-1, &previous );
      totalDistance += ( previous - vector ).GetNorm();
      }

    PointSetType::PointType point;
    point[0] = 0.0;
    pointSet->SetPoint( count, point );

    weights->InsertElement( count, l );
    count++;
    }

  RealType cumSum = 0.0;
  for ( unsigned int i = 1; i < pointSet->GetNumberOfPoints(); i++ )
    {
    VectorType vector, previous;
    pointSet->GetPointData( i, &vector );
    pointSet->GetPointData( i-1, &previous );

    cumSum += ( vector - previous ).GetNorm();
    PointSetType::PointType point;
    point[0] = cumSum / totalDistance;

    pointSet->SetPoint( i, point );
    }

  filter->SetInput( pointSet );
  filter->SetGenerateOutputImage( true );

  CurveImageType::PointType origin;
  origin.Fill( 0.0 );
  filter->SetOrigin( origin );
  CurveImageType::SpacingType spacing;
  spacing[0] = 0.001;
  if ( argc > 6 )
    {
    spacing[0] = atof( argv[6] );
    }
  filter->SetSpacing( spacing );
  CurveImageType::SizeType size;
  size[0] = static_cast<unsigned int>( 1.0 / spacing[0] + 1 );
  filter->SetSize( size );
  FilterType::ArrayType order;
  order[0] = 3;
  if ( argc > 3 )
    {
    order[0] = atoi( argv[3] );
    }
  filter->SetSplineOrder( order );
  FilterType::ArrayType ncps;
  ncps[0] = order[0] + 1;
  if ( argc > 5 )
    {
    ncps[0] = atoi( argv[5] );
    }
  filter->SetNumberOfControlPoints( ncps );
  FilterType::ArrayType nlevels;
  nlevels[0] = 10;
  if ( argc > 4 )
    {
    nlevels[0] = atoi( argv[4] );
    }
  filter->SetNumberOfLevels( nlevels );
  FilterType::ArrayType close;
  close[0] = false;
  if ( argc > 7 )
    {
    close[0] = atoi( argv[7] );
    }
  filter->SetCloseDimension( close );

  filter->Update();

  unsigned int numberOfSpans = filter->GetPhiLattice()->GetLargestPossibleRegion().GetSize()[0] - order[0];
  RealType numberOfSamplesPerSpan = static_cast<RealType>( size[0] ) / static_cast<RealType>( numberOfSpans );

  {
  std::string filename = std::string( argv[2] ) + std::string( ".txt" );
  ofstream ostr( filename.c_str() );
  ostr << "0 0 0 0" << std::endl;

  itk::ImageRegionIterator<CurveImageType> It(
    filter->GetOutput(), filter->GetOutput()->GetLargestPossibleRegion() );
  RealType sample = 0;
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    unsigned int whichSpan = 1 + static_cast<unsigned int>( sample / numberOfSamplesPerSpan );
    ostr << It.Get()[0] << " " << It.Get()[1] << " " << It.Get()[2] << " " << whichSpan << std::endl;
    sample++;
    }
  ostr << "0 0 0 0" << std::endl;
  ostr.close();
  }

  {
  std::string filename = std::string( argv[2] ) + std::string( "_cps.txt" );
  ofstream ostr( filename.c_str() );
  ostr << "0 0 0 0" << std::endl;

  itk::ImageRegionIterator<CurveImageType> It(
    filter->GetPhiLattice(), filter->GetPhiLattice()->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    ostr << It.Get()[0] << " " << It.Get()[1] << " " << It.Get()[2] << " 1" << std::endl;
    }
  ostr << "0 0 0 0" << std::endl;
  ostr.close();
  }


//   typedef itk::BSplineControlPointImageFilter<CurveImageType> BSplinerType;
//   BSplinerType::Pointer bspliner = BSplinerType::New();
//   bspliner->SetInput( filter->GetPhiLattice() );
//   bspliner->SetSpacing( filter->GetSpacing() );
//   bspliner->SetSize( filter->GetSize() );
//   bspliner->SetOrigin( filter->GetOrigin() );
//   bspliner->SetDirection( filter->GetDirection() );
//   bspliner->SetSplineOrder( filter->GetSplineOrder() );
//
//   BSplinerType::PointType params;
//   params.Fill( 0.9 );
//
//   BSplinerType::PointDataType point;
//   point[0] = 1.8;
//   point[1] = 1.0;
//   point[2] = 0.0;
//
//   std::cout << params << std::endl;
//   bspliner->CalculateParametersClosestToDataPoint( point, params );
//   std::cout << params << std::endl;
//
//   BSplinerType::PixelType data;
//
//   bspliner->EvaluateAtPoint( params, data );
//   std::cout << "Closest point: " << data << std::endl;

  return 0;
}
