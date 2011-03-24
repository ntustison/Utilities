#include "itkGaussianProbabilityDensityFunction.h"
#include "itkManifoldParzenWindowsPointSetFunction.h"
#include "itkPointSet.h"
#include "itkVariableSizeMatrix.h"

#include "itkLandmarkFileReader.h"
#include "itkLandmarkFileWriter.h"
#include "itkNumericTraits.h"

#include <stdio.h>
#include <fstream.h>

#include "global.h"

int main( int argc, char *argv[] )        
{
 
  std::cout << itk::NumericTraits<double>::min() << std::endl;
  std::cout << itk::NumericTraits<float>::min() << std::endl;
 
  typedef double RealType; 

  typedef itk::PointSet<long, ImageDimension> PointSetType;

  typedef itk::LandmarkFileReader<PointSetType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  
  for ( unsigned int i = 0; i < reader->GetOutput()->GetNumberOfPoints(); i++ )
    {
    PointSetType::PointType point;
    reader->GetOutput()->GetPoint( i, &point );
    std::cout << i << ": " << point << std::endl;
    }

  typedef itk::LandmarkFileWriter<PointSetType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( reader->GetOutput() );
  writer->Update();
  
//  if ( argc < 6 )
//    {
//    std::cout << argv[0] << " inputPointSet outputPointSet numberOfOutputPoints initializationSigma regularizationSigma numberOfNeighbors" << std::endl;
//    exit( 0 );
//    } 
//
//  typedef double RealType; 
//
//  typedef itk::VariableSizeMatrix<RealType> CovarianceMatrixType; 
//  typedef itk::PointSet<CovarianceMatrixType, ImageDimension> PointSetType;
//
//  RealType x, y, z, l;
//
//  PointSetType::Pointer points = PointSetType::New();
//  points->Initialize();
//
//  RealType initSigma = atof( argv[4] );
//
//  ifstream strF( argv[1] );
//  unsigned int count = 0;
//  while ( strF >> x >> y >> z >> l )
//    {
//    if ( x == 0 && y == 0 && z == 0 && l == 0 )
//      {
//      continue; 
//      }
//    PointSetType::PointType point;
//    point[0] = x;
//    point[1] = y;
//    if ( ImageDimension == 3 )
//      {
//      point[2] = z;
//      } 
//  
//    CovarianceMatrixType C( ImageDimension, ImageDimension );
//    C.SetIdentity();
//    C *= initSigma;
//    
//    points->SetPoint( count, point );    
//    points->SetPointData( count, C );
//    count++;      
//    }
//  strF.close();
//
//
//  typedef itk::ManifoldParzenWindowsPointSetFunction<PointSetType> ParzenFilterType;
//  ParzenFilterType::Pointer parzen = ParzenFilterType::New();
//  parzen->SetBucketSize( 4 );
//  parzen->SetSigma( atof( argv[5] ) );
//  parzen->SetKNeighborhood( atoi( argv[6] ) );
//  parzen->SetNormalize( false );
//  parzen->SetInputPointSet( points );
//
//  ofstream str( argv[2] );
//  str << "0 0 0 0" << std::endl;
//
//  for ( int n = 0; n < atoi( argv[3] ); n++ )
//    {
//    ParzenFilterType::PointType point = parzen->GenerateRandomSample();
//    str << point[0] << " " << point[1] << " ";
//    if ( ImageDimension == 3 )
//      {
//      str << point[2] << " ";
//      }
//    str << n+1 << std::endl;
//    } 
//  str << "0 0 0 0" << std::endl;
//  str.close();

  return EXIT_SUCCESS;
}     

