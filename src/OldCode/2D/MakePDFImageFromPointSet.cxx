#include "itkGaussianProbabilityDensityFunction.h"
#include "itkManifoldParzenWindowsPointSetFunction.h"
#include "itkPointSet.h"
#include "itkVariableSizeMatrix.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include <stdio.h>
#include <fstream.h>

#include "global.h"

int main( int argc, char *argv[] )        
{
  if ( argc < 7 )
    {
    std::cout << argv[0] << " inputPointSet initializationSigma "
              << "regularizationSigma numberOfNeighbors domainImage outputImage" 
              << std::endl;
    exit( 0 );
    } 

  typedef double RealType; 

  typedef itk::VariableSizeMatrix<RealType> CovarianceMatrixType; 
  typedef itk::PointSet<CovarianceMatrixType, ImageDimension> PointSetType;
  typedef PointSetType::PointType PointType;

  RealType x, y, z, l;

  PointSetType::Pointer points = PointSetType::New();
  points->Initialize();

  RealType initSigma = atof( argv[2] );

  ifstream strF( argv[1] );
  unsigned int count = 0;
  while ( strF >> x >> y >> z >> l )
    {
    if ( x == 0 && y == 0 && z == 0 && l == 0 )
      {
      continue; 
      }
    PointType point;
    point[0] = x;
    point[1] = y;
    if ( ImageDimension == 3 )
      {
      point[2] = z;
      } 
  
    CovarianceMatrixType C( ImageDimension, ImageDimension );
    C.SetIdentity();
    C *= initSigma;
    
    points->SetPoint( count, point );    
    points->SetPointData( count, C );
    count++;      
    }
  strF.close();


  typedef itk::ManifoldParzenWindowsPointSetFunction<PointSetType> ParzenFilterType;
  ParzenFilterType::Pointer parzen = ParzenFilterType::New();
  parzen->SetBucketSize( 4 );
  parzen->SetSigma( atof( argv[3] ) );
  parzen->SetKNeighborhood( atoi( argv[4] ) );
  parzen->SetNormalize( true );
  parzen->SetInputPointSet( points );




  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[5] );
  reader->Update();

  itk::ImageRegionIteratorWithIndex<RealImageType> It( reader->GetOutput(),
    reader->GetOutput()->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    RealImageType::PointType point;
    reader->GetOutput()->TransformIndexToPhysicalPoint( It.GetIndex(), point );
    PointType pt;
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      pt[d] = point[d]; 
      } 
    It.Set( parzen->Evaluate( pt ) );
    }   

  typedef itk::ImageFileWriter<RealImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[6] );
  writer->SetInput( reader->GetOutput() );
  writer->Update();


  return EXIT_SUCCESS;
}     

