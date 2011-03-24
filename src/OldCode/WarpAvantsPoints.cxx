#include "itkVectorImageFileReader.h"
#include "itkVectorLinearInterpolateImageFunction.h"

#include <stdio.h>
#include <vector>
#include <fstream.h>

#include "vnl/vnl_vector.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc != 4 )
    {
    std::cout << "Usage: " << argv[0] << "  inputLandmarkFile deformationField outputLandmarkFile " << std::endl;
    exit( 1 );
    }

  typedef float RealType;

  RealType x;
  RealType y;
  RealType z;
  int l;

  std::vector<RealType> X;
  std::vector<RealType> Y;
  std::vector<RealType> Z;
  std::vector<int> L;
  ifstream str( argv[1] );

  while ( str >> x >> y >> z >> l )
    {
    X.push_back( x );
    Y.push_back( y );
    Z.push_back( z );
    L.push_back( l );

    }
     
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;
  
  typedef itk::VectorImageFileReader<RealImageType, 
    DeformationFieldType> FieldReaderType;
  FieldReaderType::Pointer fieldreader = FieldReaderType::New();
  fieldreader->SetFileName( argv[2] );
  fieldreader->SetUseAvantsNamingConvention( true );
  fieldreader->Update();
    
  typedef itk::VectorLinearInterpolateImageFunction<DeformationFieldType, RealType> 
    DeformationFieldInterpolatorType;
  DeformationFieldInterpolatorType::Pointer interpolator = DeformationFieldInterpolatorType::New();
  interpolator->SetInputImage( fieldreader->GetOutput() );
    
  vnl_vector<RealType> Xt( X.size(), 0.0 );
  vnl_vector<RealType> Yt( Y.size(), 0.0 );
  vnl_vector<RealType> Zt( Z.size(), 0.0 );

  for ( unsigned int i = 1; i < X.size()-1; i++ )
    {
    DeformationFieldInterpolatorType::PointType point;
    point[0] = X[i]; 
    point[1] = Y[i]; 
    point[2] = Z[i]; 

    bool isOutside = false;
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      if ( point[d] < 0.0 || point[d] > fieldreader->GetOutput()->GetOrigin()[d] + 
             ( fieldreader->GetOutput()->GetLargestPossibleRegion().GetSize()[d]-1 )
               *fieldreader->GetOutput()->GetSpacing()[d] )
        {
        isOutside = true;
        break;
        }              
      }

    DeformationFieldInterpolatorType::OutputType disp;  
    disp.Fill( 0.0 );
    if ( !isOutside )
      {
      disp = interpolator->Evaluate( point );
      }  
    Xt[i] = X[i] + disp[0];  
    Yt[i] = Y[i] + disp[1];  
    Zt[i] = Z[i] + disp[2];  
    }
  
  ofstream strt( argv[3] );

  strt << "0 0 0 0" << std::endl;
  for ( unsigned int i = 1; i < Xt.size()-1; i++ )
    {
    strt << Xt[i] << " " << Yt[i] << " " << Zt[i] << " " << L[i] << std::endl;
    }
  strt << "0 0 0 0" << std::endl;

  return 0;
}
