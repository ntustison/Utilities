#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkFloatArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"

#include <stdio>
#include <vector>
#include <fstream>
#include <string>

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " fixedLandmarkFile movingLandmarkFile outputFilePrefix "
              << "[format]" << std::endl;
    std::cout << "format:  " << std::endl;
    std::cout << "   0: x y z ( 1 to 1 correspondence ) " << std::endl;
    std::cout << "   1: x y z ( many to many correpondence ) " << std::endl;
    std::cout << "   2: x y z label ( 1 to 1 correspondence ) " << std::endl;
    std::cout << "   3: x y z label ( many to many correpondence ) " << std::endl;
    exit( 1 );
    }

  unsigned int format = 0;
  if ( argc == 5 )
    {
    format = atoi( argv[4] );
    }
  typedef float RealType;

  RealType x;
  RealType y;
  RealType z;
  int l;

  std::vector<RealType> Xm;
  std::vector<RealType> Ym;
  std::vector<RealType> Zm;
  std::vector<int> Lm;
  std::ifstream strM( argv[2] );

  std::vector<RealType> Xf;
  std::vector<RealType> Yf;
  std::vector<RealType> Zf;
  std::vector<int> Lf;
  std::ifstream strF( argv[1] );

  if ( format > 1 )
    {
    while ( strM >> x >> y >> z >> l )
      {
      Xm.push_back( x );
      Ym.push_back( y );
      Zm.push_back( z );
      Lm.push_back( l );
      }

    while ( strF >> x >> y >> z >> l )
      {
      Xf.push_back( x );
      Yf.push_back( y );
      Zf.push_back( z );
      Lf.push_back( l );
      }
    }
  else
    {
    while ( strM >> x >> y >> z )
      {
      Xm.push_back( x );
      Ym.push_back( y );
      Zm.push_back( z );
      Lm.push_back( 1 );
      }

    while ( strF >> x >> y >> z )
      {
      Xf.push_back( x );
      Yf.push_back( y );
      Zf.push_back( z );
      Lf.push_back( 1 );
      }
    }

  // Write moving points to a file with the vectors pointing

  int totalPoints = Xm.size();

  vtkPolyData *moving = vtkPolyData::New();
  vtkPoints *pointsM = vtkPoints::New();
  pointsM->Allocate( totalPoints-2 );
  vtkFloatArray *scalarsM = vtkFloatArray::New();
  scalarsM->SetNumberOfTuples( totalPoints-2 );
  vtkFloatArray *vectorsM = vtkFloatArray::New();
  vectorsM->SetNumberOfComponents( 3 );
  vectorsM->SetNumberOfTuples( totalPoints );

  float xM[3], v[3];

  unsigned int begin, end;
  if ( format > 1 )
    {
    begin = 1;
    end = totalPoints - 1;
    }
  else
    {
    begin = 0;
    end = totalPoints;
    }
  for ( unsigned int i = begin; i < end; i++ )
    {
    xM[0] = Xm[i];
    xM[1] = Ym[i];
    xM[2] = Zm[i];
    pointsM->InsertPoint( i-1, xM );
    scalarsM->InsertValue( i-1, Lm[i] );

    float x, y, z;
    x = y = z = 0.0;
    float N = 0.0;

    if ( format == 1 || format == 3 )
      {
      for ( unsigned int j = 1; j < Lf.size()-1; j++ )
        {
        if ( Lf[j] == Lm[i] )
          {
          x += Xf[j];
          y += Yf[j];
          z += Zf[j];
          N++;
          }
        }
      if ( N > 0.0 )
        {
        x /= N;
        y /= N;
        z /= N;
        v[0] = x - xM[0];
        v[1] = y - xM[1];
        v[2] = z - xM[2];
        }
      else
        {
        v[0] = 0.0;
        v[1] = 0.0;
        v[2] = 0.0;
        }
      }
    else
      {
      v[0] = Xf[i] - xM[0];
      v[1] = Yf[i] - xM[1];
      v[2] = Zf[i] - xM[2];
      }
    vectorsM->InsertTuple( i-begin, v );
    }
  moving->SetPoints( pointsM );
  moving->GetPointData()->SetScalars( scalarsM );
  moving->GetPointData()->SetVectors( vectorsM );

  pointsM->Delete();
  vectorsM->Delete();
  scalarsM->Delete();

  std::string filename = std::string( argv[3] ) + std::string( "Moving.vtk" );

  if ( strcmp( argv[1], argv[2] ) != 0 )
    {
    vtkPolyDataWriter *writerM = vtkPolyDataWriter::New();
    writerM->SetInput( moving );
  //  writerM->SetFileTypeToBinary();
    writerM->SetFileName( filename.c_str() );
    writerM->Write();
    }

  // Write fixed points to a file

  float xF[3];
  totalPoints = Xf.size();

  vtkPolyData *fixed = vtkPolyData::New();
  vtkPoints *pointsF = vtkPoints::New();
  pointsF->Allocate( totalPoints-2 );
  vtkFloatArray *scalarsF = vtkFloatArray::New();
  scalarsF->SetNumberOfTuples( totalPoints-2 );

  for ( unsigned int i = 1; i < Xf.size()-1; i++ )
    {
    xF[0] = Xf[i];
    xF[1] = Yf[i];
    xF[2] = Zf[i];
    pointsF->InsertPoint( i-1, xF );
    scalarsF->InsertValue( i-1, Lf[i] );
    }
  fixed->SetPoints( pointsF );
  fixed->GetPointData()->SetScalars( scalarsF );

  pointsF->Delete();
  scalarsF->Delete();

  if ( strcmp( argv[1], argv[2] ) != 0 )
    {
    filename = std::string( argv[3] ) + std::string( "Fixed.vtk" );
    }
  else
    {
    filename = std::string( argv[3] ) + std::string( ".vtk" );
    }

  vtkPolyDataWriter *writerF = vtkPolyDataWriter::New();
  writerF->SetInput( fixed );
//  writerF->SetFileTypeToBinary();
  writerF->SetFileName( filename.c_str() );
  writerF->Write();

  return 0;
}
