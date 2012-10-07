#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkIdList.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkUnsignedIntArray.h"

#include <vector>
#include <fstream>

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputFile outputFile " << std::endl;
    exit( 1 );
    }

  typedef float RealType;

  RealType x;
  RealType y;
  RealType z;
  int l;

  int count = 0;
  int currentL = 0;

  vtkPolyData *data = vtkPolyData::New();
  data->Initialize();

  vtkPoints *points = vtkPoints::New();
  points->Initialize();

  vtkCellArray *lines = vtkCellArray::New();
  lines->Initialize();

  vtkUnsignedIntArray *scalars = vtkUnsignedIntArray::New();
  scalars->Initialize();

  vtkIdList *id = vtkIdList::New();
  id->Initialize();

  std::ifstream str( argv[1] );
  while ( str >> x >> y >> z >> l )
    {
    //std::cout << x << ' ' << y << ' ' << z << ' ' << l << std::endl;
    if ( x == 0 && y == 0 && z == 0 && l == 0 )
      {
      if ( id->GetNumberOfIds() > 0 )
        {
        lines->InsertNextCell( id );
        }
      continue;
      }

    points->InsertPoint( count, x, y, z );
    scalars->InsertNextValue( l );
    if ( l != currentL )
      {
      if ( id->GetNumberOfIds() > 0 )
        {
        lines->InsertNextCell( id );
        }
      id->Delete();
      id = vtkIdList::New();
      id->Initialize();
      id->InsertNextId( count );
      currentL = l;
      }
    else
      {
      id->InsertNextId( count );
      }
    count++;
    }
  str.close();

  data->SetPoints( points );
  data->SetLines( lines );
  data->GetPointData()->SetScalars( scalars );

  points->Delete();
  lines->Delete();
  scalars->Delete();
  id->Delete();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput( data );
//  writer->SetFileTypeToBinary();
  writer->SetFileName( argv[2] );
  writer->Write();

  return 0;
}
