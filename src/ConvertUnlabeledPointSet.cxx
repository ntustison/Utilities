#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"

#include <fstream>

#include <iomanip>


template <unsigned int ImageDimension>
int ConvertUnlabeledPointSet( unsigned int argc, char *argv[] )

{


  typedef double RealType;

  typedef itk::PointSet<long, ImageDimension> PointSetType;
  typename PointSetType::Pointer points = PointSetType::New();
  points->Initialize();

  long label = ( argc > 4 ) ? atoi( argv[4] ) : 1;

  std::fstream str( argv[2], std::ios::in );

  unsigned long count = 0;
  typename PointSetType::PointType point;
  while( str >> point )
    {
    points->SetPoint( count, point );
    points->SetPointData( count, label );
    count++;
    }
  str.close();

  std::cout << "Number of points: " << count << std::endl;

  typedef itk::LabeledPointSetFileWriter<PointSetType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( points );
  writer->Update();

  return 0;

}


int main( int argc, char *argv[] )
{
  if ( argc < 4 )

    {

    std::cout << "Usage: " << argv[0] << " imageDimension inputFile outputFile [label]" << std::endl;
    std::cout << "   Assumes the following text file format: " << std::endl;
    std::cout << "      x_0 y_0 [z_0] " << std::endl;
    std::cout << "      x_1 y_1 [z_1] " << std::endl;
    std::cout << "      .   .   .   " << std::endl;
    std::cout << "      x_N y_N [z_N] " << std::endl;

    exit( 1 );

    }


  switch( atoi( argv[1] ) )
   {
   case 2:
     ConvertUnlabeledPointSet<2>( argc, argv );
     break;
   case 3:
     ConvertUnlabeledPointSet<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
