#include "itkHilbertPath.hxx"

#include <fstream>

template<unsigned int Dimension>
int Hilbert( int argc, char *argv[] )
{
  typedef itk::HilbertPath<unsigned int, Dimension> PathType;
  typename PathType::Pointer path = PathType::New();
  path->SetHilbertOrder( atoi( argv[2] ) );
  path->Initialize();

  std::ofstream str( argv[3] );
  typedef typename PathType::IndexType IndexType;

  str << "0 0 0 0" << std::endl;
  for( unsigned int s = 0; s < path->NumberOfSteps(); s++ )
    {
    IndexType index = path->Evaluate( s );
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      str << index[d] << " ";
      }
    if( Dimension == 2 )
      {
      str << "0 ";
      }
    str << 1 << std::endl;
    }

  str << "0 0 0 0" << std::endl;
  str.close();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " dimension order output.txt" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     {
     Hilbert<2>( argc, argv );
     break;
     }
   case 3:
     {
     Hilbert<3>( argc, argv );
     break;
     }
   default:
     {
     std::cerr << "Unsupported dimension" << std::endl;
     exit( EXIT_FAILURE );
     }
   }
  return EXIT_SUCCESS;
}

