#include <fstream>
#include <iostream>

#include <vector>

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " outputPointSet inputPointSet1 inputPointSet2 ... [inputPointSetN]" << std::endl;
    exit( 0 );
    }

  typedef double RealType;
  unsigned long count = 0;
  std::vector<RealType> X;
  std::vector<RealType> Y;
  std::vector<RealType> Z;
  std::vector<long> L;


  for ( unsigned int n = 2; n < static_cast<unsigned int>( argc ); n++ )
    {
    ifstream strF( argv[n] );

    RealType x, y, z;
    long l;
    while ( strF >> x >> y >> z >> l )
      {
      if ( x == 0 && y == 0 && z == 0 && l == 0 )
        {
        continue;
        }
      X.push_back( x );
      Y.push_back( y );
      Z.push_back( z );
      L.push_back( l );
      count++;
      }
    strF.close();
    }

  std::ofstream str( argv[1] );
  str << "0 0 0 0" << std::endl;
  for ( unsigned long i = 0; i < count; i++ )
    {
    str << X[i] << " " << Y[i] << " " << Z[i] << " " << L[i] << std::endl;
    }
  str << "0 0 0 0" << std::endl;
  str.close();

  std::cout << count << " total points. " << std::endl;

  return EXIT_SUCCESS;
}

