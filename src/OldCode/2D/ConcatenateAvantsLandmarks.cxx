#include <stdio.h>
#include <fstream.h>
#include <iostream>

int main( int argc, char *argv[] )        
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " outputPointSet inputPointSet1 inputPointSet2 ... [inputPointSetN]" << std::endl;
    exit( 0 );
    } 

  typedef double RealType; 
  unsigned long count = 0;
  RealType x, y, z, l;

  ofstream str( argv[1] );
  str << "0 0 0 0" << std::endl;

  for ( unsigned int n = 2; n < static_cast<unsigned int>( argc ); n++ )
    {
    ifstream strF( argv[n] );

    while ( strF >> x >> y >> z >> l )
      {
      if ( x == 0 && y == 0 && z == 0 && l == 0 )
        {
        continue; 
        }
      str << x << " " << y << " " << z << " " << l << std::endl; 
      count++;
      } 
    strF.close();  
    } 
  str << "0 0 0 0" << std::endl;
  str.close();

  std::cout << count << " total points. " << std::endl;

  return EXIT_SUCCESS;
}     

