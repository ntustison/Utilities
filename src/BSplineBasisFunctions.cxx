#include "itkCoxDeBoorBSplineKernelFunction.h"

#include <fstream>
#include <iomanip.h>

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << argv[0] << " order ulow uhigh outputFile [delta] [whichDerivative] "
      << " [includeGaussianApproximation]" << std::endl;
    exit( 0 );
    }

  unsigned int order = static_cast<unsigned int>( atoi( argv[1] ) );

  typedef itk::CoxDeBoorBSplineKernelFunction<3> KernelType;
  KernelType::Pointer kernel = KernelType::New();
  kernel->SetSplineOrder( order );

  std::ofstream str( argv[4] );

  float delta = ( argc > 5 ) ? atof( argv[5] ) : 0.01;
  unsigned int whichDer = ( argc > 6 ) ?
    static_cast<unsigned int>( atoi( argv[6] ) ) : 0;
  bool gauss = ( argc > 7 ) ? static_cast<bool>( atoi( argv[7] ) ) : false;

  float ulow = atof( argv[2] );
  float uhigh = atof( argv[3] );

//  for( float u = ulow; u <= uhigh; u+=delta )
//    {
//
//    str << std::setprecision( 5 )
//      << std::setw( 20 ) << u << std::setw( 20 ) << kernel->Evaluate( u );
//    if( whichDer > 0 )
//      {
//      str << std::setw( 20 ) << kernel->EvaluateNthDerivative( u, whichDer );
//      }
//    if( gauss )
//      {
//      /**
//       * The variance of the n-th order B-spline is (n+1)/12.0.
//       *  See Yu-Ping Wang and S.L. Lee, "Scale-Space Derived From B-Splines",
//       *  IEEE-PAMI, 20(10):1040-1055,1998.
//       */
//
//      float variance = static_cast<float>( order + 1 ) / 12.0;
//      float gaussEval = vcl_exp( -0.5*u*u/variance ) / vcl_sqrt( 2.0*vnl_math::pi*variance );
//      str << std::setw( 20 ) << gaussEval;
//      }
//      str << std::endl;
//    }

  KernelType::MatrixType shapeFunctions = kernel->GetShapeFunctionsInZeroToOneInterval();

  str << "0 0 0 0" << std::endl;

  for( unsigned int d = 0; d < shapeFunctions.rows(); d++ )
    {
    KernelType::PolynomialType polynomial( shapeFunctions.get_row( d ) );
    for( float u = ulow; u <= uhigh; u+=delta )
      {
      float value = polynomial.evaluate( u );
      str << u << " " << value << " 0 " << d << std::endl;
      }
    }

  str << "0 0 0 0" << std::endl;

  str.close();
  kernel->Print( std::cout, 4 );


  return EXIT_SUCCESS;
}
