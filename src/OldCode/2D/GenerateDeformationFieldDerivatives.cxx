#include "itkDecomposeTensorFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorFieldGradientImageFunction.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorImageFileReader.h"
#include "itkVector.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc != 2 + ImageDimension )
    {
    if ( ImageDimension == 2 )
      {
      std::cout << "Usage: " << argv[0] << "deformationField point_x point_y" << std::endl;
      }
    else
      {
      std::cout << "Usage: " << argv[0] << "deformationField point_x point_y point_z" << std::endl;
      } 
    exit( 1 );
    }
  
  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> VectorImageType;

  /**
   * Read in vector field
   */
  typedef itk::VectorImageFileReader<RealImageType, VectorImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->SetUseAvantsNamingConvention( true );
  reader->Update();

  /**
   * Evaluate various measures at the specified point
   */ 
  typedef itk::VectorFieldGradientImageFunction<VectorImageType> FunctionType;
  FunctionType::Pointer function = FunctionType::New();
  function->SetInputImage( reader->GetOutput() );

  FunctionType::PointType point;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    point[i] = atof( argv[2+i] );
    }

  std::cout << "Evaluation at " << point << std::endl;
  std::cout << "  Deformation gradient tensor = " << std::endl 
            << function->EvaluateDeformationGradientTensor( point ) << std::endl;
  std::cout << "  Jacobian = " << std::endl << function->EvaluateJacobian( point ) << std::endl;
  std::cout << "  Jacobian Determinant = " << std::endl 
            << function->EvaluateJacobianDeterminant( point ) << std::endl;  
  std::cout << "  Lagrangian Strain Tensor = " << std::endl 
            << function->EvaluateLagrangianStrainTensor( point ) << std::endl;  
  std::cout << "  Eulerian Strain Tensor = " << std::endl 
            << function->EvaluateEulerianStrainTensor( point ) << std::endl; 
  std::cout << std::endl; 

  itk::DecomposeTensorFunction<FunctionType::MatrixType, RealType> decomposer;
  FunctionType::MatrixType E = function->EvaluateLagrangianStrainTensor( point );
  FunctionType::MatrixType L;

  std::cout << "Various decompositions of the Lagrangian strain tensor" << std::endl;
  decomposer.EvaluateCholeskyDecomposition( E, L );
  std::cout << "  Cholesky Decomposition = " << std::endl << L << std::endl;

  FunctionType::MatrixType D, V;
  decomposer.EvaluateSymmetricEigenDecomposition( E, D, V );
  std::cout << "  Symmetric EigenDecomposition (D) = " << std::endl << D << std::endl;
  std::cout << "  Symmetric EigenDecomposition (V) = " << std::endl << V << std::endl;

  std::cout << "Various decompositions of the Jacobian tensor" << std::endl;
  FunctionType::MatrixType J = function->EvaluateJacobian( point );
  FunctionType::MatrixType U;
  FunctionType::MatrixType W;
  decomposer.EvaluateSVDDecomposition( J, U, V, W );
  std::cout << "  SVD Decomposition (U) = " << std::endl << U << std::endl;
  std::cout << "  SVD Decomposition (V) = " << std::endl << V << std::endl;
  std::cout << "  SVD Decomposition (W) = " << std::endl << W << std::endl;
  decomposer.EvaluateSVDEconomyDecomposition( J, V, W );
  std::cout << "  SVD Economy Decomposition (V) = " << std::endl << V << std::endl;
  std::cout << "  SVD Economy Decomposition (W) = " << std::endl << W << std::endl;
  FunctionType::MatrixType R;
  FunctionType::MatrixType S;
  decomposer.EvaluateRightPolarDecomposition( J, R, S );
  std::cout << "  Right Polar Decomposition (R) = " << std::endl << R << std::endl;
  std::cout << "  Right Polar Decomposition (S) = " << std::endl << S << std::endl;
  decomposer.EvaluateLeftPolarDecomposition( J, S, R );
  std::cout << "  Left Polar Decomposition (S) = " << std::endl << S << std::endl;
  std::cout << "  Left Polar Decomposition (R) = " << std::endl << R << std::endl;
  decomposer.EvaluateEigenDecomposition( J, D, V );
  std::cout << "  Unsymmetric EigenDecomposition (D) = " << std::endl << D << std::endl;
  std::cout << "  Unsymmetric EigenDecomposition (V) = " << std::endl << V << std::endl;
  FunctionType::MatrixType Q;
  decomposer.EvaluateQRDecomposition( J, Q, R );
  std::cout << "  QR Decomposition (Q) = " << std::endl << Q << std::endl;
  std::cout << "  QR Decomposition (R) = " << std::endl << R << std::endl;

  
  typedef itk::ImageFileReader<RealImageType> RealImageReaderType;
  RealImageReaderType::Pointer realreader = RealImageReaderType::New();
  realreader->SetFileName( "GR9_GR21jacobian.hdr" );
  realreader->Update();
  
  RealImageType::Pointer diffImage = RealImageType::New(); 
  diffImage->SetRegions( realreader->GetOutput()->GetLargestPossibleRegion() );
  diffImage->Allocate();

  itk::ImageRegionIteratorWithIndex<RealImageType> It( realreader->GetOutput(),
    realreader->GetOutput()->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    RealType j = function->EvaluateJacobianDeterminantAtIndex( It.GetIndex() );
    diffImage->SetPixel( It.GetIndex(), fabs( j - It.Get() ) );
    }

  typedef itk::ImageFileWriter<RealImageType> RealImageWriterType;
  RealImageWriterType::Pointer realwriter = RealImageWriterType::New();
  realwriter->SetFileName( "diff.hdr" ); 
  realwriter->SetInput( diffImage );
  realwriter->Update(); 
  return 0;
}
