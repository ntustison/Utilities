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
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " deformationField outputImage [UseAvantsNamingConvention]" << std::endl;
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
  if ( argc == 4 )
    {
    reader->SetUseAvantsNamingConvention( atoi( argv[3] ) );
    } 
  reader->Update();

  RealImageType::Pointer jacobian = RealImageType::New();
  jacobian->SetOrigin( reader->GetOutput()->GetOrigin() );
  jacobian->SetSpacing( reader->GetOutput()->GetSpacing() ); 
  jacobian->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  jacobian->Allocate();

  /**
   * Evaluate various measures at the specified point
   */ 
  typedef itk::VectorFieldGradientImageFunction<VectorImageType> FunctionType;
  FunctionType::Pointer function = FunctionType::New();
  function->SetInputImage( reader->GetOutput() );

  itk::ImageRegionIteratorWithIndex<VectorImageType> It
    ( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    jacobian->SetPixel( It.GetIndex(), function->EvaluateJacobianDeterminantAtIndex( It.GetIndex() ) ); 
    }    

  typedef itk::ImageFileWriter<RealImageType> RealImageWriterType;
  RealImageWriterType::Pointer realwriter = RealImageWriterType::New();
  realwriter->SetFileName( argv[2] ); 
  realwriter->SetInput( jacobian );
  realwriter->Update(); 
  return 0;
}
