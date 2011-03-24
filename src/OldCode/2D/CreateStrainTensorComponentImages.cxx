#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorFieldGradientImageFunction.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorImageFileReader.h"
#include "itkVectorImageFileWriter.h"
#include "itkVector.h"

#include <string>

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " deformationField outputImageFile strainType [maskImage] " << std::endl;
    std::cout << "       Strain Types " << std::endl;
    std::cout << "         0. Lagrangian strain tensor " << std::endl;
    std::cout << "         1. Eulerian strain tensor " << std::endl;
    std::cout << "         2. Right Cauchy-Green deformation tensor " << std::endl;
    std::cout << "         3. Left Cauchy-Green deformation tensor " << std::endl;
    std::cout << "         4. Right stretch tensor " << std::endl;
    std::cout << "         5. Left stretch tensor " << std::endl;
    exit( 1 );
    }

  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Image<int, ImageDimension> MaskImageType;
  typedef itk::Vector<RealType, 6> VectorType;
  typedef itk::Image<VectorType, ImageDimension> VectorImageType;

  typedef itk::Vector<RealType, 3> FieldVectorType;
  typedef itk::Image<FieldVectorType, ImageDimension> DeformationFieldType;

  /**
   * Read in vector field
   */
  typedef itk::VectorImageFileReader<RealImageType, DeformationFieldType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  MaskImageType::Pointer mask = MaskImageType::New();
  if ( argc < 5 )
    {
    mask->SetOrigin( reader->GetOutput()->GetOrigin() );
    mask->SetSpacing( reader->GetOutput()->GetSpacing() ); 
    mask->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
    mask->Allocate();
    mask->FillBuffer( 1 );
    }
  else
    {
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    MaskReaderType::Pointer maskreader = MaskReaderType::New();
    maskreader->SetFileName( argv[4] );
    maskreader->Update();
    mask = maskreader->GetOutput();
    }

  VectorImageType::Pointer strain = VectorImageType::New();
  strain->SetOrigin( reader->GetOutput()->GetOrigin() );
  strain->SetSpacing( reader->GetOutput()->GetSpacing() ); 
  strain->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  strain->Allocate();
  VectorType V;
  V.Fill( 0 );
  strain->FillBuffer( V );

  /**
   * Evaluate various measures at the specified point
   */ 
  typedef itk::VectorFieldGradientImageFunction<DeformationFieldType> FunctionType;
  FunctionType::Pointer function = FunctionType::New();
  function->SetInputImage( reader->GetOutput() );

  typedef itk::DecomposeTensorFunction<FunctionType::MatrixType> DecomposerType;
  DecomposerType::Pointer decomposer = DecomposerType::New();

  itk::ImageRegionIteratorWithIndex<MaskImageType> ItM
    ( mask, mask->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<VectorImageType> ItS
    ( strain, strain->GetLargestPossibleRegion() );

  ItM.GoToBegin();
  ItS.GoToBegin();
  while ( !ItM.IsAtEnd() )
    {
    if ( ItM.Get() == 0 )
      {
      ++ItM;
      ++ItS;
      continue;
      }  
    FunctionType::MatrixType M, F;
    if ( atoi( argv[3] ) == 0 )
      {
      M = function->EvaluateLagrangianStrainTensorAtIndex( ItM.GetIndex() );
      }
    else if ( atoi( argv[3] ) == 1 )
      {
      M = function->EvaluateEulerianStrainTensorAtIndex( ItM.GetIndex() );
      }
    else if ( atoi( argv[3] ) == 2 )
      {
      M = function->EvaluateRightCauchyGreenDeformationTensorAtIndex( ItM.GetIndex() );
      }
    else if ( atoi( argv[3] ) == 3 )
      {
      M = function->EvaluateLeftCauchyGreenDeformationTensorAtIndex( ItM.GetIndex() );
      }
    else if ( atoi( argv[3] ) == 4 )
      {
      M = function->EvaluateRightStretchTensorAtIndex( ItM.GetIndex() );
      }
    else if ( atoi( argv[3] ) == 5 )
      {
      M = function->EvaluateLeftStretchTensorAtIndex( ItM.GetIndex() );
      }

    VectorType V;
    V[0] = M[0][0];
    V[1] = M[1][0];
    V[2] = M[1][1];
    V[3] = M[2][0];
    V[4] = M[2][1];
    V[5] = M[2][2];

    ItS.Set( V ); 

    ++ItM;
    ++ItS;
    }    

  typedef itk::VectorImageFileWriter<VectorImageType, RealImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetUseZhangNamingConvention( true );
  writer->SetFileName( argv[2] ); 
  writer->SetInput( strain );
  writer->Update(); 

  return 0;
}
