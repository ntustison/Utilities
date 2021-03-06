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
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " deformationField outputImagePrefix [maskImage]" << std::endl;
    exit( 1 );
    }

  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Image<int, ImageDimension> MaskImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> VectorImageType;

  /**
   * Read in vector field
   */
  typedef itk::VectorImageFileReader<RealImageType, VectorImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  MaskImageType::Pointer mask = MaskImageType::New();
  if ( argc < 4 )
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
    maskreader->SetFileName( argv[3] );
    maskreader->Update();
    mask = maskreader->GetOutput();
    }

  VectorImageType::Pointer strain1 = VectorImageType::New();
  strain1->SetOrigin( reader->GetOutput()->GetOrigin() );
  strain1->SetSpacing( reader->GetOutput()->GetSpacing() ); 
  strain1->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  strain1->Allocate();
  VectorImageType::Pointer strain2 = VectorImageType::New();
  strain2->SetOrigin( reader->GetOutput()->GetOrigin() );
  strain2->SetSpacing( reader->GetOutput()->GetSpacing() ); 
  strain2->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  strain2->Allocate();
  VectorImageType::Pointer strain3 = VectorImageType::New();
  strain3->SetOrigin( reader->GetOutput()->GetOrigin() );
  strain3->SetSpacing( reader->GetOutput()->GetSpacing() ); 
  strain3->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  strain3->Allocate();

  VectorType V;
  V.Fill( 0 );
  strain1->FillBuffer( V );
  strain2->FillBuffer( V );
  strain3->FillBuffer( V );

  /**
   * Evaluate various measures at the specified point
   */ 
  typedef itk::VectorFieldGradientImageFunction<VectorImageType> FunctionType;
  FunctionType::Pointer function = FunctionType::New();
  function->SetInputImage( reader->GetOutput() );

  typedef itk::DecomposeTensorFunction<FunctionType::MatrixType> DecomposerType;
  DecomposerType::Pointer decomposer = DecomposerType::New();

  itk::ImageRegionIteratorWithIndex<MaskImageType> ItM
    ( mask, mask->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<VectorImageType> ItS1
    ( strain1, strain1->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<VectorImageType> ItS2
    ( strain2, strain2->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<VectorImageType> ItS3
    ( strain3, strain3->GetLargestPossibleRegion() );

  RealType N = 0.0;
  RealType mean1 = 0.0;
  RealType var1 = 0.0;
  RealType mean2 = 0.0;
  RealType var2 = 0.0;
  RealType mean3 = 0.0;
  RealType var3 = 0.0;

  ItM.GoToBegin();
  ItS1.GoToBegin();
  ItS2.GoToBegin();
  ItS3.GoToBegin();

  while ( !ItM.IsAtEnd() )
    {
    if ( ItM.Get() == 0 )
      {
      ++ItM;
      ++ItS1;
      ++ItS2;
      ++ItS3;
      continue;
      }  
    FunctionType::MatrixType E = function->EvaluateLagrangianStrainTensorAtIndex( ItM.GetIndex() );
    FunctionType::MatrixType D;
    FunctionType::MatrixType V; 
    decomposer->EvaluateSymmetricEigenDecomposition( E, D, V );    
    VectorType P1;  P1.Fill( 0 );
    VectorType P2; 
    VectorType P3; 
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      if ( ImageDimension == 3 )
        {
        P1[i] = D[2][2]*V[i][2];  
        } 
      P2[i] = D[1][1]*V[i][1];
      P3[i] = D[0][0]*V[i][0];  
      }
    ItS1.Set( P1 ); 
    ItS2.Set( P2 ); 
    ItS3.Set( P3 ); 

    N += 1.0;
    if ( ImageDimension == 3 )
      { 
      mean1 = mean1*( N - 1.0 )/N + D[2][2]/N;
      } 
    mean2 = mean2*( N - 1.0 )/N + D[1][1]/N;
    mean3 = mean3*( N - 1.0 )/N + D[0][0]/N;
    if ( N > 1.0 )
      {
      if ( ImageDimension == 3 )
        { 
        var1 = var1*( N - 1.0 )/N + ( D[2][2] - mean1 )*( D[2][2] - mean1 )/( N - 1.0 );
        } 
      var2 = var2*( N - 1.0 )/N + ( D[1][1] - mean2 )*( D[1][1] - mean2 )/( N - 1.0 );
      var3 = var3*( N - 1.0 )/N + ( D[0][0] - mean3 )*( D[0][0] - mean3 )/( N - 1.0 );
      }
    ++ItM;
    ++ItS1;
    ++ItS2;
    ++ItS3;
    }    
  typedef itk::VectorImageFileWriter<VectorImageType, RealImageType> WriterType;

  if ( ImageDimension == 3 )
    { 
    std::cout << argv[1] << ": mean1 = " << mean1 << ", std1 = " << sqrt( var1 ) << std::endl;
    std::string file1 = std::string( argv[2] ) + std::string( "1.nii.gz" );
    WriterType::Pointer writer1 = WriterType::New();
    writer1->SetFileName( file1 ); 
    writer1->SetInput( strain1 );
    writer1->Update(); 
    }   
  std::cout << argv[1] << ": mean2 = " << mean2 << ", std1 = " << sqrt( var2 ) << std::endl;
  std::cout << argv[1] << ": mean3 = " << mean3 << ", std1 = " << sqrt( var3 ) << std::endl;

  std::string file2 = std::string( argv[2] ) + std::string( "2.nii.gz" );
  std::string file3 = std::string( argv[2] ) + std::string( "3.nii.gz" );

  WriterType::Pointer writer2 = WriterType::New();
  writer2->SetFileName( file2 ); 
  writer2->SetInput( strain2 );
  writer2->Update(); 
  WriterType::Pointer writer3 = WriterType::New();
  writer3->SetFileName( file3 ); 
  writer3->SetInput( strain3 );
  writer3->Update(); 

  return 0;
}
