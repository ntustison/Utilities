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
    std::cout << "Usage: " << argv[0] << " deformationField normalizedDirectionalField outputImagePrefix [maskImage]" << std::endl;
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

  ReaderType::Pointer directions = ReaderType::New();
  directions->SetFileName( argv[2] );
  directions->Update();

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

  VectorImageType::Pointer lagrangian = VectorImageType::New();
  lagrangian->SetOrigin( reader->GetOutput()->GetOrigin() );
  lagrangian->SetSpacing( reader->GetOutput()->GetSpacing() ); 
  lagrangian->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  lagrangian->Allocate();
  VectorImageType::Pointer eulerian = VectorImageType::New();
  eulerian->SetOrigin( reader->GetOutput()->GetOrigin() );
  eulerian->SetSpacing( reader->GetOutput()->GetSpacing() ); 
  eulerian->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  eulerian->Allocate();

  VectorType V;
  V.Fill( 0 );
  lagrangian->FillBuffer( V );
  eulerian->FillBuffer( V );

  /**
   * Evaluate various measures at the specified point
   */ 
  typedef itk::VectorFieldGradientImageFunction<VectorImageType> FunctionType;
  FunctionType::Pointer function = FunctionType::New();
  function->SetInputImage( reader->GetOutput() );

  itk::DecomposeTensorFunction<FunctionType::MatrixType> decomposer;

  itk::ImageRegionIteratorWithIndex<MaskImageType> ItM
    ( mask, mask->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<VectorImageType> ItF
    ( directions->GetOutput(), directions->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<VectorImageType> ItL
    ( lagrangian, lagrangian->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<VectorImageType> ItE
    ( eulerian, eulerian->GetLargestPossibleRegion() );

  RealType N = 0.0;
  RealType mean1 = 0.0;
  RealType var1 = 0.0;
  RealType mean2 = 0.0;
  RealType var2 = 0.0;

  ItM.GoToBegin();
  ItL.GoToBegin();
  ItE.GoToBegin();

  while ( !ItM.IsAtEnd() )
    {
    if ( ItM.Get() == 0 )
      {
      ++ItM;
      ++ItL;
      ++ItE;
      ++ItF;
      continue;
      }  
    RealType l = function->EvaluateLagrangianDirectionalStrainAtIndex( ItM.GetIndex(), ItF.Get() );
    RealType e = function->EvaluateEulerianDirectionalStrainAtIndex( ItM.GetIndex(), ItF.Get() );

    ItL.Set( l*ItF.Get() ); 
    ItE.Set( e*ItF.Get() ); 

    N += 1.0;
    mean1 = mean1*( N - 1.0 )/N + l/N;
    mean2 = mean2*( N - 1.0 )/N + e/N;
    if ( N > 1.0 )
      {
      var1 = var1*( N - 1.0 )/N + ( l - mean1 )*( l - mean1 )/( N - 1.0 );
      var2 = var2*( N - 1.0 )/N + ( e - mean2 )*( e - mean2 )/( N - 1.0 );
      }
    ++ItM;
    ++ItL;
    ++ItE;
    ++ItF;
    }    
  std::cout << argv[1] << ": lagrangian mean = " << mean1 << ", std1 = " << sqrt( var1 ) << std::endl;
  std::cout << argv[1] << ": eulerian mean   = " << mean2 << ", std1 = " << sqrt( var2 ) << std::endl;

  std::string file1 = std::string( argv[3] ) + std::string( "Lagrangian.hdr" );
  std::string file2 = std::string( argv[3] ) + std::string( "Eulerian.hdr" );

  typedef itk::VectorImageFileWriter<VectorImageType, RealImageType> WriterType;
  WriterType::Pointer writer1 = WriterType::New();
  writer1->SetFileName( file1 ); 
  writer1->SetInput( lagrangian );
  writer1->Update(); 
  WriterType::Pointer writer2 = WriterType::New();
  writer2->SetFileName( file2 ); 
  writer2->SetInput( eulerian );
  writer2->Update(); 

  return 0;
}
