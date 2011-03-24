#include "itkDecomposeTensorFunction.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPoint.h"
#include "itkMatrix.h"
#include "itkVariableSizeMatrix.h"
#include "itkVector.h"
#include "itkVectorContainer.h"
#include "itkVectorImageFileWriter.h"

#include <fstream.h>
#include <iomanip.h>

template <unsigned int ImageDimension>
int GenerateTPSDeformationField( unsigned int argc, char *argv[] )
{

  typedef float RealType;
  typedef itk::Point<RealType, ImageDimension> PointType;
  typedef itk::Matrix<RealType, ImageDimension+1, ImageDimension> 
    AffineMatrixType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::VectorContainer<unsigned, VectorType> CoefficientContainerType;  
  
  // Read in affine matrix coefficients
  ifstream tpsAffine( argv[4] );
  
  AffineMatrixType A;
  A.Fill( 1 );
  for( unsigned int d = 0; d < ImageDimension+1; d++ )
    {
    for( unsigned int e = 0; e < ImageDimension; e++ )
      {
      tpsAffine >> A( d, e );
      }
    }
  tpsAffine.close();  

  // Read in tps coefficients
  ifstream tpsDeform( argv[5] );
  
  typename CoefficientContainerType::Pointer C = CoefficientContainerType::New();  
  C->Initialize();
  unsigned long count = 0;
  while( !tpsDeform.eof() )
    {
    bool isEOF = false; 
    VectorType coeff;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      tpsDeform >> coeff[d];
      if( tpsDeform.eof() )
        {
        isEOF = true;
        break;
        }
      } 
    if( isEOF )
      {
      break;
      }  
    C->InsertElement( count++, coeff );     
    }
  tpsDeform.close();  
  
  // Join tps and affine parameters into single matrix
  typedef itk::VariableSizeMatrix<RealType> MatrixType;
  MatrixType params;
  params.SetSize( ImageDimension + 1 + C->Size(), ImageDimension );
  for( unsigned int n = 0; n < ImageDimension + 1; n++ )
    {
    for( unsigned int d = 0; d < ImageDimension; d++ )
      { 
      params( n, d ) = A( n, d );
      }
    }
  for( unsigned int n = ImageDimension + 1; n < params.Rows(); n++ )
    {
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      VectorType V = C->GetElement( n - ImageDimension - 1 ); 
      params( n, d ) = V[d];
      }
    }    

  // Read in control points

  ifstream cpFile( argv[3] );
  typename CoefficientContainerType::Pointer 
    L = CoefficientContainerType::New();
  L->Initialize();
  count = 0;
  while( !cpFile.eof() )
    {
    bool isEOF = false; 
    VectorType coeff;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      cpFile >> coeff[d];
      if( cpFile.eof() )
        {
        isEOF = true;
        break;
        }
      } 
    if( isEOF )
      {
      break;
      }  
    L->InsertElement( count++, coeff );     
    }
  cpFile.close();  

  MatrixType Pn;
  Pn.SetSize( L->Size(), ImageDimension + 1 ); 
  Pn.Fill( 1 );
  for( unsigned int n = 0; n < L->Size(); n++ )
    {
    VectorType coeff = L->GetElement( n ); 
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      Pn(n, d + 1) = coeff[d];
      }
    }

  typedef itk::DecomposeTensorFunction<MatrixType> DecomposerType;
  typename DecomposerType::Pointer decomposer = DecomposerType::New();  
  
  MatrixType Q;
  MatrixType R;
  decomposer->EvaluateQRDecomposition( Pn, Q, R );

  vnl_matrix<RealType> PP = ( Q.GetVnlMatrix() ).extract( 
    L->Size(), L->Size() - ImageDimension - 1, 0, ImageDimension + 1 );
    
  // Read in image
  typedef itk::Image<RealType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::Image<VectorType, ImageDimension> VectorImageType;
  typename VectorImageType::Pointer output = VectorImageType::New();
  output->SetOrigin( reader->GetOutput()->GetOrigin() );
  output->SetSpacing( reader->GetOutput()->GetSpacing() );
  output->SetRegions( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  output->SetDirection( reader->GetOutput()->GetDirection() );
  output->Allocate();

  itk::ImageRegionIteratorWithIndex<ImageType> ItR( reader->GetOutput(),
    reader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<VectorImageType> ItV( output,
    output->GetLargestPossibleRegion() );
  for( ItR.GoToBegin(), ItV.GoToBegin(); !ItR.IsAtEnd(); ++ItR, ++ItV )
    {
    PointType point;
    output->TransformIndexToPhysicalPoint( ItV.GetIndex(), point );

    vnl_matrix<RealType> U;
    U.set_size( 1, L->Size() );
    U.fill( 0 );
    for( unsigned int n = 0; n < U.cols(); n++ )
      {
      RealType r = 0;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        r += vnl_math_sqr( Pn(n, d+1) - point[d] );
        }    
      r = vcl_sqrt( r );  
      if( r > 0 )
        {
        if( ImageDimension == 2 )
          {
          U(0, n) = r * r * vcl_log( r );
          }
        else
          {
          U(0, n) = -r;
          }
        }
      }
    vnl_matrix<RealType> basis;
    basis.set_size( 1, L->Size() );
    basis(0, 0) = 1;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      basis(0, 1 + d) = point[d]; 
      }
      
    basis.update( U * PP, static_cast<unsigned int>( 0 ), ImageDimension+1 );

    vnl_matrix<RealType> warpedPoint = basis * params.GetVnlMatrix();

    VectorType V;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      V[d] = warpedPoint(0, d) - point[d];
      }
    ItV.Set( V );
    }  
    
  typedef itk::Image<RealType, ImageDimension> RealImageType;  
  
  typedef itk::VectorImageFileWriter<VectorImageType, RealImageType> 
    WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( output );
  writer->SetFileName( argv[6] );
  writer->Update();
  
  return EXIT_SUCCESS;
  
}

int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension domainImage "
      << "controlPointFile tpsAffineFile tpsDeformFile outputField" << std::endl; 
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GenerateTPSDeformationField<2>( argc, argv );
     break;
   case 3:
     GenerateTPSDeformationField<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
