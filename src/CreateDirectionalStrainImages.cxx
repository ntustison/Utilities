#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorFieldGradientImageFunction.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVector.h"

#include <string>
#include <vector>
#include "Common.h"

template <unsigned int ImageDimension>
int CreateDirectionalStrainImages( int argc, char *argv[] )
{

  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Image<int, ImageDimension> MaskImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> VectorImageType;

  /**
   * Read in vector field
   */
  typedef itk::ImageFileReader<VectorImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename VectorImageType::Pointer directions = NULL;
  std::vector<unsigned int> whichComponent;
  try
    {
    typename ReaderType::Pointer readerD = ReaderType::New();
    readerD->SetFileName( argv[4] );
    readerD->Update();
    directions = readerD->GetOutput();
    }
  catch(...)
    {
    whichComponent = ConvertVector<unsigned int>( std::string( argv[4] ) );
    if( whichComponent.size() != 2 )
      {
      std::cerr << "Incorrect component specification." << std::endl;
      return EXIT_FAILURE;
      }
    }

  typename MaskImageType::Pointer mask = NULL;
  if ( argc > 5 )
    {
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    typename MaskReaderType::Pointer maskreader = MaskReaderType::New();
    maskreader->SetFileName( argv[5] );
    maskreader->Update();
    mask = maskreader->GetOutput();
    }

  typename RealImageType::Pointer lagrangian = RealImageType::New();
  lagrangian->SetOrigin( reader->GetOutput()->GetOrigin() );
  lagrangian->SetSpacing( reader->GetOutput()->GetSpacing() );
  lagrangian->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  lagrangian->Allocate();
  lagrangian->FillBuffer( 0 );

  typename RealImageType::Pointer eulerian = RealImageType::New();
  eulerian->SetOrigin( reader->GetOutput()->GetOrigin() );
  eulerian->SetSpacing( reader->GetOutput()->GetSpacing() );
  eulerian->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  eulerian->Allocate();
  eulerian->FillBuffer( 0 );



  /**
   * Evaluate various measures at the specified point
   */
  typedef itk::VectorFieldGradientImageFunction<VectorImageType> FunctionType;
  typename FunctionType::Pointer function = FunctionType::New();
  function->SetInputImage( reader->GetOutput() );

  itk::ImageRegionIteratorWithIndex<RealImageType> ItL
    ( lagrangian, lagrangian->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<RealImageType> ItE
    ( eulerian, eulerian->GetLargestPossibleRegion() );

  RealType N = 0.0;
  RealType mean1 = 0.0;
  RealType var1 = 0.0;
  RealType mean2 = 0.0;
  RealType var2 = 0.0;

  ItL.GoToBegin();
  ItE.GoToBegin();

  while ( !ItL.IsAtEnd() )
    {
    if ( mask && mask->GetPixel( ItL.GetIndex() ) == 0 )
      {
      ++ItL;
      ++ItE;
      continue;
      }

    RealType l, e;

    if( directions )
      {
      l = function->EvaluateLagrangianDirectionalStrainAtIndex( ItL.GetIndex(),
        directions->GetPixel( ItL.GetIndex() ) );
      e = function->EvaluateEulerianDirectionalStrainAtIndex( ItL.GetIndex(),
        directions->GetPixel( ItL.GetIndex() ) );
      }
    else
      {
      typename FunctionType::MatrixType L
        = function->EvaluateLagrangianStrainTensorAtIndex( ItL.GetIndex() );
      typename FunctionType::MatrixType E
        = function->EvaluateEulerianStrainTensorAtIndex( ItL.GetIndex() );
      l = L( whichComponent[0], whichComponent[1] );
      e = E( whichComponent[0], whichComponent[1] );
      }

    ItL.Set( l );
    ItE.Set( e );

    N += 1.0;
    mean1 = mean1*( N - 1.0 )/N + l/N;
    mean2 = mean2*( N - 1.0 )/N + e/N;
    if ( N > 1.0 )
      {
      var1 = var1*( N - 1.0 )/N + ( l - mean1 )*( l - mean1 )/( N - 1.0 );
      var2 = var2*( N - 1.0 )/N + ( e - mean2 )*( e - mean2 )/( N - 1.0 );
      }
    ++ItL;
    ++ItE;
    }

  std::cout << argv[2] << ": lagrangian mean = " << mean1 << ", std1 = " << sqrt( var1 ) << std::endl;
  std::cout << argv[2] << ": eulerian mean   = " << mean2 << ", std1 = " << sqrt( var2 ) << std::endl;

  std::string file1 = std::string( argv[3] ) + std::string( "Lagrangian.nii.gz" );
  std::string file2 = std::string( argv[3] ) + std::string( "Eulerian.nii.gz" );
  typedef itk::ImageFileWriter<RealImageType> WriterType;
  typename WriterType::Pointer writer1 = WriterType::New();
  writer1->SetFileName( file1 );
  writer1->SetInput( lagrangian );
  writer1->Update();
  typename WriterType::Pointer writer2 = WriterType::New();
  writer2->SetFileName( file2 );
  writer2->SetInput( eulerian );
  writer2->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension deformationField"
       << " outputImagePrefix " <<
       " \'whichComponent,e.g. e_xy=0x1 or normalizedDirectionalField\' [maskImage]" << std::endl;
    std::cout << "   Note:  if normalizedDirectionalField is not specified, "
      << " the component strains are saved." << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     CreateDirectionalStrainImages<2>( argc, argv );
     break;
   case 3:
     CreateDirectionalStrainImages<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

