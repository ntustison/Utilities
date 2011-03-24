#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVector.h"

template <unsigned int ImageDimension>
int GenerateDeformationFieldMagnitudeImage( int argc, char *argv[] )
{

  typedef float RealType;

  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  /**
   * Read in vector field
   */
  typedef itk::ImageFileReader<DeformationFieldType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename RealImageType::Pointer image = RealImageType::New();
  image->SetOrigin( reader->GetOutput()->GetOrigin() );
  image->SetSpacing( reader->GetOutput()->GetSpacing() );
  image->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  image->Allocate();
  image->FillBuffer( 0 );

  itk::ImageRegionIteratorWithIndex<DeformationFieldType>
    It1( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<RealImageType>
    It2( image, image->GetLargestPossibleRegion() );
  for ( It1.GoToBegin(), It2.GoToBegin(); !It1.IsAtEnd(); ++It1, ++It2 )
    {
    It2.Set( It1.Get().GetNorm() );
    }

  typedef itk::ImageFileWriter<RealImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( image );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension inputDeformationField "
              << "outputImage" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GenerateDeformationFieldMagnitudeImage<2>( argc, argv );
     break;
   case 3:
     GenerateDeformationFieldMagnitudeImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

