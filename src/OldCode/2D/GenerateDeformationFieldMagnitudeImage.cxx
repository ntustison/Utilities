#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVector.h"
#include "itkVectorImageFileReader.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] 
              << " inputDeformationField outputImage [maskImage]" << std::endl;
    exit( 1 );
    } 

  typedef float RealType;

  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Image<unsigned int, ImageDimension> MaskImageType;

  MaskImageType::Pointer mask = MaskImageType::New();
  mask = NULL;

  if ( argc <= 4 )
    {
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName( argv[3] );
    maskReader->Update();
    mask = maskReader->GetOutput();
    }

  /**
   * Read in vector field
   */
  typedef itk::VectorImageFileReader<RealImageType, DeformationFieldType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  RealImageType::Pointer image = RealImageType::New();
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
    if ( mask && mask->GetPixel( It1.GetIndex() ) > 0 )
      {
      It2.Set( It1.Get().GetNorm() );
      } 
    }   
 
  typedef itk::ImageFileWriter<RealImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( image );
  writer->Update();

  return 0;
}
