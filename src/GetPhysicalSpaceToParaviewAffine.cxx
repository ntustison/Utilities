#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTransformFileWriter.h"

template<unsigned int ImageDimension>
int GetAffineTransform( int argc, char *argv[] )
{
  typedef itk::Image<char, ImageDimension> ImageType;

  typedef itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension>
    TransformType;
  typename TransformType::OffsetType offset;
  typename TransformType::MatrixType matrix;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  matrix = ( reader->GetOutput()->GetDirection() ).GetInverse();

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    offset[i] = -reader->GetOutput()->GetOrigin()[i];
    }

  typename TransformType::Pointer transform = TransformType::New();
  transform->SetMatrix( matrix );
  transform->SetOffset( offset );

  typedef itk::TransformFileWriter WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( transform );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension "
      << "inputImage outputTransform" << std::endl;
    exit( 0 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GetAffineTransform<2>( argc, argv );
     break;
   case 3:
     GetAffineTransform<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

