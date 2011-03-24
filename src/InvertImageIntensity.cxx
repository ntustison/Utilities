#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIterator.h"

template <unsigned int ImageDimension>
int Invert( unsigned int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  PixelType minimumValue = itk::NumericTraits<PixelType>::max();
  PixelType maximumValue = itk::NumericTraits<PixelType>::NonpositiveMin();

  itk::ImageRegionIterator<ImageType> It( reader->GetOutput(),
    reader->GetOutput()->GetRequestedRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( It.Get() > maximumValue )
      {
      maximumValue = It.Get();
      }
    if( It.Get() < minimumValue )
      {
      minimumValue = It.Get();
      }
    }

  PixelType slope = ( minimumValue - maximumValue ) / ( maximumValue - minimumValue );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    It.Set( slope * ( It.Get() - maximumValue ) + minimumValue );
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( reader->GetOutput() );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension inputImage outputImage" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     Invert<2>( argc, argv );
     break;
   case 3:
     Invert<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

