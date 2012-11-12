
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

template <unsigned int ImageDimension>
int AlignImageCenters( int argc, char *argv[] )
{
  typedef float PixelType;

  typedef float RealType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  typename ReaderType::Pointer fixed = ReaderType::New();
  fixed->SetFileName( argv[2] );
  fixed->Update();

  typename ReaderType::Pointer moving = ReaderType::New();
  moving->SetFileName( argv[3] );
  moving->Update();

  typename ImageType::PointType movingCenter;
  typename ImageType::PointType newOrigin;
  typename ImageType::PointType fixedCenter;
  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    movingCenter[d] = moving->GetOutput()->GetOrigin()[d]
      + 0.5 * ( moving->GetOutput()->GetLargestPossibleRegion().GetSize()[d] - 1 )
        * moving->GetOutput()->GetSpacing()[d] ;
    fixedCenter[d] = fixed->GetOutput()->GetOrigin()[d]
      + 0.5 * ( fixed->GetOutput()->GetLargestPossibleRegion().GetSize()[d] - 1 )
        * fixed->GetOutput()->GetSpacing()[d] ;
    newOrigin[d] =  moving->GetOutput()->GetOrigin()[d]+ fixedCenter[d] - movingCenter[d];
    }

  moving->GetOutput()->SetOrigin( newOrigin );


  typedef itk::ImageFileWriter<ImageType> WriterType;

  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( moving->GetOutput() );
  writer->SetFileName( argv[4] );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << "Usage: " << argv[0]
              << " imageDimension fixedImage movingImage outputImage" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     AlignImageCenters<2>( argc, argv );
     break;
   case 3:
     AlignImageCenters<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}



