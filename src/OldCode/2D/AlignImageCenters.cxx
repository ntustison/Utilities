#include <stdio.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "global.h"

int main( int argc, char *argv[] )
{

  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0]
              << " fixedImage movingImage outputImage" << std::endl;
    exit( 1 );
    }

  typedef float RealType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer fixed = ReaderType::New();
  fixed->SetFileName( argv[1] );
  fixed->Update();

  ReaderType::Pointer moving = ReaderType::New();
  moving->SetFileName( argv[2] );
  moving->Update();

  ImageType::PointType movingCenter;
  ImageType::PointType newOrigin;
  ImageType::PointType fixedCenter;
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

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( moving->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();

  return 0;
}
