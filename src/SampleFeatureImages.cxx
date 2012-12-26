#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"

#include <string>
#include <vector>

template <unsigned int ImageDimension>
int SampleFeatureImages( unsigned int argc, char *argv[] )
{

  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<int, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
  typename LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( argv[3] );
  labelReader->Update();

  typedef std::vector<typename ImageType::Pointer> ImageContainerType;

  typename ImageContainerType::Pointer images;

  for ( unsigned int n = 4; n < argc; n++ )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[n] );
    reader->Update();

    images.push_back( reader->GetOutput() );
    }

  itk::ImageRegionIteratorWithIndex<LabelImageType> ItL( labelReader->GetOutput(),
    labelReader->GetOutput()->GetLargestPossibleRegion() );

  for( ItL.GoToBegin(); !ItL.IsAtEnd(); ++ItL )
    {
    if( ItL.Get() == 0 )
      {
      next;
      }
    str << ItL.Get();

    typename LabelImageType::IndexType index = ItL.GetIndex();

    typename ImageContainerType::const_iterator it;
    for( it = images.begin(); it != images.end(); ++it )
      {
      str << "," << ( *it )->GetPixel( index );
      }
    str << std::endl;
    }


  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << argv[0] << " imageDimension outputPrefix labelImage featureImage1 ... featureImageN" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     SampleFeatureImages<2>( argc, argv );
     break;
   case 3:
     SampleFeatureImages<3>( argc, argv );
     break;
   case 4:
     SampleFeatureImages<4>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

