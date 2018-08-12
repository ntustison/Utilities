#include "itkArray.h"
#include "itkImageDuplicator.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"

#include <string>
#include <vector>
#include <fstream>

template <unsigned int ImageDimension>
int GetJointSamples( int argc, char * argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType>  ReaderType;

  typedef unsigned int LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;
  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;

  typename LabelImageType::Pointer maskImage = ITK_NULLPTR;
  try
    {
    typename LabelReaderType::Pointer maskReader = LabelReaderType::New();
    maskReader->SetFileName( argv[3] );
    maskReader->Update();
    maskImage = maskReader->GetOutput();
    }
  catch(...)
    {
    }

  /**
   * list the files
   */
  std::vector<typename ImageType::Pointer> images;

  std::cout << "Using the following files: " << std::endl;
  for( unsigned int n = 4; n < static_cast<unsigned int>( argc ); n++ )
    {
    std::cout << "   " << n-3 << ": " << argv[n] << std::endl;

    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[n] );
    reader->Update();

    images.push_back( reader->GetOutput() );
    }

  std::ofstream os( argv[2] );

  itk::ImageRegionIteratorWithIndex<ImageType> It( images[0],
    images[0]->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( !maskImage || maskImage->GetPixel( It.GetIndex() ) != 0 )
      {
      for( unsigned int n = 0; n < images.size(); n++ )
        {
        os << images[n]->GetPixel( It.GetIndex() ) << ' ';
        }
      os << std::endl;
      }
    }
  os.close();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " imageDimension outputFile "
      << "maskImage imageList" << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GetJointSamples<2>( argc, argv );
     break;
   case 3:
     GetJointSamples<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}



