#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include <fstream>
#include <vector>

template <unsigned int ImageDimension>
int GenerateImage( int argc, char *argv[] )
{

  typedef int PixelType;
  typedef float RealType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[2] );
  typename ImageType::Pointer referenceImage = imageReader->GetOutput();
  referenceImage->Update();
  referenceImage->DisconnectPipeline();

  typename RealImageType::Pointer outputImage = RealImageType::New();
  outputImage->CopyInformation( referenceImage );
  outputImage->SetRegions( referenceImage->GetLargestPossibleRegion() );
  outputImage->Allocate();
  outputImage->FillBuffer( 0 );

  typedef typename ImageType::IndexType IndexType;
  std::vector<IndexType> indices;
  std::ifstream indexFile( argv[4] );

  unsigned int count = 0;
  while( indexFile.good() )
    {
    std::string line;
    std::getline( indexFile, line );
    char * pch = std::strtok( const_cast<char *>( line.c_str() ), "," );

    IndexType index;
    unsigned int d = 0;
    while( pch != NULL )
      {
      index[d++] = atoi( pch );
      pch = std::strtok( NULL, "," );
      }
    indices.push_back( index );
    }

  std::vector<RealType> pixelValues;
  if( argc > 5 )
    {
    std::ifstream valuesFile( argv[5] );
    while( valuesFile.good() )
      {
      std::string line;
      std::getline( valuesFile, line );
      pixelValues.push_back( atof( line.c_str() ) );
      }
    if( pixelValues.size() != indices.size() )
      {
      std::cerr << "The size of the indices (=" << indices.size()
        << ") does not match the number of pixel values (=" << pixelValues.size() << ")." << std::endl;
      return EXIT_FAILURE;
      }
    }

  typename std::vector<IndexType>::const_iterator it;
  for( it = indices.begin(); it != indices.end(); ++it )
    {
    if( pixelValues.size() == indices.size() )
      {
      outputImage->SetPixel( *it, pixelValues[it - indices.begin()] );
      }
    else
      {
      outputImage->SetPixel( *it, 1 );
      }
    }

  typedef itk::ImageFileWriter<RealImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( outputImage );
  writer->SetFileName( argv[3] );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension referenceImage outputImage indexFile [sampleFile]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GenerateImage<2>( argc, argv );
     break;
   case 3:
     GenerateImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

