#include <stdio.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNeighborhoodFirstOrderStatisticsImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

#include <fstream>

#include <string>
#include <vector>

#include "Common.h"

template <unsigned int ImageDimension>
int CalculateStatisticsImage( int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::VectorImage<PixelType, ImageDimension> VectorImageType;
  typedef itk::FlatStructuringElement<ImageDimension> KernelType;
  typedef itk::NeighborhoodFirstOrderStatisticsImageFilter<ImageType, VectorImageType, KernelType> TextureFilterType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename KernelType::SizeType radius;
  std::vector<PixelType> rad = ConvertVector<PixelType>( std::string( argv[5] ) );
  if( rad.size() != ImageDimension )
    {
    radius.Fill( rad[0] );
    }
  else
    {
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      radius[d] = rad[d];
      }
    }
  radius.Fill( atoi( argv[5] ) );
  KernelType kernel = KernelType::Box( radius );

  typename TextureFilterType::Pointer filter = TextureFilterType::New();
  filter->SetKernel( kernel );
  filter->SetInput( reader->GetOutput() );
  filter->Update();

  typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ImageType> IndexSelectionType;
  typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
  indexSelectionFilter->SetInput( filter->GetOutput() );

  int operation = atoi( argv[4] );

  switch( operation )
    {
    case 0:
      {
      indexSelectionFilter->SetIndex( 0 );
      break;
      }
    case 1:
      {
      indexSelectionFilter->SetIndex( 1 );
      break;
      }
    case 2:
      {
      indexSelectionFilter->SetIndex( 2 );
      break;
      }
    case 3:
      {
      indexSelectionFilter->SetIndex( 3 );
      break;
      }
    case 4:
      {
      indexSelectionFilter->SetIndex( 4 );
      break;
      }
    case 5:
      {
      indexSelectionFilter->SetIndex( 5 );
      break;
      }
    case 6:
      {
      indexSelectionFilter->SetIndex( 6 );
      break;
      }
    case 7:
      {
      indexSelectionFilter->SetIndex( 7 );
      break;
      }
    default:
      {
      std::cerr << "Unrecognized option: " << operation << std::endl;
      return EXIT_FAILURE;
      }
    }

  indexSelectionFilter->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( indexSelectionFilter->GetOutput() );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension inputImage "
              << "outputImage operation radius" << std::endl;
    std::cerr << "  operation: " << std::endl;
    std::cerr << "    0. mean " << std::endl;
    std::cerr << "    1. min " << std::endl;
    std::cerr << "    2. max " << std::endl;
    std::cerr << "    3. variance " << std::endl;
    std::cerr << "    4. sigma " << std::endl;
    std::cerr << "    5. skewness " << std::endl;
    std::cerr << "    6. kurtosis " << std::endl;
    std::cerr << "    7. entropy " << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     CalculateStatisticsImage<2>( argc, argv );
     break;
   case 3:
     CalculateStatisticsImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
