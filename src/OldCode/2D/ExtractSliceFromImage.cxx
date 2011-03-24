#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc != 5 )
    {
    std::cout << "Usage: " << argv[0] << " input3DImage output2DImage directions(i.e. 0, 1, 2) slice_number" << std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<PixelType, ImageDimension-1> SliceType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  ImageType::RegionType region;
  ImageType::RegionType::SizeType size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  size[atoi( argv[3] )] = 0;
  ImageType::IndexType index;
  index.Fill( 0 );
  index[atoi( argv[3] )] = atoi( argv[4] );
  region.SetIndex( index );
  region.SetSize( size );

  typedef itk::ExtractImageFilter<ImageType, SliceType> ExtracterType;
  ExtracterType::Pointer extracter = ExtracterType::New();
  extracter->SetInput( reader->GetOutput() );
  extracter->SetExtractionRegion( region );
  extracter->Update();

  typedef itk::ImageFileWriter<SliceType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( extracter->GetOutput() );
  writer->Update();


  return 0;
}
