#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"

int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cout
      << argv[0] << " Aimage Bimage T1image inversionTime outputImage"
      << std::endl;
    exit( 1 );
    }

  typedef float PixelType;
  typedef itk::Image<PixelType, 2> ImageType;

  float inversionTime = atof( argv[4] );

  ImageType::Pointer A;
  {
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  A = reader->GetOutput();
  A->Update();
  A->DisconnectPipeline();
  }

  ImageType::Pointer B;
  {
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );

  B = reader->GetOutput();
  B->Update();
  B->DisconnectPipeline();
  }

  ImageType::Pointer T1;
  {
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[3] );

  T1 = reader->GetOutput();
  T1->Update();
  T1->DisconnectPipeline();
  }

  ImageType::Pointer outputImage = ImageType::New();
  outputImage->CopyInformation( A );
  outputImage->SetRegions( A->GetRequestedRegion() );
  outputImage->Allocate();
  outputImage->FillBuffer( 0.0 );

  itk::ImageRegionIteratorWithIndex<ImageType> It( outputImage, outputImage->GetRequestedRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    ImageType::IndexType index = It.GetIndex();

    float Axy = A->GetPixel( index );
    float Bxy = B->GetPixel( index );
    float T1xy = T1->GetPixel( index );

    PixelType value = Axy - Bxy * vcl_exp( -inversionTime / T1xy );
    It.Set( value );
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[5] );
  writer->SetInput( outputImage );
  writer->Update();

  return 0;
}

