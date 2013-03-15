#include "itkGDCMImageIO.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageFileReader.h"
#include "itkImageSeriesWriter.h"
#include "itkMetaDataDictionary.h"
#include "itkShiftScaleImageFilter.h"

int copy( int argc, char* argv[] )
{
  typedef signed short                          PixelType;
  typedef itk::Image< PixelType, 2>             ImageType;
  typedef itk::ImageFileReader<ImageType>       ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::GDCMImageIO                      DicomImageIOType;

  DicomImageIOType::Pointer dicomImageIO = DicomImageIOType::New();

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  reader2->SetImageIO( dicomImageIO );
  reader2->Update();

  itk::MetaDataDictionary & dicomImageDictionary = reader2->GetOutput()->GetMetaDataDictionary();
  reader->GetOutput()->SetMetaDataDictionary( dicomImageDictionary );

  typedef itk::ImageFileWriter<ImageType>       WriterType;

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( reader->GetOutput() );
  writer->SetImageIO( dicomImageIO );
  writer->SetFileName( argv[3] );
  writer->Update();

  return EXIT_SUCCESS;
}


int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " inputImageSlice";
    std::cerr << " referenceDicomSlice outputDicomSlice";

    return EXIT_FAILURE;
    }

  copy( argc, argv );
}
