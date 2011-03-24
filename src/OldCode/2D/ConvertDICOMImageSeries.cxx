#include "itkDICOMImageIO2Factory.h"
#include "itkDICOMImageIO2.h"
#include "itkImageSeriesReader.h"
#include "itkDICOMSeriesFileNames.h"
#include "itkImageFileWriter.h"

int main( int ac, char* av[] )
{

  if ( ac != 4 )
    {
    std::cerr << "Usage: " << av[0] << " DicomDirectory ReverseOrder(0/1) OutputImage\n";
    return EXIT_FAILURE;
    }

  typedef itk::Image<int, 3> ImageNDType;
  typedef itk::ImageSeriesReader<ImageNDType> ReaderType;

  itk::DICOMImageIO2::Pointer io = itk::DICOMImageIO2::New();

  // Get the DICOM filenames from the directory
  itk::DICOMSeriesFileNames::Pointer names = itk::DICOMSeriesFileNames::New();
  names->SetDirectory( av[1] );
  
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileNames( names->GetFileNames() );
  reader->SetImageIO( io );
  std::cout << names;

  try
    {
    if ( atoi( av[2] ) )
      {
      reader->ReverseOrderOn();
      }
    reader->Update();
    reader->GetOutput()->Print( std::cout );
    }
  catch ( itk::ExceptionObject &ex )
    {
    std::cout << ex;
    return EXIT_FAILURE;
    }

  typedef itk::ImageFileWriter<ImageNDType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( av[3] );
  writer->SetInput( reader->GetOutput() );
  writer->Update();


  return EXIT_SUCCESS;

}
