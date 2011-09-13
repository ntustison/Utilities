#include "itkGDCMImageIOFactory.h"
#include "itkGDCMImageIO.h"
#include "itkImageSeriesReader.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageFileWriter.h"

template<unsigned int ImageDimension>
int dicom( int argc, char* argv[] )
{
  typedef itk::Image<int, ImageDimension> ImageNDType;
  typedef itk::ImageSeriesReader<ImageNDType> ReaderType;

  itk::GDCMImageIO::Pointer io = itk::GDCMImageIO::New();

  // Get the GDCM filenames from the directory
  itk::GDCMSeriesFileNames::Pointer names = itk::GDCMSeriesFileNames::New();
  names->SetDirectory( argv[2] );

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileNames( names->GetInputFileNames() );
  reader->SetImageIO( io );
  std::cout << names;

  try
    {
    reader->Update();
    reader->GetOutput()->Print( std::cout );
    }
  catch ( itk::ExceptionObject &ex )
    {
    std::cout << ex;
    return EXIT_FAILURE;
    }

  typedef itk::ImageFileWriter<ImageNDType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( reader->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}


int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " ImageDimension DicomDirectory OutputImage\n";
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     dicom<2>( argc, argv );
     break;
   case 3:
     dicom<3>( argc, argv );
     break;
   case 4:
     dicom<4>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
