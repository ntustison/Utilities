#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericSeriesFileNames.h"

#include "global.h"

int main( int ac, char* av[] )
{

  if ( ac != 5 )
    {
    std::cerr << "Usage: " << av[0] << " directory/prefix outputImage startIndex endIndex" << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<PixelType, 3> ImageNDType;
  typedef itk::ImageSeriesReader<ImageNDType> ReaderType;

  // Get the filenames from the directory
  itk::NumericSeriesFileNames::Pointer names 
    = itk::NumericSeriesFileNames::New();
  names->SetSeriesFormat( av[1] );
  names->SetStartIndex( atoi( av[3] ) );
  names->SetEndIndex( atoi( av[4] ) );
  names->SetIncrementIndex( 1 );
  
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileNames( names->GetFileNames() );
  std::cout << names << std::endl;

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
  writer->SetFileName( av[2] );
  writer->SetInput( reader->GetOutput() );
  writer->Update();


  return EXIT_SUCCESS;

}
