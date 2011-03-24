#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericSeriesFileNames.h"

#include <string>

int main( int ac, char* av[] )
{

  if ( ac < 4 )
    {
    std::cerr << "Usage: " << av[0] << " directory format 3Dimage [startIndex] [endIndex] [readSeries=true]" << std::endl;
    return EXIT_FAILURE;
    }

  if( ac == 7 && atoi( av[6] ) )
    {
    typedef float PixelType;
    typedef itk::Image<PixelType, 3> ImageNDType;
    typedef itk::ImageSeriesReader<ImageNDType> ReaderType;
  
    std::string format = std::string( av[1] ) + std::string( "/" )
      + std::string( av[2] );
  
    // Get the filenames from the directory
    itk::NumericSeriesFileNames::Pointer names 
      = itk::NumericSeriesFileNames::New();
    names->SetSeriesFormat( format.c_str() );
    names->SetStartIndex( atoi( av[4] ) );
    names->SetEndIndex( atoi( av[5] ) );
    names->SetIncrementIndex( 1 );
    
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileNames( names->GetFileNames() );
//    std::cout << names << std::endl;
  
    try
      {
      if ( atoi( av[2] ) )
        {
        reader->ReverseOrderOn();
        }
      reader->Update();
//      reader->GetOutput()->Print( std::cout );
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
    }
  else
    {
    const char * outputDirectory = av[1];
    
    itksys::SystemTools::MakeDirectory( outputDirectory );
    
    typedef itk::Image<float, 3> ImageType;
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( av[3] );
    reader->Update();
    
    std::string format = std::string( av[1] ) + std::string( "/" )
      + std::string( av[2] );
  
    // Get the filenames from the directory
    itk::NumericSeriesFileNames::Pointer names 
      = itk::NumericSeriesFileNames::New();
    names->SetSeriesFormat( format.c_str() );
    names->SetStartIndex( 0 );
    names->SetEndIndex( 
      reader->GetOutput()->GetLargestPossibleRegion().GetSize()[2]-1 );
    names->SetIncrementIndex( 1 );
        
    typedef signed short OutputPixelType;
    const unsigned int OutputDimension = 2;

    typedef itk::Image<OutputPixelType, OutputDimension> Image2DType;
    
    typedef itk::ImageSeriesWriter<ImageType, Image2DType> SeriesWriterType;
    SeriesWriterType::Pointer writer = SeriesWriterType::New();
    writer->SetInput( reader->GetOutput() );
    writer->SetFileNames( names->GetFileNames() );
//    writer->SetMetaDataDictionaryArray( reader->GetMetaDataDictionaryArray() );
     
    try
      {
      writer->Update(); 
      }
    catch ( itk::ExceptionObject &ex )
      {
      std::cout << "Exception thrown hile writing the series" << std::endl; 
      std::cout << ex << std::endl;;
      return EXIT_FAILURE;
      }
    }
    


  return EXIT_SUCCESS;

}
