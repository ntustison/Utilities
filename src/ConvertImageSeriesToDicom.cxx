#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"

#include <string>

int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cerr << "Usage: " << argv[0]
      << " inputDirectory inputSeriesFormat startIndex endIndex outputDirectory"
      << "[referenceDirectory] [referenceDicomSeriesFormat]... " << std::endl;
    exit( EXIT_FAILURE );
    }


  typedef signed short    PixelType;
  const unsigned int      Dimension = 3;

  typedef itk::Image< PixelType, Dimension >      ImageType;
  typedef itk::ImageSeriesReader< ImageType >     ReaderType;

  typedef itk::GDCMImageIO                        ImageIOType;

  typedef itk::GDCMSeriesFileNames                NamesGeneratorType;
  NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();

  ImageIOType::Pointer gdcmIO = ImageIOType::New();

  std::string format = std::string( argv[1] ) + std::string( "/" )
    + std::string( argv[2] );

  // Get the filenames from the directory
  itk::NumericSeriesFileNames::Pointer names =
    itk::NumericSeriesFileNames::New();
  names->SetSeriesFormat( format.c_str() );
  names->SetStartIndex( atoi( argv[3] ) );
  names->SetEndIndex( atoi( argv[4] ) );
  names->SetIncrementIndex( 1 );

  const ReaderType::FileNamesContainer & filenames = names->GetFileNames();
  unsigned int numberOfFilenames =  filenames.size();
  std::cout << numberOfFilenames << std::endl;
  for(unsigned int fni = 1; fni<=numberOfFilenames; fni++)
    {
    std::cout << "filename # " << fni << " = ";
    std::cout << filenames[fni-1] << std::endl;
    }

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileNames( filenames );
  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject &excp)
    {
    std::cerr << "Exception thrown while writing the image" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  itksys::SystemTools::MakeDirectory( argv[5] );
  typedef signed short    OutputPixelType;
  const unsigned int      OutputDimension = 2;

  typedef itk::Image< OutputPixelType, OutputDimension >    Image2DType;

  typedef itk::ImageSeriesWriter<
                             ImageType, Image2DType >  SeriesWriterType;
  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();

  seriesWriter->SetInput( reader->GetOutput() );
  seriesWriter->SetImageIO( gdcmIO );

  std::string outformat = std::string( argv[5] ) + std::string( "/" )
    + std::string( argv[2] );
  itk::NumericSeriesFileNames::Pointer outnames = itk::NumericSeriesFileNames::New();
  outnames->SetSeriesFormat( outformat.c_str() );
  outnames->SetStartIndex( atoi( argv[3] ) );
  outnames->SetEndIndex( atoi( argv[4] ) );
  outnames->SetIncrementIndex( 1 );

  seriesWriter->SetFileNames( outnames->GetFileNames() );

  // Replace meta data dictionary values if reference dicom series is specified
  // by user
  if( argc > 6 )
    {
    std::string refFormat = std::string( argv[6] ) + std::string( "/" )
      + std::string( argv[7] );

    // Get the filenames from the directory
    itk::NumericSeriesFileNames::Pointer refNames =
      itk::NumericSeriesFileNames::New();
    refNames->SetSeriesFormat( refFormat.c_str() );
    refNames->SetStartIndex( atoi( argv[3] ) );
    refNames->SetEndIndex( atoi( argv[4] ) );
    refNames->SetIncrementIndex( 1 );

    const ReaderType::FileNamesContainer & refFilenames =
      refNames->GetFileNames();

    for(unsigned int fni = 1; fni<=refFilenames.size(); fni++)
      {
      std::cout << "reffilename # " << fni << " = ";
      std::cout << refFilenames[fni-1] << std::endl;
      }

    ReaderType::Pointer refReader = ReaderType::New();
    refReader->SetFileNames( refFilenames );
    refReader->SetImageIO( gdcmIO );
    refReader->Update();

//    typedef itk::MetaDataDictionary DictionaryType;
//    DictionaryType *dictionary = (*(refReader->GetMetaDataDictionaryArray()))[0];
//
//    ReaderType::DictionaryArrayType outputArray;
//
//    for( unsigned int i = 0; i < refFilenames.size(); i++ )
//      {
//      dictionary = (*(refReader->GetMetaDataDictionaryArray()))[i];
//      outputArray.push_back( dictionary );
//      }
//    seriesWriter->SetMetaDataDictionaryArray( &outputArray );

    seriesWriter->SetMetaDataDictionaryArray(
      refReader->GetMetaDataDictionaryArray() );
    try
      {
      seriesWriter->DebugOn();
      seriesWriter->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Exception thrown while writing the series " << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
      }

    }
  else
    {
    seriesWriter->SetMetaDataDictionaryArray(
      reader->GetMetaDataDictionaryArray() );

    try
      {
      seriesWriter->DebugOn();
      seriesWriter->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Exception thrown while writing the series " << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
      }

    }

  return EXIT_SUCCESS;

}
