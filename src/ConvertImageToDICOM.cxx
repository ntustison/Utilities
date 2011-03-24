#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageSeriesWriter.h"
#include "itkGDCMImageIO.h"
#include "itkNumericSeriesFileNames.h"
#include "itkRGBPixel.h"

#include <string>

template<typename PixelType>
int ConvertImageToDICOM( int argc, char *argv[] )
{
  typedef itk::Image<PixelType, 3> ImageType;
  typedef itk::Image<PixelType, 2> Image2DType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  itksys::SystemTools::MakeDirectory( argv[2] );
  
  std::string filePrefix = std::string( argv[2] ) 
    + std::string( "/" ) + std::string( argv[3] )
    + std::string( "%03d" );

  typedef itk::GDCMImageIO ImageIOType; 
  typedef itk::NumericSeriesFileNames NamesGeneratorType; 
 
  typename ImageIOType::Pointer gdcmIO = ImageIOType::New(); 
  typename NamesGeneratorType::Pointer namesGenerator 
    = NamesGeneratorType::New(); 
  namesGenerator->SetSeriesFormat( filePrefix.c_str() );
  namesGenerator->SetStartIndex( 0 );
  namesGenerator->SetEndIndex( 
    reader->GetOutput()->GetLargestPossibleRegion().GetSize()[2] - 1 );
       

  typedef itk::ImageSeriesWriter<ImageType, Image2DType> SeriesWriterType;
  typename SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
  seriesWriter->SetInput( reader->GetOutput() );
  seriesWriter->SetImageIO( gdcmIO );
  seriesWriter->SetFileNames( namesGenerator->GetFileNames() );
  seriesWriter->SetMetaDataDictionary( 
    reader->GetMetaDataDictionary() );

  try 
    { 
    seriesWriter->Update(); 
    } 
  catch( itk::ExceptionObject & excp ) 
    { 
    std::cerr << "Exception thrown while writing the series " << std::endl; 
    std::cerr << excp << std::endl; 
    return EXIT_FAILURE; 
    } 

  return 0;
}

int main( int argc, char *argv[] )        
{
  if ( argc < 4  )
    {
    std::cout << "Usage: " << argv[0] 
      << " inputImage outputDirectory prefix [isRGB]" << std::endl;
    exit( 1 );
    }

  typedef short PixelType;
  typedef itk::RGBPixel<char> RGBPixelType;
    
  bool isRGB = ( argc > 4 ) ? static_cast<bool>( atoi( argv[4] ) ) : false;    

  if( isRGB )
    {
    ConvertImageToDICOM<RGBPixelType>( argc, argv ); 
    }
  else
    {
    ConvertImageToDICOM<PixelType>( argc, argv ); 
    } 
}      

