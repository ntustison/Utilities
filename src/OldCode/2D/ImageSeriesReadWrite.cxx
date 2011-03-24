#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkImage.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericSeriesFileNames.h"
#include "itkAnalyzeImageIO.h"
#include "itkNiftiImageIO.h"
#include "itkPNGImageIO.h"
#include "global.h"

int main( int argc, char ** argv )
{
  if( argc < 6 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " format firstSliceValue lastSliceValue spacing_1 spacing_2 spacing_3 outputImageFile " << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<PixelType, 3>           ImageType;
  typedef itk::ImageSeriesReader<ImageType>  ReaderType;
  typedef itk::ImageFileWriter<ImageType>  WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  typedef itk::NumericSeriesFileNames    NameGeneratorType;
  NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
  nameGenerator->SetSeriesFormat( argv[1] );
  nameGenerator->SetStartIndex( atoi( argv[2] ) );
  nameGenerator->SetEndIndex( atoi( argv[3] ) );
  nameGenerator->SetIncrementIndex( 1 );

  reader->SetImageIO( itk::PNGImageIO::New() );
  reader->SetFileNames( nameGenerator->GetFileNames() );
  reader->Update();

  ImageType::SpacingType spacing;
  spacing[0] = atof( argv[4] );
  spacing[1] = atof( argv[5] );
  spacing[2] = atof( argv[6] );

  ImageType::PointType origin;
  origin.Fill( 0.0 );

  reader->GetOutput()->SetSpacing( spacing ); 
  reader->GetOutput()->SetOrigin( origin );

  try 
    {
    writer->SetInput( reader->GetOutput() );
    writer->SetFileName( argv[7] );
    writer->Update(); 
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    } 

  return EXIT_SUCCESS;
}



