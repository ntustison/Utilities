#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRegionalMaximaImageFilter.h"
#include "itkRegionalMinimaImageFilter.h"

template <unsigned int ImageDimension>
int RegionalMinMax( int argc, char *argv[] )
{
 
  typedef float PixelType;
   
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  if( atoi( argv[4] ) )
    {
    typedef itk::RegionalMaximaImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
  
    filter->SetInput( reader->GetOutput() );
    filter->SetBackgroundValue( 0 );
    filter->SetForegroundValue( 1 );
  
    filter->Update();

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( filter->GetOutput() );
    writer->SetFileName( argv[3] );                                          
    writer->Update();
    }
  else
    {
    typedef itk::RegionalMinimaImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
  
    filter->SetInput( reader->GetOutput() );
    filter->SetBackgroundValue( 0 );
    filter->SetForegroundValue( 1 );
  
    filter->Update();

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( filter->GetOutput() );
    writer->SetFileName( argv[3] );                                          
    writer->Update();
    }
    
  
    
  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage outputImage doMinimum" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     RegionalMinMax<2>( argc, argv );
     break;
   case 3:
     RegionalMinMax<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

