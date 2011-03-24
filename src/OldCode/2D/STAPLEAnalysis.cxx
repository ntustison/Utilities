#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSTAPLEImageFilter.h"
#include "itkTimeProbe.h"

#include <string>
#include <fstream.h>

#include "global.h"

int main( unsigned int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " outputImageFile label segmentationImage1 ... segmentationImageN" << std::endl;     
    exit( 0 );
    }   

  itk::TimeProbe timer;
  timer.Start();

  typedef float RealType;

  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;

  typedef itk::STAPLEImageFilter<LabelImageType, RealImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetForegroundValue( static_cast<LabelImageType::PixelType>( atoi( argv[2] ) ) );

  for ( unsigned int i = 3; i < argc; i++ )
    {
    typedef itk::ImageFileReader<LabelImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[i] );
    reader->Update();

    filter->SetInput( i-3, reader->GetOutput() );
    }
  filter->Update();
 
  for ( unsigned int i = 3; i < argc; i++ )
    {
    std::cout << argv[i] << std::endl;
    std::cout << "   Specificity: " << filter->GetSpecificity( i-3 ) << std::endl;
    std::cout << "   Sensitivity: " << filter->GetSensitivity( i-3 ) << std::endl;
    }

  typedef itk::ImageFileWriter<RealImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[1] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();
  
  timer.Stop();
}
