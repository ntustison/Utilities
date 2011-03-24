#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkConstantPadImageFilter.h"

#include "global.h"
#include <string>

int main( unsigned int argc, char *argv[] )
{
  if ( argc < 2 + ImageDimension*2 )
    {
    std::cout << "Usage: " << argv[0] 
              << " inputImage outputImage " 
              << "padSizeBefore[0] padSizeAfter[0] ... "
              << "padSizeBefore[n-1] padSizeAfter[n-1] [value]" << std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  
  unsigned long lowerBound[ImageDimension];
  unsigned long upperBound[ImageDimension];
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    lowerBound[i] = atoi( argv[3 + 2*i] );
    upperBound[i] = atoi( argv[3 + (2*i+1)] );
    } 

  typedef itk::ConstantPadImageFilter<ImageType, ImageType> PadderType;
  PadderType::Pointer padder = PadderType::New();
  padder->SetInput( reader->GetOutput() );
  padder->SetPadLowerBound( lowerBound );
  padder->SetPadUpperBound( upperBound );
  if ( argc > 3 + 2*ImageDimension )
    {
    padder->SetConstant( static_cast<PixelType>( atof( argv[3 + 2*ImageDimension] ) ) );
    }  
  padder->Update();

  ImageType::PointType origin = padder->GetOutput()->GetOrigin();
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin[i] -= ( padder->GetOutput()->GetSpacing()[i] 
      * static_cast<double>( lowerBound[i] ) ); 
    }
  padder->GetOutput()->SetOrigin( origin );
  
  
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( padder->GetOutput() );
  writer->SetFileName( argv[2] );                                          
  writer->Update();
    
  return 0;
}
