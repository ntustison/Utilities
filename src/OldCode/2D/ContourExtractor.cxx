#include "itkImageFileReader.h"
#include "itkContourExtractor2DImageFilter.h"

#include "fstream.h"

#include "global.h"

int main( int argc, char *argv[] ) 
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage contourValue outputFile.txt " << std::endl; 
    exit( 1 );
    }  

  const unsigned int ImageDimension = 2;

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  typedef itk::ContourExtractor2DImageFilter<ImageType> ExtractorType;
  ExtractorType::Pointer extractor = ExtractorType::New();
  extractor->SetInput( reader->GetOutput() );
  extractor->SetContourValue( atoi( argv[2] ) );
  extractor->VertexConnectHighPixelsOff();
  extractor->ReverseContourOrientationOff();
  extractor->Update();

  ofstream str( argv[3] );  

  str << "0 0 0 0" << std::endl;
  for ( unsigned int i = 0; i < extractor->GetNumberOfOutputs(); i++ )
    { 
    for ( unsigned int j = 0; j < extractor->GetOutput( i )->GetVertexList()->Size(); j++ )
      {
      ExtractorType::VertexType V = extractor->GetOutput( i )->GetVertexList()->GetElement( j );
      str << V[0] << ' ' << V[1] << " 0 " << i+1 << std::endl;
      }   
    }
  str << "0 0 0 0" << std::endl;
  str.close();

  return 0;
}
