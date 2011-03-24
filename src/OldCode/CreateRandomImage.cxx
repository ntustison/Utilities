#include <stdio.h>

#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

#include <vector>

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc != 3 + 3*ImageDimension )
    {
    std::cout << "Usage: " << argv[0] << " outputImage pixelValueSupremum " 
              << " origin[0] ... origin[n]" 
              << " spacing[0] ... spacing[n]" 
              << " size[0] ... size[n]" << std::endl;
    std::cout << "   Returns an image with pixel values in the range [0, pixelValueSupremum]. " << std::endl;
    exit( 1 );
    }
  
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  ImageType::Pointer image = ImageType::New();

  ImageType::SpacingType spacing;
  ImageType::PointType origin;
  ImageType::RegionType::SizeType size;

  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    size[i] = atoi( argv[3+2*ImageDimension + i] );
    spacing[i] = atof( argv[3+ImageDimension + i] ); 
    origin[i] = atof( argv[3+i] );
    }

  image->SetRegions( size );
  image->SetOrigin( origin );
  image->SetSpacing( spacing );
  image->Allocate();
  
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
  GeneratorType::Pointer generator = GeneratorType::New();
  
  PixelType a;
  int t_int;
  unsigned int t_uint;
  long t_long;
  unsigned long t_ulong;
  short t_sh;
  unsigned short t_ush;
  char t_ch;
  unsigned char t_uch;
  
  std::vector<unsigned int> count( atoi( argv[2] )+1 );
  for ( unsigned int i = 0; i < count.size(); i++ )
    {
    count[i] = 0;  
    }

  itk::ImageRegionIterator<ImageType> It( image, image->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( typeid( a ) == typeid( t_int ) ||
         typeid( a ) == typeid( t_uint ) ||  
         typeid( a ) == typeid( t_long ) ||  
         typeid( a ) == typeid( t_ulong ) ||  
         typeid( a ) == typeid( t_sh ) ||  
         typeid( a ) == typeid( t_ush ) ||  
         typeid( a ) == typeid( t_ch ) ||  
         typeid( a ) == typeid( t_uch ) )
      {     
      It.Set( generator->GetIntegerVariate( atoi( argv[2] ) ) );
      count[It.Get()]++;
      }
    else
      {
      double lim = static_cast<double>( atof( argv[2] ) );
      PixelType val = static_cast<PixelType>( 
        generator->GetVariateWithClosedRange( lim ) );
      It.Set( val );
      }
    }

  for ( unsigned int i = 0; i < count.size(); i++ )
    {
    std::cout << "count[" << i << "] = " << count[i] << std::endl;  
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[1] );
  writer->SetInput( image );
  writer->Update();

 return 0;
}
