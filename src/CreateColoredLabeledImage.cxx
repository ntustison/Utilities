#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkRGBPixel.h"
#include "itkRGBAPixel.h"

#include "itkRedColormapFunction.h"
#include "itkGreenColormapFunction.h"
#include "itkBlueColormapFunction.h"
#include "itkGreyColormapFunction.h"
#include "itkHotColormapFunction.h"
#include "itkCoolColormapFunction.h"
#include "itkSpringColormapFunction.h"
#include "itkSummerColormapFunction.h"
#include "itkAutumnColormapFunction.h"
#include "itkWinterColormapFunction.h"
#include "itkCopperColormapFunction.h"
#include "itkHSVColormapFunction.h"
#include "itkJetColormapFunction.h"
#include "itkCustomColormapFunction.h"
#include "itkOverUnderColormapFunction.h"

#include "itkScalarToRGBColormapImageFilter.h"

#include <fstream>
#include <sstream>
#include <string>

template <unsigned int ImageDimension>
int ConvertScalarImageToRGB( int argc, char *argv[] )
{
  typedef unsigned char PixelType;
  typedef itk::RGBPixel<unsigned char> RGBPixelType;
//  typedef itk::RGBAPixel<unsigned char> RGBPixelType;

  typedef float RealType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<float, ImageDimension> RealImageType;
  typedef itk::Image<RGBPixelType, ImageDimension> RGBImageType;

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::RescaleIntensityImageFilter<RealImageType, ImageType> RescalerType;
  typename RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetInput( reader->GetOutput() );
  rescaler->SetOutputMinimum( static_cast<PixelType>( 0 ) );
  rescaler->SetOutputMaximum( static_cast<PixelType>( 255 ) );
  rescaler->Update();

  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;
  typename LabelImageType::Pointer maskImage = LabelImageType::New();
  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
  typename LabelReaderType::Pointer labelImageReader = LabelReaderType::New();
  labelImageReader->SetFileName( argv[3] );
  labelImageReader->Update();

  std::string colormapString( argv[5] );

  typedef itk::ScalarToRGBColormapImageFilter<LabelImageType,
    RGBImageType> RGBFilterType;
  typename RGBFilterType::Pointer rgbfilter = RGBFilterType::New();
  rgbfilter->SetInput( labelImageReader->GetOutput() );

  if ( colormapString == "red" )
    {
    rgbfilter->SetColormap( RGBFilterType::Red );
    }
  else if ( colormapString == "green"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Green );
    }
  else if ( colormapString == "blue"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Blue );
    }
  else if ( colormapString == "grey"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Grey );
    }
  else if ( colormapString == "cool"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Cool );
    }
  else if ( colormapString == "hot"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Hot );
    }
  else if ( colormapString == "spring"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Spring );
    }
  else if ( colormapString == "autumn"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Autumn );
    }
  else if ( colormapString == "winter"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Winter );
    }
  else if ( colormapString == "copper"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Copper );
    }
  else if ( colormapString == "summer"  )
    {
    rgbfilter->SetColormap( RGBFilterType::Summer );
    }
  else if ( colormapString == "jet"  )
    {
//    rgbfilter->SetColormap( RGBFilterType::Jet );
    typedef itk::Function::JetColormapFunction<typename LabelImageType::PixelType,
      typename RGBImageType::PixelType> ColormapType;
    typename ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if ( colormapString == "hsv"  )
    {
//    rgbfilter->SetColormap( RGBFilterType::HSV );
    typedef itk::Function::HSVColormapFunction<typename LabelImageType::PixelType,
      typename RGBImageType::PixelType> ColormapType;
    typename ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if ( colormapString == "overunder"  )
    {
    rgbfilter->SetColormap( RGBFilterType::OverUnder );
    }
  else if ( colormapString == "custom"  )
    {
    typedef itk::Function::CustomColormapFunction<typename LabelImageType::PixelType,
      typename RGBImageType::PixelType> ColormapType;
    typename ColormapType::Pointer colormap = ColormapType::New();

    std::ifstream str( argv[6] );
    std::string line;

    // Get red values
    {
    std::getline( str, line );
    std::istringstream iss( line );
    float value;
    typename ColormapType::ChannelType channel;
    while ( iss >> value )
      {
      channel.push_back( value );
      }
    colormap->SetRedChannel( channel );
    }

    // Get green values
    {
    std::getline( str, line );
    std::istringstream iss( line );
    float value;
    typename ColormapType::ChannelType channel;
    while ( iss >> value )
      {
      channel.push_back( value );
      }
    colormap->SetGreenChannel( channel );
    }
    // Get blue values
    {
    std::getline( str, line );
    std::istringstream iss( line );
    float value;
    typename ColormapType::ChannelType channel;
    while ( iss >> value )
      {
      channel.push_back( value );
      }
    colormap->SetBlueChannel( channel );
    }
    rgbfilter->SetColormap( colormap );
    }

  try
    {
    rgbfilter->Update();
    }
  catch (...)
    {
    return EXIT_FAILURE;
    }

  itk::ImageRegionIterator<LabelImageType> ItM( labelImageReader->GetOutput(),
    labelImageReader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<RGBImageType> ItC( rgbfilter->GetOutput(),
    rgbfilter->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> ItS( rescaler->GetOutput(),
    rescaler->GetOutput()->GetLargestPossibleRegion() );

  ItM.GoToBegin();
  ItC.GoToBegin();
  ItS.GoToBegin();

  while( !ItM.IsAtEnd() )
    {
    if( ItM.Get() == 0 )
      {
      RGBPixelType rgbpixel;
      rgbpixel.Fill( static_cast<typename RGBPixelType::ComponentType>( ItS.Get() ) );
      ItC.Set( rgbpixel );
      }
    ++ItM;
    ++ItC;
    ++ItS;
    }

  typedef itk::ImageFileWriter<RGBImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( rgbfilter->GetOutput() );
  writer->SetFileName( argv[4] );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension inputImage labelImage outputImage "
      << "colormap [customColormapFile] " << std::endl;
    std::cerr << "  Note, output images are automatically rescaled to [0, 255]" << std::endl;
    std::cerr << "  Possible colormaps: grey, red, green, blue, copper, jet, hsv, ";
    std::cerr << "spring, summer, autumn, winter, hot, cool, overunder, custom" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     ConvertScalarImageToRGB<2>( argc, argv );
     break;
   case 3:
     ConvertScalarImageToRGB<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

