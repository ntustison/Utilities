#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkScalarToRGBColormapImageFilter.h"

#include "itkScalarToRGBRedColormapFunctor.h"
#include "itkScalarToRGBGreenColormapFunctor.h"
#include "itkScalarToRGBBlueColormapFunctor.h"
#include "itkScalarToRGBGreyColormapFunctor.h"
#include "itkScalarToRGBHotColormapFunctor.h"
#include "itkScalarToRGBCoolColormapFunctor.h"
#include "itkScalarToRGBSpringColormapFunctor.h"
#include "itkScalarToRGBSummerColormapFunctor.h"
#include "itkScalarToRGBAutumnColormapFunctor.h"
#include "itkScalarToRGBWinterColormapFunctor.h"
#include "itkScalarToRGBCopperColormapFunctor.h"
#include "itkScalarToRGBJetColormapFunctor.h"
#include "itkScalarToRGBHSVColormapFunctor.h"
#include "itkScalarToRGBTustiColormapFunctor.h"
#include "itkScalarToRGBCustomColormapFunctor.h"

#include "itkStatisticsImageFilter.h"

#include <fstream.h>
#include <sstream>
#include <string>

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage outputImage colormap [customColormapFile]" << std::endl;
    std::cout << "  Possible colormaps: grey, red, green, blue, copper, jet, hsv, ";
    std::cout << "spring, summer, autumn, winter, hot, cool, tusti, custom" << std::endl;
    exit( 1 );
    }
 
  typedef itk::RGBPixel<float> RGBPixelType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<float, ImageDimension> RealImageType;
  typedef itk::Image<RGBPixelType, ImageDimension> RGBImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::StatisticsImageFilter<ImageType> StatsFilterType;
  StatsFilterType::Pointer stats = StatsFilterType::New();
  stats->SetInput( reader->GetOutput() );
  stats->Update();

  std::string colormap( argv[3] );

  typedef itk::ScalarToRGBColormapImageFilter<ImageType, 
    RGBImageType> RGBFilterType;
  RGBFilterType::Pointer rgbfilter = RGBFilterType::New();
  rgbfilter->SetInput( reader->GetOutput() );
  
  if ( colormap == "red" )
    {
    typedef itk::Functor::ScalarToRGBRedColormapFunctor<ImageType::PixelType, 
      RGBImageType::PixelType::ValueType> ColormapType;
    ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if ( colormap == "green"  )  
    {
    typedef itk::Functor::ScalarToRGBGreenColormapFunctor<ImageType::PixelType, 
      RGBImageType::PixelType::ValueType> ColormapType;
    ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if ( colormap == "blue"  )  
    {
    typedef itk::Functor::ScalarToRGBBlueColormapFunctor<ImageType::PixelType, 
      RGBImageType::PixelType::ValueType> ColormapType;
    ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if ( colormap == "grey"  )  
    {
    typedef itk::Functor::ScalarToRGBGreyColormapFunctor<ImageType::PixelType, 
      RGBImageType::PixelType::ValueType> ColormapType;
    ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if ( colormap == "cool"  )  
    {
    typedef itk::Functor::ScalarToRGBCoolColormapFunctor<ImageType::PixelType, 
      RGBImageType::PixelType::ValueType> ColormapType;
    ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if ( colormap == "hot"  )  
    {
    typedef itk::Functor::ScalarToRGBHotColormapFunctor<ImageType::PixelType, 
      RGBImageType::PixelType::ValueType> ColormapType;
    ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if ( colormap == "spring"  )  
    {
    typedef itk::Functor::ScalarToRGBSpringColormapFunctor<ImageType::PixelType, 
      RGBImageType::PixelType::ValueType> ColormapType;
    ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if ( colormap == "autumn"  )  
    {
    typedef itk::Functor::ScalarToRGBAutumnColormapFunctor<ImageType::PixelType, 
      RGBImageType::PixelType::ValueType> ColormapType;
    ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if ( colormap == "winter"  )  
    {
    typedef itk::Functor::ScalarToRGBWinterColormapFunctor<ImageType::PixelType, 
      RGBImageType::PixelType::ValueType> ColormapType;
    ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if ( colormap == "copper"  )  
    {
    typedef itk::Functor::ScalarToRGBCopperColormapFunctor<ImageType::PixelType, 
      RGBImageType::PixelType::ValueType> ColormapType;
    ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if ( colormap == "summer"  )  
    {
    typedef itk::Functor::ScalarToRGBSummerColormapFunctor<ImageType::PixelType, 
      RGBImageType::PixelType::ValueType> ColormapType;
    ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if ( colormap == "jet"  )  
    {
    typedef itk::Functor::ScalarToRGBJetColormapFunctor<ImageType::PixelType, 
      RGBImageType::PixelType::ValueType> ColormapType;
    ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if ( colormap == "hsv"  )  
    {
    typedef itk::Functor::ScalarToRGBHSVColormapFunctor<ImageType::PixelType, 
      RGBImageType::PixelType::ValueType> ColormapType;
    ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if ( colormap == "tusti"  )  
    {
    typedef itk::Functor::ScalarToRGBTustiColormapFunctor<ImageType::PixelType, 
      RGBImageType::PixelType::ValueType> ColormapType;
    ColormapType::Pointer colormap = ColormapType::New();
    rgbfilter->SetColormap( colormap );
    }
  else if ( colormap == "custom"  )  
    {
    typedef itk::Functor::ScalarToRGBCustomColormapFunctor<ImageType::PixelType, 
      RGBImageType::PixelType::ValueType> ColormapType;
    ColormapType::Pointer colormap = ColormapType::New();

    ifstream str( argv[4] );
    std::string line;
    
    // Get red values 
    {
    std::getline( str, line );
    std::istringstream iss( line );
    float value;
    ColormapType::ChannelType channel;
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
    ColormapType::ChannelType channel;
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
    ColormapType::ChannelType channel;
    while ( iss >> value )
      {
      channel.push_back( value );
      }
    colormap->SetBlueChannel( channel );      
    }   
    rgbfilter->SetColormap( colormap );
    }
    
  rgbfilter->GetColormap()->SetMinimumRGBComponentValue( 0 );
  rgbfilter->GetColormap()->SetMaximumRGBComponentValue( 255 );
  rgbfilter->GetColormap()->SetMinimumInputValue( stats->GetMinimum() );
  rgbfilter->GetColormap()->SetMaximumInputValue( stats->GetMaximum() );
  rgbfilter->Update();

  typedef itk::ImageFileWriter<RGBImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( rgbfilter->GetOutput() );
  writer->SetFileName( argv[2] );
  writer->Update();

//  if ( argc > 4 )
//    {
//    typedef itk::Image<float, 2> ColormapImageType; 
//    ColormapImageType::Pointer colormapImage = ColormapImageType::New();
//    ColormapImageType::SizeType size;
//    size[0] = 255;
//    size[1] = 30;
//    colormapImage->SetRegions( size );
//    colormapImage->Allocate();
//    itk::ImageRegionIteratorWithIndex<ColormapImageType> It( colormapImage, 
//      colormapImage->GetLargestPossibleRegion() );
//    for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
//      {
//      It.Set( It.GetIndex()[0] ); 
//      }   
//    typedef itk::RescaleIntensityImageFilter
//      <ColormapImageType, ColormapImageType> RescalerType;
//    RescalerType::Pointer rescaler = RescalerType::New();
//    rescaler->SetInput( colormapImage );
//    rescaler->SetOutputMinimum( stats->GetMinimum() );
//    rescaler->SetOutputMaximum( stats->GetMaximum() );
//    rescaler->Update();
//      
//    typedef itk::Image<RGBPixelType, 2> RGBImageType;
//    typedef itk::ScalarToRGBColormapImageFilter<ColormapImageType, 
//      RGBImageType> RGBFilterType;
//    RGBFilterType::Pointer rgbfilter = RGBFilterType::New();
//    rgbfilter->SetInput( rescaler->GetOutput() );
//      
//    if ( colormap == "red" )
//      {
//      typedef itk::Functor::ScalarToRGBRedColormapFunctor<ColormapImageType::PixelType, 
//        RGBImageType::PixelType::ValueType> ColormapType;
//      ColormapType::Pointer colormap = ColormapType::New();
//      rgbfilter->SetColormap( colormap );
//      }
//    else if ( colormap == "green"  )  
//      {
//      typedef itk::Functor::ScalarToRGBGreenColormapFunctor<ColormapImageType::PixelType, 
//        RGBImageType::PixelType::ValueType> ColormapType;
//      ColormapType::Pointer colormap = ColormapType::New();
//      rgbfilter->SetColormap( colormap );
//      }
//    else if ( colormap == "blue"  )  
//      {
//      typedef itk::Functor::ScalarToRGBBlueColormapFunctor<ColormapImageType::PixelType, 
//        RGBImageType::PixelType::ValueType> ColormapType;
//      ColormapType::Pointer colormap = ColormapType::New();
//      rgbfilter->SetColormap( colormap );
//      }
//    else if ( colormap == "grey"  )  
//      {
//      typedef itk::Functor::ScalarToRGBGreyColormapFunctor<ColormapImageType::PixelType, 
//        RGBImageType::PixelType::ValueType> ColormapType;
//      ColormapType::Pointer colormap = ColormapType::New();
//      rgbfilter->SetColormap( colormap );
//      }
//    else if ( colormap == "cool"  )  
//      {
//      typedef itk::Functor::ScalarToRGBCoolColormapFunctor<ColormapImageType::PixelType, 
//        RGBImageType::PixelType::ValueType> ColormapType;
//      ColormapType::Pointer colormap = ColormapType::New();
//      rgbfilter->SetColormap( colormap );
//      }
//    else if ( colormap == "hot"  )  
//      {
//      typedef itk::Functor::ScalarToRGBHotColormapFunctor<ColormapImageType::PixelType, 
//        RGBImageType::PixelType::ValueType> ColormapType;
//      ColormapType::Pointer colormap = ColormapType::New();
//      rgbfilter->SetColormap( colormap );
//      }
//    else if ( colormap == "spring"  )  
//      {
//      typedef itk::Functor::ScalarToRGBSpringColormapFunctor<ColormapImageType::PixelType, 
//        RGBImageType::PixelType::ValueType> ColormapType;
//      ColormapType::Pointer colormap = ColormapType::New();
//      rgbfilter->SetColormap( colormap );
//      }
//    else if ( colormap == "autumn"  )  
//      {
//      typedef itk::Functor::ScalarToRGBAutumnColormapFunctor<ColormapImageType::PixelType, 
//        RGBImageType::PixelType::ValueType> ColormapType;
//      ColormapType::Pointer colormap = ColormapType::New();
//      rgbfilter->SetColormap( colormap );
//      }
//    else if ( colormap == "winter"  )  
//      {
//      typedef itk::Functor::ScalarToRGBWinterColormapFunctor<ColormapImageType::PixelType, 
//        RGBImageType::PixelType::ValueType> ColormapType;
//      ColormapType::Pointer colormap = ColormapType::New();
//      rgbfilter->SetColormap( colormap );
//      }
//    else if ( colormap == "copper"  )  
//      {
//      typedef itk::Functor::ScalarToRGBCopperColormapFunctor<ColormapImageType::PixelType, 
//        RGBImageType::PixelType::ValueType> ColormapType;
//      ColormapType::Pointer colormap = ColormapType::New();
//      rgbfilter->SetColormap( colormap );
//      }
//    else if ( colormap == "summer"  )  
//      {
//      typedef itk::Functor::ScalarToRGBSummerColormapFunctor<ColormapImageType::PixelType, 
//        RGBImageType::PixelType::ValueType> ColormapType;
//      ColormapType::Pointer colormap = ColormapType::New();
//      rgbfilter->SetColormap( colormap );
//      }
//    else if ( colormap == "jet"  )  
//      {
//      typedef itk::Functor::ScalarToRGBJetColormapFunctor<ColormapImageType::PixelType, 
//        RGBImageType::PixelType::ValueType> ColormapType;
//      ColormapType::Pointer colormap = ColormapType::New();
//      rgbfilter->SetColormap( colormap );
//      }
//    else if ( colormap == "hsv"  )  
//      {
//      typedef itk::Functor::ScalarToRGBHSVColormapFunctor<ColormapImageType::PixelType, 
//        RGBImageType::PixelType::ValueType> ColormapType;
//      ColormapType::Pointer colormap = ColormapType::New();
//      rgbfilter->SetColormap( colormap );
//      }
//    else if ( colormap == "tusti"  )  
//      {
//      typedef itk::Functor::ScalarToRGBTustiColormapFunctor<ColormapImageType::PixelType, 
//        RGBImageType::PixelType::ValueType> ColormapType;
//      ColormapType::Pointer colormap = ColormapType::New();
//      rgbfilter->SetColormap( colormap );
//      }
//    else if ( colormap == "custom"  )  
//      {
//      typedef itk::Functor::ScalarToRGBCustomColormapFunctor<ColormapImageType::PixelType, 
//        RGBImageType::PixelType::ValueType> ColormapType;
//      ColormapType::Pointer colormap = ColormapType::New();
//      rgbfilter->SetColormap( colormap );
//      }
//
//    rgbfilter->GetColormap()->SetMinimumRGBComponentValue( 0 );
//    rgbfilter->GetColormap()->SetMaximumRGBComponentValue( 255 );
//    rgbfilter->GetColormap()->SetMinimumInputValue( stats->GetMinimum() );
//    rgbfilter->GetColormap()->SetMaximumInputValue( stats->GetMaximum() );
//    rgbfilter->Update();
//  
//    typedef itk::ImageFileWriter<RGBImageType> WriterType;
//    WriterType::Pointer writer = WriterType::New();
//    writer->SetInput( rgbfilter->GetOutput() );
//    writer->SetFileName( argv[4] );
//    writer->Update();
//    }



  return 0;
}
