/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ImageReadDicomSeriesWrite.cxx,v $
  Language:  C++
  Date:      $Date: 2009-03-17 20:36:50 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkGDCMImageIO.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageFileReader.h"
#include "itkImageSeriesWriter.h"
#include "itkMetaDataObject.h"
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



#include "gdcmUIDGenerator.h"

#include <vector>
#include <string>
#include <itksys/SystemTools.hxx>


int main( int argc, char* argv[] )
{

  if( argc < 11 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage ";
    std::cerr << "mask ";
    std::cerr << "colormap customColormapFile ";
    std::cerr << "minimumInput maximumInput ";
    std::cerr << "minimumRGBOutput maximumRGBOutput ";
    std::cerr << "outputDicomDirectory outputPrefix <uidPrefix>";
    std::cerr << " <TagID_1,Value_1> ... <TagID_1,Value_1> " << std::endl;
    std::cerr << "  Note:  instead of tags being specified 0008|0020, ";
    std::cerr << "  one needs to put 0008\\|00020." << std::endl;
    std::cerr << "  Set uid prefix to  0 if one is to be generated in its entirety. "<< std::endl;
    std::cerr << "  ALso some common tagID,value pairs: "<< std::endl;
    std::cerr << "    0008|0008 -> ImageType (DERIVED\\SECONDARY)" << std::endl;
    std::cerr << "    0008|0020 -> StudyDate " << std::endl;
    std::cerr << "    0008|0021 -> SeriesDate " << std::endl;
    std::cerr << "    0008|0022 -> AcquisitionDate " << std::endl;
    std::cerr << "    0008|0030 -> StudyTime " << std::endl;
    std::cerr << "    0008|0031 -> SeriesTime " << std::endl;
    std::cerr << "    0008|0032 -> AcquisitionTime " << std::endl;
    std::cerr << "    0008|0050 -> AccessionNumber " << std::endl;
    std::cerr << "    0008|1030 -> StudyDescription " << std::endl;
    std::cerr << "    0010|0010 -> PatientName " << std::endl;
    std::cerr << "    0010|0020 -> PatientID " << std::endl;
    std::cerr << "    0020|000d -> StudyID " << std::endl;


    return EXIT_FAILURE;
    }

  const unsigned int      ImageDimension = 3;

  typedef signed short    PixelType;
  typedef itk::RGBPixel<PixelType> RGBPixelType;

  typedef float RealType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<float, ImageDimension> RealImageType;
  typedef itk::Image<RGBPixelType, ImageDimension> RGBImageType;

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::Image<unsigned char, ImageDimension> MaskImageType;
  MaskImageType::Pointer maskImage = MaskImageType::New();
  typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
  MaskReaderType::Pointer maskreader = MaskReaderType::New();
  maskreader->SetFileName( argv[2] );
  try
    {
    maskreader->Update();
    maskImage = maskreader->GetOutput();
    }
  catch(...)
    {
    maskImage = NULL;
    };


  std::string colormapString( argv[3] );

  typedef itk::ScalarToRGBColormapImageFilter<RealImageType,
    RGBImageType> RGBFilterType;
  RGBFilterType::Pointer rgbfilter = RGBFilterType::New();
  rgbfilter->SetInput( reader->GetOutput() );

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
    rgbfilter->SetColormap( RGBFilterType::Jet );
    }
  else if ( colormapString == "hsv"  )
    {
    rgbfilter->SetColormap( RGBFilterType::HSV );
    }
  else if ( colormapString == "overunder"  )
    {
    rgbfilter->SetColormap( RGBFilterType::OverUnder );
    }
  else if ( colormapString == "custom"  )
    {
    typedef itk::Function::CustomColormapFunction<RealImageType::PixelType,
      RGBImageType::PixelType> ColormapType;
    ColormapType::Pointer colormap = ColormapType::New();

    std::ifstream str( argv[4] );
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

  if( maskImage )
    {
    RealType maskMinimumValue = itk::NumericTraits<RealType>::max();
    RealType maskMaximumValue = itk::NumericTraits<RealType>::NonpositiveMin();

    itk::ImageRegionIterator<MaskImageType> ItM( maskImage,
      maskImage->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<RealImageType> ItS( reader->GetOutput(),
      reader->GetOutput()->GetLargestPossibleRegion() );
    for( ItM.GoToBegin(), ItS.GoToBegin(); !ItM.IsAtEnd(); ++ItM, ++ItS )
      {
      if( ItM.Get() != 0 )
        {
        if( maskMinimumValue > ItS.Get() )
          {
          maskMinimumValue = ItS.Get();
          }
        if( maskMaximumValue < ItS.Get() )
          {
          maskMaximumValue = ItS.Get();
          }
        }
      }

    rgbfilter->SetUseInputImageExtremaForScaling( false );
    rgbfilter->GetColormap()->SetMinimumInputValue( maskMinimumValue );
    rgbfilter->GetColormap()->SetMaximumInputValue( maskMaximumValue );
    }

  rgbfilter->SetUseInputImageExtremaForScaling( false );
  rgbfilter->GetColormap()->SetMinimumInputValue(
    static_cast<RealType>( atof( argv[5] ) ) );
  rgbfilter->GetColormap()->SetMaximumInputValue(
    static_cast<RealType>( atof( argv[6] ) ) );
  rgbfilter->GetColormap()->SetMinimumRGBComponentValue(
    static_cast<RGBPixelType::ComponentType>( atof( argv[7] ) ) );
  rgbfilter->GetColormap()->SetMaximumRGBComponentValue(
    static_cast<RGBPixelType::ComponentType>( atof( argv[8] ) ) );

  try
    {
    rgbfilter->Update();
    }
  catch (...)
    {
    return EXIT_FAILURE;
    }

  if( maskImage )
    {
    itk::ImageRegionIterator<MaskImageType> ItM( maskImage,
      maskImage->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<RGBImageType> ItC( rgbfilter->GetOutput(),
      rgbfilter->GetOutput()->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<RealImageType> ItS( reader->GetOutput(),
      reader->GetOutput()->GetLargestPossibleRegion() );

    ItM.GoToBegin();
    ItC.GoToBegin();
    ItS.GoToBegin();

    while( !ItM.IsAtEnd() )
      {
      if( ItM.Get() == 0 )
        {
        RGBPixelType rgbpixel;

        rgbpixel.Fill( itk::NumericTraits<RGBPixelType::ComponentType>::Zero );

        ItC.Set( rgbpixel );
        }
      ++ItM;
      ++ItC;
      ++ItS;
      }
    }

  RGBImageType::Pointer rgbImage = rgbfilter->GetOutput();



  typedef itk::GDCMImageIO                        ImageIOType;
  typedef itk::NumericSeriesFileNames             NamesGeneratorType;

  ImageIOType::Pointer gdcmIO = ImageIOType::New();
  gdcmIO->SetKeepOriginalUID( false );

  const char * outputDirectory = argv[9];

  itksys::SystemTools::MakeDirectory( outputDirectory );


  const unsigned int      OutputDimension = 2;

  typedef itk::Image< RGBPixelType, OutputDimension >    Image2DType;

  typedef itk::ImageSeriesWriter<
                         RGBImageType, Image2DType >  SeriesWriterType;
  SeriesWriterType::DictionaryArrayType dictionaryArray;

  NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();

//  tagkey = "0008|0060"; // Modality
//  value = "MR";
//  itk::EncapsulateMetaData<std::string>(dict, tagkey, value );
//  tagkey = "0008|0008"; // Image Type
//  value = "DERIVED\\SECONDARY";
//  itk::EncapsulateMetaData<std::string>(dict, tagkey, value);
//  tagkey = "0008|0064"; // Conversion Type
//  value = "DV";
//  itk::EncapsulateMetaData<std::string>(dict, tagkey, value);

  gdcm::UIDGenerator suid;
  if( argc > 11 && atoi( argv[11] ) != 0 )
    {
    suid.SetRoot( argv[11] );
    }
  std::string seriesUID = suid.Generate();

  gdcm::UIDGenerator fuid;
  std::string frameOfReferenceUID = fuid.Generate();

  ImageType::RegionType region =
     reader->GetOutput()->GetLargestPossibleRegion();

  ImageType::IndexType start = region.GetIndex();
  ImageType::SizeType  size  = region.GetSize();

  for( unsigned int s = 0; s < size[2]; s++ )
    {
    SeriesWriterType::DictionaryRawPointer dict =
      new SeriesWriterType::DictionaryType;

    std::string tagkey, value;

    for( int n = 12; n < argc; n++ )
      {
      std::cout << argv[n] << std::endl;

      std::string pair = std::string( argv[n] );
      std::string::size_type crosspos = pair.find( ',', 0 );

      tagkey = pair.substr( 0, crosspos );
      value = pair.substr( crosspos + 1, pair.length() );

      itk::EncapsulateMetaData<std::string>( *dict, tagkey, value );
      }

    itk::EncapsulateMetaData<std::string>( *dict,"0020|000e", seriesUID );
    itk::EncapsulateMetaData<std::string>( *dict,"0020|0052", frameOfReferenceUID );

    gdcm::UIDGenerator sopuid;
    std::string sopInstanceUID = sopuid.Generate();
    itk::EncapsulateMetaData<std::string>( *dict, "0008|0018", sopInstanceUID );
    itk::EncapsulateMetaData<std::string>( *dict, "0002|0003", sopInstanceUID );

    // Slice number
    itksys_ios::ostringstream value2;
    value2.str( "" );
    value2 << s + 1;
    itk::EncapsulateMetaData<std::string>( *dict, "0020|0013", value2.str() );

   // Image Position Patient: This is calculated by computing the
    // physical coordinate of the first pixel in each slice.
    ImageType::PointType position;
    ImageType::IndexType index;
    index[0] = 0;
    index[1] = 0;
    index[2] = s;
    reader->GetOutput()->TransformIndexToPhysicalPoint( index, position );

    value2.str("");
    value2 << position[0] << "\\" << position[1] << "\\" << position[2];
    itk::EncapsulateMetaData<std::string>( *dict,"0020|0032", value2.str() );
    // Slice Location: For now, we store the z component of the Image
    // Position Patient.
    value2.str( "" );
    value2 << position[2];

    itk::EncapsulateMetaData<std::string>( *dict,"0020|1041", value2.str() );

    dictionaryArray.push_back( dict );
    }

  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
  seriesWriter->SetInput( rgbImage );
  seriesWriter->SetImageIO( gdcmIO );
  seriesWriter->SetMetaDataDictionaryArray( &dictionaryArray );

  std::string format = std::string( argv[9] )
    + std::string( "/" ) + std::string( argv[10] )
    + std::string( "%04d.dcm" );

  namesGenerator->SetSeriesFormat( format.c_str() );

  namesGenerator->SetStartIndex( start[2] );
  namesGenerator->SetEndIndex( start[2] + size[2] - 1 );
  namesGenerator->SetIncrementIndex( 1 );


  seriesWriter->SetFileNames( namesGenerator->GetFileNames() );


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


  return EXIT_SUCCESS;
}


