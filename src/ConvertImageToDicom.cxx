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
#include "itkShiftScaleImageFilter.h"

#include "gdcmUIDGenerator.h"

#include <vector>
#include <string>
#include <itksys/SystemTools.hxx>


int main( int argc, char* argv[] )
{

  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputDicomDirectory outputPrefix <uidPrefix>";
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

  typedef signed short    PixelType;
  const unsigned int      Dimension = 3;

  typedef itk::Image< PixelType, Dimension >      ImageType;
  typedef itk::ImageFileReader< ImageType >       ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject &excp)
    {
    std::cerr << "Exception thrown while writing the image" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::GDCMImageIO                        ImageIOType;
  typedef itk::NumericSeriesFileNames             NamesGeneratorType;

  ImageIOType::Pointer gdcmIO = ImageIOType::New();
  gdcmIO->SetKeepOriginalUID( false );

  const char * outputDirectory = argv[2];

  itksys::SystemTools::MakeDirectory( outputDirectory );


  typedef signed short    OutputPixelType;
  const unsigned int      OutputDimension = 2;

  typedef itk::Image< OutputPixelType, OutputDimension >    Image2DType;

  typedef itk::ImageSeriesWriter<
                         ImageType, Image2DType >  SeriesWriterType;
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
  if( argc > 4 && atoi( argv[4] ) != '0' )
    {
    suid.SetRoot( argv[4] );
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
    SeriesWriterType::DictionaryRawPointer dict = new SeriesWriterType::DictionaryType;

    std::string tagkey, value;

    for( int n = 5; n < argc; n++ )
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

////////////////////////////////////////////////
// 4) Shift data to undo the effect of a rescale intercept by the
//    DICOM reader
  std::string interceptTag("0028|1052");
  typedef itk::MetaDataObject< std::string > MetaDataStringType;
  itk::MetaDataObjectBase::Pointer entry = (reader->GetOutput()->GetMetaDataDictionary())[interceptTag];

  MetaDataStringType::ConstPointer interceptValue =
    dynamic_cast<const MetaDataStringType *>( entry.GetPointer() ) ;

  int interceptShift = 0;
  if( interceptValue )
    {
    std::string tagValue = interceptValue->GetMetaDataObjectValue();
    interceptShift = -atoi ( tagValue.c_str() );
    }

  typedef itk::ShiftScaleImageFilter<ImageType, ImageType> ShiftScaleType;
  ShiftScaleType::Pointer shiftScale = ShiftScaleType::New();
  shiftScale->SetInput( reader->GetOutput());
  shiftScale->SetShift( interceptShift );

  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();

  seriesWriter->SetInput( shiftScale->GetOutput() );
  seriesWriter->SetImageIO( gdcmIO );
  seriesWriter->SetMetaDataDictionaryArray( &dictionaryArray );




  std::string format = std::string( argv[2] )
    + std::string( "/" ) + std::string( argv[3] )
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



//#include "itkImage.h"
//#include "itkImageFileReader.h"
//#include "itkImageSeriesReader.h"
//#include "itkImageSeriesWriter.h"
//#include "itkGDCMImageIO.h"
//#include "itkGDCMSeriesFileNames.h"
//#include "itkNumericSeriesFileNames.h"
//#include "itkRGBPixel.h"
//
//#include <string>
//
//template<typename PixelType>
//int ConvertImageToDICOM( int argc, char *argv[] )
//{
//  typedef itk::Image<PixelType, 4> Image4DType;
//  typedef itk::Image<PixelType, 3> ImageType;
//  typedef itk::Image<PixelType, 2> Image2DType;
//
//  typedef itk::ImageFileReader<ImageType> ReaderType;
//  typename ReaderType::Pointer reader = ReaderType::New();
//  reader->SetFileName( argv[1] );
//  reader->Update();
//
//  // Read the input dicom series from which to copy the header information
//  // Cf example DicomSeriesReadImageWrite2.cxx in ITK/Examples/
////  typedef itk::ImageSeriesReader<Image4DType> SeriesReaderType;
////  typename SeriesReaderType::Pointer seriesReader = SeriesReaderType::New();
//  typedef itk::GDCMImageIO ImageIOType;
////  typename ImageIOType::Pointer dicomIO = ImageIOType::New();
////  seriesReader->SetImageIO( dicomIO );
////
////  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
////  typename NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
////  nameGenerator->SetUseSeriesDetails( true );
////  nameGenerator->AddSeriesRestriction("0008|0021" );
////  nameGenerator->SetDirectory( argv[2] );
////
////  typedef std::vector< std::string >    SeriesIdContainer;
////  const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
////
////  // Use the first series found
////  std::string seriesIdentifier = seriesUID.begin()->c_str();
////  seriesReader->SetFileNames( nameGenerator->GetFileNames( seriesIdentifier ) );
////  seriesReader->Update();
//
//  // Now get the meta data dictionary to attach to the new dicom image series.
//
////  itk::MetaDataDictionary &dict = dicomIO->GetMetaDataDictionary();
////  itk::MetaDataDictionary &dictNew = reader->GetOutput()->GetMetaDataDictionary();
////
////  std::string tagkey, value;
////
////  tagkey = "0008|0008";
////  value ="DERIVED\\SECONDARY";
////  itk::EncapsulateMetaData<std::string>( dictNew, tagkey, value );
////  tagkey = "0020|4000";
////  value = std::string( argv[5] );
////  itk::EncapsulateMetaData<std::string>( dictNew, tagkey, value );
//
//
//  // Now write out the resulting image
//  itksys::SystemTools::MakeDirectory( argv[3] );
//
//  std::string filePrefix = std::string( argv[3] )
//    + std::string( "/" ) + std::string( argv[4] )
//    + std::string( "%06d" );
//
//  typedef itk::NumericSeriesFileNames NamesGeneratorType2;
//  typename NamesGeneratorType2::Pointer namesGenerator2
//    = NamesGeneratorType2::New();
//  namesGenerator2->SetSeriesFormat( filePrefix.c_str() );
//  namesGenerator2->SetStartIndex( 0  );
//  namesGenerator2->SetEndIndex(
//    reader->GetOutput()->GetLargestPossibleRegion().GetSize()[2] - 1 );
//  namesGenerator2->SetIncrementIndex( 1 );
//
//
//  typename ImageIOType::Pointer dicomIO2 = ImageIOType::New();
//
//  typedef itk::ImageSeriesWriter<ImageType, Image2DType> SeriesWriterType;
//  typename SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
//  seriesWriter->SetInput( reader->GetOutput() );
//  seriesWriter->SetImageIO( dicomIO2 );
//  seriesWriter->SetFileNames( namesGenerator2->GetFileNames() );
////  seriesWriter->SetMetaDataDictionary( reader->GetOutput()->GetMetaDataDictionary() );
//
//  try
//    {
//    seriesWriter->Update();
//    }
//  catch( itk::ExceptionObject & excp )
//    {
//    std::cerr << "Exception thrown while writing the series " << std::endl;
//    std::cerr << excp << std::endl;
//    return EXIT_FAILURE;
//    }
//
//  return 0;
//}
//
//int main( int argc, char *argv[] )
//{
//  if ( argc < 4  )
//    {
//    std::cout << "Usage: " << argv[0]
//      << " inputImage referenceDicomSeriesDirectory outputDirectory prefix [Comments] [isRGB] "
//      << std::endl;
//
//    exit( 1 );
//    }
//
//  typedef short PixelType;
//  typedef itk::RGBPixel<char> RGBPixelType;
//
//  bool isRGB = ( argc > 4 ) ? static_cast<bool>( atoi( argv[4] ) ) : false;
//
//  if( isRGB )
//    {
//    ConvertImageToDICOM<RGBPixelType>( argc, argv );
//    }
//  else
//    {
//    ConvertImageToDICOM<PixelType>( argc, argv );
//    }
//}
//
