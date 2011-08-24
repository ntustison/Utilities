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
#include <fstream>

template <unsigned int Dimension>
int convert( int argc, char* argv[] )
{
  typedef signed short    PixelType;

  typedef itk::Image< PixelType, Dimension >      ImageType;
  typedef itk::ImageFileReader< ImageType >       ReaderType;

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );

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

  typename ImageIOType::Pointer gdcmIO = ImageIOType::New();
  gdcmIO->SetKeepOriginalUID( false );

  const char * outputDirectory = argv[3];

  itksys::SystemTools::MakeDirectory( outputDirectory );


  typedef signed short    OutputPixelType;
  const unsigned int      OutputDimension = 2;

  typedef itk::Image< OutputPixelType, OutputDimension >    Image2DType;

  typedef itk::ImageSeriesWriter<
                         ImageType, Image2DType >  SeriesWriterType;
  typename SeriesWriterType::DictionaryArrayType dictionaryArray;

  typename NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();

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
  if( argc > 5 && atoi( argv[5] ) != 0 )
    {
    suid.SetRoot( argv[5] );
    }
  std::string seriesUID = suid.Generate();

  gdcm::UIDGenerator fuid;
  std::string frameOfReferenceUID = fuid.Generate();

  typename ImageType::RegionType region =
     reader->GetOutput()->GetLargestPossibleRegion();

  typename ImageType::IndexType start = region.GetIndex();
  typename ImageType::SizeType  size  = region.GetSize();

  unsigned int numberOfTimePoints = 1;
  unsigned int numberOfSlices = 1;

  if( Dimension == 4 )
    {
    numberOfTimePoints = size[3];
    }
  if( Dimension >= 3 )
    {
    numberOfSlices = size[2];
    }

  std::vector<std::vector<std::string> > values;
  std::vector<std::string> tagkeys;

  for( int n = 6; n < argc; n++ )
    {
    std::cout << argv[n] << std::endl;

    std::string pair = std::string( argv[n] );
    std::string::size_type crosspos = pair.find( ',', 0 );

    std::string tagkey = pair.substr( 0, crosspos );
//    value = pair.substr( crosspos + 1, pair.length() );

    tagkeys.push_back( tagkey );
    }

  for( int n = 6; n < argc; n++ )
    {
    std::string pair = std::string( argv[n] );
    std::string::size_type crosspos = pair.find( ',', 0 );

    std::string tagkey = pair.substr( 0, crosspos );
    std::string value = pair.substr( crosspos + 1, pair.length() );

    std::vector<std::string> tagkeyvalues;
    if( value.find( ".txt" ) == std::string::npos )
      {
      for( unsigned int t = 0; t < numberOfTimePoints; t++ )
        {
        for( unsigned int s = 0; s < numberOfSlices; s++ )
          {
          tagkeyvalues.push_back( value );
          }
        }
      values.push_back( tagkeyvalues );
      }
    else
      {
      std::string file = value;

      std::fstream str( file.c_str() );

      while( str >> value )
        {
        tagkeyvalues.push_back( value );
        }
      if( tagkeyvalues.size() != numberOfTimePoints * numberOfSlices )
        {
        std::cerr << "Number of values in file " << file << " do no match "
          << "the number of time points x the number of slices ("
          << tagkeyvalues.size() << " != " << numberOfTimePoints * numberOfSlices
          << ")" << std::endl;
        return EXIT_FAILURE;
        }
      else
        {
        values.push_back( tagkeyvalues );
        }
      }
    }

  for( unsigned int t = 0; t < numberOfTimePoints; t++ )
    {
    for( unsigned int s = 0; s < numberOfSlices; s++ )
      {
      typename SeriesWriterType::DictionaryRawPointer dict =
        new typename SeriesWriterType::DictionaryType;

      std::string tagkey, value;

      std::cout << "t = " << t << ", s = " << s << std::endl;
      for( unsigned int n = 0; n < tagkeys.size(); n++ )
        {
        std::cout << "  " << tagkeys[n] << " -> " << values[n][t*numberOfSlices + s]
          << std::endl;
        itk::EncapsulateMetaData<std::string>( *dict, tagkeys[n],
          values[n][t*numberOfSlices + s] );
        }

      itk::EncapsulateMetaData<std::string>( *dict,"0020|000e", seriesUID );
      itk::EncapsulateMetaData<std::string>( *dict,"0020|0052", frameOfReferenceUID );

      gdcm::UIDGenerator sopuid;
      std::string sopInstanceUID = sopuid.Generate();
      itk::EncapsulateMetaData<std::string>( *dict, "0008|0018", sopInstanceUID );
      itk::EncapsulateMetaData<std::string>( *dict, "0002|0003", sopInstanceUID );

      // instance number
      typename itksys_ios::ostringstream value2;
      value2.str( "" );
      value2 << ( t * numberOfSlices + s + 1 );
      itk::EncapsulateMetaData<std::string>( *dict, "0020|0013", value2.str() );

     // Image Position Patient: This is calculated by computing the
      // physical coordinate of the first pixel in each slice.
      typename ImageType::PointType position;
      typename ImageType::IndexType index;
      index[0] = start[0];
      index[1] = start[1];
      if( Dimension >= 3)
        {
        index[2] = start[2] + s;
        }
      if( Dimension == 4 )
        {
        index[3] = start[3] + t;
        }
      reader->GetOutput()->TransformIndexToPhysicalPoint( index, position );

      value2.str( "" );
      if( Dimension >= 3 )
        {
        value2 << position[0] << "\\" << position[1] << "\\" << position[2];
        }
      else
        {
        value2 << position[0] << "\\" << position[1] << "\\" << 0.0;
        }

      itk::EncapsulateMetaData<std::string>( *dict,"0020|0032", value2.str() );
      // Slice Location: For now, we store the z component of the Image
      // Position Patient.
      value2.str( "" );
      if( Dimension >= 3 )
        {
        value2 << position[2];
        }
      else
        {
        value2 << 0.0;
        }

      itk::EncapsulateMetaData<std::string>( *dict,"0020|1041", value2.str() );

      dictionaryArray.push_back( dict );
      }
    }

////////////////////////////////////////////////
// 4) Shift data to undo the effect of a rescale intercept by the
//    DICOM reader
  std::string interceptTag("0028|1052");
  typedef itk::MetaDataObject< std::string > MetaDataStringType;
  itk::MetaDataObjectBase::Pointer entry = (
    reader->GetOutput()->GetMetaDataDictionary() )[interceptTag];

  typename MetaDataStringType::ConstPointer interceptValue =
    dynamic_cast<const MetaDataStringType *>( entry.GetPointer() ) ;

  int interceptShift = 0;
  if( interceptValue )
    {
    std::string tagValue = interceptValue->GetMetaDataObjectValue();
    interceptShift = -atoi ( tagValue.c_str() );
    }

  typedef itk::ShiftScaleImageFilter<ImageType, ImageType> ShiftScaleType;
  typename ShiftScaleType::Pointer shiftScale = ShiftScaleType::New();
  shiftScale->SetInput( reader->GetOutput());
  shiftScale->SetShift( interceptShift );

  typename SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();

  seriesWriter->SetInput( shiftScale->GetOutput() );
  seriesWriter->SetImageIO( gdcmIO );
  seriesWriter->SetMetaDataDictionaryArray( &dictionaryArray );

  std::string format = std::string( argv[3] )
    + std::string( "/" ) + std::string( argv[4] )
    + std::string( "%04d.dcm" );

  namesGenerator->SetSeriesFormat( format.c_str() );

  namesGenerator->SetStartIndex( 0 );
  namesGenerator->SetEndIndex( numberOfTimePoints * numberOfSlices - 1 );
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


int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension";
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

  switch( atoi( argv[1] ) )
   {
   case 2:
     convert<2>( argc, argv );
     break;
   case 3:
     convert<3>( argc, argv );
     break;
   case 4:
     convert<4>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
