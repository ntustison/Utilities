#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkGDCMImageIO.h"

#include <string>

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0]
      << " inputImage dicomImage outputImage" << std::endl;
    exit( EXIT_FAILURE );
    }

  typedef signed short InputPixelType;
  const unsigned int Dimension = 2;

  typedef itk::Image<InputPixelType, Dimension> InputImageType;

  typedef itk::ImageFileReader<InputImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  ReaderType::Pointer dicomReader = ReaderType::New();
  typedef itk::GDCMImageIO ImageIOType;
  ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
  dicomReader->SetImageIO( gdcmImageIO );
  dicomReader->SetFileName( argv[2] );
  dicomReader->Update();

  typedef itk::MetaDataDictionary DictionaryType;
  DictionaryType dictionary;
  DictionaryType srcDictionary = dicomReader->GetMetaDataDictionary();

  DictionaryType::ConstIterator its = srcDictionary.Begin();

  while( its != srcDictionary.End() )
    {
    itk::MetaDataObjectBase::Pointer entry = its->second;

    itk::MetaDataObject<std::string>::Pointer entryvalue =
      dynamic_cast<itk::MetaDataObject<std::string> *>( entry.GetPointer() );

    std::string tagkey = its->first;

    if( std::strcmp( tagkey.c_str(), "0008|0060" ) == 0 )
      {
      tagkey = "0008|0060";  // Modality
      std::string value = "CT";
      itk::EncapsulateMetaData<std::string>( dictionary, tagkey, value );
      }
    else
      {
      itk::EncapsulateMetaData( dictionary, tagkey,
        entryvalue->GetMetaDataObjectValue() );
      }

    ++its;
    }


//  std::string tagkey, value;
//  tagkey = "0020|000e";  // Series Instance UID
//  value = "1.2.826.0.1.3680043.2.1125.1.18358515732600720380145022665620682";
//  itk::EncapsulateMetaData<std::string>( dictionary, tagkey, value );

//  tagkey = "0008|0060";  // Modality
//  value = "OT";
//  itk::EncapsulateMetaData<std::string>( dictionary, tagkey, value );
//
  InputImageType::Pointer inputImage = reader->GetOutput();
  inputImage->SetMetaDataDictionary( dictionary );

  typedef itk::ImageFileWriter<InputImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( inputImage );
  writer->SetFileName( argv[3] );
  writer->SetImageIO( gdcmImageIO );
  writer->Update();

  return 0;
}
