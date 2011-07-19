/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDeformationFieldReader.hxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:21:04 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkDeformationFieldReader_hxx
#define _itkDeformationFieldReader_hxx
#include "itkDeformationFieldReader.h"

#include "itkObjectFactory.h"
#include "itkImageIOFactory.h"
#include "itkConvertPixelBuffer.h"
#include "itkImageRegion.h"
#include "itkPixelTraits.h"
#include "itkVectorImage.h"
#include "itkImageRegionIterator.h"

#include <itksys/SystemTools.hxx>
#include <fstream>

namespace itk
{

template <class TImage, class TDeformationField, class ConvertPixelTraits>
DeformationFieldReader<TImage, TDeformationField, ConvertPixelTraits>
::DeformationFieldReader()
{
  m_ImageIO = 0;
  m_FileName = "";
  m_UserSpecifiedImageIO = false;
  m_UseAvantsNamingConvention = false;
  
  this->m_Image = TImage::New();
}

template <class TImage, class TDeformationField, class ConvertPixelTraits>
DeformationFieldReader<TImage, TDeformationField, ConvertPixelTraits>
::~DeformationFieldReader()
{
}

template <class TImage, class TDeformationField, class ConvertPixelTraits>
void DeformationFieldReader<TImage, TDeformationField, ConvertPixelTraits>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  if (m_ImageIO)
    {
    os << indent << "ImageIO: \n";
    m_ImageIO->Print(os, indent.GetNextIndent());
    }
  else
    {
    os << indent << "ImageIO: (null)" << "\n";
    }

  os << indent << "UserSpecifiedImageIO flag: " << m_UserSpecifiedImageIO << "\n";
  os << indent << "m_FileName: " << m_FileName << "\n";
}


template <class TImage, class TDeformationField, class ConvertPixelTraits>
void 
DeformationFieldReader<TImage, TDeformationField, ConvertPixelTraits>
::SetImageIO( ImageIOBase * imageIO)
{
  itkDebugMacro("setting ImageIO to " << imageIO ); 
  if (this->m_ImageIO != imageIO ) 
    {
    this->m_ImageIO = imageIO;
    this->Modified(); 
    } 
  m_UserSpecifiedImageIO = true;
}


template <class TImage, class TDeformationField, class ConvertPixelTraits>
void 
DeformationFieldReader<TImage, TDeformationField, ConvertPixelTraits>
::GenerateOutputInformation(void)
{

  typename TDeformationField::Pointer output = this->GetOutput();

  itkDebugMacro(<<"Reading file for GenerateOutputInformation()" << m_FileName);
  
  // Check to see if we can read the file given the name or prefix
  //
  if ( m_FileName == "" )
    {
    throw DeformationFieldReaderException(__FILE__, __LINE__, "FileName must be specified", ITK_LOCATION);
    }

  // Test if the files exist and if it can be open.
  // and exception will be thrown otherwise.
  //

  this->TestFileExistanceAndReadability();
  
  unsigned int dimension = itk::GetVectorDimension
     <DeformationFieldPixelType>::VectorDimension;

  std::string filename = this->m_FileName; 
  std::string::size_type Pos = this->m_FileName.rfind( "." );
  std::string extension( this->m_FileName, Pos, this->m_FileName.length()-1 );

  for ( unsigned int i = 0; i < dimension; i++ )
    {
    this->m_FileName = std::string( filename, 0, Pos );

    if ( this->m_UseAvantsNamingConvention )
      {
      switch ( i )
        {
        case 0:
          this->m_FileName += std::string( "xvec" );
          break;
        case 1:
          this->m_FileName += std::string( "yvec" );
          break;
        case 2:
          this->m_FileName += std::string( "zvec" );
          break;
        default:
          this->m_FileName += std::string( "you_are_screwed_vec" );
          break;
        }  
      }
    else
      {
      itk::OStringStream buf;
      buf << i;
      this->m_FileName += ( std::string( "." )  + std::string( buf.str().c_str() ) );
      }
    this->m_FileName += extension;

    if ( i == 0 )
      {

      itkDebugMacro( << "Generating output information from the file " << this->m_FileName );
    
      if ( m_UserSpecifiedImageIO == false ) //try creating via factory
        {
        m_ImageIO = ImageIOFactory::CreateImageIO( m_FileName.c_str(), ImageIOFactory::ReadMode );
        }
      
      if ( m_ImageIO.IsNull() )
        {
        OStringStream msg;
        msg << " Could not create IO object for file "
            << m_FileName.c_str() << std::endl;
        msg << "  Tried to create one of the following:" << std::endl;
        std::list<LightObject::Pointer> allobjects = 
          ObjectFactoryBase::CreateAllInstance("itkImageIOBase");
        for(std::list<LightObject::Pointer>::iterator i = allobjects.begin();
            i != allobjects.end(); ++i)
          {
          ImageIOBase* io = dynamic_cast<ImageIOBase*>(i->GetPointer());
          msg << "    " << io->GetNameOfClass() << std::endl; 
          }
        msg << "  You probably failed to set a file suffix, or" << std::endl;
        msg << "    set the suffix to an unsupported type." << std::endl;
        DeformationFieldReaderException e(__FILE__, __LINE__, msg.str().c_str(), ITK_LOCATION);
        throw e;
        return;
        }    

      // Got to allocate space for the image. Determine the characteristics of
      // the image.
      //
      m_ImageIO->SetFileName(m_FileName.c_str());
      m_ImageIO->ReadImageInformation();

      typename TDeformationField::SizeType dimSize;
      double spacing[ TDeformationField::ImageDimension ];
      double origin[ TDeformationField::ImageDimension ];
      typename TDeformationField::DirectionType direction;
      std::vector<double> axis;

      for(unsigned int i=0; i<TImage::ImageDimension; i++)
        {
        if ( i < m_ImageIO->GetNumberOfDimensions() )
          {
          dimSize[i] = m_ImageIO->GetDimensions(i);
          spacing[i] = m_ImageIO->GetSpacing(i);
          origin[i]  = m_ImageIO->GetOrigin(i);
    // Please note: direction cosines are stored as columns of the
    // direction matrix
          axis = m_ImageIO->GetDirection(i);
          for (unsigned j=0; j<TImage::ImageDimension; j++)
            {
            if (j < m_ImageIO->GetNumberOfDimensions())
              {
              direction[j][i] = axis[j];
              }
            else
              {
              direction[j][i] = 0.0;
              }
            }
          }
        else
          {
          // Number of dimensions in the output is more than number of dimensions
          // in the ImageIO object (the file).  Use default values for the size,
          // spacing, origin and direction for the final (degenerate) dimensions.
          dimSize[i] = 1;  
          spacing[i] = 1.0;
          origin[i] = 0.0;
          for (unsigned j = 0; j < TImage::ImageDimension; j++)
            {
            if (i == j)
              {
              direction[j][i] = 1.0;
              }
            else
              {
              direction[j][i] = 0.0;
              }
            }
          }
        }

      output->SetSpacing( spacing );     // Set the image spacing
      output->SetOrigin( origin );       // Set the image origin
      output->SetDirection( direction ); // Set the image direction cosines

      //Copy MetaDataDictionary from instantiated reader to output image.
      output->SetMetaDataDictionary(m_ImageIO->GetMetaDataDictionary());
      this->SetMetaDataDictionary(m_ImageIO->GetMetaDataDictionary());
      
      this->m_Image->SetSpacing( spacing );
      this->m_Image->SetOrigin( origin );
      this->m_Image->SetDirection( direction );
      this->m_Image->SetMetaDataDictionary(m_ImageIO->GetMetaDataDictionary());

      typedef typename TDeformationField::IndexType   IndexType;

      IndexType start;
      start.Fill(0);

      DeformationFieldRegionType region;
      region.SetSize(dimSize);
      region.SetIndex(start);
      
      ImageRegionType imageregion;
      imageregion.SetSize(dimSize);
      imageregion.SetIndex(start);

     
      // If a VectorImage, this requires us to set the 
      // VectorLength before allocate
      //if( strcmp( output->GetNameOfClass(), "VectorImage" ) == 0 ) 
      //  {
      //  typedef typename TImage::AccessorFunctorType AccessorFunctorType;
      //  AccessorFunctorType::SetVectorLength( output, m_ImageIO->GetNumberOfComponents() );
      //  }
      
      output->SetLargestPossibleRegion( region );
      this->m_Image->SetLargestPossibleRegion( imageregion );
      }
    }
  this->m_FileName = filename;  
}



template <class TImage, class TDeformationField, class ConvertPixelTraits>
void
DeformationFieldReader<TImage, TDeformationField, ConvertPixelTraits>
::TestFileExistanceAndReadability()
{
  unsigned int dimension = itk::GetVectorDimension
     <DeformationFieldPixelType>::VectorDimension;

  std::string filename = this->m_FileName; 
  std::string::size_type Pos = this->m_FileName.rfind( "." );
  std::string extension( this->m_FileName, Pos, this->m_FileName.length()-1 );

  for ( unsigned int i = 0; i < dimension; i++ )
    {
    this->m_FileName = std::string( filename, 0, Pos );

    if ( this->m_UseAvantsNamingConvention )
      {
      switch ( i )
        {
        case 0:
          this->m_FileName += std::string( "xvec" );
          break;
        case 1:
          this->m_FileName += std::string( "yvec" );
          break;
        case 2:
          this->m_FileName += std::string( "zvec" );
          break;
        default:
          this->m_FileName += std::string( "you_are_screwed_vec" );
          break;
        }  
      }
    else
      {
      itk::OStringStream buf;
      buf << i;
      this->m_FileName += ( std::string( "." )  + std::string( buf.str().c_str() ) );
      }
    this->m_FileName += extension;
    
    itkDebugMacro( << "Checking for the file " << this->m_FileName );
    
    // Test if the file exists.
    if( ! itksys::SystemTools::FileExists( m_FileName.c_str() ) )
      {
      DeformationFieldReaderException e(__FILE__, __LINE__);
      OStringStream msg;
      msg <<"The file doesn't exists. "
          << std::endl << "Filename = " << m_FileName
          << std::endl;
      e.SetDescription(msg.str().c_str());
      throw e;
      return;
      }

    // Test if the file can be open for reading access.
    std::ifstream readTester;
    readTester.open( m_FileName.c_str() );
    if( readTester.fail() )
      {
      readTester.close();
      OStringStream msg;
      msg <<"The file couldn't be opened for reading. "
          << std::endl << "Filename: " << m_FileName
          << std::endl;
      DeformationFieldReaderException e(__FILE__, __LINE__,msg.str().c_str(),ITK_LOCATION);
      throw e;
      return;

      }
    readTester.close();
    }
  this->m_FileName = filename;  
}




template <class TImage, class TDeformationField, class ConvertPixelTraits>
void
DeformationFieldReader<TImage, TDeformationField, ConvertPixelTraits>
::EnlargeOutputRequestedRegion(DataObject *output)
{
  typename TDeformationField::Pointer out = dynamic_cast<TDeformationField*>(output);

  // the ImageIO object cannot stream, then set the RequestedRegion to the
  // LargestPossibleRegion
  if (!m_ImageIO->CanStreamRead())
    {
    if (out)
      {
      out->SetRequestedRegion( out->GetLargestPossibleRegion() );
      }
    else
      {
      throw DeformationFieldReaderException(__FILE__, __LINE__,
                                     "Invalid output object type");
      }
    }
}


template <class TImage, class TDeformationField, class ConvertPixelTraits>
void DeformationFieldReader<TImage, TDeformationField, ConvertPixelTraits>
::GenerateData()
{
  typename TDeformationField::Pointer output = this->GetOutput();

  // allocate the output buffer
  output->SetBufferedRegion( output->GetRequestedRegion() );
  output->Allocate();

  this->m_Image->SetBufferedRegion( output->GetRequestedRegion() );
  this->m_Image->Allocate();

  // Test if the file exist and if it can be open.
  // and exception will be thrown otherwise.
  this->TestFileExistanceAndReadability();

  unsigned int dimension = itk::GetVectorDimension
     <DeformationFieldPixelType>::VectorDimension;

  std::string filename = this->m_FileName; 
  std::string::size_type Pos = this->m_FileName.rfind( "." );
  std::string extension( this->m_FileName, Pos, this->m_FileName.length()-1 );

  for ( unsigned int i = 0; i < dimension; i++ )
    {
    this->m_FileName = std::string( filename, 0, Pos );

    if ( this->m_UseAvantsNamingConvention )
      {
      switch ( i )
        {
        case 0:
          this->m_FileName += std::string( "xvec" );
          break;
        case 1:
          this->m_FileName += std::string( "yvec" );
          break;
        case 2:
          this->m_FileName += std::string( "zvec" );
          break;
        default:
          this->m_FileName += std::string( "you_are_screwed_vec" );
          break;
        }  
      }
    else
      {
      itk::OStringStream buf;
      buf << i;
      this->m_FileName += ( std::string( "." )  + std::string( buf.str().c_str() ) );
      }
    this->m_FileName += extension;

    itkDebugMacro( << "Reading image buffer from the file " << this->m_FileName );
    
    // Tell the ImageIO to read the file
    //
    ImagePixelType *buffer = 
      this->m_Image->GetPixelContainer()->GetBufferPointer();
    m_ImageIO->SetFileName(m_FileName.c_str());

    ImageIORegion ioRegion(TImage::ImageDimension);
    
    ImageIORegion::SizeType ioSize = ioRegion.GetSize();
    ImageIORegion::IndexType ioStart = ioRegion.GetIndex();

    typename TImage::SizeType dimSize;
    for(unsigned int j=0; j<TImage::ImageDimension; j++)
      {
      if (j < m_ImageIO->GetNumberOfDimensions())
        {
        dimSize[j] = m_ImageIO->GetDimensions(j);
        }
      else
        {
        // Number of dimensions in the output is more than number of dimensions
        // in the ImageIO object (the file).  Use default values for the size,
        // spacing, and origin for the final (degenerate) dimensions.
        dimSize[j] = 1;
        }
      }

    for(unsigned int j = 0; j < dimSize.GetSizeDimension(); ++j)
      {
      ioSize[j] = dimSize[j];
      }

    typedef typename TImage::IndexType   IndexType;
    IndexType start;
    start.Fill(0);
    for(unsigned int j = 0; j < start.GetIndexDimension(); ++j)
      {
      ioStart[j] = start[j];
      }

    ioRegion.SetSize(ioSize);
    ioRegion.SetIndex(ioStart);

    itkDebugMacro (<< "ioRegion: " << ioRegion);
 
    m_ImageIO->SetIORegion(ioRegion);

    if ( m_ImageIO->GetComponentTypeInfo()
         == typeid(ITK_TYPENAME ConvertPixelTraits::ComponentType)
         && (m_ImageIO->GetNumberOfComponents()
             == ConvertPixelTraits::GetNumberOfComponents()))
      {
      itkDebugMacro(<< "No buffer conversion required.");
      // allocate a buffer and have the ImageIO read directly into it
      m_ImageIO->Read(buffer);
//      return;
      }
    else // a type conversion is necessary
      {
      itkDebugMacro(<< "Buffer conversion required.");
      // note: char is used here because the buffer is read in bytes
      // regardles of the actual type of the pixels.
      ImageRegionType region = this->m_Image->GetBufferedRegion();
      char * loadBuffer = 
        new char[m_ImageIO->GetImageSizeInBytes()];

      m_ImageIO->Read(loadBuffer);
      
      itkDebugMacro(<< "Buffer conversion required from: "
                    << m_ImageIO->GetComponentTypeInfo().name()
                    << " to: "
                    << typeid(ITK_TYPENAME ConvertPixelTraits::ComponentType).name());

      this->DoConvertBuffer(loadBuffer, region.GetNumberOfPixels());
      delete [] loadBuffer;
      }

    ImageRegionIterator<TDeformationField> Id( this->GetOutput(),
      this->GetOutput()->GetLargestPossibleRegion() );
    ImageRegionIterator<TImage> It( this->m_Image,
      this->m_Image->GetLargestPossibleRegion() );

    Id.GoToBegin();
    It.GoToBegin();
    while ( !Id.IsAtEnd() || !It.IsAtEnd() )
      {
      DeformationFieldPixelType V = Id.Get();
      V[i] = static_cast< typename DeformationFieldPixelType::ValueType>( It.Get() );
      Id.Set( V );     
      ++Id;
      ++It;
      }
    }  
}




template <class TImage, class TDeformationField, class ConvertPixelTraits>
void 
DeformationFieldReader<TImage, TDeformationField, ConvertPixelTraits>
::DoConvertBuffer(void* inputData,
                  unsigned long numberOfPixels)
{
  // get the pointer to the destination buffer
  ImagePixelType *imageData =
    this->m_Image->GetPixelContainer()->GetBufferPointer();


  // TODO:
  // Pass down the PixelType (RGB, VECTOR, etc.) so that any vector to
  // scalar conversion be type specific. i.e. RGB to scalar would use
  // a formula to convert to luminance, VECTOR to scalar would use
  // vector magnitude.
  
  
// Create a macro as this code is a bit lengthy and repetitive
// if the ImageIO pixel type is typeid(type) then use the ConvertPixelBuffer
// class to convert the data block to TImage's pixel type
// see DefaultConvertPixelTraits and ConvertPixelBuffer

 
// The first else if block applies only to images of type itk::VectorImage  
// VectorImage needs to copy out the buffer differently.. The buffer is of
// type InternalPixelType, but each pixel is really 'k' consecutive pixels.

#define ITK_CONVERT_BUFFER_IF_BLOCK(type)               \
 else if( m_ImageIO->GetComponentTypeInfo() == typeid(type) )   \
   {                                                   \
   if( strcmp( this->GetOutput()->GetNameOfClass(), "VectorImage" ) == 0 ) \
     { \
     ConvertPixelBuffer<                                 \
      type,                                             \
      ImagePixelType,                             \
      ConvertPixelTraits                                \
      >                                                 \
      ::ConvertVectorImage(                             \
        static_cast<type*>(inputData),                  \
        m_ImageIO->GetNumberOfComponents(),             \
        imageData,                                     \
        numberOfPixels);                                \
     } \
   else \
     { \
     ConvertPixelBuffer<                                 \
      type,                                             \
      ImagePixelType,                             \
      ConvertPixelTraits                                \
      >                                                 \
      ::Convert(                                        \
        static_cast<type*>(inputData),                  \
        m_ImageIO->GetNumberOfComponents(),             \
        imageData,                                     \
        numberOfPixels);                                \
      } \
    }
  if(0)
    {
    }
  ITK_CONVERT_BUFFER_IF_BLOCK(unsigned char)
    ITK_CONVERT_BUFFER_IF_BLOCK(char)
    ITK_CONVERT_BUFFER_IF_BLOCK(unsigned short)
    ITK_CONVERT_BUFFER_IF_BLOCK( short)
    ITK_CONVERT_BUFFER_IF_BLOCK(unsigned int)
    ITK_CONVERT_BUFFER_IF_BLOCK( int)
    ITK_CONVERT_BUFFER_IF_BLOCK(unsigned long)
    ITK_CONVERT_BUFFER_IF_BLOCK( long)
    ITK_CONVERT_BUFFER_IF_BLOCK(float)
    ITK_CONVERT_BUFFER_IF_BLOCK( double)
    else
      {
      DeformationFieldReaderException e(__FILE__, __LINE__);
      OStringStream msg;
      msg <<"Couldn't convert component type: "
          << std::endl << "    "
          << m_ImageIO->GetComponentTypeAsString(m_ImageIO->GetComponentType())
          << std::endl << "to one of: "
          << std::endl << "    " << typeid(unsigned char).name()
          << std::endl << "    " << typeid(char).name()
          << std::endl << "    " << typeid(unsigned short).name()
          << std::endl << "    " << typeid(short).name()
          << std::endl << "    " << typeid(unsigned int).name()
          << std::endl << "    " << typeid(int).name()
          << std::endl << "    " << typeid(unsigned long).name()
          << std::endl << "    " << typeid(long).name()
          << std::endl << "    " << typeid(float).name()
          << std::endl << "    " << typeid(double).name()
          << std::endl;
      e.SetDescription(msg.str().c_str());
      e.SetLocation(ITK_LOCATION);
      throw e;
      return;
      }
#undef ITK_CONVERT_BUFFER_IF_BLOCK
}


} //namespace ITK

#endif
