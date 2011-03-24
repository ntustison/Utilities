/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDeformationFieldWriter.txx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:21:04 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkDeformationFieldWriter_txx
#define _itkDeformationFieldWriter_txx

#include "itkDeformationFieldWriter.h"
#include "itkImageFileWriter.h"
#include "itkDataObject.h"
#include "itkObjectFactoryBase.h"
#include "itkImageIOFactory.h"
#include "itkCommand.h"
#include "vnl/vnl_vector.h"
#include "itkVectorImage.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

namespace itk
{

//---------------------------------------------------------
template <class TDeformationField, class TImage>
DeformationFieldWriter<TDeformationField, TImage>
::DeformationFieldWriter():
  m_FileName(""),
  m_ImageIO(0), m_UserSpecifiedImageIO(false),
  m_UserSpecifiedIORegion(false)
{
  m_UseCompression = false;
  m_UseInputMetaDataDictionary = true;
  m_FactorySpecifiedImageIO = false;
  m_UseAvantsNamingConvention = false;
}


//---------------------------------------------------------
template <class TDeformationField, class TImage>
DeformationFieldWriter<TDeformationField, TImage>
::~DeformationFieldWriter()
{
}

//---------------------------------------------------------
template <class TDeformationField, class TImage>
void 
DeformationFieldWriter<TDeformationField, TImage>
::SetInput(const DeformationFieldType *input)
{
  this->ProcessObject::SetNthInput(0, 
                                   const_cast<TDeformationField *>(input ) );
}


//---------------------------------------------------------
template <class TDeformationField, class TImage>
const typename DeformationFieldWriter<TDeformationField, TImage>::DeformationFieldType *
DeformationFieldWriter<TDeformationField, TImage>
::GetInput(void)
{
  if (this->GetNumberOfInputs() < 1)
    {
    return 0;
    }
  
  return static_cast<TDeformationField*>
    (this->ProcessObject::GetInput(0));
}
  
//---------------------------------------------------------
template <class TDeformationField, class TImage>
const typename DeformationFieldWriter<TDeformationField, TImage>::DeformationFieldType *
DeformationFieldWriter<TDeformationField, TImage>
::GetInput(unsigned int idx)
{
  return static_cast<TDeformationField*>
    (this->ProcessObject::GetInput(idx));
}

//---------------------------------------------------------
template <class TDeformationField, class TImage>
void 
DeformationFieldWriter<TDeformationField, TImage>
::SetIORegion (const ImageIORegion& region) 
{
  itkDebugMacro("setting IORegion to " << region );
  if ( m_IORegion != region)
    {
    m_IORegion = region;
    this->Modified();
    m_UserSpecifiedIORegion = true;
    }
} 

//---------------------------------------------------------
template <class TDeformationField, class TImage>
void 
DeformationFieldWriter<TDeformationField, TImage>
::GenerateData(void)
{

  itkDebugMacro(<<"Writing file: " << m_ComponentImageFileName);
  
  // Make sure that the image is the right type and no more than 
  // four components.
  typedef typename ImageType::PixelType ScalarType;

  if( strcmp( m_Image->GetNameOfClass(), "VectorImage" ) == 0 ) 
    {
    typedef typename ImageType::InternalPixelType VectorImageScalarType;
    m_ImageIO->SetPixelTypeInfo( typeid(VectorImageScalarType) );
    
    typedef typename ImageType::AccessorFunctorType AccessorFunctorType;
    m_ImageIO->SetNumberOfComponents( AccessorFunctorType::GetVectorLength( m_Image ) );
    }
  else
    {
    // Set the pixel and component type; the number of components.
    m_ImageIO->SetPixelTypeInfo(typeid(ScalarType));  
    }

  // Setup the image IO for writing.
  //
  m_ImageIO->SetFileName(m_ComponentImageFileName.c_str());

  //okay, now extract the data as a raw buffer pointer
  const void* dataPtr = (const void*) this->m_Image->GetBufferPointer();

  m_ImageIO->Write(dataPtr);

}
//---------------------------------------------------------
template <class TDeformationField, class TImage>
void 
DeformationFieldWriter<TDeformationField, TImage>
::Write()
{
  std::string::size_type Pos = this->m_FileName.rfind( "." );
  std::string extension( this->m_FileName, Pos, this->m_FileName.length()-1 );

  typedef VectorIndexSelectionCastImageFilter
       <DeformationFieldType, ImageType> SelectorType;
  typename SelectorType::Pointer selector = SelectorType::New();
  selector->SetInput( this->GetInput() );

  unsigned int dimension = itk::GetVectorDimension
      <typename DeformationFieldType::PixelType>::VectorDimension;
  for ( unsigned int i = 0; i < dimension; i++ )
    {
    selector->SetIndex( i );
    selector->Update();

    std::string filename( this->m_FileName, 0, Pos );    

    if ( this->m_UseAvantsNamingConvention )
      {
      switch ( i )
        {
        case 0:
          filename += std::string( "xvec" );
          break;
        case 1:
          filename += std::string( "yvec" );
          break;
        case 2:
          filename += std::string( "zvec" );
          break;
        default:
          filename += std::string( "you_are_screwed_vec" );
          break;
        }  
      }
    else
      {
      itk::OStringStream buf;
      buf << i;
      filename += ( std::string( "." )  + std::string( buf.str().c_str() ) );
      }
    filename += extension;
    
    m_ComponentImageFileName = filename;

    ImageType *input = selector->GetOutput();
    this->m_Image = selector->GetOutput();

    itkDebugMacro( <<"Writing an image file" );

    // Make sure input is available
    if ( input == 0 )
      {
      itkExceptionMacro(<< "No input to writer!");
      }

    // Make sure that we can write the file given the name
    //
    if ( filename == "" )
      {
      itkExceptionMacro(<<"No filename was specified");
      }

    if ( m_ImageIO.IsNull() ) //try creating via factory
      {
      itkDebugMacro(<<"Attempting factory creation of ImageIO for file: " 
                    << filename);
      m_ImageIO = ImageIOFactory::CreateImageIO( filename.c_str(), 
                                                 ImageIOFactory::WriteMode );
      m_FactorySpecifiedImageIO = true;
      }
    else
      {
      if( m_FactorySpecifiedImageIO && !m_ImageIO->CanWriteFile( filename.c_str() ) )
        {
        itkDebugMacro(<<"ImageIO exists but doesn't know how to write file:" 
                      << m_FileName );
        itkDebugMacro(<<"Attempting creation of ImageIO with a factory for file:"
                      << m_FileName);
        m_ImageIO = ImageIOFactory::CreateImageIO( filename.c_str(), 
                                                   ImageIOFactory::WriteMode );
        m_FactorySpecifiedImageIO = true;
        }
      }

    if ( m_ImageIO.IsNull() )
      {
      ImageFileWriterException e(__FILE__, __LINE__);
      OStringStream msg;
      msg << " Could not create IO object for file "
          << filename.c_str() << std::endl;
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
      e.SetDescription(msg.str().c_str());
      e.SetLocation(ITK_LOCATION);
      throw e;
      }

    // NOTE: this const_cast<> is due to the lack of const-correctness
    // of the ProcessObject.
    ImageType *nonConstImage = input;

    typedef typename TImage::RegionType   RegionType;

    if ( ! m_UserSpecifiedIORegion )
      {
      // Make sure the data is up-to-date.
      if( nonConstImage->GetSource() )
        {
        nonConstImage->GetSource()->UpdateLargestPossibleRegion();
        }
      // Write the whole image
      ImageIORegion ioRegion(TImage::ImageDimension);
      RegionType region = this->m_Image->GetLargestPossibleRegion();

      for(unsigned int i=0; i<TDeformationField::ImageDimension; i++)
        {
        ioRegion.SetSize(i,region.GetSize(i));
        ioRegion.SetIndex(i,region.GetIndex(i));
        }
      m_IORegion = ioRegion; //used by GenerateData
      }
    else
      {
      nonConstImage->Update();
      }

    // Setup the ImageIO
    //
    m_ImageIO->SetNumberOfDimensions(TImage::ImageDimension);
    RegionType region = this->m_Image->GetLargestPossibleRegion();
    const typename TImage::SpacingType& spacing = this->m_Image->GetSpacing();
    const typename TImage::PointType& origin = this->m_Image->GetOrigin();
    const typename TImage::DirectionType& direction = this->m_Image->GetDirection();

    for(unsigned int i=0; i<TDeformationField::ImageDimension; i++)
      {
      m_ImageIO->SetDimensions(i,region.GetSize(i));
      m_ImageIO->SetSpacing(i,spacing[i]);
      m_ImageIO->SetOrigin(i,origin[i]);
      vnl_vector< double > axisDirection(TDeformationField::ImageDimension);
  // Please note: direction cosines are stored as columns of the
  // direction matrix
      for(unsigned int j=0; j<TImage::ImageDimension; j++)
        {
        axisDirection[j] = direction[j][i];
        }
      m_ImageIO->SetDirection( i, axisDirection );
      }

    m_ImageIO->SetUseCompression(m_UseCompression);
    m_ImageIO->SetIORegion(m_IORegion);
    if( m_UseInputMetaDataDictionary )
      {
      m_ImageIO->SetMetaDataDictionary(input->GetMetaDataDictionary());
      }
    // Notify start event observers
    this->InvokeEvent( StartEvent() );

    // Actually do something
    this->GenerateData();
    
    // Notify end event observers
    this->InvokeEvent( EndEvent() );

    // Release upstream data if requested
    if ( input->ShouldIReleaseData() )
      {
      nonConstImage->ReleaseData();
      }
    }  
}


//---------------------------------------------------------
template <class TDeformationField, class TImage>
void 
DeformationFieldWriter<TDeformationField, TImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "File Name: " 
     << (m_FileName.data() ? m_FileName.data() : "(none)") << std::endl;

  os << "Number of vector components (i.e. number of images created): " << 
    itk::GetVectorDimension<typename DeformationFieldType::PixelType>::VectorDimension
    << std::endl;

  os << indent << "Image IO: ";
  if ( m_ImageIO.IsNull() )
    {
    os << "(none)\n";
    }
  else
    {
    os << m_ImageIO << "\n";
    }

  os << indent << "IO Region: " << m_IORegion << "\n";


  if (m_UseCompression)
    {
    os << indent << "Compression: On\n";
    }
  else
    {
    os << indent << "Compression: Off\n";
    }

  if (m_UseInputMetaDataDictionary)
    {
    os << indent << "UseInputMetaDataDictionary: On\n";
    }
  else
    {
    os << indent << "UseInputMetaDataDictionary: Off\n";
    }

  if (m_FactorySpecifiedImageIO)
    {
    os << indent << "FactorySpecifiedmageIO: On\n";
    }
  else
    {
    os << indent << "FactorySpecifiedmageIO: Off\n";
    }
}

} // end namespace itk

#endif
