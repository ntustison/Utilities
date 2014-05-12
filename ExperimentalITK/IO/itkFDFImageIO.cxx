/*  Copyright (C) 2004 Glenn Pierce.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
#include "itkFDFImageIO.h"
#include "itkFDFCommonImageIO.h"

#include "itkByteSwapper.h"
#include "itkRGBPixel.h"
#include "itkRGBAPixel.h"
#include <stdio.h>
#include <fstream>

namespace itk
{

bool FDFImageIO::CanReadFile(const char* file)
{
  this->SetFileName(file);

  resolutions = new float[3];

  // First check the extension
  std::string filename = file;
  if(  filename == "" )
    {
    itkDebugMacro(<<"No filename specified.");
    return false;
    }

  bool extensionFound = false;
  std::string::size_type FDFPos = filename.rfind(".fdf");
  if ((FDFPos != std::string::npos)
      && (FDFPos == filename.length() - 4))
    {
    extensionFound = true;
    }

  FDFPos = filename.rfind(".FDF");
  if ((FDFPos != std::string::npos)
      && (FDFPos == filename.length() - 4))
    {
    extensionFound = true;
    }

  if( !extensionFound )
    {
    itkDebugMacro(<<"The filename extension is not recognized");
    return false;
    }

  std::ifstream inFile;
  inFile.open(m_FileName.c_str(), std::ios::in | std::ios::binary );
  if( !inFile )
    {
    ExceptionObject exception(__FILE__, __LINE__);
    std::string msg = "File \"" + m_FileName + "\" cannot be read.";
    exception.SetDescription(msg.c_str());
    throw exception;
    }

  // Check for a neccessary header variable
  // HERE

  return true;
}

void FDFImageIO::ReadImageInformation()
{
  if(!this->CanReadFile(m_FileName.c_str()))
    RAISE_EXCEPTION();

  std::string line;
  std::vector<std::string> tokens;
  std::string type, name, value;
  int matrix_size = 1;

  std::ifstream inFile(m_FileName.c_str(), std::ios::in | std::ios::binary);

  // Check if there was an error opening the file
  if (!inFile)
  {
    std::cout << "Unable to open the file\n";
    RAISE_EXCEPTION();
  }

  while (getline(inFile, line, '\n')) {

    if ( line == "\0" ) {
      break;
    }

    // Formats the lines in the FDF header such as removing whitespace between {}
    line = ParseLine(line);

    Tokenize(line, tokens, " ;");

    if(tokens.size() == 4) {

            type = tokens[0];
            name = tokens[1];
            value = tokens[3];

            if (name == "spatial_rank") {
                this->spatial_rank = value;
            }

            if (name == "matrix") {
                StringToVector(value, matrix);

                // Set the number of dimensions
                this->SetNumberOfDimensions(matrix.size());

                for(int i=0; i<matrix.size(); i++) {
                    matrix_size *= this->matrix[i];
                    // Set the size of each dimension
                    this->SetDimensions(i,this->matrix[i]);
                }
            }

            if (name == "span") {
                StringToVector(value, span);
            }

            if (name == "roi") {
                StringToVector(value, roi);
            }

            if (name == "location") {
                StringToVector(value, location);
            }

            // Get the binary data type
            if( name == "storage" )
              {
              this->storage = value;

              this->SetPixelType( SCALAR );

              if( value == "double" )
                {
                this->SetComponentType( DOUBLE );
                }
              else if( value == "float" )
                {
                this->SetComponentType( FLOAT );
                }
              else if( value == "long" )
                {
                this->SetComponentType( LONG );
                }
              else if( value == "unsigned long" )
                {
                this->SetComponentType( ULONG );
                }
              else if( value == "int" )
                {
                this->SetComponentType( INT );
                }
              else if( value == "unsigned int" )
                {
                this->SetComponentType( UINT );
                }
              else if( value == "short" )
                {
                this->SetComponentType( SHORT );
                }
              else if( value == "unsigned short" )
                {
                this->SetComponentType( USHORT );
                }
              else if( value == "char" )
                {
                this->SetComponentType( CHAR );
                }
              else if( value == "unsigned char" )
                {
                this->SetComponentType( UCHAR );
                }
              else
                {
                itkExceptionMacro( "Unknown component type: " << value );
                }
              }

            // Get the bits
            if (name == "bits") {
                ConvertFromString (value, this->bits);
            }

            // Get the checksum
            if (name == "checksum") {
                ConvertFromString (value, this->checksum);
            }

        }

        tokens.clear();
    }

    inFile.seekg (0, std::ios::end);
    long int fileSize = inFile.tellg();
    this->m_InputPosition = fileSize - (matrix_size * 4);
    this->resolutions[0] = (this->roi[0] * 10 ) / this->matrix[0];
    this->resolutions[1] = (this->roi[1] * 10 ) / this->matrix[1];
    this->resolutions[2] = (this->roi[2] * 10 ) / this->matrix[2];
}


void FDFImageIO::ReadVolume(void*)
{

}

// const std::type_info& FDFImageIO::GetPixelType() const
// {
//   switch(m_PixelType)
//     {
//     case UCHAR:
//       return typeid(unsigned char);
//     case USHORT:
//       return typeid(unsigned short);
//     case CHAR:
//       return typeid(char);
//     case SHORT:
//       return typeid(short);
//     case UINT:
//       return typeid(unsigned int);
//     case INT:
//       return typeid(int);
//     case ULONG:
//       return typeid(unsigned long);
//     case LONG:
//       return typeid(long);
//     case FLOAT:
//       return typeid(float);
//     case DOUBLE:
//       return typeid(double);
//     case RGB:
//       return typeid(RGBPixel<unsigned char>);
//     case RGBA:
//       return typeid(RGBAPixel<unsigned char>);
//     default:
//     {
//     itkExceptionMacro ("Invalid type: " << m_PixelType << ", only unsigned char, unsigned short, RGB<unsigned char> are allowed.");
//     return this->ConvertToTypeInfo(m_PixelType);
//     }
//     case UNKNOWN:
//       itkExceptionMacro ("Unknown pixel type: " << m_PixelType);
//     }
//   return typeid(ImageIOBase::UnknownType);
// }

// unsigned int FDFImageIO::GetComponentSize() const
// {
//     switch(m_PixelType)
//     {
//     case UCHAR:
//       return sizeof(unsigned char);
//     case USHORT:
//       return sizeof(unsigned short);
//     case CHAR:
//       return sizeof(char);
//     case SHORT:
//       return sizeof(short);
//     case UINT:
//       return sizeof(unsigned int);
//     case INT:
//       return sizeof(int);
//     case ULONG:
//       return sizeof(unsigned long);
//     case LONG:
//       return sizeof(long);
//     case FLOAT:
//       return sizeof(float);
//     case DOUBLE:
//       return sizeof(double);
//     case RGB:
//       return sizeof(unsigned char);
//     case RGBA:
//       return sizeof(unsigned char);
//     case UNKNOWNPIXELTYPE:
//     default:
//     {
//     itkExceptionMacro ("Invalid type: " << m_PixelType
//                        << ", only unsigned char and unsigned short are allowed.");
//     return 0;
//     }
//     }
//
//     return 1;
// }

void FDFImageIO::Read(void* buffer)
{
  unsigned int dimensions = this->GetNumberOfDimensions();
  unsigned int numberOfPixels = 1;

  for( unsigned int dim=0; dim< dimensions; dim++ )
  {
    numberOfPixels *= m_Dimensions[ dim ];
  }

  std::ifstream inFile(m_FileName.c_str(), std::ios::in | std::ios::binary);

  // Check if there was an error opening the file
  if (!inFile)
  {
    RAISE_EXCEPTION();
  }

  inFile.seekg( m_InputPosition );

  if (!inFile)
  {
    RAISE_EXCEPTION();
  }

  char * p = static_cast<char *>(buffer);

  inFile.read( p, this->GetImageSizeInBytes() );

  bool success = !inFile.bad();
  inFile.close();
  if( !success )
  {
    itkExceptionMacro("Error reading image data.");
  }

  SwapBytesIfNecessary( buffer, numberOfPixels );

}


FDFImageIO::FDFImageIO()
{
  this->SetNumberOfDimensions(2);
  this->m_FileType = Binary;
  this->m_ByteOrder = BigEndian;
}

FDFImageIO::~FDFImageIO() {}

void FDFImageIO::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "PixelType " << m_PixelType << "\n";
  os << indent << "Start of image in bytes from start of file " << this->m_InputPosition << "\n";
  os << indent << "Number of pixels in image: " << this->GetImageSizeInPixels() << "\n";
  os << indent << "Image size in bytes: " << this->GetImageSizeInBytes() << "\n";
  os << indent << "Checksum: " << this->checksum << "\n";
  os << indent << "Storage: " << this->storage << "\n";
  os << indent << "Spatial Rank: " << this->spatial_rank << "\n";
  os << indent << "Bits: " << this->bits << "\n";
  os << indent; PrintVector(os, "Matrix", this->matrix);
  os << indent; PrintVector(os, "Location", this->location);
  os << indent; PrintVector(os, "ROI", this->roi);
  os << indent; PrintVector(os, "Span", this->span);
  os << indent << "resolutions: " << this->resolutions[0] << ", " << this->resolutions[1]
                                                          << ", " << this->resolutions[2] << "\n";
}

bool FDFImageIO::CanWriteFile( const char * name )
{
  //not possible to write a fdf file
  return false;
}

void FDFImageIO::SwapBytesIfNecessary( void* buffer, unsigned long numberOfPixels )
{
  switch(m_PixelType)
    {
    case CHAR:
    {
    if ( m_ByteOrder == LittleEndian )
      {
      ByteSwapper<int>::SwapRangeFromSystemToLittleEndian(
        (int*)buffer, numberOfPixels );
      }
    else if ( m_ByteOrder == BigEndian )
      {
      ByteSwapper<int>::SwapRangeFromSystemToBigEndian(
        (int *)buffer, numberOfPixels );
      }
    break;
    }
    case FLOAT:
    {
    if ( m_ByteOrder == LittleEndian )
      {
      ByteSwapper<float>::SwapRangeFromSystemToLittleEndian(
        (float *)buffer, numberOfPixels );
      }
    else if ( m_ByteOrder == BigEndian )
      {
      ByteSwapper<float>::SwapRangeFromSystemToBigEndian(
        (float *)buffer, numberOfPixels );
      }
    break;
    }
    case UCHAR:
    {
    if ( m_ByteOrder == LittleEndian )
      {
      ByteSwapper<unsigned int>::SwapRangeFromSystemToLittleEndian(
        (unsigned int*)buffer, numberOfPixels );
      }
    else if ( m_ByteOrder == BigEndian )
      {
      ByteSwapper<unsigned int>::SwapRangeFromSystemToBigEndian(
        (unsigned int *)buffer, numberOfPixels );
      }
    break;
    }
    case SHORT:
    {
    if ( m_ByteOrder == LittleEndian )
      {
      ByteSwapper<short>::SwapRangeFromSystemToLittleEndian(
        (short*)buffer, numberOfPixels );
      }
    else if ( m_ByteOrder == BigEndian )
      {
      ByteSwapper<short>::SwapRangeFromSystemToBigEndian(
        (short *)buffer, numberOfPixels );
      }
    break;
    }
    case USHORT:
    {
    if ( m_ByteOrder == LittleEndian )
      {
      ByteSwapper<unsigned short>::SwapRangeFromSystemToLittleEndian(
        (unsigned short*)buffer, numberOfPixels );
      }
    else if ( m_ByteOrder == BigEndian )
      {
      ByteSwapper<unsigned short>::SwapRangeFromSystemToBigEndian(
        (unsigned short *)buffer, numberOfPixels );
      }
    break;
    }
    default:
      ExceptionObject exception(__FILE__, __LINE__);
      exception.SetDescription("Pixel Type Unknown");
      throw exception;
    }
}

void FDFImageIO::WriteImageInformation(void)
{
  //not possible to write a fdf file
}

void FDFImageIO::Write(const void* buffer)
{
  //not possible to write a fdf file
}

} // end namespace itk

