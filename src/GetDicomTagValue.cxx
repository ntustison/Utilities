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

#include <vector>
#include <itksys/SystemTools.hxx>

#include <string>

int main( int argc, char* argv[] )
{

  if( argc < 2 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage tagID" << std::endl;


    return EXIT_FAILURE;
    }

  typedef signed short    PixelType;
  const unsigned int      Dimension = 2;

  typedef itk::Image< PixelType, Dimension >      ImageType;
  typedef itk::ImageFileReader< ImageType >       ReaderType;
  typedef itk::GDCMImageIO                        ImageIOType;
  ImageIOType::Pointer gdcmIO = ImageIOType::New();

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetImageIO( gdcmIO );
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

  const itk::MetaDataDictionary & dictionary = reader->GetOutput()->GetMetaDataDictionary();
  std::string value;

  bool exists = itk::ExposeMetaData<std::string>( dictionary, argv[2], value );

  std::string output = std::string( argv[2] ) + std::string( "," ) + value;

  if( exists )
    {
    std::cout << output << std::endl;
    }
  else
    {
    std::cout << argv[2] << std::endl;
    }

  return EXIT_SUCCESS;
}
