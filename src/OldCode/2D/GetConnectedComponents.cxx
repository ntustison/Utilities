/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: GetConnectedComponents.cxx,v $
  Language:  C++
  Date:      $Date: 2008/05/03 01:47:04 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include "global.h"

int main(int argc, char* argv[] )
{
  if ( argc != 3 )
    {
    std::cout << "Usage: GetConnectedComponents binary_image output_image" << std::endl;
    exit( 0 );
    }
    
  typedef itk::Image<PixelType, ImageDimension> ImageType;
    
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();  
  reader->SetFileName( argv[1] );
        
  typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentType;
  ConnectedComponentType::Pointer filter = ConnectedComponentType::New();
  filter->SetInput( reader->GetOutput() );

  typedef itk::RelabelComponentImageFilter<ImageType, ImageType> RelabelerType;
  RelabelerType::Pointer relabeler = RelabelerType::New();
  relabeler->SetInput( filter->GetOutput() );
  relabeler->Update();
  
  std::cout << "NumberOfObjects: " << relabeler->GetNumberOfObjects() << std::endl;
  for ( unsigned int i = 0; i < relabeler->GetNumberOfObjects(); i++ )
    {
    std::cout << "  Object[" << i+1 << "] consists of " << relabeler->GetSizeOfObjectsInPixels()[i] << " pixels." << std::endl;
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( relabeler->GetOutput() );
  writer->Update();
   
  return 0;
}
