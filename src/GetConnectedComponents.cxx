/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: GetConnectedComponents.cxx,v $
  Language:  C++
  Date:      $Date: 2010/02/22 15:34:48 $
  Version:   $Revision: 1.4 $

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
#include "itkImageRegionIterator.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include <vector>
#include <algorithm>

template <unsigned int ImageDimension>
int GetConnectedComponents(int argc, char* argv[] )
{
  typedef int PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::RelabelComponentImageFilter<ImageType, ImageType> RelabelerType;
  typename RelabelerType::Pointer relabeler = RelabelerType::New();

  if( argc > 4 && atoi( argv[4] ) == 1 )
    {
    std::vector<PixelType> labels;

    itk::ImageRegionIterator<ImageType> It( reader->GetOutput(),
      reader->GetOutput()->GetLargestPossibleRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if( It.Get() != 0 &&
        std::find( labels.begin(), labels.end(), It.Get() ) == labels.end() )
        {
        labels.push_back( It.Get() );
        }
      }
    std::sort( labels.begin(), labels.end() );

    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if( It.Get() != 0 )
        {
        typename std::vector<PixelType>::iterator it =
          std::find( labels.begin(), labels.end(), It.Get() );
        It.Set( static_cast<PixelType>( it - labels.begin() ) + 1 );
        }
      }

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[3] );
    writer->SetInput( reader->GetOutput() );
    writer->Update();

    return EXIT_SUCCESS;
    }
  else
    {
    typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentType;
    typename ConnectedComponentType::Pointer filter = ConnectedComponentType::New();
    filter->SetInput( reader->GetOutput() );
    filter->Update();

    relabeler->SetInput( filter->GetOutput() );
    relabeler->Update();
    }


  unsigned int numberOfObjects = relabeler->GetNumberOfObjects();
  std::cout << "NumberOfObjects: " << numberOfObjects << std::endl;
  for ( unsigned int i = 1; i <= numberOfObjects; i++ )
    {
    std::cout << "  Object[" << i << "] consists of "
              << relabeler->GetSizeOfObjectsInPixels()[i-1]
              << std::endl;
    }

  if( argc > 5 && atoi( argv[5] ) > 1 )
    {
    itk::ImageRegionIterator<ImageType> It( relabeler->GetOutput(),
      relabeler->GetOutput()->GetLargestPossibleRegion() );
    for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      typename ImageType::PixelType label = It.Get();
      if ( label != 0 && relabeler->GetSizeOfObjectsInPixels()[label-1] < static_cast<unsigned int>( atoi( argv[5] ) ) )
        {
        It.Set( 0 );
        }
      }
    }

  if( argc > 6 )
    {
    typename ReaderType::Pointer maskReader = ReaderType::New();
    maskReader->SetFileName( argv[6] );
    maskReader->Update();

    std::vector<typename ImageType::PixelType> maskedLabels;

    itk::ImageRegionIterator<ImageType> ItM( maskReader->GetOutput(),
      maskReader->GetOutput()->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<ImageType> It( relabeler->GetOutput(),
      relabeler->GetOutput()->GetLargestPossibleRegion() );
    for( ItM.GoToBegin(), It.GoToBegin(); !ItM.IsAtEnd(); ++ItM, ++It )
      {
      typename ImageType::PixelType label = It.Get();
      if( ItM.Get() > 0 && label > 0 )
        {
        if( std::find( maskedLabels.begin(), maskedLabels.end(), label )
          == maskedLabels.end() )
          {
          maskedLabels.push_back( label );
          }
        }
      }
    for( ItM.GoToBegin(), It.GoToBegin(); !ItM.IsAtEnd(); ++ItM, ++It )
      {
      typename ImageType::PixelType label = It.Get();
      if( label > 0 )
        {
        if( std::find( maskedLabels.begin(), maskedLabels.end(), label )
          != maskedLabels.end() )
          {
          It.Set( 0 );
          }
        }
      }
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( relabeler->GetOutput() );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension "
              << "inputImage outputImage [relabelOnly] "
              << "[smallestAllowableSize] [maskImage]" << std::endl;
    std::cerr << "The [maskImage] option functions as follows.  After calculating "
              << "the connected component image, any label in the masked region is "
              << "removed." << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GetConnectedComponents<2>( argc, argv );
     break;
   case 3:
     GetConnectedComponents<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
