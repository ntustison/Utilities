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

#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include <vector>
#include <algorithm>

template <unsigned int ImageDimension>
int GetConnectedComponents(int argc, char* argv[] )
{
  typedef unsigned int PixelType;
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
    }
  else
    {
    typename ImageType::Pointer output = ImageType::New();
    output->CopyInformation( reader->GetOutput() );
    output->SetRegions( reader->GetOutput()->GetRequestedRegion() );
    output->Allocate();
    output->FillBuffer( 0 );

    typedef itk::RelabelComponentImageFilter<ImageType, ImageType> RelabelerType;
    typename RelabelerType::Pointer relabeler = RelabelerType::New();
    relabeler->SetInput( reader->GetOutput() );
    relabeler->Update();

    unsigned int count = 0;

    for ( unsigned int i = 1; i <= relabeler->GetNumberOfObjects(); i++ )
      {
      typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholderType;
      typename ThresholderType::Pointer thresholder = ThresholderType::New();
      thresholder->SetInput( relabeler->GetOutput() );
      thresholder->SetLowerThreshold( i );
      thresholder->SetUpperThreshold( i );
      thresholder->SetOutsideValue( 0 );
      thresholder->SetInsideValue( 1 );
      thresholder->Update();

      typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentType;
      typename ConnectedComponentType::Pointer filter = ConnectedComponentType::New();
      filter->SetInput( thresholder->GetOutput() );
      filter->Update();

      typename RelabelerType::Pointer relabeler2 = RelabelerType::New();
      relabeler2->SetInput( filter->GetOutput() );
      relabeler2->Update();

      itk::ImageRegionIterator<ImageType> It2( relabeler2->GetOutput(),
        relabeler2->GetOutput()->GetRequestedRegion() );

      itk::ImageRegionIterator<ImageType> ItO( output, output->GetRequestedRegion() );

      for( It2.GoToBegin(), ItO.GoToBegin(); !It2.IsAtEnd(); ++It2, ++ItO )
        {
        unsigned int label = It2.Get();
        if( label != 0 )
          {
          ItO.Set( label + count );
          }
        }
      count += relabeler2->GetNumberOfObjects();
      }

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[3] );
    writer->SetInput( output );
    writer->Update();
    }

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension "
              << "inputImage outputImage [relabelOnly] "
              << std::endl;
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
