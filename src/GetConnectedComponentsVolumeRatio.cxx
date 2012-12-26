/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: GetConnectedComponentsVolumeRatio.cxx,v $
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

template <unsigned int ImageDimension>
int GetConnectedComponents(int argc, char* argv[] )
{
  typedef int PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<float, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename RealImageType::Pointer output = RealImageType::New();
  output->CopyInformation( reader->GetOutput() );
  output->SetRegions( reader->GetOutput()->GetRequestedRegion() );
  output->Allocate();
  output->FillBuffer( 0.0 );

  typedef itk::RelabelComponentImageFilter<ImageType, ImageType> RelabelerType;
  typename RelabelerType::Pointer relabeler = RelabelerType::New();
  relabeler->SetInput( reader->GetOutput() );
  relabeler->Update();

  for ( unsigned int i = 1; i <= relabeler->GetNumberOfObjects(); i++ )
    {
    typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput( relabeler->GetOutput() );
    thresholder->SetLowerThreshold( i );
    thresholder->SetUpperThreshold( i );
    thresholder->SetInsideValue( 1 );
    thresholder->SetOutsideValue( 0 );
    thresholder->Update();

    typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentType;
    typename ConnectedComponentType::Pointer filter = ConnectedComponentType::New();
    filter->SetInput( thresholder->GetOutput() );
    filter->Update();

    typename RelabelerType::Pointer relabeler2 = RelabelerType::New();
    relabeler2->SetInput( filter->GetOutput() );
    relabeler2->Update();

    itk::ImageRegionIterator<ImageType> It( relabeler->GetOutput(),
      relabeler->GetOutput()->GetRequestedRegion() );
    itk::ImageRegionIterator<ImageType> It2( relabeler2->GetOutput(),
      relabeler2->GetOutput()->GetRequestedRegion() );
    itk::ImageRegionIterator<RealImageType> ItO( output, output->GetRequestedRegion() );

    for( It.GoToBegin(), It2.GoToBegin(), ItO.GoToBegin(); !It.IsAtEnd(); ++It, ++It2, ++ItO )
      {
      int label = It2.Get();
      if( label != 0 )
        {
        ItO.Set( relabeler2->GetSizeOfObjectInPhysicalUnits( label ) / relabeler->GetSizeOfObjectInPhysicalUnits( i ) );
        }
      }
    }

  typedef itk::ImageFileWriter<RealImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( output );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension "
              << "inputImage outputImage" << std::endl;
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
