#include <stdio.h>

#include "itkBinaryThresholdImageFilter.h"
#include "itkDerivativeImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkZeroCrossingImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " labelImage outputImage direction" << std::endl;
    exit( 1 );
    }

  typedef itk::Image<float, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  ImageType::Pointer output = ImageType::New();
  output->SetOrigin( reader->GetOutput()->GetOrigin() );
  output->SetSpacing( reader->GetOutput()->GetSpacing() );
  output->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  output->Allocate();
  output->FillBuffer( 0 );

  for ( unsigned int label = 1; label <= 2; label++ )
    {
    typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> BinaryType;
    BinaryType::Pointer binary = BinaryType::New();
    binary->SetInput( reader->GetOutput() );
    binary->SetUpperThreshold( label );
    binary->SetLowerThreshold( label );
    binary->SetInsideValue( 1 );
    binary->SetOutsideValue( 0 );
    binary->Update();

    typedef itk::SignedMaurerDistanceMapImageFilter<ImageType, ImageType> DistanceType;
    DistanceType::Pointer distance = DistanceType::New();
    distance->SetInput( binary->GetOutput() );
    distance->Update();

    typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
    NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 2 );

    unsigned int direction = atoi( argv[3] );

    NeighborhoodIteratorType ItD( radius, distance->GetOutput(),
      distance->GetOutput()->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<ImageType> It( output,
      output->GetLargestPossibleRegion() );
    for ( It.GoToBegin(), ItD.GoToBegin(); !It.IsAtEnd(); ++It, ++ItD )
      {
      ImageType::PixelType centerPixel = ItD.GetCenterPixel();
      ImageType::PixelType prevPixel = ItD.GetNext( direction );
      ImageType::PixelType nextPixel = ItD.GetPrevious( direction );
      if ( ( ( ( centerPixel < nextPixel && centerPixel < prevPixel ) ||
               ( centerPixel == prevPixel && centerPixel < nextPixel &&
                 centerPixel < ItD.GetPrevious( direction, 2 ) ) ) &&
             centerPixel <= 1 ) ) 
        {
        It.Set( label );
        }
      } 
    } 

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( output );
  writer->SetFileName( argv[2] );                                          
  writer->Update();
    
  return 0;
}
