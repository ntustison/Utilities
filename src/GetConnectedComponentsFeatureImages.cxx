#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkLabelPerimeterEstimationCalculator.h"

#include <vector>

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

  // Output images:
  // [0] = volume (in physical coordinates)
  // [1] = volume / surface area
  // [2] = eccentricity
  // [3] = elongation
  // [4] = orientation

  std::vector<typename RealImageType::Pointer> outputImages;
  for( unsigned int n = 0; n < 5; n++ )
    {
    typename RealImageType::Pointer output = RealImageType::New();
    output->CopyInformation( reader->GetOutput() );
    output->SetRegions( reader->GetOutput()->GetRequestedRegion() );
    output->Allocate();
    output->FillBuffer( 0.0 );

    outputImages.push_back( output );
    }

  typename ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();
  float prefactor = 1.0;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    prefactor *= spacing[d];
    }

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

    typedef itk::LabelGeometryImageFilter<ImageType, RealImageType> GeometryFilterType;
    typename GeometryFilterType::Pointer geometry = GeometryFilterType::New();
    geometry->SetInput( relabeler2->GetOutput() );
    geometry->CalculatePixelIndicesOff();
    geometry->CalculateOrientedBoundingBoxOff();
    geometry->CalculateOrientedLabelRegionsOff();
    geometry->Update();

    typedef itk::LabelPerimeterEstimationCalculator<ImageType> AreaFilterType;
    typename AreaFilterType::Pointer area = AreaFilterType::New();
    area->SetImage( relabeler2->GetOutput() );
    area->Compute();

    itk::ImageRegionIteratorWithIndex<ImageType> It( relabeler->GetOutput(),
      relabeler->GetOutput()->GetRequestedRegion() );
    itk::ImageRegionIterator<ImageType> It2( relabeler2->GetOutput(),
      relabeler2->GetOutput()->GetRequestedRegion() );

    for( It.GoToBegin(), It2.GoToBegin(); !It.IsAtEnd(); ++It, ++It2 )
      {
      int label = It2.Get();
      if( label != 0 )
        {
        typename ImageType::IndexType index = It.GetIndex();

        // Output images:
        // [0] = volume (in physical coordinates)
        // [1] = volume / surface area
        // [2] = eccentricity
        // [3] = elongation
        // [4] = orientation

        float volume = prefactor * static_cast<float>( geometry->GetVolume( label ) );

        outputImages[0]->SetPixel( index, volume );
        outputImages[1]->SetPixel( index, area->GetPerimeter( label ) / volume );
        outputImages[2]->SetPixel( index, geometry->GetEccentricity( label ) );
        outputImages[3]->SetPixel( index, geometry->GetElongation( label ) );
        outputImages[4]->SetPixel( index, geometry->GetOrientation( label ) );
        }
      }
    }

  typedef itk::ImageFileWriter<RealImageType> WriterType;

  {
  std::string filename = std::string( argv[3] ) + std::string( "PHYSICAL_VOLUME.nii.gz" );
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput( outputImages[0] );
  writer->Update();
  }

  {
  std::string filename = std::string( argv[3] ) + std::string( "VOLUME_TO_SURFACE_AREA_RATIO.nii.gz" );
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput( outputImages[1] );
  writer->Update();
  }

  {
  std::string filename = std::string( argv[3] ) + std::string( "ECCENTRICITY.nii.gz" );
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput( outputImages[2] );
  writer->Update();
  }

  {
  std::string filename = std::string( argv[3] ) + std::string( "ELONGATION.nii.gz" );
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput( outputImages[3] );
  writer->Update();
  }

  {
  std::string filename = std::string( argv[3] ) + std::string( "ORIENTATION.nii.gz" );
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput( outputImages[4] );
  writer->Update();
  }


  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension "
              << "inputSegmentationImage outputImagePrefix" << std::endl;
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
