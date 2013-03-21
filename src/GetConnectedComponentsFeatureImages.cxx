#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkLabelPerimeterEstimationCalculator.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

#include <string>
#include <vector>

template<class TValue>
TValue Convert( std::string optionString )
{
  TValue value;
  std::istringstream iss( optionString );
  iss >> value;
  return value;
}

template<class TValue>
std::vector<TValue> ConvertVector( std::string optionString )
{
  std::vector<TValue> values;
  std::string::size_type crosspos = optionString.find( 'x', 0 );

  if ( crosspos == std::string::npos )
    {
    values.push_back( Convert<TValue>( optionString ) );
    }
  else
    {
    std::string element = optionString.substr( 0, crosspos ) ;
    TValue value;
    std::istringstream iss( element );
    iss >> value;
    values.push_back( value );
    while ( crosspos != std::string::npos )
      {
      std::string::size_type crossposfrom = crosspos;
      crosspos = optionString.find( 'x', crossposfrom + 1 );
      if ( crosspos == std::string::npos )
        {
        element = optionString.substr( crossposfrom + 1, optionString.length() );
        }
      else
        {
        element = optionString.substr( crossposfrom + 1, crosspos ) ;
        }
      std::istringstream iss2( element );
      iss2 >> value;
      values.push_back( value );
      }
    }
  return values;
}

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

  std::vector<typename RealImageType::Pointer> outputImages;
  for( unsigned int n = 0; n < 4; n++ )
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

//   std::vector<unsigned int> tumorLabels = ConvertVector<unsigned int>( std::string( argv[4] ) );

  for( unsigned int i = 1; i <= relabeler->GetNumberOfObjects(); i++ )
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

        float volume = prefactor * static_cast<float>( geometry->GetVolume( label ) );

        outputImages[0]->SetPixel( index, volume );
        outputImages[1]->SetPixel( index, area->GetPerimeter( label ) / volume );
        outputImages[2]->SetPixel( index, geometry->GetEccentricity( label ) );
        outputImages[3]->SetPixel( index, geometry->GetElongation( label ) );
        }
      }
    }

//   typename ImageType::Pointer relabeledImage = ImageType::New();
//   relabeledImage->CopyInformation( relabeler->GetOutput() );
//   relabeledImage->SetRegions( relabeler->GetOutput()->GetRequestedRegion() );
//   relabeledImage->Allocate();
//   relabeledImage->FillBuffer( 0 );
//
//   itk::ImageRegionIteratorWithIndex<ImageType> ItR( relabeler->GetOutput(),
//     relabeler->GetOutput()->GetRequestedRegion() );
//   for( ItR.GoToBegin(); !ItR.IsAtEnd(); ++ItR )
//     {
//     PixelType label = ItR.Get();
//     if( label != 0 )
//       {
//       relabeledImage->SetPixel( ItR.GetIndex(), label );
// //       for( unsigned int n = 1; n < tumorLabels.size(); n++ )
// //         {
// //         if( label == tumorLabels[n] )
// //           {
// //           relabeledImage->SetPixel( ItR.GetIndex(), tumorLabels[0] );
// //           }
// //         }
//       }
//     }
//
//   typename RelabelerType::Pointer relabeler3 = RelabelerType::New();
//   relabeler3->SetInput( relabeledImage );
//   relabeler3->Update();
//
//   for( unsigned int i = 1; i <= relabeler3->GetNumberOfObjects(); i++ )
//     {
//     typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholderType;
//     typename ThresholderType::Pointer thresholder = ThresholderType::New();
//     thresholder->SetInput( relabeler3->GetOutput() );
//     thresholder->SetLowerThreshold( i );
//     thresholder->SetUpperThreshold( i );
//     thresholder->SetInsideValue( 1 );
//     thresholder->SetOutsideValue( 0 );
//     thresholder->Update();
//
//     typedef itk::SignedMaurerDistanceMapImageFilter<ImageType, RealImageType> FilterType;
//     typename FilterType::Pointer filter = FilterType::New();
//     filter->SetInput( thresholder->GetOutput() );
//     filter->SetSquaredDistance( false );
//     filter->SetUseImageSpacing( true );
//     filter->SetInsideIsPositive( true );
//     filter->Update();
//
//     typedef itk::LabelStatisticsImageFilter<RealImageType, ImageType> HistogramGeneratorType;
//     typename HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
//     stats->SetInput( filter->GetOutput() );
//     stats->SetLabelInput( thresholder->GetOutput() );
//     stats->Update();
//
//     float maxDistance = stats->GetMaximum( 1 );
//
//     typedef itk::MultiplyImageFilter<RealImageType, RealImageType, RealImageType> MultiplierType;
//     typename MultiplierType::Pointer multiplier = MultiplierType::New();
//     multiplier->SetInput( filter->GetOutput() );
//     multiplier->SetConstant( 1.0 /* / maxDistance */ );
//     multiplier->Update();
//
//     itk::ImageRegionIterator<ImageType> ItT( thresholder->GetOutput(),
//       thresholder->GetOutput()->GetRequestedRegion() );
//     itk::ImageRegionIteratorWithIndex<RealImageType> ItM( multiplier->GetOutput(),
//       multiplier->GetOutput()->GetRequestedRegion() );
//     for( ItT.GoToBegin(), ItM.GoToBegin(); !ItT.IsAtEnd(); ++ItT, ++ItM )
//       {
//       int label = ItT.Get();
//       if( label != 0 )
//         {
//         outputImages[2]->SetPixel( ItM.GetIndex(), ItM.Get() );
//         }
//       }
//     }

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
