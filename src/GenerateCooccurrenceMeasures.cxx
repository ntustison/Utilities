#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkScalarImageToTextureFeaturesFilter.h"
#include "itkDenseFrequencyContainer2.h"

template <unsigned int ImageDimension>
int GenerateCooccurrenceMeasures( int argc, char *argv[] )
{

  typedef float PixelType;
  typedef float RealType;


  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[2] );
  imageReader->Update();

  // we have to rescale the input image since the cooccurrence filter
  // adds 1 to the upper bound of the joint histogram and small ranges
  // will be greatly affected by this.

  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescalerType;
  typename RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetInput( imageReader->GetOutput() );
  rescaler->SetOutputMinimum( 0 );
  rescaler->SetOutputMaximum( 10000 );
  rescaler->Update();

  typedef itk::Statistics::DenseFrequencyContainer2 HistogramFrequencyContainerType;

  typedef itk::Statistics::ScalarImageToTextureFeaturesFilter
    <RealImageType, HistogramFrequencyContainerType> TextureFilterType;
  typename TextureFilterType::Pointer textureFilter = TextureFilterType::New();
  textureFilter->SetInput( rescaler->GetOutput() );

  typename ImageType::Pointer mask = ITK_NULLPTR;
  PixelType label = itk::NumericTraits<PixelType>::One;
  if ( argc > 4 )
    {
    typename ReaderType::Pointer labelImageReader = ReaderType::New();
    labelImageReader->SetFileName( argv[4] );
    labelImageReader->Update();
    mask = labelImageReader->GetOutput();
    textureFilter->SetMaskImage( mask );

    if ( argc > 5 )
      {
      label = static_cast<PixelType>( atoi( argv[5] ) );
      }
    textureFilter->SetInsidePixelValue( label );
    }


  unsigned int numberOfBins = 256;
  if ( argc > 3 )
    {
    numberOfBins = static_cast<PixelType>( atoi( argv[3] ) );
    }
  textureFilter->SetNumberOfBinsPerAxis( numberOfBins );


  itk::ImageRegionIteratorWithIndex<ImageType> ItI( rescaler->GetOutput(),
    rescaler->GetOutput()->GetLargestPossibleRegion() );

  PixelType maxValue = itk::NumericTraits<PixelType>::NonpositiveMin();
  PixelType minValue = itk::NumericTraits<PixelType>::max();

  for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
    {
    if ( !mask || ( mask->GetPixel( ItI.GetIndex() ) == label ) )
      {
      if ( ItI.Get() < minValue )
        {
        minValue = ItI.Get();
        }
      else if ( ItI.Get() > maxValue )
        {
        maxValue = ItI.Get();
        }
      }
    }

  textureFilter->SetPixelValueMinMax( minValue, maxValue );
  textureFilter->SetNumberOfBinsPerAxis( numberOfBins );
  textureFilter->FastCalculationsOff();

  /**
   * Second order measurements
   * These include:
   *   1. energy (f1) *cth, *amfm
   *   2. entropy (f2) *cth, *amfm
   *   3. correlation (f3) *amfm
   *   4. inverse difference moment (f4) *cth, *amfm
   *   5. inertia (f5) *cth, *amfm
   *   6. cluster shade (f6) *cth
   *   7. cluster prominence (f7) *cth
   *   8. haralick's correlation (f8)
   */

  textureFilter->Update();

  typename TextureFilterType::FeatureValueVectorPointer means =
    textureFilter->GetFeatureMeans();
  const typename TextureFilterType::FeatureNameVector* names =
    textureFilter->GetRequestedFeatures();

  typename TextureFilterType::FeatureValueVector::ConstIterator mIt =
    means->Begin();
  typename TextureFilterType::FeatureNameVector::ConstIterator nIt =
    names->Begin();


  std::cout << "Energy,Entropy,InverseDifferenceMoment,Inertia,ClusterShade,ClusterProminence" << std::endl;
  while( mIt != means->End() )
    {
//    std::cout << nIt.Value() << ": " << mIt.Value() << std::endl;
    std::cout << mIt.Value() << " ";
    ++mIt;
    ++nIt;
    }
  std::cout << std::endl;

//  RealType entropy = textureFilter->GetFeatureMeansOutput()[1];
//  RealType inverseDifferenceMoment = textureFilter->GetFeatureMeansOutput()[2];
//  RealType inertia = textureFilter->GetFeatureMeansOutput()[3];
//  RealType clusterShade = textureFilter->GetFeatureMeansOutput()[4];
//  RealType clusterProminence = textureFilter->GetFeatureMeansOutput()[5];
//
//  std::cout << energy << " "
//            << entropy << " "
//            << inverseDifferenceMoment << " "
//            << inertia << " "
//            << clusterShade << " "
//            << clusterProminence << std::endl;

/*
  std::cout << "energy             : " << energy << std::endl;
  std::cout << "entropy            : " << entropy << std::endl;
  std::cout << "inverse diff moment: " << inverseDifferenceMoment << std::endl;
  std::cout << "inertia            : " << inertia << std::endl;
  std::cout << "cluster Shade      : " << clusterShade << std::endl;
  std::cout << "cluster prominence : " << clusterProminence << std::endl;
*/

  return 0;
}


int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage "
     << "[numberOfBinsPerAxis=256] [maskImage] [maskLabel=1]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GenerateCooccurrenceMeasures<2>( argc, argv );
     break;
   case 3:
     GenerateCooccurrenceMeasures<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

