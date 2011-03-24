#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include "itkLabelStatisticsImageFilter.h"

template <unsigned int ImageDimension>
int GeneratePercentileAttenuationMask( int argc, char *argv[] )
{
  typedef int PixelType;
  typedef float RealType;

  unsigned int numberOfBins = 200;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[2] );
  imageReader->Update();

  PixelType outputLabel = itk::NumericTraits<PixelType>::One;
  if ( argc > 6 )
    {
    outputLabel = static_cast<PixelType>( atoi( argv[6] ) );
    }

  typename ImageType::Pointer maskImage = ImageType::New();
  if ( argc > 7 )
    {
    typename ReaderType::Pointer labelImageReader = ReaderType::New();
    labelImageReader->SetFileName( argv[7] );
    labelImageReader->Update();
    maskImage = labelImageReader->GetOutput();
    }
  else
    {
    maskImage->SetOrigin( imageReader->GetOutput()->GetOrigin() );
    maskImage->SetSpacing( imageReader->GetOutput()->GetSpacing() );
    maskImage->SetRegions( imageReader->GetOutput()->GetLargestPossibleRegion() );
    maskImage->Allocate();
    maskImage->FillBuffer( itk::NumericTraits<PixelType>::One );
    }
  PixelType label = itk::NumericTraits<PixelType>::One;
  if ( argc > 8 )
    {
    label = static_cast<PixelType>( atoi( argv[8] ) );
    }
  itk::ImageRegionIterator<ImageType> ItI( imageReader->GetOutput(),
    imageReader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> ItM( maskImage,
    maskImage->GetLargestPossibleRegion() );

  PixelType maxValue = itk::NumericTraits<PixelType>::min();
  PixelType minValue = itk::NumericTraits<PixelType>::max();

  for ( ItM.GoToBegin(), ItI.GoToBegin(); !ItM.IsAtEnd(); ++ItM, ++ItI )
    {
    if ( ItM.Get() == label )
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

  typedef itk::LabelStatisticsImageFilter<ImageType, ImageType> HistogramGeneratorType;
  typename HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
  stats->SetInput( imageReader->GetOutput() );
  stats->SetLabelInput( maskImage );
  stats->UseHistogramsOn();
  stats->SetHistogramParameters( numberOfBins, minValue, maxValue );
  stats->Update();

  typedef typename HistogramGeneratorType::HistogramType  HistogramType;
  const HistogramType *histogram = stats->GetHistogram( label );

  double minPercentileValue = histogram->Quantile( 0, atof( argv[4] ) );
  double maxPercentileValue = histogram->Quantile( 0, atof( argv[5] ) );

  std::cout << "Min percentile: " << minPercentileValue << std::endl;
  std::cout << "Max percentile: " << maxPercentileValue << std::endl;

  for ( ItM.GoToBegin(), ItI.GoToBegin(); !ItM.IsAtEnd(); ++ItM, ++ItI )
    {
    RealType value = ItI.Get();
    if ( ItM.Get() == label &&
         value >= minPercentileValue && value <= maxPercentileValue )
      {
      ItM.Set( outputLabel );
      }
    else
      {
      ItM.Set( 0 );
      }
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( maskImage );
  writer->SetFileName( argv[3] );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage outputImage minPercentile maxPercentile ";
    std::cout << " [outputLabel] [labelImage] [label]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GeneratePercentileAttenuationMask<2>( argc, argv );
     break;
   case 3:
     GeneratePercentileAttenuationMask<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


