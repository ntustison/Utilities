#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include "itkLabelStatisticsImageFilter.h"

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << " image outputImage minPercentile maxPercentile ";
    std::cout << " [outputLabel] [labelImage] [label]" << std::endl;
    exit( 1 );
    }

  typedef int PixelType;
  const unsigned int ImageDimension = 3;
  typedef float RealType;

  unsigned int numberOfBins = 100;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[1] );
  imageReader->Update();

  PixelType outputLabel = itk::NumericTraits<PixelType>::One;
  if ( argc > 5 )
    {
    outputLabel = static_cast<PixelType>( atoi( argv[5] ) );
    }
  
  ImageType::Pointer maskImage = ImageType::New();
  if ( argc > 6 )
    {
    ReaderType::Pointer labelImageReader = ReaderType::New();
    labelImageReader->SetFileName( argv[6] );
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
  if ( argc > 7 )
    {
    label = static_cast<PixelType>( atoi( argv[7] ) );
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
  HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
  stats->SetInput( imageReader->GetOutput() );
  stats->SetLabelInput( maskImage );
  stats->UseHistogramsOn();
  stats->SetHistogramParameters( numberOfBins, minValue, maxValue );
  stats->Update();

  typedef HistogramGeneratorType::HistogramType  HistogramType;
  const HistogramType *histogram = stats->GetHistogram( label );

  double minPercentileValue = histogram->Quantile( 0, atof( argv[3] ) );
  double maxPercentileValue = histogram->Quantile( 0, atof( argv[4] ) );

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
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( maskImage );
  writer->SetFileName( argv[2] );
  writer->Update();

  return 0;
}
