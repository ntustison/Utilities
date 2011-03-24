#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkMaskImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkScalarImageToHistogramGenerator.h"

#include <fstream.h>

int main( int argc, char *argv[] )
{
  if ( argc != 4 )
    {
    std::cout << "Usage: " << argv[0] << "inputImage boundingBoxFile number_of_bins" << std::endl;
    exit( 1 );
    }
  typedef int PixelType;
  const unsigned int ImageDimension = 3;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  ImageType::Pointer mask = ImageType::New();
  mask->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  mask->Allocate();
  mask->FillBuffer( itk::NumericTraits<PixelType>::Zero );
  
  int indices[2*ImageDimension];
  ifstream str( argv[2], std::ios::in );
  for ( unsigned int i = 0; i < 2*ImageDimension; i++ )
    {
    str >> indices[i];
    } 
    
  itk::ImageRegionIteratorWithIndex<ImageType> It( mask,
    mask->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    bool inside = true;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      { 
      if ( It.GetIndex()[i] < indices[i] || 
           It.GetIndex()[i] > indices[i+ImageDimension] )
        {
        inside = false;
        }          
      }
    if ( inside )
      {
      It.Set( itk::NumericTraits<PixelType>::One );
      }
    }

  typedef itk::MaskImageFilter<ImageType, ImageType, ImageType> MaskerType;
  MaskerType::Pointer masker = MaskerType::New();
  masker->SetInput1( reader->GetOutput() );
  masker->SetInput2( mask );
  masker->Update();
  
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( "mask.hdr" );
  writer->SetInput( masker->GetOutput() );
  writer->Update();
  
  
  typedef itk::LabelStatisticsImageFilter<ImageType, ImageType> HistogramGeneratorType;
  HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
  stats->SetInput( masker->GetOutput() );
  stats->SetLabelInput( mask );
  stats->UseHistogramsOn();
  stats->SetHistogramParameters( atoi( argv[3] ), 0.5, 1500.5 );
  stats->Update();
  
  typedef HistogramGeneratorType::HistogramType  HistogramType;
  const HistogramType *histogram = stats->GetHistogram( 1 );

  ofstream str0( "histogram.dat" );

  unsigned int i;
  for( i = 0; i < histogram->Size(); i++ )
    {
    str0 << histogram->GetMeasurement( i, 0 ) << " "
        << histogram->GetFrequency( i, 0 ) << " "
        << histogram->GetTotalFrequency() << std::endl;
    }
    
  ofstream str1( "stats.dat" );  
  str1 << "Mean  = " << stats->GetMean( 1 ) << std::endl;
  str1 << "Sigma = " << stats->GetSigma( 1 ) << std::endl;

      
  return 0;
}
