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
  if ( argc != 5 )
    {
    std::cout << "Usage: " << argv[0] << "image_filename segmented_image_filename number_of_bins label" << std::endl;
    exit( 1 );
    }
  typedef int PixelType;
  const unsigned int ImageDimension = 3;


  typedef itk::Image<PixelType, ImageDimension> ImageType;
  
  typedef itk::MaskImageFilter<ImageType, ImageType, ImageType> MaskerType;
  MaskerType::Pointer masker = MaskerType::New();
  
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[1] );
  reader1->Update();

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  reader2->Update();
  
  if ( argc == 5 )
    {
    PixelType label = static_cast<PixelType>( atoi( argv[4] ) );
    
    itk::ImageRegionIteratorWithIndex<ImageType> It( reader2->GetOutput(),
      reader2->GetOutput()->GetLargestPossibleRegion() );
    
    for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if ( It.Get() != label )
        {
        It.Set( itk::NumericTraits<PixelType>::Zero );
        }
      else
        {
        It.Set( itk::NumericTraits<PixelType>::One );
        }
      }
    }


  masker->SetInput1( reader1->GetOutput() );
  masker->SetInput2( reader2->GetOutput() );
  masker->Update();
  
  typedef itk::LabelStatisticsImageFilter<ImageType, ImageType> HistogramGeneratorType;
  HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
  stats->SetInput( masker->GetOutput() );
  stats->SetLabelInput( reader2->GetOutput() );
  stats->UseHistogramsOn();
  stats->SetHistogramParameters( atoi( argv[3] ), 0.5, 1500.5 );
  stats->Update();
  
  typedef HistogramGeneratorType::HistogramType  HistogramType;
  const HistogramType *histogram = stats->GetHistogram( 1 );

  std::string histogram_file = std::string( "histogram_" )
    + std::string( argv[4] ) + std::string( ".dat" );
  ofstream str( histogram_file.c_str() );

  unsigned int i;
  for( i = 0; i < histogram->Size(); i++ )
    {
    str << histogram->GetMeasurement( i, 0 ) << " "
        << histogram->GetFrequency( i, 0 ) << " "
        << histogram->GetTotalFrequency() << std::endl;
    }
    
  std::string stats_file = std::string( "stats_" )
    + std::string( argv[4] ) + std::string( ".dat" );
  ofstream str1( stats_file.c_str() );

  str1 << "Mean  = " << stats->GetMean( 1 ) << std::endl;
  str1 << "Sigma = " << stats->GetSigma( 1 ) << std::endl;

      
  return 0;
}
