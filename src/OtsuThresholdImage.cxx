#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkLabelStatisticsImageFilter.h"
#include "itkOtsuMultipleThresholdsCalculator.h"

template <unsigned int ImageDimension>
int OtsuThresholdImage( int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<int, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  unsigned int numberOfThresholds = 2;
  if ( argc > 4 )
    {
    numberOfThresholds = static_cast<unsigned int>( atoi( argv[4] ) );
    }
  unsigned int numberOfBins = 100;
  if ( argc > 5 )
    {
    numberOfBins = static_cast<unsigned int>( atoi( argv[5] ) );
    }

  int maskLabel = 1;
  if ( argc > 7 )
    {
    maskLabel = atoi( argv[7] );
    }

  typedef itk::Image<int, ImageDimension> MaskImageType;
  typename MaskImageType::Pointer maskImage = MaskImageType::New();
  if ( argc > 6 )
    {
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    typename MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName( argv[6] );
    maskReader->Update();
    maskImage = maskReader->GetOutput();
    }
  else
    {
    maskImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
    maskImage->SetOrigin( reader->GetOutput()->GetOrigin() );
    maskImage->SetSpacing( reader->GetOutput()->GetSpacing() );
    maskImage->SetDirection( reader->GetOutput()->GetDirection() );
    maskImage->Allocate();
    maskImage->FillBuffer( maskLabel );
    }

  itk::ImageRegionIterator<ImageType> ItI( reader->GetOutput(),
    reader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<MaskImageType> ItM( maskImage,
    maskImage->GetLargestPossibleRegion() );
  PixelType maxValue = itk::NumericTraits<PixelType>::min();
  PixelType minValue = itk::NumericTraits<PixelType>::max();
  for ( ItM.GoToBegin(), ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItM, ++ItI )
    {
    if ( ItM.Get() == maskLabel )
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

  typedef itk::LabelStatisticsImageFilter<ImageType, MaskImageType> StatsType;
  typename StatsType::Pointer stats = StatsType::New();
  stats->SetInput( reader->GetOutput() );
  stats->SetLabelInput( maskImage );
  stats->UseHistogramsOn();
  stats->SetHistogramParameters( numberOfBins, minValue, maxValue );
  stats->Update();

  typedef itk::OtsuMultipleThresholdsCalculator<typename StatsType::HistogramType>
    OtsuType;
  typename OtsuType::Pointer otsu = OtsuType::New();
  otsu->SetInputHistogram( stats->GetHistogram( maskLabel ) );
  otsu->SetNumberOfThresholds( numberOfThresholds );
  otsu->Update();

  typename OtsuType::OutputType thresholds = otsu->GetOutput();

  typename MaskImageType::Pointer output = MaskImageType::New();
  output->SetRegions( maskImage->GetLargestPossibleRegion() );
  output->SetOrigin( maskImage->GetOrigin() );
  output->SetSpacing( maskImage->GetSpacing() );
  output->SetDirection( maskImage->GetDirection() );
  output->Allocate();
  output->FillBuffer( 0 );

  itk::ImageRegionIterator<MaskImageType> ItO( output,
    output->GetLargestPossibleRegion() );
  for ( unsigned int i = 0; i < thresholds.size(); i++ )
    {

    ItI.GoToBegin();
    ItM.GoToBegin();
    ItO.GoToBegin();
    while ( !ItM.IsAtEnd() )
      {
      if ( ItM.Get() == maskLabel )
        {
        if ( ItO.Get() == 0 && ItI.Get() < thresholds[i] )
          {
          ItO.Set( i+1 );
          }
        }
      ++ItI;
      ++ItM;
      ++ItO;
      }
    }

  ItI.GoToBegin();
  ItM.GoToBegin();
  ItO.GoToBegin();
  while ( !ItM.IsAtEnd() )
    {
    if ( ItM.Get() == maskLabel && ItO.Get() == 0 )
      {
      ItO.Set( thresholds.size()+1 );
      }
    ++ItI;
    ++ItM;
    ++ItO;
    }

  typedef itk::ImageFileWriter<MaskImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( output );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage "
      << "outputImage [numberOfThresholds] [numberOfBins] [maskImage] "
      << "[maskLabel]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     OtsuThresholdImage<2>( argc, argv );
     break;
   case 3:
     OtsuThresholdImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

