#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkLabelStatisticsImageFilter.h"

template <unsigned int ImageDimension>
int TruncateImageIntensity( unsigned int argc, char *argv[] )
{
  typedef int PixelType;
  typedef float RealType;

  unsigned int numberOfBins = 200;
  if ( argc > 8 )
    {
    numberOfBins = atoi( argv[8] );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  typename ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[2] );
  imageReader->Update();

  typename ImageType::Pointer mask = ImageType::New();
  if ( argc > 4 )
    {
    try
      {
      typedef itk::ImageFileReader<ImageType> ReaderType;
      typename ReaderType::Pointer labelImageReader = ReaderType::New();
      labelImageReader->SetFileName( argv[4] );
      labelImageReader->Update();
      mask = labelImageReader->GetOutput();
      }
    catch(...)
      {
      mask->SetOrigin( imageReader->GetOutput()->GetOrigin() );
      mask->SetSpacing( imageReader->GetOutput()->GetSpacing() );
      mask->SetRegions( imageReader->GetOutput()->GetLargestPossibleRegion() );
      mask->SetDirection( imageReader->GetOutput()->GetDirection() );
      mask->Allocate();
      mask->FillBuffer( itk::NumericTraits<PixelType>::One );
      };
    }
  if( !mask )
    {
    mask->SetOrigin( imageReader->GetOutput()->GetOrigin() );
    mask->SetSpacing( imageReader->GetOutput()->GetSpacing() );
    mask->SetRegions( imageReader->GetOutput()->GetLargestPossibleRegion() );
    mask->SetDirection( imageReader->GetOutput()->GetDirection() );
    mask->Allocate();
    mask->FillBuffer( itk::NumericTraits<PixelType>::One );
    }
  PixelType label = itk::NumericTraits<PixelType>::One;
  if ( argc > 5 )
    {
    label = static_cast<PixelType>( atoi( argv[5] ) );
    }

  itk::ImageRegionIterator<RealImageType> ItI( imageReader->GetOutput(),
    imageReader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> ItM( mask,
    mask->GetLargestPossibleRegion() );

  RealType maxValue = itk::NumericTraits<RealType>::NonpositiveMin();
  RealType minValue = itk::NumericTraits<RealType>::max();

  for ( ItM.GoToBegin(), ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItM, ++ItI )
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
      ItM.Set( itk::NumericTraits<PixelType>::One );
      }
    else
      {
      ItM.Set( itk::NumericTraits<PixelType>::Zero );
      }
    if( std::isnan( ItI.Get() ) || std::isinf( ItI.Get() ) )
      {
      ItM.Set( itk::NumericTraits<PixelType>::Zero );
      }
    }

  typedef itk::LabelStatisticsImageFilter<RealImageType, ImageType> HistogramGeneratorType;
  typename HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
  stats->SetInput( imageReader->GetOutput() );
  stats->SetLabelInput( mask );
  stats->SetUseHistograms( true );
  stats->SetHistogramParameters( numberOfBins, minValue, maxValue );
  stats->Update();
  typedef typename HistogramGeneratorType::HistogramType  HistogramType;
  const HistogramType *histogram = stats->GetHistogram( 1 );

  double lowerValue = 0.025;
  if( argc > 6 )
    {
    lowerValue = atof( argv[6] );
    }
		double lowerQuantile = histogram->Quantile( 0, lowerValue );
  double upperValue = 0.975;
  if( argc > 7 )
    {
    upperValue = atof( argv[7] );
    }
		double	upperQuantile = histogram->Quantile( 0, upperValue );

  std::cout << "Lower quantile: " << lowerQuantile << std::endl;
  std::cout << "Upper quantile: " << upperQuantile << std::endl;


  for ( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
    {
    if ( ItI.Get() < lowerQuantile )
      {
      ItI.Set( lowerQuantile );
      }
    else if( ItI.Get() > upperQuantile )
      {
      ItI.Set( upperQuantile );
      }
    }

  typedef itk::ImageFileWriter<RealImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( imageReader->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension inputImage outputImage "
      << "[maskImage] [maskLabel=1] [lowerQuantile=0.025] [upperQuantile=0.975] "
      << "[numberOfBins=200]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     TruncateImageIntensity<2>( argc, argv );
     break;
   case 3:
     TruncateImageIntensity<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

