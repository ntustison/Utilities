#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkThresholdImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include <itkGradientMagnitudeImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include "itkBinaryBallStructuringElement.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include <list>
#include "itkCastImageFilter.h"
#include "vnl/vnl_math.h"
#include "itkHoughTransformRadialVotingImageFilter.h"

#include "time.h"


template <unsigned int ImageDimension>
int HoughTransform( unsigned int argc, char *argv[] )
{
	clock_t beginT, endT;
  beginT = clock();

  const    unsigned int    Dimension = ImageDimension;
  typedef  unsigned int   InputPixelType;
  typedef  float           InternalPixelType;
  typedef  unsigned int   OutputPixelType;

  typedef itk::Image< InputPixelType, Dimension >  InputImageType;
  typedef itk::Image< InternalPixelType, Dimension >    InternalImageType;
  typedef itk::Image< OutputPixelType, Dimension >  OutputImageType;

  typename InputImageType::IndexType localIndex;
  typename InputImageType::SpacingType spacing;

  typedef  itk::ImageFileReader< InputImageType > ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1+1] );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }
  typename InputImageType::Pointer localImage = reader->GetOutput();
  spacing = localImage->GetSpacing();
  std::cout << "Computing Hough Map" << std::endl;

  typedef itk::HoughTransformRadialVotingImageFilter< InputImageType,
               InternalImageType > HoughTransformFilterType;
  typename HoughTransformFilterType::Pointer houghFilter = HoughTransformFilterType::New();
  houghFilter->SetInput( reader->GetOutput() );
  houghFilter->SetNumberOfSpheres( atoi(argv[1+4]) );
  houghFilter->SetMinimumRadius(   atof(argv[1+5]) );
  houghFilter->SetMaximumRadius(   atof(argv[1+6]) );

  if( argc > 8 )
    {
    houghFilter->SetSigmaGradient( atof(argv[1+7]) );
    }
  if( argc > 9 )
    {
    houghFilter->SetVariance( atof(argv[1+8]) );
    }
  if( argc > 10 )
    {
    houghFilter->SetSphereRadiusRatio( atof(argv[1+9]) );
    }
  if( argc > 11 )
    {
  houghFilter->SetVotingRadiusRatio( atof(argv[1+10]) );
    }
  if( argc > 12 )
    {
    houghFilter->SetThreshold( atof(argv[1+11]) );
    }
  if( argc > 13 )
    {
    houghFilter->SetOutputThreshold( atof(argv[1+12]) );
    }
  if( argc > 14 )
    {
    houghFilter->SetGradientThreshold( atof(argv[1+13]) );
    }
  if( argc > 15 )
    {
  houghFilter->SetNbOfThreads( atoi(argv[1+14]) );
    }
  if( argc > 16 )
    {
  houghFilter->SetSamplingRatio( atof(argv[1+15]) );
    }

  std::cout << "Try updating." << std::endl;
  houghFilter->Update();
  std::cout << "Done updating." << std::endl;


  typename InternalImageType::Pointer localAccumulator = houghFilter->GetOutput();

  typename HoughTransformFilterType::SpheresListType circles;
  circles = houghFilter->GetSpheres( );

  endT = clock();
  std::cout << ( static_cast< double >( endT - beginT )/static_cast< double >( CLOCKS_PER_SEC ) ) << std::endl;

  std::cout << "Found " << circles.size() << " circle(s)." << std::endl;


  // Computing the circles output
  typename OutputImageType::Pointer  localOutputImage = OutputImageType::New();

  typename OutputImageType::RegionType region;
  region.SetSize( localImage->GetLargestPossibleRegion().GetSize() );
  region.SetIndex( localImage->GetLargestPossibleRegion().GetIndex() );
  localOutputImage->SetRegions( region );
  localOutputImage->SetOrigin(localImage->GetOrigin());
  localOutputImage->SetSpacing(localImage->GetSpacing());
  localOutputImage->Allocate();
  localOutputImage->FillBuffer(0);

  typedef typename HoughTransformFilterType::SpheresListType SpheresListType;
  typename SpheresListType::const_iterator itSpheres = circles.begin();

  unsigned int count = 1;
  while( itSpheres != circles.end() )
    {
    std::cout << "Center: ";
    std::cout << (*itSpheres)->GetObjectToParentTransform()->GetOffset()
              << std::endl;
    std::cout << "Radius: " << (*itSpheres)->GetRadius()[0] << std::endl;

    itk::ImageRegionIteratorWithIndex<OutputImageType> It( localImage,
      localImage->GetLargestPossibleRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      typename OutputImageType::IndexType index = It.GetIndex();
      float sum = 0.0;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        sum += vnl_math_sqr( index[d] - (*itSpheres)->GetObjectToParentTransform()->GetOffset()[d] );
        }
      if( sum <= vnl_math_sqr( (*itSpheres)->GetRadius()[0] ) )
        {
        localOutputImage->SetPixel( index, count );
        }
      }
    itSpheres++;
    count++;
    }

//  int radius = 2;
//  typedef itk::BinaryBallStructuringElement< OutputPixelType, Dimension >
//    SEType;
//  SEType sE;
//  sE.SetRadius ( radius );
//  sE.CreateStructuringElement();
//
//  typedef itk::GrayscaleDilateImageFilter< OutputImageType, OutputImageType, SEType >
//    DilateFilterType;
//  typename DilateFilterType::Pointer grayscaleDilate = DilateFilterType::New();
//  grayscaleDilate->SetKernel ( sE );
//  grayscaleDilate->SetInput ( localOutputImage );
//  grayscaleDilate->Update();

  typedef itk::ImageFileWriter< OutputImageType > CirclesWriterType;
  typename CirclesWriterType::Pointer cwriter = CirclesWriterType::New();
  cwriter->SetInput( localOutputImage );
  cwriter->SetFileName( argv[1+2] );
  cwriter->Update();

  try
    {
    cwriter->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

  typedef  itk::ImageFileWriter< InternalImageType  > InputWriterType;
  typename InputWriterType::Pointer writer2 = InputWriterType::New();
  writer2->SetFileName( argv[1+3] );
  writer2->SetInput( localAccumulator );

  try
    {
    writer2->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

  return 0;
}

int main( int argc, char *argv[] )
{
  if( argc < 8 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0] << std::endl;
    std::cerr << " imageDimension" << std::endl;
    std::cerr << " inputImage " << std::endl;
    std::cerr << " outputImage" << std::endl;
    std::cerr << " accumulatorImage" << std::endl;
    std::cerr << " numberOfSpheres " << std::endl;
    std::cerr << " radius Min " << std::endl;
    std::cerr << " radius Max " << std::endl;
    std::cerr << " SigmaGradient (default = 1) " << std::endl;
    std::cerr << " variance of the accumulator blurring (default = 1) " << std::endl;
    std::cerr << " radius ratio of the disk to remove from the accumulator (default = 1) "<< std::endl;
    std::cerr << " voting radius ratio (default = 0.5) "<< std::endl;
    std::cerr << " input threshold "<< std::endl;
    std::cerr << " output threshold "<< std::endl;
    std::cerr << " gradient threshold "<< std::endl;
    std::cerr << " number of threads "<< std::endl;
    std::cerr << " sampling ratio "<< std::endl;
    return 1;
    }


  switch( atoi( argv[1] ) )
   {
   case 2:
     HoughTransform<2>( argc, argv );
     break;
    case 3:
      HoughTransform<3>( argc, argv );
      break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

