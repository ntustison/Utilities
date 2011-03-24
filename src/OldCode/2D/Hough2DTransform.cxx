#include "itkGradientMagnitudeImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkHoughTransform2DCirclesImageFilter.h"

#include "global.h"

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImage outputImage [minimumRadius] [maximumRadius] "
              << "[threshold] [sigmaGradient] [numberOfCircles] [sweepAngle] [variance] "
              << "[discRadiusRatio]" << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType> GradientFilterType;
  GradientFilterType::Pointer gradFilter =  GradientFilterType::New();
  gradFilter->SetInput( reader->GetOutput() );
  gradFilter->Update();

  typedef itk::HoughTransform2DCirclesImageFilter<
    ImageType::PixelType, ImageType::PixelType> HoughFilterType;
  HoughFilterType::Pointer hough = HoughFilterType::New();
  hough->SetInput( gradFilter->GetOutput() );     
  if ( argc > 3 )
    {
    hough->SetMinimumRadius( atof( argv[3] ) );  
    } 
  if ( argc > 4 )
    {
    hough->SetMaximumRadius( atof( argv[4] ) );  
    } 
  if ( argc > 5 )
    {
    hough->SetThreshold( atof( argv[5] ) );  
    } 
  if ( argc > 6 )
    {
    hough->SetSigmaGradient( atof( argv[6] ) );  
    } 
  if ( argc > 7 )
    {
    hough->SetNumberOfCircles( atoi( argv[7] ) );  
    } 
  if ( argc > 8 )
    {
    hough->SetSweepAngle( atof( argv[8] ) );  
    } 
  if ( argc > 9 )
    {
    hough->SetVariance( atof( argv[9] ) );  
    } 
  if ( argc > 10 )
    {
    hough->SetDiscRadiusRatio( atof( argv[10] ) );  
    } 
  hough->Update();

  HoughFilterType::CirclesListType circles = hough->GetCircles( hough->GetNumberOfCircles() );
  HoughFilterType::CirclesListType::iterator iter;
  for ( iter = circles.begin(); iter != circles.end(); ++iter )
    {  
    std::cout << ": " << (*iter)->GetObjectToParentTransform()->GetOffset()
      << ", " << (*iter)->GetRadius() << std::endl;
    } 

  typedef itk::ImageFileWriter<ImageType>  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( hough->GetRadiusImage() );
  writer->SetFileName( argv[2] );
  writer->Update();

  return EXIT_SUCCESS;
}

