#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkAddImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageDuplicator.h"
#include "itkMultiplyImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"

#include <vector>
#include <iomanip>
#include <string>

typedef float RealType;

float GetCorrelation( std::vector<float> X, std::vector<float> Y )
{
  if( X.size() != Y.size() )
    {
    std::cerr << "Vectors are not the same size" << std::endl;
    exit( 1 );
    }

  RealType N = X.size();

  std::vector<RealType>::const_iterator itX;
  std::vector<RealType>::const_iterator itY;

  RealType sumX = 0.0;
  RealType sumY = 0.0;
  RealType sumX2 = 0.0;
  RealType sumY2 = 0.0;
  RealType sumXY = 0.0;

  for( itX = X.begin(), itY = Y.begin(); itX != X.end(); ++itX, ++itY )
    {
    sumX  += (*itX);
    sumY  += (*itY);
    sumXY += (*itX) * (*itY);
    sumX2 += (*itX) * (*itX);
    sumY2 += (*itY) * (*itY);
    }

  float numerator = N * sumXY - sumX * sumY;
  float denominator = vcl_sqrt( N * sumX2 - vnl_math_sqr( sumX ) ) *
    vcl_sqrt( N * sumY2 - vnl_math_sqr( sumY ) );

  float r = numerator / denominator;

  return r;
}

template<class ImageType>
typename ImageType::Pointer
CreateCorrelationMap( std::vector<typename ImageType::Pointer> images )
{
  typename ImageType::Pointer w = ImageType::New();
  w->CopyInformation( images[0] );
  w->SetRegions( images[0]->GetLargestPossibleRegion() );
  w->Allocate();
  w->FillBuffer( 0.0 );

  typedef itk::ConstNeighborhoodIterator<ImageType> NeighborhoodIteratorType;
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );

  NeighborhoodIteratorType It( radius, images[0], images[0]->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( !It.InBounds() )
      {
      continue;
      }

    std::vector<float> X;
    std::vector<float> Y0;
    std::vector<float> Y1;
    std::vector<float> Y2;
    std::vector<float> Y3;
    for( unsigned int n = 0; n < images.size(); n++ )
      {
      X.push_back( images[n]->GetPixel( It.GetIndex() ) );
      Y0.push_back( images[n]->GetPixel( It.GetIndex( 1 ) ) );
      Y1.push_back( images[n]->GetPixel( It.GetIndex( 3 ) ) );
      Y2.push_back( images[n]->GetPixel( It.GetIndex( 5 ) ) );
      Y3.push_back( images[n]->GetPixel( It.GetIndex( 7 ) ) );
      }

    float r0 = GetCorrelation( X, Y0 );
    float r1 = GetCorrelation( X, Y1 );
    float r2 = GetCorrelation( X, Y2 );
    float r3 = GetCorrelation( X, Y3 );

    float sum = ( r0 + r1 + r2 + r3 ) / 4.0;
    if( sum > 0.75 )
      {
      w->SetPixel( It.GetIndex(), sum );
      }
    }

  return w;
}

template<class ImageType>
float
CalculateEnergy( std::vector<typename ImageType::Pointer> images,
  std::vector<typename ImageType::Pointer> signalImages,
  std::vector<typename ImageType::Pointer> syntheticImages,
  typename ImageType::Pointer correlationMap, float alpha, float beta )
{
  std::vector<typename ImageType::Pointer> gradImages;

  for( unsigned int n = 0; n < images.size(); n++ )
    {
    typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType, ImageType> GradientFilterType;
    typename GradientFilterType::Pointer gradFilter = GradientFilterType::New();
    gradFilter->SetSigma( 1.2 );
    gradFilter->SetInput( images[n] );
    gradFilter->SetNormalizeAcrossScale( false );
    gradFilter->Update();

    gradImages.push_back( gradFilter->GetOutput() );
    }

  float summation = 0.0;
  float N = 0.0;

  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
  IteratorType It( images[0], images[0]->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    typename ImageType::IndexType index = It.GetIndex();

    float w = correlationMap->GetPixel( index );

    for( unsigned int n = 0; n < images.size(); n++ )
      {
      float I = images[n]->GetPixel( index );
      float M = syntheticImages[n]->GetPixel( index );
      float mag = gradImages[n]->GetPixel( index );
      float absS = vnl_math_abs( signalImages[n]->GetPixel( index ) );

      summation += ( vnl_math_sqr( I - M ) + alpha * w * mag + beta * vnl_math_sqr( absS -  M ) );
      N+=1.0;
      }
    }
  summation /= N;

  typename ImageType::SpacingType spacing = images[0]->GetSpacing();

  for( unsigned int d = 0; d < spacing.Size(); d++ )
    {
    summation *= spacing[d];
    }

  return summation;
}

int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cout
      << argv[0] << " syntheticImagePrefix deltaK alpha beta inputImage1 signalImage1 "
      << "inputImage2 signalImage2 ... inputImageN signalImageN "
      << std::endl;
    exit( 1 );
    }

  typedef float PixelType;
  typedef itk::Image<PixelType, 2> ImageType;

  std::vector<ImageType::Pointer> inputImages;
  std::vector<ImageType::Pointer> signalImages;
  std::vector<ImageType::Pointer> syntheticImages;
  std::vector<ImageType::Pointer> newSyntheticImages;

  float deltaK = atof( argv[2] );
  float alpha = atof( argv[3] );
  float beta = atof( argv[4] );

  for( unsigned int n = 5; n < argc; n+=2 )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[n] );
    reader->Update();

    ReaderType::Pointer reader2 = ReaderType::New();
    reader2->SetFileName( argv[n+1] );
    reader2->Update();

    inputImages.push_back( reader->GetOutput() );
    signalImages.push_back( reader2->GetOutput() );

    typedef itk::ImageDuplicator<ImageType> DuplicatorType;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage( reader2->GetOutput() );
    duplicator->Update();

    syntheticImages.push_back( duplicator->GetOutput() );
    }

  ImageType::Pointer correlationMap = CreateCorrelationMap<ImageType>( inputImages );

  typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> MultiplierType;

  MultiplierType::Pointer multiplier = MultiplierType::New();
  multiplier->SetInput( correlationMap );
  multiplier->SetConstant( alpha );

  ImageType::Pointer alphaW = multiplier->GetOutput();
  alphaW->Update();
  alphaW->DisconnectPipeline();

  float energy = itk::NumericTraits<float>::max();
  float newEnergy = itk::NumericTraits<float>::max();
  float epsilon = itk::NumericTraits<float>::max();

  unsigned int numberOfIterations = 0;

  while( epsilon > 1e-3 || numberOfIterations++ > 30 )
    {
    energy = newEnergy;

    newSyntheticImages.clear();
    for( unsigned int n = 0; n < inputImages.size(); n++ )
      {
      typedef itk::LaplacianRecursiveGaussianImageFilter<ImageType> LogFilterType;
      LogFilterType::Pointer logFilter = LogFilterType::New();
      logFilter->SetSigma( 1.2 );
      logFilter->SetInput( syntheticImages[n] );
      logFilter->SetNormalizeAcrossScale( true );
      logFilter->Update();

      MultiplierType::Pointer term1 = MultiplierType::New();
      term1->SetInput1( alphaW );
      term1->SetInput2( logFilter->GetOutput() );
      term1->Update();

      MultiplierType::Pointer term2 = MultiplierType::New();
      term2->SetInput( syntheticImages[n] );
      term2->SetConstant( -(1 + beta) );
      term2->Update();

      MultiplierType::Pointer term4 = MultiplierType::New();
      term4->SetInput( signalImages[n] );
      term4->SetConstant( beta );
      term4->Update();

      typedef itk::AddImageFilter<ImageType, ImageType, ImageType> AdderType;

      AdderType::Pointer adder1 = AdderType::New();
      adder1->SetInput1( term1->GetOutput() );
      adder1->SetInput2( term2->GetOutput() );

      AdderType::Pointer adder2 = AdderType::New();
      adder2->SetInput1( adder1->GetOutput() );
      adder2->SetInput2( inputImages[n] );

      AdderType::Pointer adder3 = AdderType::New();
      adder3->SetInput1( adder2->GetOutput() );
      adder3->SetInput2( term4->GetOutput() );

      MultiplierType::Pointer multiplier2 = MultiplierType::New();
      multiplier2->SetInput( adder3->GetOutput() );
      multiplier2->SetConstant( deltaK );

      AdderType::Pointer adder4 = AdderType::New();
      adder4->SetInput1( multiplier2->GetOutput() );
      adder4->SetInput2( syntheticImages[n] );
      adder4->Update();

      newSyntheticImages.push_back( adder4->GetOutput() );
      }

    syntheticImages = newSyntheticImages;

    newEnergy = CalculateEnergy<ImageType>( inputImages, signalImages, syntheticImages, correlationMap, alpha, beta );

    epsilon = energy - newEnergy;
    }

  for( unsigned int n = 0; n < inputImages.size(); n++ )
    {
    std::ostringstream iss;
    iss << std::setw( 5 ) << std::setfill( '0' ) << n;

    std::string filename = std::string( argv[1] ) + iss.str() + std::string( ".nii.gz" );

    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( filename.c_str() );
    writer->SetInput( newSyntheticImages[n] );
    writer->Update();
    }

  return 0;
}

