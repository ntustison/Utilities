#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkCastImageFilter.h"
#include "itkGreyLevelCooccurrenceMatrixTextureCoefficientsCalculator.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMaskedScalarImageToGreyLevelCooccurrenceMatrixGenerator.h"

int main( int argc, char *argv[] )
{
  if ( argc < 2 )
    {
    std::cout << "Usage: " << argv[0] << " image [labelImage] [label]" << std::endl;
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

  ImageType::Pointer mask = ImageType::New();
  if ( argc > 2 )
    {
    ReaderType::Pointer labelImageReader = ReaderType::New();
    labelImageReader->SetFileName( argv[2] );
    labelImageReader->Update();
    mask = labelImageReader->GetOutput();
    } 
  else
    {
    mask->SetOrigin( imageReader->GetOutput()->GetOrigin() );
    mask->SetSpacing( imageReader->GetOutput()->GetSpacing() );
    mask->SetRegions( imageReader->GetOutput()->GetLargestPossibleRegion() );
    mask->Allocate();
    mask->FillBuffer( itk::NumericTraits<PixelType>::Zero );
    } 
  PixelType label = itk::NumericTraits<PixelType>::One;
  if ( argc > 3 )
    {
    label = static_cast<PixelType>( atoi( argv[3] ) );
    } 

  itk::ImageRegionIterator<ImageType> ItI( imageReader->GetOutput(),
    imageReader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> ItM( mask,
    mask->GetLargestPossibleRegion() );

  PixelType maxValue = itk::NumericTraits<PixelType>::min();
  PixelType minValue = itk::NumericTraits<PixelType>::max();

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
    }   
     
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
  typedef itk::Statistics::MaskedScalarImageToGreyLevelCooccurrenceMatrixGenerator<ImageType>
    CooccurrenceMatrixGeneratorType;
  CooccurrenceMatrixGeneratorType::Pointer generator = CooccurrenceMatrixGeneratorType::New();
  CooccurrenceMatrixGeneratorType::OffsetType offset;
  offset.Fill( -1 );
  generator->SetOffset( offset );
      


  generator->SetInput( imageReader->GetOutput() );
  generator->SetImageMask( mask );
  generator->SetNumberOfBinsPerAxis( numberOfBins );
  generator->SetPixelValueMinMax( minValue, maxValue );
  generator->Compute();   

  typedef itk::Statistics::GreyLevelCooccurrenceMatrixTextureCoefficientsCalculator
    <CooccurrenceMatrixGeneratorType::HistogramType> CalculatorType;
  CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetHistogram( generator->GetOutput() );
  calculator->Compute();

  RealType energy = calculator->GetEnergy();
  RealType entropy = calculator->GetEntropy();
  RealType correlation = calculator->GetCorrelation();
  RealType inverseDifferenceMoment = calculator->GetInverseDifferenceMoment();
  RealType inertia = calculator->GetInertia();
  RealType clusterShade = calculator->GetClusterShade();
  RealType clusterProminence = calculator->GetClusterProminence();
  RealType haralickCorrelation = calculator->GetHaralickCorrelation();

  std::cout << "[" << argv[0] << "]" << std::endl;
  std::cout << energy << " "
            << entropy << " " 
            << correlation << " "
            << inverseDifferenceMoment << " "
            << inertia << " "
            << clusterShade << " "
            << clusterProminence << " "
            << haralickCorrelation << std::endl;

/*
  std::cout << "energy             : " << energy << std::endl;
  std::cout << "entropy            : " << entropy << std::endl;
  std::cout << "correlation        : " << correlation << std::endl;
  std::cout << "inverse diff moment: " << inverseDifferenceMoment << std::endl;
  std::cout << "inertia            : " << inertia << std::endl;
  std::cout << "cluster Shade      : " << clusterShade << std::endl;
  std::cout << "cluster prominence : " << clusterProminence << std::endl;
  std::cout << "haralickCorrelation: " << haralickCorrelation << std::endl;
*/

  return 0;
}
