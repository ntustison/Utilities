#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkBoundingBox.h"
#include "itkCastImageFilter.h"
#include "itkGreyLevelRunLengthMatrixTextureCoefficientsCalculator.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMaskedScalarImageToGreyLevelRunLengthMatrixGenerator.h"

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

  unsigned long numberOfPixelsInMask = 0;
  for ( ItM.GoToBegin(), ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItM, ++ItI )
    {
    if ( ItM.Get() == label )
      {
      numberOfPixelsInMask++;
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

  // Generate bounding box of the masked image to determine

  typedef itk::BoundingBox<unsigned long,
       ImageDimension, RealType> BoundingBoxType;
  BoundingBoxType::Pointer bbox = BoundingBoxType::New();
  BoundingBoxType::PointsContainerPointer Points
       = BoundingBoxType::PointsContainer::New();
  itk::Point<RealType, ImageDimension> point;

  int idx = 0;
  itk::ImageRegionIteratorWithIndex<ImageType> It(
    mask, mask->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( It.Get() > 0 )
      {
      mask->TransformIndexToPhysicalPoint( It.GetIndex(), point );
      Points->InsertElement( idx++, point );
      }
    }
  bbox->SetPoints( Points );
  bbox->ComputeBoundingBox();

  BoundingBoxType::PointType pointMin = bbox->GetMinimum();
  BoundingBoxType::PointType pointMax = bbox->GetMaximum();

  typedef itk::Statistics::MaskedScalarImageToGreyLevelRunLengthMatrixGenerator<ImageType>
    RunLengthMatrixGeneratorType;
  RunLengthMatrixGeneratorType::Pointer generator = RunLengthMatrixGeneratorType::New();

  generator->SetInput( imageReader->GetOutput() );
  generator->SetImageMask( mask );
  generator->SetNumberOfBinsPerAxis( numberOfBins );
  generator->SetPixelValueMinMax( minValue, maxValue );
  generator->SetDistanceValueMinMax( 0, pointMin.EuclideanDistanceTo( pointMax ) );
  generator->SetInsidePixelValue( 1 );
  generator->Compute();

  typedef itk::Statistics::GreyLevelRunLengthMatrixTextureCoefficientsCalculator
    <RunLengthMatrixGeneratorType::HistogramType> CalculatorType;
  CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetHistogram( generator->GetOutput() );
  calculator->Compute();

  RealType sre = calculator->GetShortRunEmphasis();
  RealType lre = calculator->GetLongRunEmphasis();
  RealType gln = calculator->GetGreyLevelNonuniformity();
  RealType rln = calculator->GetRunLengthNonuniformity();
  RealType rp = static_cast<RealType>( calculator->GetTotalNumberOfRuns() )
    / static_cast<RealType>( numberOfPixelsInMask );
  RealType lgre = calculator->GetLowGreyLevelRunEmphasis();
  RealType hgre = calculator->GetHighGreyLevelRunEmphasis();
  RealType srlge = calculator->GetShortRunLowGreyLevelEmphasis();
  RealType srhge = calculator->GetShortRunHighGreyLevelEmphasis();
  RealType lrlge = calculator->GetLongRunLowGreyLevelEmphasis();
  RealType lrhge = calculator->GetLongRunHighGreyLevelEmphasis();

  std::cout << "[" << argv[0] << "]" << std::endl;
  std::cout << sre << " "
            << lre << " "
            << gln << " "
            << rln << " "
            << rp << " "
            << lgre << " "
            << hgre << " "
            << srlge << " "
            << srhge << " "
            << lrlge << " "
            << lrhge << std::endl;

/*
  std::cout << "SRE: " << sre << std::endl;
  std::cout << "LRE: " << lre << std::endl;
  std::cout << "GLN: " << gln << std::endl;
  std::cout << "RLN: " << rln << std::endl;
  std::cout << "RP : " << rp << std::endl;
  std::cout << "LGRE: " << lgre << std::endl;
  std::cout << "HGRE: " << hgre << std::endl;
  std::cout << "SRLGE: " << srlge << std::endl;
  std::cout << "SRHGE: " << srhge << std::endl;
  std::cout << "LRLGE: " << lrlge << std::endl;
  std::cout << "LRHGE: " << lrhge << std::endl;
*/

  return 0;
}
