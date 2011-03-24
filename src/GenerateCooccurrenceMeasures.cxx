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

#include <string>
#include <vector>

template<class TValue>
TValue Convert( std::string optionString )
			{
			TValue value;
			std::istringstream iss( optionString );
			iss >> value;
			return value;
			}

template<class TValue>
std::vector<TValue> ConvertVector( std::string optionString )
			{
			std::vector<TValue> values;
			std::string::size_type crosspos = optionString.find( 'x', 0 );

			if ( crosspos == std::string::npos )
					{
					values.push_back( Convert<TValue>( optionString ) );
					}
			else
					{
					std::string element = optionString.substr( 0, crosspos ) ;
					TValue value;
					std::istringstream iss( element );
					iss >> value;
					values.push_back( value );
					while ( crosspos != std::string::npos )
							{
							std::string::size_type crossposfrom = crosspos;
							crosspos = optionString.find( 'x', crossposfrom + 1 );
							if ( crosspos == std::string::npos )
									{
									element = optionString.substr( crossposfrom + 1, optionString.length() );
									}
							else
									{
									element = optionString.substr( crossposfrom + 1, crosspos ) ;
									}
							std::istringstream iss( element );
							iss >> value;
							values.push_back( value );
							}
					}
			return values;
			}

template <unsigned int ImageDimension>
int GenerateCooccurrenceMeasures( int argc, char *argv[] )
{

  typedef float PixelType;
  typedef float RealType;

  unsigned int numberOfBins = 200;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[2] );
  imageReader->Update();

  typedef itk::Statistics::MaskedScalarImageToGreyLevelCooccurrenceMatrixGenerator<ImageType>
    CooccurrenceMatrixGeneratorType;
  typename CooccurrenceMatrixGeneratorType::Pointer generator
    = CooccurrenceMatrixGeneratorType::New();

  typename CooccurrenceMatrixGeneratorType::OffsetVectorPointer offsets
    = CooccurrenceMatrixGeneratorType::OffsetVector::New();

  unsigned int nOffsets = atoi( argv[3] );
  for( unsigned int n = 0; n < nOffsets; n++ )
    {
    typename CooccurrenceMatrixGeneratorType::OffsetType offset;
    std::vector<int> vector
      = ConvertVector<int>( std::string( argv[4 + n] ) );
    if( vector.size() > 0 && vector.size() != ImageDimension )
      {
      std::cerr << "Error:  offset size does not equal image dimension." << std::endl;
      return EXIT_FAILURE;
      }
    else
      {
      for ( unsigned int d = 0; d < ImageDimension; d++ )
        {
        offset[d] = vector[d];
        }
      }
    offsets->push_back( offset );
    }
  generator->SetOffsets( offsets );


  typename ImageType::Pointer mask = ImageType::New();
  if ( argc > static_cast<int>( 4 + nOffsets ) )
    {
    typename ReaderType::Pointer labelImageReader = ReaderType::New();
    labelImageReader->SetFileName( argv[4 + nOffsets] );
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
  if ( argc >static_cast<int>( 5 + nOffsets ) )
    {
    label = static_cast<PixelType>( atoi( argv[5 + nOffsets] ) );
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

  generator->SetInput( imageReader->GetOutput() );
  generator->SetImageMask( mask );
  generator->SetNumberOfBinsPerAxis( numberOfBins );
  generator->SetPixelValueMinMax( minValue, maxValue );
  generator->Compute();

  typedef itk::Statistics::GreyLevelCooccurrenceMatrixTextureCoefficientsCalculator
    <typename CooccurrenceMatrixGeneratorType::HistogramType> CalculatorType;
  typename CalculatorType::Pointer calculator = CalculatorType::New();
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


int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage "
     << "nOffsets [offset1] [offset2] ... [offsetn]  [labelImage] [label]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GenerateCooccurrenceMeasures<2>( argc, argv );
     break;
   case 3:
     GenerateCooccurrenceMeasures<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

