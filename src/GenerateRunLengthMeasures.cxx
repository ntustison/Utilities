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
int GenerateRunLengthMeasures( int argc, char *argv[] )
{
  typedef int PixelType;
  typedef float RealType;

  unsigned int numberOfBins = 200;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[2] );
  imageReader->Update();

  typedef itk::Statistics::MaskedScalarImageToGreyLevelRunLengthMatrixGenerator<ImageType>
    RunLengthMatrixGeneratorType;
  typename RunLengthMatrixGeneratorType::Pointer generator = RunLengthMatrixGeneratorType::New();
  typename RunLengthMatrixGeneratorType::OffsetVectorPointer offsets
    = RunLengthMatrixGeneratorType::OffsetVector::New();

  unsigned int nOffsets = atoi( argv[3] );
  for( unsigned int n = 0; n < nOffsets; n++ )
    {
    typename RunLengthMatrixGeneratorType::OffsetType offset;
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
  if ( argc > static_cast<int>( 5 + nOffsets ) )
    {
    label = static_cast<PixelType>( atoi( argv[5 + nOffsets] ) );
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
  typename BoundingBoxType::Pointer bbox = BoundingBoxType::New();
  typename BoundingBoxType::PointsContainerPointer Points
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

  typename BoundingBoxType::PointType pointMin = bbox->GetMinimum();
  typename BoundingBoxType::PointType pointMax = bbox->GetMaximum();

  generator->SetInput( imageReader->GetOutput() );
  generator->SetImageMask( mask );
  generator->SetNumberOfBinsPerAxis( numberOfBins );
  generator->SetPixelValueMinMax( minValue, maxValue );
  generator->SetDistanceValueMinMax( 0, pointMin.EuclideanDistanceTo( pointMax ) );
  generator->SetInsidePixelValue( 1 );
  generator->Compute();

  typedef itk::Statistics::GreyLevelRunLengthMatrixTextureCoefficientsCalculator
    <typename RunLengthMatrixGeneratorType::HistogramType> CalculatorType;
  typename CalculatorType::Pointer calculator = CalculatorType::New();
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

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage "
     << "[nOffsets] [offset1] [offset2] ... [offsetn]  [labelImage] [label]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     GenerateRunLengthMeasures<2>( argc, argv );
     break;
   case 3:
     GenerateRunLengthMeasures<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

