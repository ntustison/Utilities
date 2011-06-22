#include <stdio.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkBinaryThresholdImageFunction.h"
#include "itkShapedFloodFilledImageFunctionConditionalConstIterator.h"

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
int FloodFill( int argc, char *argv[] )
{
  typedef float PixelType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef unsigned int LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;
  typename LabelImageType::Pointer labelImage = LabelImageType::New();
  labelImage->SetOrigin( reader->GetOutput()->GetOrigin() );
  labelImage->SetSpacing( reader->GetOutput()->GetSpacing() );
  labelImage->SetDirection( reader->GetOutput()->GetDirection() );
  labelImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  labelImage->Allocate();
  labelImage->FillBuffer( 0 );

  PixelType lowerThreshold = static_cast<PixelType>( atof( argv[4] ) );
  PixelType upperThreshold = static_cast<PixelType>( atof( argv[5] ) );

  typedef itk::BinaryThresholdImageFunction<ImageType> FunctionType;
  typename FunctionType::Pointer function = FunctionType::New();

  function->SetInputImage ( reader->GetOutput() );
  function->ThresholdBetween( lowerThreshold, upperThreshold );

  std::vector<typename ImageType::IndexType> seedList;
  for( unsigned int n = 7; n < static_cast<unsigned int>( argc ); n++ )
    {
    std::vector<int> idx = ConvertVector<int>( std::string( argv[n] ) );
    typename ImageType::IndexType index;
    if( idx.size() != ImageDimension )
      {
      std::cerr << "Incorrect index specification." << std::endl;
      return EXIT_FAILURE;
      }
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      index[d] = idx[d];
      }
    seedList.push_back( index );
    }

  typedef itk::ShapedFloodFilledImageFunctionConditionalConstIterator<
    ImageType, FunctionType> ShapedFloodFilledIteratorType;

  ShapedFloodFilledIteratorType fIt( reader->GetOutput(), function, seedList );
  fIt.SetFullyConnected( static_cast<bool>( atoi( argv[6] ) ) );
  for( fIt.GoToBegin(); !fIt.IsAtEnd(); ++fIt )
    {
    labelImage->SetPixel( fIt.GetIndex(), itk::NumericTraits<LabelType>::One );
    }

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( labelImage );
  writer->SetFileName( argv[3] );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 8 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension "
      << "inputImage outputImage lowerThreshold upperThreshold useFullConnectedness "
      << "seedIndex1 [seedIndex2] ... [seedIndexN]"<< std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     FloodFill<2>( argc, argv );
     break;
   case 3:
     FloodFill<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

