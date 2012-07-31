#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPathIterator.h"
#include "itkMinimalPathImageFunction.h"

#include <string>
#include <vector>

#include <fstream>

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

template<unsigned int ImageDimension>
int MinimalPath( unsigned int argc, char *argv[] )
{
  typedef float PixelType;
  typedef unsigned char PathPixelType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<PathPixelType, ImageDimension> PathImageType;
  typedef typename ImageType::IndexType IndexType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename PathImageType::Pointer pathImage = PathImageType::New();
  pathImage->SetDirection( reader->GetOutput()->GetDirection() );
  pathImage->SetOrigin( reader->GetOutput()->GetOrigin() );
  pathImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  pathImage->SetSpacing( reader->GetOutput()->GetSpacing() );
  pathImage->Allocate();
  pathImage->FillBuffer( itk::NumericTraits<PathPixelType>::Zero );

  std::vector<unsigned int> anchor
    = ConvertVector<unsigned int>( std::string( argv[4] ) );
  std::vector<unsigned int> free
    = ConvertVector<unsigned int>( std::string( argv[5] ) );

  IndexType anchorIndex;
  IndexType freeIndex;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    anchorIndex[d] = anchor[d];
    freeIndex[d] = free[d];
    }

  typedef itk::MinimalPathImageFunction<ImageType> FunctionType;
  typename FunctionType::Pointer function = FunctionType::New();
  function->SetUseFaceConnectedness( argc > 6 ?
    static_cast<bool>( atoi( argv[6] ) ) : true );
  function->SetUseImageSpacing( true );
  function->SetInputImage( reader->GetOutput() );
  function->SetAnchorSeed( anchorIndex );

  std::string txtFileName = std::string( argv[3] ) + std::string( ".txt" );
  std::ofstream str( txtFileName.c_str() );
  str << "0 0 0 0" << std::endl;

  typedef itk::PathIterator<ImageType,
    typename FunctionType::OutputType> IteratorType;

  typename FunctionType::OutputType::Pointer path
    = function->EvaluateAtIndex( freeIndex );
  IteratorType It( reader->GetOutput(), path );
  It.GoToBegin();
  while ( !It.IsAtEnd() )
    {
    typename PathImageType::PointType point;
    pathImage->TransformIndexToPhysicalPoint( It.GetIndex(), point );
    str << point[0] << " " << point[1];
    if( ImageDimension == 3 )
      {
      str << " " << point[2];
      }
    else
      {
      str << " 0";
      }
    str << " 1" << std::endl;

    pathImage->SetPixel( It.GetIndex(),
      itk::NumericTraits<PathPixelType>::One );
    ++It;
    }
  str << "0 0 0 0" << std::endl;
  str.close();


  std::string imageFileName = std::string( argv[3] ) + std::string( ".nii.gz" );
  typedef itk::ImageFileWriter<PathImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( imageFileName.c_str() );
  writer->SetInput( pathImage );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 6 )

    {

    std::cout << "Usage: "<< argv[0]
      << " imageDimension inputImage outputPrefix anchorIndex"
      << " freeIndex [useFaceConnectedness]" << std::endl;

    exit( 1 );

    }


  switch( atoi( argv[1] ) )
   {
   case 2:
     MinimalPath<2>( argc, argv );
     break;
   case 3:
     MinimalPath<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

