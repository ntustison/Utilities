#include "itkBresenhamLine.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

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
int DrawLines( int argc, char *argv[] )
{
  typedef unsigned int PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename ImageType::IndexType targetIndex;

  if( argc > 4 )
    {
    std::vector<int> point = ConvertVector<int>( std::string( argv[4] ) );
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      targetIndex[d] = point[d];
      }
    }
  else
    {
    vnl_vector<float> centerOfMass( ImageDimension );
    centerOfMass.fill( 0.0 );
    float N = 0.0;

    itk::ImageRegionIteratorWithIndex<ImageType> It(
      reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if( It.Get() == 1 )
        {
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          centerOfMass[d] += It.GetIndex()[d];
          }
        N++;
        }
      }
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      centerOfMass[d] /= N;
      }
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      targetIndex[d] = static_cast<int>( centerOfMass[d] );
      }
    }

  std::cout << "Target index = " << targetIndex << std::endl;

  typedef itk::BresenhamLine<ImageDimension> LinerType;
  LinerType liner;

  itk::ImageRegionIteratorWithIndex<ImageType> It(
    reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( It.Get() == 1 )
      {
      typename ImageType::IndexType startIndex = It.GetIndex();

      typename LinerType::LType direction;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        direction[d] = targetIndex[d] - startIndex[d];
        }
      unsigned int length = static_cast<unsigned int>( direction.GetNorm() );
      typename LinerType::OffsetArray offsets = liner.BuildLine( direction, length );

      typename LinerType::OffsetArray::const_iterator it;
      for( it = offsets.begin(); it != offsets.end(); it++ )
        {
        typename ImageType::IndexType currentIndex = startIndex + *it;
        if( reader->GetOutput()->GetPixel( currentIndex ) == 0 )
          {
          reader->GetOutput()->SetPixel( currentIndex, 2 );
          }
        }
      }
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( reader->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 2 )
    {
    std::cout << argv[0] << " imageDimension mask1 mask2 centerPoint outputFile" << std::endl;
    std::cout << "     Note:  inputMasks are assumed to be 1/0 with mask label = 1." << std::endl;
    exit( 0 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     DrawLines<2>( argc, argv );
     break;
   case 3:
     DrawLines<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

