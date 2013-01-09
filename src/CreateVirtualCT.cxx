#include "itkBresenhamLine.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLabelContourImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkVector.h"

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
  typedef float PixelType;
  typedef unsigned int LabelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
  typename LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( argv[3] );
  labelReader->Update();

  typename ImageType::IndexType targetIndex;

  vnl_vector<float> centerOfMass( ImageDimension );
  centerOfMass.fill( 0.0 );
  float N = 0.0;

  itk::ImageRegionIteratorWithIndex<ImageType> It(
    reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    typename ImageType::IndexType index = It.GetIndex();
    PixelType weight = It.Get();

    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      centerOfMass[d] += ( weight * index[d] );
      }
    N += weight;
    }
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    centerOfMass[d] /= N;
    }
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    targetIndex[d] = static_cast<int>( centerOfMass[d] );
    }
  std::cout << "Target index = " << targetIndex << std::endl;


  typename ImageType::Pointer output = ImageType::New();
  output->CopyInformation( reader->GetOutput() );
  output->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  output->Allocate();
  output->FillBuffer( 0 );

  typedef itk::BresenhamLine<ImageDimension> LinerType;
  LinerType liner;
  typedef typename LinerType::LType VectorType;
  typedef typename LinerType::OffsetType OffsetType;
  typedef typename LinerType::IndexType IndexType;

  typename ImageType::SizeType size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();

  unsigned long maxLength = 0;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    maxLength += size[d] * size[d];
    }
  maxLength = static_cast<unsigned long>( vcl_sqrt( maxLength ) );

  itk::ImageRegionIteratorWithIndex<LabelImageType> ItL(
    labelReader->GetOutput(), labelReader->GetOutput()->GetLargestPossibleRegion() );
  for( ItL.GoToBegin(); !ItL.IsAtEnd(); ++ItL )
    {
    if( ItL.Get() == 1 )
      {
      typename ImageType::IndexType startIndex = ItL.GetIndex();

      VectorType direction;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        direction[d] = startIndex[d] - targetIndex[d];
        }
      direction.Normalize();

      typename LinerType::OffsetArray offsets = liner.BuildLine( direction, maxLength );

      IndexType currentIndex = targetIndex;

      bool isFound = false;
      unsigned int begOffsetIndex = 0;
      unsigned int endOffsetIndex = 0;
      IndexType begIndex;
      IndexType endIndex;

      typename LinerType::OffsetArray::const_iterator it;

      for( it = offsets.begin(); it != offsets.end(); ++it )
        {
        if( begOffsetIndex == 0 && labelReader->GetOutput()->GetLargestPossibleRegion().IsInside( currentIndex ) == 3 )
          {
          begOffsetIndex = it - offsets.begin();
          begIndex = currentIndex;
          }
        if( !labelReader->GetOutput()->GetLargestPossibleRegion().IsInside( currentIndex ) )
          {
          break;
          }
        if( startIndex == currentIndex && begOffsetIndex != 0 )
          {
          if( labelReader->GetOutput()->GetPixel( currentIndex + offsets[it - offsets.begin() + 1] ) == 0 )
            {
            isFound = true;
            endOffsetIndex = it - offsets.begin();
            endIndex = currentIndex;
            break;
            }
          }
        currentIndex = targetIndex + *it;
        }

      // Debug
      if( isFound )
        {
        output->SetPixel( begIndex, 1 );
        output->SetPixel( endIndex, 2 );
        }

      // check to see if we're on the outer boundary

//         if( offsetIndex > offsets.size()



      }
    }








  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( output );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << argv[0] << " imageDimension inputT1 labelMask outputVirtualCT" << std::endl;
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

