#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkCollidingFrontsImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"

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
							std::istringstream iss2( element );
							iss2 >> value;
							values.push_back( value );
							}
					}
			return values;
			}

template <unsigned int ImageDimension>
int CollidingFronts( unsigned int argc, char *argv[] )
{

  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::ImageFileReader<LabelImageType> SeedsReaderType;
  typename SeedsReaderType::Pointer seedsReader = SeedsReaderType::New();
  seedsReader->SetFileName( argv[3] );
  seedsReader->Update();

  typedef itk::CollidingFrontsImageFilter<ImageType, ImageType> FilterType;


  typedef typename FilterType::NodeContainer NodeContainerType;
  typedef typename FilterType::NodeType NodeType;
  typename NodeContainerType::Pointer seeds1 = NodeContainerType::New();
  typename NodeContainerType::Pointer seeds2 = NodeContainerType::New();

  unsigned long seeds1Counter = 0;
  unsigned long seeds2Counter = 0;

  itk::ImageRegionConstIteratorWithIndex<LabelImageType> It( seedsReader->GetOutput(),
    seedsReader->GetOutput()->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( It.Get() == 1 )
      {
      typename LabelImageType::IndexType position = It.GetIndex();

      NodeType node;
      const double value = 0.0;

      node.SetValue( value );
      node.SetIndex( position );

      seeds1->InsertElement( seeds1Counter++, node );
      }
    else if( It.Get() == 2 )
      {
      typename LabelImageType::IndexType position = It.GetIndex();

      NodeType node;
      const double value = 0.0;

      node.SetValue( value );
      node.SetIndex( position );

      seeds2->InsertElement( seeds2Counter++, node );
      }
    }

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetSeedPoints1( seeds1 );
  filter->SetSeedPoints2( seeds2 );
//   filter->ApplyConnectivityOn();
  filter->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension inputImage seedImage outputImage" << std::endl;
    std::cout << "Note:  The seed image needs to have two labels (1 and 2) designating the two sets of seed points." << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     CollidingFronts<2>( argc, argv );
     break;
   case 3:
     CollidingFronts<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

