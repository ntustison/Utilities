#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkFlipImageFilter.h"

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
int FlipImage( int argc, char *argv[] )
{

  typedef float PixelType; 
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();
  
  typedef itk::FlipImageFilter<ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  typename FilterType::FlipAxesArrayType array;

  std::vector<bool> which = ConvertVector<bool>( std::string( argv[4] ) );
  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    array[d] = which[d]; 
    }
  filter->SetInput( reader->GetOutput() );
  filter->SetFlipAxes( array ); 
  filter->Update();
  
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << argv[0] << " imageDimension inputImage outputImage whichAxes[e.g.1x0x0]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     FlipImage<2>( argc, argv );
     break;
   case 3:
     FlipImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

