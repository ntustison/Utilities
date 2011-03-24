#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkConstantPadImageFilter.h"

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
int PadImage( unsigned int argc, char *argv[] )
{

  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();
  
  unsigned long lowerBound[ImageDimension];
  unsigned long upperBound[ImageDimension];

  std::vector<unsigned int> before 
    = ConvertVector<unsigned int>( std::string( argv[4] ) );
  std::vector<unsigned int> after 
    = ConvertVector<unsigned int>( std::string( argv[5] ) );

  if ( before.size() == 1 )
    {
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      lowerBound[d] = before[0]; 
      }
    }
  else if ( before.size() == ImageDimension )
    {
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      lowerBound[d] = before[d]; 
      }
    }
  else
    {
    std::cerr << "Invalid padding." << std::endl; 
    }

  if ( after.size() == 1 )
    {
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      upperBound[d] = after[0]; 
      }
    }
  else if ( after.size() == ImageDimension )
    {
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      upperBound[d] = after[d]; 
      }
    }
  else
    {
    std::cerr << "Invalid padding." << std::endl; 
    }

  typedef itk::ConstantPadImageFilter<ImageType, ImageType> PadderType;
  typename PadderType::Pointer padder = PadderType::New();
  padder->SetInput( reader->GetOutput() );
  padder->SetPadLowerBound( lowerBound );
  padder->SetPadUpperBound( upperBound );
  if ( argc > 6 )
    {
    padder->SetConstant( static_cast<PixelType>( atof( argv[6] ) ) );
    }  
  padder->Update();

  typename ImageType::PointType origin = padder->GetOutput()->GetOrigin();
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin[i] -= ( padder->GetOutput()->GetSpacing()[i] 
      * static_cast<double>( lowerBound[i] ) ); 
    }
  padder->GetOutput()->SetOrigin( origin );
  
  
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( padder->GetOutput() );
  writer->SetFileName( argv[3] );                                          
  writer->Update();
    
  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage "
      << "outputImage lower[0]xlower[1]... upper[0]xupper[1]... [constant]" << std::endl;
    exit( 1 );
    }
  
  switch( atoi( argv[1] ) ) 
   {
   case 2:
     PadImage<2>( argc, argv );
     break;
   case 3:
     PadImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


