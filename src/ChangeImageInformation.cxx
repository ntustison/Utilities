#include <stdio.h>

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
int ChangeImageInformation( int argc, char *argv[] )
{
  typedef float PixelType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  switch( atoi( argv[4] ) )
    {
    case 0:
      {
      std::vector<float> what = ConvertVector<float>( std::string( argv[5] ) );
      if( what.size() != ImageDimension )
        {
        std::cerr << "Origin/spacing dimension does not equal image dimension."
          << std::endl;
        return EXIT_FAILURE;
        }
      typename ImageType::PointType origin;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        origin[d] = what[d];
        }
      reader->GetOutput()->SetOrigin( origin );
      break;
      }
    case 1:
      {
      std::vector<float> what = ConvertVector<float>( std::string( argv[5] ) );
      if( what.size() != ImageDimension )
        {
        std::cerr << "Origin/spacing dimension does not equal image dimension."
          << std::endl;
        return EXIT_FAILURE;
        }
      typename ImageType::SpacingType spacing;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        spacing[d] = what[d];
        }
      reader->GetOutput()->SetSpacing( spacing );
      break;
      }
    case 2:
      {
      typename ImageType::DirectionType direction;
      direction.SetIdentity();
      reader->GetOutput()->SetDirection( direction );
      break;
      }
    case 3:
      {
      std::vector<float> what = ConvertVector<float>( std::string( argv[5] ) );
      if( what.size() != ImageDimension*ImageDimension )
        {
        std::cerr << "matrix dimension does not equal image dimension*imagedimension."
          << std::endl;
        return EXIT_FAILURE;
        }
      typename ImageType::DirectionType direction;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        for( unsigned int j = 0; j < ImageDimension; j++ )
          {
          direction[i][j] = what[i*ImageDimension+j];
          }
        }
      reader->GetOutput()->SetDirection( direction );
      break;
      }
    case 4:
      {
						typename ReaderType::Pointer refReader = ReaderType::New();
						refReader->SetFileName( argv[5] );
      try
        {
  						refReader->Update();
        }
      catch( ... )
        {
        std::cerr << "Unable to read image file." << std::endl;
        return EXIT_FAILURE;
        }

      reader->GetOutput()->SetOrigin( refReader->GetOutput()->GetOrigin() );
      reader->GetOutput()->SetSpacing( refReader->GetOutput()->GetSpacing() );
      reader->GetOutput()->SetDirection( refReader->GetOutput()->GetDirection() );
      break;
      }
    case 5:
      {
      std::vector<int> what = ConvertVector<int>( std::string( argv[5] ) );
      if( what.size() != ImageDimension )
        {
        std::cerr << "Origin/spacing dimension does not equal image dimension."
          << std::endl;
        return EXIT_FAILURE;
        }
      typename ImageType::IndexType index;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        index[d] = what[d];
        }
      typename ImageType::RegionType region;
      region.SetIndex( index );
      region.SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );

      reader->GetOutput()->SetRegions( region );
      break;
      }
    default:
      {
						std::cerr << "Unrecognized option." << std::endl;
						return EXIT_FAILURE;
      }
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( reader->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension "
      << "inputImage outputImage what" << std::endl
      << " what = 0: origin" << std::endl
      << " what = 1: spacing" << std::endl
      << " what = 2: set direction to identity" << std::endl
      << " what = 3: set direction to vector m_{11}xm_{12}xm_{13}xm_{21}xm_{22}xm_{23}xm_{31}xm_{32}xm_{33}" << std::endl
      << " what = 4: imageFile  copies image header information" << std::endl
      << " what = 5: index" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     ChangeImageInformation<2>( argc, argv );
     break;
   case 3:
     ChangeImageInformation<3>( argc, argv );
     break;
   case 4:
     ChangeImageInformation<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

