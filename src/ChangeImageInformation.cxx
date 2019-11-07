#include <stdio.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkChangeInformationImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkVector.h"

#include <string>
#include <vector>

#include "Common.h"

template <unsigned int ImageDimension, class PixelType>
int ChangeImageInformation( int argc, char *argv[] )
{
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  using FilterType = itk::ChangeInformationImageFilter<ImageType>;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );

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
      filter->SetOutputOrigin( origin );
      filter->ChangeOriginOn();
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
      filter->SetOutputSpacing( spacing );
      filter->ChangeSpacingOn();
      break;
      }
    case 2:
      {
      typename ImageType::DirectionType direction;
      direction.SetIdentity();
      filter->SetOutputDirection( direction );
      filter->ChangeDirectionOn();
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
      filter->SetOutputDirection( direction );
      filter->ChangeDirectionOn();
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

      filter->SetOutputDirection( refReader->GetOutput()->GetDirection() );
      filter->ChangeDirectionOn();
      filter->SetOutputSpacing( refReader->GetOutput()->GetSpacing() );
      filter->ChangeSpacingOn();
      filter->SetOutputOrigin( refReader->GetOutput()->GetOrigin() );
      filter->ChangeOriginOn();

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
  writer->SetInput( filter->GetOutput() );
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
      << " what = 4: imageFile (copies_image_header_information)" << std::endl
      << " what = 5: index" << std::endl;
    exit( 1 );
    }

  std::string type = std::string( argv[1] );

  char pixelType = 'f';

  char dimension = type[0];
  if( type.length() > 1 )
    {
    pixelType = type[1];
    }

  switch( dimension )
   {
   case '2':
     {
     typedef itk::SymmetricSecondRankTensor<float, 2> TensorType;
     typedef itk::Vector<float, 2> VectorType;
     switch( pixelType )
       {
       case 'f':  default:
         ChangeImageInformation<2, float>( argc, argv );
         break;
       case 'v':
         ChangeImageInformation<2, VectorType>( argc, argv );
         break;
       case 't':
         ChangeImageInformation<2, TensorType>( argc, argv );
         break;
       }
     break;
     }
   case '3':
     {
     typedef itk::SymmetricSecondRankTensor<float, 3> TensorType;
     typedef itk::Vector<float, 3> VectorType;
     switch( pixelType )
       {
       case 'f':  default:
         ChangeImageInformation<3, float>( argc, argv );
         break;
       case 'v':
         ChangeImageInformation<3, VectorType>( argc, argv );
         break;
       case 't':
         ChangeImageInformation<3, TensorType>( argc, argv );
         break;
       }
     break;
     }
   case '4':
     {
     typedef itk::SymmetricSecondRankTensor<float, 4> TensorType;
     typedef itk::Vector<float, 4> VectorType;
     switch( pixelType )
       {
       case 'f':  default:
         ChangeImageInformation<4, float>( argc, argv );
         break;
       case 'v':
         ChangeImageInformation<4, VectorType>( argc, argv );
         break;
       case 't':
         ChangeImageInformation<4, TensorType>( argc, argv );
         break;
       }
     break;
     }
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

