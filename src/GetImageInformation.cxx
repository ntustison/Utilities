#include "itkImage.h"
#include "itkImageFileReader.h"

#include "Common.h"

template <unsigned int ImageDimension>
int GetImageInformation( int argc, char *argv[] )
{
  typedef float PixelType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  if( argc > 3 )
    {
    switch( atoi( argv[3] ) )
      {
      case 0:
        {
        for( unsigned int d = 0; d < ImageDimension-1; d++ )
          {
          std::cout << reader->GetOutput()->GetOrigin()[d] << 'x';
          }
        std::cout << reader->GetOutput()->GetOrigin()[ImageDimension-1] << std::endl;
        break;
        }
      case 1:
        {
        for( unsigned int d = 0; d < ImageDimension-1; d++ )
          {
          std::cout << reader->GetOutput()->GetSpacing()[d] << 'x';
          }
        std::cout << reader->GetOutput()->GetSpacing()[ImageDimension-1] << std::endl;
        break;
        }
      case 2:
        {
        for( unsigned int d = 0; d < ImageDimension-1; d++ )
          {
          std::cout << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[d] << 'x';
          }
        std::cout << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[ImageDimension-1] << std::endl;
        break;
        }
      case 3:
        {
        for( unsigned int d = 0; d < ImageDimension-1; d++ )
          {
          std::cout << reader->GetOutput()->GetLargestPossibleRegion().GetIndex()[d] << 'x';
          }
        std::cout << reader->GetOutput()->GetLargestPossibleRegion().GetIndex()[ImageDimension-1] << std::endl;
        break;
        }
      case 4:
        {
        for( unsigned int di = 0; di < ImageDimension; di++ )
          {
          for( unsigned int dj = 0; dj < ImageDimension; dj++ )
            {
            std::cout << reader->GetOutput()->GetDirection()[di][dj];
            if( di == dj && di == ImageDimension-1 )
              {
              std::cout << std::endl;
              }
            else
              {
              std::cout << 'x';
              }
            }
          }
        break;
        }
      case 5:
        {
        typename ImageType::IndexType index;
        std::vector<unsigned int> idx = ConvertVector<unsigned int>(
          std::string( argv[4] ) );
        for ( unsigned int d = 0; d < ImageDimension; d++ )
          {
          index[d] = idx[d];
          }
        std::cout << reader->GetOutput()->GetPixel( index ) << std::endl;
        }
      }
    }
  else
    {
    typename ImageType::PointType point;
    typename ImageType::PointType center;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      point[d] = static_cast<float>(
        reader->GetOutput()->GetLargestPossibleRegion().GetSize()[d] - 1 )
        * reader->GetOutput()->GetSpacing()[d];
      center[d] = reader->GetOutput()->GetOrigin()[d] + 0.5*static_cast<float>(
        reader->GetOutput()->GetLargestPossibleRegion().GetSize()[d] - 1 )
        * reader->GetOutput()->GetSpacing()[d];
      }

    std::cout << "Image information" << std::endl;
    std::cout << "  Size:          " << reader->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout << "  Origin:        " << reader->GetOutput()->GetOrigin() << std::endl;
    std::cout << "  SpatialExtent: " << point << std::endl;
    std::cout << "  Center:        " << center << std::endl;
    std::cout << "  Spacing:       " << reader->GetOutput()->GetSpacing() << std::endl;
    std::cout << "  Index:         " << reader->GetOutput()->GetLargestPossibleRegion().GetIndex() << std::endl;
    std::cout << "  Direction:     " << std::endl << reader->GetOutput()->GetDirection() << std::endl;
    }

  return 0;
}


int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension image [what]" << std::endl;
    std::cerr << "   what = 0: origin" << std::endl;
    std::cerr << "   what = 1: spacing" << std::endl;
    std::cerr << "   what = 2: size" << std::endl;
    std::cerr << "   what = 3: index" << std::endl;
    std::cerr << "   what = 4: direction" << std::endl;
    std::cerr << "   what = 5: voxel value [index]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GetImageInformation<2>( argc, argv );
     break;
   case 3:
     GetImageInformation<3>( argc, argv );
     break;
   case 4:
     GetImageInformation<4>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
