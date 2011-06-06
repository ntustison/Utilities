#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkConstNeighborhoodIterator.h"


#include "topological_numbers.h"

template <unsigned int ImageDimension>
int CreateTopologicalNumberFromBinaryImage( int argc, char *argv[] )
{

  typedef unsigned char PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  /**
   * Read in input
   */
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  unsigned int connectivity = ( argc > 4 ) ? atoi( argv[4] ) : 1; 
  bool calculateInverse = ( argc > 5 ) ? 
    static_cast<bool>( atoi( argv[5] ) ) : false; 

  typename ImageType::Pointer map = ImageType::New();
  map->SetOrigin( reader->GetOutput()->GetOrigin() );
  map->SetDirection( reader->GetOutput()->GetDirection() );
  map->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  map->SetSpacing( reader->GetOutput()->GetSpacing() );
  map->Allocate();
  map->FillBuffer( 0 );

  typedef itk::ConstNeighborhoodIterator<ImageType> NeighborhoodIteratorType;
  typedef typename NeighborhoodIteratorType::RadiusType   RadiusType;
  RadiusType radius;
  radius.Fill( 1 );

  typedef typename itk::NeighborhoodAlgorithm
    ::ImageBoundaryFacesCalculator<ImageType> FaceCalculatorType;
  FaceCalculatorType faceCalculator;

  typename FaceCalculatorType::FaceListType faceList = faceCalculator(
    reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion(), radius );
  typename FaceCalculatorType::FaceListType::iterator fit = faceList.begin();

//  for( fit = faceList.begin(); fit != faceList.end(); ++fit )
//    {
    NeighborhoodIteratorType inIt( radius, reader->GetOutput(), *fit );
    itk::ImageRegionIterator<ImageType> outIt( map, *fit );

    for( inIt.GoToBegin(), outIt.GoToBegin(); !inIt.IsAtEnd();
      ++inIt, ++outIt )
      {
      NBH neighborhood;
      if( ImageDimension == 2 )
        {
        for( unsigned int i = 0; i < 3; i++ )
          {
          for( unsigned int j = 0; j < 3; j++ )
            {
            neighborhood[i][0][j] = neighborhood[i][2][j] = 0;
            neighborhood[i][1][j] = static_cast<unsigned char>( 
              inIt.GetPixel( 3*i+j ) != 0 );
            }
          } 
        }
      else if( ImageDimension == 3 )
        {
        for( unsigned int i = 0; i < 3; i++ )
          {
          for( unsigned int j = 0; j < 3; j++ )
            {
            for( unsigned int k = 0; k < 3; k++ )
              {
              neighborhood[i][j][k] = static_cast<unsigned char>( 
                inIt.GetPixel( 9*i + 3*j + k ) != 0 );
              }
            }
          } 
        }
      NBH dnbh;
      if( calculateInverse )
        {
        for( unsigned int i = 0; i < 3; i++ )
          {
          for( unsigned int j = 0; j < 3; j++ )
            {
            for( unsigned int k = 0; k < 3; k++ )
              {
              neighborhood[i][j][k] = !neighborhood[i][j][k];
              }
            }
          } 
        outIt.Set( checkTn( &neighborhood, &dnbh, 
          associatedConnectivity( connectivity ) ) );
        }
      else
        {  
        outIt.Set( checkTn( &neighborhood, &dnbh, connectivity ) );
        }
      }
//    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( map );
  writer->Update();


  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputBinaryImage "
      << "outputTopologyNumberMap [connectivity] [calculateInverse]" << std::endl;
    std::cout << "  connectivity: " << std::endl;
    std::cout << "    1.  (6,18) -> default" << std::endl; 
    std::cout << "    2.  (18,6)" << std::endl; 
    std::cout << "    3.  (6,26)" << std::endl; 
    std::cout << "    4.  (26,6)" << std::endl; 
    exit( 1 );
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     CreateTopologicalNumberFromBinaryImage<2>( argc, argv );
     break;
   case 3:
     CreateTopologicalNumberFromBinaryImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


