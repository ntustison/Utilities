#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryReinhardtMorphologicalImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

template <unsigned int ImageDimension>
int ShapeMorphology( int argc, char *argv[] )
{

  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
  typename LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( argv[2] );
  labelReader->Update();

  typedef itk::BinaryBallStructuringElement< 
                      typename LabelImageType::PixelType,
                      ImageDimension>             StructuringElementType;

  typedef itk::BinaryReinhardtMorphologicalImageFilter<
    LabelImageType, LabelImageType, StructuringElementType>  FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( labelReader->GetOutput() );

  if ( argc > 4 )
    {   
    filter->SetForegroundValue( atoi( argv[4] ) );
    }

  filter->SetEmploySaltAndPepperRepair( false );
  if ( argc > 5 )
    {
    filter->SetEmploySaltAndPepperRepair( atoi( argv[5] ) );
    }
  if ( argc > 6 )
    {
    filter->SetSaltAndPepperMinimumSizeInPixels( atoi( argv[6] ) );
    }

  filter->SetEmployMinimumDiameterFilter( false );
  if ( argc > 7 )
    {
    filter->SetEmployMinimumDiameterFilter( atoi( argv[7] ) );
    }
  if ( argc > 8 )
    {
    filter->SetMinimumDiameterStructuringElementRadius( atoi( argv[8] ) );
    } 
   
  filter->SetEmployUnwantedCavityDeletion( false );
  if ( argc > 9 )
    {
    filter->SetEmployUnwantedCavityDeletion( atoi( argv[9] ) );
    }

  filter->SetEmployMinimumSizeFilter( false );
  if ( argc > 10 )
    {
    filter->SetEmployMinimumSizeFilter( atoi( argv[10] ) );
    }
  if ( argc > 11 )
    {
    filter->SetMinimumSizeStructuringElementRadius( atoi( argv[11] ) );
    }

  filter->SetEmployMaximumDiameterFilter( false );
  if ( argc > 12 )
    {
    filter->SetEmployMaximumDiameterFilter( atoi( argv[12] ) );
    }
  if ( argc > 13 )
    {
    filter->SetMaximumDiameterStructuringElementRadius( atoi( argv[13] ) );
    }

  filter->SetEmployConnectivityFilter( false );
  if ( argc > 14 )
    {
    filter->SetEmployConnectivityFilter( atoi( argv[14] ) );
    }
  if ( argc > 15 )
    {
    filter->SetNumberOfConnectedComponents( atoi( argv[15] ) );
    } 

  filter->SetEmployBoundarySmoother( false );
  if ( argc > 16 )
    {
    filter->SetEmployBoundarySmoother( atoi( argv[16] ) );
    }
  if ( argc > 17 )
    {
    filter->SetBoundarySmootherStructuringElementRadius( atoi( argv[17] ) );
    }
   
  filter->SetEmployUnclassifiedPixelProcessing( false );
  if ( argc > 18 )
    {
    filter->SetEmployUnclassifiedPixelProcessing( atoi( argv[18] ) );
    }

  filter->Update();

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();
 
  return EXIT_SUCCESS;

}


int main( int argc, char *argv[] )
{

  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension labelImage "
              <<  " outputImage foregroundValue" << std::endl;
    std::cout << "   [SaltAndPepperRepair] [MinimumSize]" << std::endl;
    std::cout << "   [MinimumDiameter] [SE radius]" << std::endl;
    std::cout << "   [UnwantedCavityDeletion] " << std::endl;
    std::cout << "   [MinimumSize] [SE radius]" << std::endl;
    std::cout << "   [MaximumDiameter] [SE radius]" << std::endl;
    std::cout << "   [Connectivity] [NumberOfConnectedComponents]" << std::endl;
    std::cout << "   [BoundarySmoothing] [SE radius]" << std::endl;
    std::cout << "   [UnclassifiedPixelProcessing] " << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     ShapeMorphology<2>( argc, argv );
     break;
   case 3:
     ShapeMorphology<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}