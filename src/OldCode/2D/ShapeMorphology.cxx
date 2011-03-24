#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryReinhardtMorphologicalImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " labelImage outputImage foregroundValue" << std::endl;
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

  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
  LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( argv[1] );
  labelReader->Update();

  typedef itk::BinaryBallStructuringElement< 
                      LabelImageType::PixelType,
                      ImageDimension>             StructuringElementType;

  typedef itk::BinaryReinhardtMorphologicalImageFilter<
    LabelImageType, LabelImageType, StructuringElementType>  FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( labelReader->GetOutput() );

  if ( argc > 3 )
    {   
    filter->SetForegroundValue( atoi( argv[3] ) );
    }

  filter->SetEmploySaltAndPepperRepair( false );
  if ( argc > 4 )
    {
    filter->SetEmploySaltAndPepperRepair( atoi( argv[4] ) );
    }
  if ( argc > 5 )
    {
    filter->SetSaltAndPepperMinimumSizeInPixels( atoi( argv[5] ) );
    }

  filter->SetEmployMinimumDiameterFilter( false );
  if ( argc > 6 )
    {
    filter->SetEmployMinimumDiameterFilter( atoi( argv[6] ) );
    }
  if ( argc > 7 )
    {
    filter->SetMinimumDiameterStructuringElementRadius( atoi( argv[7] ) );
    } 
   
  filter->SetEmployUnwantedCavityDeletion( false );
  if ( argc > 8 )
    {
    filter->SetEmployUnwantedCavityDeletion( atoi( argv[8] ) );
    }

  filter->SetEmployMinimumSizeFilter( false );
  if ( argc > 9 )
    {
    filter->SetEmployMinimumSizeFilter( atoi( argv[9] ) );
    }
  if ( argc > 10 )
    {
    filter->SetMinimumSizeStructuringElementRadius( atoi( argv[10] ) );
    }

  filter->SetEmployMaximumDiameterFilter( false );
  if ( argc > 11 )
    {
    filter->SetEmployMaximumDiameterFilter( atoi( argv[11] ) );
    }
  if ( argc > 12 )
    {
    filter->SetMaximumDiameterStructuringElementRadius( atoi( argv[12] ) );
    }

  filter->SetEmployConnectivityFilter( false );
  if ( argc > 13 )
    {
    filter->SetEmployConnectivityFilter( atoi( argv[13] ) );
    }
  if ( argc > 14 )
    {
    filter->SetNumberOfConnectedComponents( atoi( argv[14] ) );
    } 

  filter->SetEmployBoundarySmoother( false );
  if ( argc > 15 )
    {
    filter->SetEmployBoundarySmoother( atoi( argv[15] ) );
    }
  if ( argc > 16 )
    {
    filter->SetBoundarySmootherStructuringElementRadius( atoi( argv[16] ) );
    }
   
  filter->SetEmployUnclassifiedPixelProcessing( false );
  if ( argc > 17 )
    {
    filter->SetEmployUnclassifiedPixelProcessing( atoi( argv[17] ) );
    }

  filter->Update();

  std::cout << filter << std::endl;

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();
 
  return EXIT_SUCCESS;

}