#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGeodesicActiveContourLevelSetImageFilter.h"

template <unsigned int ImageDimension>
int LevelSetSegmentation( unsigned int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::GeodesicActiveContourLevelSetImageFilter<ImageType, ImageType, PixelType> LevelSetFilterType;
  typename LevelSetFilterType::Pointer levelSetFilter = LevelSetFilterType::New();

  typedef  itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[2] );

  typename ImageType::Pointer initialLevelSet = reader1->GetOutput();
  initialLevelSet->Update();
  initialLevelSet->DisconnectPipeline();

  levelSetFilter->SetInitialImage( initialLevelSet );

  typename ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[3] );

  typename ImageType::Pointer edgePotentialMap = reader2->GetOutput();
  edgePotentialMap->Update();
  edgePotentialMap->DisconnectPipeline();

  levelSetFilter->SetFeatureImage( edgePotentialMap );

  float propagationScaling = 1.0;
  if( argc > 5 )
    {
    propagationScaling = atof( argv[5] );
    }

  float advectionScaling = 0.5;
  if( argc > 6 )
    {
    advectionScaling = atof( argv[6] );
    }

  float curvatureScaling = 0.1;
  if( argc > 7 )
    {
    curvatureScaling = atof( argv[7] );
    }

  unsigned int numberOfIterations = 200;
  if( argc > 8 )
    {
    numberOfIterations = atoi( argv[8] );
    }
  float maxRMSError = 0.03;
  if( argc > 9 )
    {
    maxRMSError = atof( argv[9] );
    }

  levelSetFilter->SetPropagationScaling( propagationScaling );
  levelSetFilter->SetAdvectionScaling( advectionScaling );
  levelSetFilter->SetCurvatureScaling( curvatureScaling );
  levelSetFilter->SetNumberOfIterations( numberOfIterations );
  levelSetFilter->SetMaximumRMSError( maxRMSError );

  try
    {
    levelSetFilter->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

		typedef  itk::ImageFileWriter<ImageType> WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetInput( levelSetFilter->GetOutput() );
		writer->SetFileName( argv[4] );
		writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension initialLevelSet";
    std::cerr << " edgePotentialMap outputImage";
    std::cerr << " [propagatingScaling=1] [advectionScaling=0.5] [curvatureScaling=0.1]";
    std::cerr << " [numberOfIterations=200] [maximumRMSError=0.03]" << std::endl;
    return 1;
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     LevelSetSegmentation<2>( argc, argv );
     break;
   case 3:
     LevelSetSegmentation<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

