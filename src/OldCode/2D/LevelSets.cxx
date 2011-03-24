#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGeodesicActiveContourLevelSetImageFilter.h"

#include "global.h"

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " initialLevelSetImage featureImage outputImage "
              << "[isoSurfaceValue] [advectionScaling] " 
              << "[curvatureScaling] [propagationScaling] {maxRMSerror] " 
              << "[numberOfIterations] " << std::endl;
    return EXIT_FAILURE;
    }
 
  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[1] );
  reader1->Update();

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  reader2->Update();

  typedef itk::GeodesicActiveContourLevelSetImageFilter
    <ImageType, ImageType, RealType> LevelSetFilterType;
  LevelSetFilterType::Pointer levelsets = LevelSetFilterType::New();
  levelsets->SetInitialImage( reader1->GetOutput() );
  levelsets->SetFeatureImage( reader2->GetOutput() );
  levelsets->SetIsoSurfaceValue( 0.0 );
  levelsets->SetAdvectionScaling( 0.5 );
  levelsets->SetCurvatureScaling( 0.1 );
  levelsets->SetPropagationScaling( 1.0 );
  levelsets->SetMaximumRMSError( 0.03 );
  levelsets->SetNumberOfIterations( 200 );

  if ( argc > 4 )
    {
    levelsets->SetIsoSurfaceValue( atof( argv[4] ) );
    }  
  if ( argc > 5 )
    {
    levelsets->SetAdvectionScaling( atof( argv[5] ) );
    }  
  if ( argc > 6 )
    {
    levelsets->SetCurvatureScaling( atof( argv[6] ) );
    }  
  if ( argc > 7 )
    {
    levelsets->SetPropagationScaling( atof( argv[7] ) );
    }  
  if ( argc > 8 )
    {
    levelsets->SetMaximumRMSError( atof( argv[8] ) );
    }  
  if ( argc > 9 )
    {
    levelsets->SetNumberOfIterations( atoi( argv[9] ) );
    }  
  levelsets->Update();

  typedef itk::ImageFileWriter<ImageType>  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( levelsets->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();

  return EXIT_SUCCESS;
}

