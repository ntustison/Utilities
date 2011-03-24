
#include "itkGaborFilterBankImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkTimeProbe.h"

#include <string>
#include <fstream.h>

template <unsigned int ImageDimension>
int GaborFilterBank( int argc, char *argv[] )
{
  itk::TimeProbe timer;
  timer.Start();

  typedef float RealType;

  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::GaborFilterBankImageFilter<RealImageType,RealImageType> GaborFilterType;
  typename GaborFilterType::Pointer gaborFilter = GaborFilterType::New();

  gaborFilter->SetInput( reader->GetOutput() );

  /**
    * Rotation[0, 1, 2] = [theta, psi, phi] = [rot_x, rot_y, rot_z]
    * where imaging plane is parallel to x-y plane
    */  

  typename GaborFilterType::ArrayType beginRotation;
  beginRotation.Fill( 0.0 );
  typename GaborFilterType::ArrayType endRotation;
  endRotation.Fill( 0.0 );
  typename GaborFilterType::UnsignedIntArrayType samples;
  samples.Fill( 0 );
  unsigned int numberOfSamples = 3;  
  if ( argc > 5 )
    {
    numberOfSamples = atoi( argv[4] );
    }  

  beginRotation[0] = 0.0;
  endRotation[0] = vnl_math::pi;
  beginRotation[1] = 0.0;
  endRotation[1] = vnl_math::pi;
  beginRotation[2] = 0.0;
  endRotation[2] = 0.0;
  samples[0] = numberOfSamples;
  samples[1] = numberOfSamples;
  samples[2] = numberOfSamples;

  gaborFilter->SetRotationAngleMinimum( beginRotation );
  gaborFilter->SetRotationAngleMaximum( endRotation );
  gaborFilter->SetNumberOfRotationAngleSteps( samples );

  RealType gaborSpacingFactor = 0.25;
  if ( argc > 6 )
    {
    gaborSpacingFactor = atof( argv[6] );
    } 

  gaborFilter->SetGaborSpacingMinimum( ( 1.0 - gaborSpacingFactor ) * atof( argv[5] ) );
  gaborFilter->SetGaborSpacingMaximum( ( 1.0 + gaborSpacingFactor ) * atof( argv[5] ) );

  if ( argc > 7 )
    {
    gaborFilter->SetNumberOfGaborSpacingSteps( atoi( argv[7] ) );
    }  
  else
    {
    gaborFilter->SetNumberOfGaborSpacingSteps( 3 );
    }  


  gaborFilter->Update();
  timer.Stop();
  std::cout << "Done. (" << timer.GetMeanTime() << ")" << std::endl;

  typedef itk::ImageFileWriter<RealImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( gaborFilter->GetOutput() );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage outputImage "
      << "[numberOfAngleSteps] [gaborSpacing] [gaborSpacingFactor] "
      << "[numberOfGaborSpacingSteps]" << std::endl;     
    exit( 0 );
    }   

  switch( atoi( argv[1] ) ) 
   {
//   case 2:
//     GaborFilterBank<2>( argc, argv );
//     break;
   case 3:
     GaborFilterBank<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

