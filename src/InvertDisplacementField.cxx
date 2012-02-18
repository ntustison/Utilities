#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkInvertDisplacementFieldImageFilter.h"
#include "itkVector.h"

template <unsigned int ImageDimension>
int Invert( int argc, char *argv[] )
{
  typedef itk::Vector<float, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;

  typedef itk::ImageFileReader<DisplacementFieldType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  unsigned int numberOfIterations = 20;
  float meanTolerance = 0.1;
  float maxTolerance = 0.001;
  if( argc > 4 )
    {
    numberOfIterations = atoi( argv[4] );
    }
  if( argc > 5 )
    {
    meanTolerance = atof( argv[5] );
    }
  if( argc > 6 )
    {
    maxTolerance = atoi( argv[6] );
    }

  typedef itk::InvertDisplacementFieldImageFilter<DisplacementFieldType> InverterType;
  typename InverterType::Pointer inverter = InverterType::New();
  inverter->SetInput( reader->GetOutput() );
  inverter->SetMaximumNumberOfIterations( numberOfIterations );
  inverter->SetMeanErrorToleranceThreshold( meanTolerance );
  inverter->SetMaxErrorToleranceThreshold( maxTolerance );
  inverter->Update();

  typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( inverter->GetOutput() );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension inputField outputField [numberOfIterations=20] [meanTolerance=0.1] [maxTolerance=0.001]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     Invert<2>( argc, argv );
     break;
   case 3:
     Invert<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

