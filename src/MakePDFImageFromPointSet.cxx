#include "itkGaussianProbabilityDensityFunction.h"
#include "itkManifoldParzenWindowsPointSetFunction.h"
#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"
#include "itkPointSet.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"

template <unsigned int ImageDimension>
int MakePDFImageFromPointSet( int argc, char *argv[] )
{

  typedef float RealType;
  typedef float PixelType;

  typedef itk::PointSet<long, ImageDimension> PointSetType;
  typedef typename PointSetType::PointType PointType;

  typedef itk::LabeledPointSetFileReader<PointSetType> PointSetReaderType;
  typename PointSetReaderType::Pointer psreader = PointSetReaderType::New();
  psreader->SetFileName( argv[2] );
  psreader->Update();

  typedef itk::ManifoldParzenWindowsPointSetFunction<PointSetType> ParzenFilterType;
  typename ParzenFilterType::Pointer parzen = ParzenFilterType::New();
  parzen->SetBucketSize( 4 );
  parzen->SetRegularizationSigma( atof( argv[3] ) );
  parzen->SetKernelSigma( atof( argv[4] ) );
  parzen->SetEvaluationKNeighborhood( atoi( argv[5] ) );
  parzen->SetCovarianceKNeighborhood( atoi( argv[6] ) );
  parzen->SetUseAnisotropicCovariances( atoi( argv[7] ) );
  parzen->SetInputPointSet( psreader->GetOutput() );

  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[8] );
  reader->Update();

  itk::ImageRegionIteratorWithIndex<RealImageType> It( reader->GetOutput(),
    reader->GetOutput()->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    typename RealImageType::PointType point;
    reader->GetOutput()->TransformIndexToPhysicalPoint( It.GetIndex(), point );
    PointType pt;
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      pt[d] = point[d];
      }
    It.Set( parzen->Evaluate( pt ) );
    }

  typedef itk::ImageFileWriter<RealImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[9] );
  writer->SetInput( reader->GetOutput() );
  writer->Update();

  if( argc > 10 )
    {
    typename PointSetType::Pointer samples = PointSetType::New();
    samples->Initialize();

    unsigned long numberOfSamples = 100;
    if( argc > 11 )
      {
      std::string av = std::string( argv[11] );
      if( av[av.length()-1] == '%' )
        {
        av = av.substr( 0, av.length()-1 );
        float percentage = atof( av.c_str() );
        numberOfSamples = static_cast<unsigned long>( percentage *
          psreader->GetOutput()->GetNumberOfPoints() / 100.0 );
        }
      else
        {
        numberOfSamples = atoi( argv[11] );
        }
      }

    for( unsigned int n = 0; n < numberOfSamples; n++ )
      {
      samples->SetPoint( n, parzen->GenerateRandomSample() );
      samples->SetPointData( n, 1 );
      }

    typedef itk::LabeledPointSetFileWriter<PointSetType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[10] );
    writer->SetInput( samples );
    writer->Update();

    }

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 10 )
    {
    std::cout << argv[0] << " imageDimension inputPointSet "
              << "pointSetSigma kernelSigma evaluationKNeighborhood "
              << "covarianceKNeighborhood "
              << "useAnisotropicCovariances domainImage outputImage "
              << "[samplePointSetFile] [numberOfSamples=100 or in percentage, e.g. 20%]"
              << std::endl;
    exit( 0 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     MakePDFImageFromPointSet<2>( argc, argv );
     break;
   case 3:
     MakePDFImageFromPointSet<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


