#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodIterator.h"
#include "itkTimeProbe.h"
#include "itkTranslationTransform.h"
#include "itkVectorImageFileWriter.h"

#include "itkOptMeanSquaresImageToImageMetric.h"
#include "itkMeanSquareRegistrationFunction.h"

template <unsigned int ImageDimension>
int Test( unsigned int argc, char *argv[] )
{
  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> ImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;

  typedef itk::ImageFileReader<ImageType> ReaderType;

  typename ReaderType::Pointer fixedReader = ReaderType::New();
  fixedReader->SetFileName( argv[2] );
  fixedReader->Update();

  typename ReaderType::Pointer movingReader = ReaderType::New();
  movingReader->SetFileName( argv[3] );
  movingReader->Update();

  typename ImageType::Pointer msqImage = ImageType::New();
  msqImage->SetOrigin( fixedReader->GetOutput()->GetOrigin() );
  msqImage->SetSpacing( fixedReader->GetOutput()->GetSpacing() );
  msqImage->SetRegions( fixedReader->GetOutput()->GetLargestPossibleRegion() );
  msqImage->SetDirection( fixedReader->GetOutput()->GetDirection() );
  msqImage->Allocate();
  msqImage->FillBuffer( 0.0 );

  VectorType vector( 0.0 );
  typename DeformationFieldType::Pointer gradientField =
    DeformationFieldType::New();
  gradientField->SetOrigin( fixedReader->GetOutput()->GetOrigin() );
  gradientField->SetSpacing( fixedReader->GetOutput()->GetSpacing() );
  gradientField->SetRegions( fixedReader->GetOutput()->GetLargestPossibleRegion() );
  gradientField->SetDirection( fixedReader->GetOutput()->GetDirection() );
  gradientField->Allocate();
  gradientField->FillBuffer( vector );

  /**
   * Test the registration function
   */
//  {
//  itk::TimeProbe timer;
//
//  timer.Start();
//
//  typename DeformationFieldType::Pointer deformationField =
//    DeformationFieldType::New();
//  deformationField->SetOrigin( fixedReader->GetOutput()->GetOrigin() );
//  deformationField->SetSpacing( fixedReader->GetOutput()->GetSpacing() );
//  deformationField->SetRegions( fixedReader->GetOutput()->GetLargestPossibleRegion() );
//  deformationField->SetDirection( fixedReader->GetOutput()->GetDirection() );
//  deformationField->Allocate();
//  deformationField->FillBuffer( vector );
//
//  typedef itk::MeanSquareRegistrationFunction<ImageType, ImageType,
//    DeformationFieldType> RegistrationFunctionType;
//  typename RegistrationFunctionType::Pointer msq = RegistrationFunctionType::New();
//  msq->SetDeformationField( deformationField );
//  msq->SetFixedImage( fixedReader->GetOutput() );
//  msq->SetMovingImage( movingReader->GetOutput() );
//  msq->InitializeIteration();
//
//  typedef typename itk::NeighborhoodAlgorithm
//    ::ImageBoundaryFacesCalculator<DeformationFieldType> FaceCalculatorType;
//  FaceCalculatorType faceCalculator;
//  typedef itk::NeighborhoodIterator<DeformationFieldType> NeighborhoodIteratorType;
//
//  typename FaceCalculatorType::FaceListType faceList = faceCalculator(
//    deformationField, deformationField->GetLargestPossibleRegion(),
//    msq->GetRadius() );
//  typename FaceCalculatorType::FaceListType::iterator fit;
//
//  for( fit = faceList.begin(); fit != faceList.end(); ++fit )
//    {
//    NeighborhoodIteratorType It( msq->GetRadius(),
//      deformationField, *fit );
//    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
//      {
//      msq->SetEnergy( 0.0 );
//      VectorType grad = msq->ComputeUpdate( It, NULL );
//      msqImage->SetPixel( It.GetIndex(), msq->GetEnergy() );
//      gradientField->SetPixel( It.GetIndex(), grad );
//      }
//    }
//
//  timer.Stop();
//  std::cout << "Elapsed time: " << timer.GetMeanTime() << std::endl;
//  }

  /**
   * Test the image to image msq
   */
  {
  itk::TimeProbe timer;

  timer.Start();

  typedef itk::TranslationTransform<double, ImageDimension> TransformType;
  typename TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( movingReader->GetOutput() );

  typedef itk::MeanSquaresImageToImageMetric<ImageType, ImageType> MetricType;
  typename MetricType::Pointer msq = MetricType::New();
  msq->SetFixedImage( fixedReader->GetOutput() );
  msq->SetMovingImage( movingReader->GetOutput() );
  msq->SetTransform( transform.GetPointer() );
  msq->SetInterpolator( interpolator.GetPointer() );
  msq->SetUseFixedImageIndexes( true );
  msq->SetFixedImageRegion( fixedReader->GetOutput()->GetLargestPossibleRegion() );
  msq->SetNumberOfThreads( 1 );

  typedef typename itk::NeighborhoodAlgorithm
    ::ImageBoundaryFacesCalculator<ImageType> FaceCalculatorType;
  FaceCalculatorType faceCalculator;
  typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );

  typename FaceCalculatorType::FaceListType faceList = faceCalculator(
    msqImage, msqImage->GetLargestPossibleRegion(), radius );
  typename FaceCalculatorType::FaceListType::iterator fit;

  for( fit = faceList.begin(); fit != faceList.end(); ++fit )
    {
    itk::ImageRegionIteratorWithIndex<ImageType> It( msqImage, *fit );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      typename MetricType::DerivativeType derivative( ImageDimension );
      typename MetricType::FixedImageIndexContainer container;
      container.push_back( It.GetIndex() );
      msq->SetFixedImageIndexes( container );
      msq->Initialize();
      typename MetricType::MeasureType value;// = msq->GetValue( transform->GetParameters() );
      msq->GetValueAndDerivative( transform->GetParameters(), value, derivative );
      msqImage->SetPixel( It.GetIndex(), value );

      VectorType grad;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        grad[d] = derivative[d];
        }
      gradientField->SetPixel( It.GetIndex(), grad );
      }
    }

  timer.Stop();
  std::cout << "Elapsed time: " << timer.GetMeanTime() << std::endl;
  }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( msqImage );
  writer->SetFileName( argv[4] );
  writer->Update();

  typedef itk::VectorImageFileWriter<DeformationFieldType, ImageType>
    FieldWriterType;
  typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
  fieldWriter->SetInput( gradientField );
  fieldWriter->SetFileName( argv[4] );
  fieldWriter->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
   if ( argc < 4 )
     {
     std::cout << argv[0] << " imageDimension fixedImage movingImage outputImage" << std::endl;
     exit( 1 );
     }

   switch( atoi( argv[1] ) )
    {
    case 2:
      Test<2>( argc, argv );
      break;
    case 3:
      Test<3>( argc, argv );
      break;
    default:
       std::cerr << "Unsupported dimension" << std::endl;
       exit( EXIT_FAILURE );
   }
}

