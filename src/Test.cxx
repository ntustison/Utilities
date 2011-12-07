#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTimeVaryingBSplineVelocityFieldIntegrationImageFilter.h"
#include "itkTimeVaryingBSplineVelocityFieldTransform.h"
#include "itkTimeVaryingBSplineVelocityFieldTransformParametersAdaptor.h"
#include "itkVector.h"

template <unsigned int ImageDimension>
int Test( unsigned int argc, char *argv[] )
{
  typedef itk::Vector<double, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;
  typedef itk::Image<VectorType, ImageDimension+1> TimeVaryingVelocityFieldType;

  typedef itk::ImageFileReader<TimeVaryingVelocityFieldType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::TimeVaryingBSplineVelocityFieldIntegrationImageFilter
    <TimeVaryingVelocityFieldType, DeformationFieldType> IntegratorType;
  typename IntegratorType::Pointer integrator = IntegratorType::New();
  integrator->SetInput( reader->GetOutput() );
  integrator->SetSplineOrder( 3 );
  integrator->SetLowerTimeBound( atof( argv[5] ) );
  integrator->SetUpperTimeBound( atof( argv[6] ) );

//  integrator->SetNumberOfIntegrationSteps( atoi( argv[4] ) );
//  integrator->Update();
//
//  typedef itk::ImageFileWriter<DeformationFieldType> WriterType;
//  typename WriterType::Pointer writer = WriterType::New();
//  writer->SetFileName( argv[3] );
//  writer->SetInput( integrator->GetOutput() );
//  writer->Update();
//
//
//  TimeVaryingVelocityFieldType::PointType origin;
//  origin.Fill( 0.0 );
//  TimeVaryingVelocityFieldType::SpacingType spacing;
//  spacing.Fill( 2.0 );
//  TimeVaryingVelocityFieldType::SizeType size;
//  size.Fill( 25 );
//  VectorType velocity;
//  velocity.Fill( 0.1 );
//
//  TimeVaryingVelocityFieldType::Pointer timeVaryingVelocityField =
//    TimeVaryingVelocityFieldType::New();
//  timeVaryingVelocityField->SetOrigin( origin );
//  timeVaryingVelocityField->SetSpacing( spacing );
//  timeVaryingVelocityField->SetRegions( size );
//  timeVaryingVelocityField->Allocate();
//  timeVaryingVelocityField->FillBuffer( velocity );
//
//  typedef itk::TimeVaryingVelocityFieldIntegrationImageFilter
//    <TimeVaryingVelocityFieldType, DeformationFieldType> IntegratorType;
//
//  IntegratorType::Pointer integrator = IntegratorType::New();
//  integrator->SetInput( timeVaryingVelocityField );
//  integrator->SetLowerTimeBound( 0.3 );
//  integrator->SetUpperTimeBound( 0.75 );
//  integrator->SetNumberOfIntegrationSteps( 10 );
//  integrator->Update();
//
//  integrator->Print( std::cout, 3 );
//
//  DeformationFieldType::IndexType index;
//  index.Fill( 0 );
//  VectorType displacement;
//
//  // This integration should result in a constant image of value
//  // 0.75 * 0.1 - 0.3 * 0.1 = 0.045 with ~epsilon deviation
//  // due to numerical computations
//  displacement = integrator->GetOutput()->GetPixel( index );
//  if( vnl_math_abs( displacement[0] - 0.045 ) > 0.0001 )
//    {
//    std::cerr << "Failed to produce the correct forward integration."
//      << std::endl;
//    return EXIT_FAILURE;
//    }
//
//  IntegratorType::Pointer inverseIntegrator = IntegratorType::New();
//  inverseIntegrator->SetInput( timeVaryingVelocityField );
//  inverseIntegrator->SetLowerTimeBound( 1.0 );
//  inverseIntegrator->SetUpperTimeBound( 0.0 );
//  inverseIntegrator->SetNumberOfIntegrationSteps( 10 );
//  inverseIntegrator->Update();
//
//  // This integration should result in a constant image of value
//  // -( 0.1 * 1.0 - ( 0.1 * 0.0 ) ) = -0.1 with ~epsilon deviation
//  // due to numerical computations
//  displacement = inverseIntegrator->GetOutput()->GetPixel( index );
//  if( vnl_math_abs( displacement[0] + 0.1 ) > 0.0001 )
//    {
//    std::cerr << "Failed to produce the correct inverse integration."
//      << std::endl;
//    return EXIT_FAILURE;
//    }
//
//  // Now test the transform
//
//  typedef itk::TimeVaryingVelocityFieldTransform<double, 3> TransformType;
//  TransformType::Pointer transform = TransformType::New() ;
//  transform->SetLowerTimeBound( 0.0 );
//  transform->SetUpperTimeBound( 1.0 );
//  transform->SetTimeVaryingVelocityField( timeVaryingVelocityField );
//
//  TransformType::InputPointType point;
//  point.Fill( 1.3 );
//
//  TransformType::OutputPointType transformedPoint =
//    transform->TransformPoint( point );
//
//  point += velocity;
//  if( point.EuclideanDistanceTo( transformedPoint ) > 0.001 )
//    {
//    std::cerr << "Failed to produce the expected transformed point."
//      << std::endl;
//    return EXIT_FAILURE;
//    }
//  point -= velocity;
//
//  TransformType::InputPointType point2;
//  point2.CastFrom( transformedPoint );
//  transformedPoint = transform->GetInverseTransform()->TransformPoint( point2 );
//
//  if( point.EuclideanDistanceTo( transformedPoint ) > 0.001 )
//    {
//    std::cerr << "Failed to produce the expected inverse transformed point."
//      << std::endl;
//    return EXIT_FAILURE;
//    }
//
//  transform->Print( std::cout, 3 );
//
  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension velocityField integratedField numberOfIntegrationPoints" << std::endl;
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
   case 4:
     Test<4>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


