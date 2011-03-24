#include "itkBSplineKernelTransform.h"
#include "itkLandmarkFileReader.h"
#include "itkPointSet.h"
#include "itkVector.h"
#include "itkVectorImageFileWriter.h"

#include "itkImageRegionIteratorWithIndex.h"

#include "global.h"

int main( unsigned int argc, char *argv[] )
{
  if ( argc < 3*ImageDimension + 7 )
    {
    std::cout << "Usage: " << argv[0] << "  inputMovingLandmarkFile inputFixedLandmarkFile outputField "
              << "origin[0] ... origin[dimension-1] " << std::endl
              << "spacing[0] ... spacing[dimension-1] " << std::endl
              << "size[0] ... size[dimension-1] " << std::endl
              << "nlevels bsplineOrder ncps" << std::endl
              << "(optional: outputControlPointImage)" << std::endl
              << std::endl;
    exit( 1 );
    }

//  typedef float RealType;
//
//  typedef itk::Vector<RealType, ImageDimension> VectorType;
//  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;
//
//
//  typedef itk::BSplineKernelTransform<RealType, ImageDimension> TransformType;
//  TransformType::Pointer transform = TransformType::New();
//
//  typedef TransformType::PointSetType PointSetType;
//  typedef itk::LandmarkFileReader<PointSetType> ReaderType;
//
//  ReaderType::Pointer movingReader = ReaderType::New();
//  movingReader->SetFileName( argv[1] );
//  movingReader->Update();
//
//  ReaderType::Pointer fixedReader = ReaderType::New();
//  fixedReader->SetFileName( argv[2] );
//  fixedReader->Update();
//
//  transform->SetTargetLandmarks( fixedReader->GetOutput() );
//  transform->SetSourceLandmarks( movingReader->GetOutput() );
//
//  TransformType::ArrayType ncps;
//  ncps.Fill( atoi( argv[4+ImageDimension*3+2] ) );
//
//  transform->SetNumberOfLevels( atoi( argv[4+ImageDimension*3] ) );
//  transform->SetSplineOrder( atoi( argv[4+ImageDimension*3+1] ) );
//  transform->SetNumberOfControlPoints( ncps );
//
//  DeformationFieldType::PointType origin;
//  DeformationFieldType::SpacingType spacing;
//  DeformationFieldType::SizeType size;
//
//  for ( unsigned int i = 0; i < ImageDimension; i++ )
//    {
//    origin[i] = atof( argv[4+i] ); 
//    spacing[i] = atof( argv[4+ImageDimension+i] );
//    size[i] = atoi( argv[4+2*ImageDimension+i] );
//    } 
//  
//  DeformationFieldType::Pointer field = DeformationFieldType::New();
//  field->SetOrigin( origin );
//  field->SetRegions( size );
//  field->SetSpacing( spacing );
//  field->Update();
//
//  transform->SetOrigin( origin );
//  transform->SetSize( size );
//  transform->SetSpacing( spacing );
//
//  itk::ImageRegionIteratorWithIndex<DeformationFieldType> It(
//    field, field->GetLargestPossibleRegion() );
//  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
//    {
//    DeformationFieldType::PointType point;
//    field->TransformIndexToPhysicalPoint( It.GetIndex(), point );
//    TransformType::InputPointType input;
//    input.CastFrom( point );
//    TransformType::OutputPointType output = transform->TransformPoint( input );
//    It.Set( output - input );
//    }
//
//  typedef itk::Image<RealType, ImageDimension> RealImageType;
//  typedef itk::VectorImageFileWriter<DeformationFieldType, RealImageType> WriterType;
//  WriterType::Pointer writer = WriterType::New();
//  writer->SetFileName( argv[3] );
//  writer->SetUseAvantsNamingConvention( false );
//  writer->SetInput( field );
//  writer->Update();

  return 0;
}
