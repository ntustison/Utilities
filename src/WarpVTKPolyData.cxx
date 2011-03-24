#include "itkVTKPolyDataReader.h"
#include "itkVTKPolyDataWriter.h"
#include "itkMesh.h"
#include "itkImageFileReader.h"
#include "itkVectorLinearInterpolateImageFunction.h"

template <unsigned int ImageDimension>
int WarpVTKPolyData( int argc, char *argv[] )
{
  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;

  typedef itk::Mesh<RealType, ImageDimension> MeshType;
  typedef itk::VTKPolyDataReader<MeshType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::ImageFileReader<DeformationFieldType> FieldReaderType;
  typename FieldReaderType::Pointer fieldreader = FieldReaderType::New();
  fieldreader->SetFileName( argv[3] );
  fieldreader->Update();

  typedef itk::VectorLinearInterpolateImageFunction
    <DeformationFieldType, RealType> DeformationFieldInterpolatorType;
  typename DeformationFieldInterpolatorType::Pointer interpolator
    = DeformationFieldInterpolatorType::New();
  interpolator->SetInputImage( fieldreader->GetOutput() );

  typename MeshType::PointsContainerIterator It
    = reader->GetOutput()->GetPoints()->Begin();

  while( It != reader->GetOutput()->GetPoints()->End() )
    {
    typename DeformationFieldInterpolatorType::PointType point;
    point.CastFrom( It.Value() );

    if( interpolator->IsInsideBuffer( point ) )
      {
      typename DeformationFieldInterpolatorType::OutputType vec
        = interpolator->Evaluate( point );
      typename MeshType::PointType newPoint = It.Value() + vec;
      reader->GetOutput()->SetPoint( It.Index(), newPoint );
      }
    ++It;
    }

  typedef itk::VTKPolyDataWriter<MeshType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( reader->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << argv[0] << " ImageDimension inputPointSet "
     << "deformationField outputPointSet" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     WarpVTKPolyData<2>( argc, argv );
     break;
   case 3:
     WarpVTKPolyData<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
