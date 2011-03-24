#include "itkContinuousIndex.h"
#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"
#include "itkPointSet.h"
#include "itkImageFileReader.h"
#include "itkVectorImageFileReader.h"
#include "itkVectorLinearInterpolateImageFunction.h"

template <unsigned int ImageDimension>
int WarpPoints( int argc, char *argv[] )
{
  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;

  typedef itk::PointSet<unsigned long, ImageDimension> PointSetType;
  typedef typename PointSetType::PointType PointType;
  typedef itk::LabeledPointSetFileReader<PointSetType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  try
    {
    typedef itk::VectorImageFileReader<RealImageType,
      DeformationFieldType> FieldReaderType;
    typename FieldReaderType::Pointer fieldreader = FieldReaderType::New();
    fieldreader->SetFileName( argv[3] );
    fieldreader->SetUseAvantsNamingConvention( true );
    fieldreader->Update();

    typedef itk::VectorLinearInterpolateImageFunction
      <DeformationFieldType, RealType> DeformationFieldInterpolatorType;
    typename DeformationFieldInterpolatorType::Pointer interpolator
      = DeformationFieldInterpolatorType::New();
    interpolator->SetInputImage( fieldreader->GetOutput() );

    typename PointSetType::PointsContainerIterator It
      = reader->GetOutput()->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator ItL
      = reader->GetOutput()->GetPointData()->Begin();

    typename PointSetType::Pointer output = PointSetType::New();
    output->Initialize();

    unsigned long idx = 0;
    while( It != reader->GetOutput()->GetPoints()->End() )
      {
      typename DeformationFieldInterpolatorType::PointType point;
      point.CastFrom( It.Value() );

      if( interpolator->IsInsideBuffer( point ) )
        {
        typename DeformationFieldInterpolatorType::OutputType vec
          = interpolator->Evaluate( point );
        PointType newPoint = It.Value() + vec;
        output->SetPoint( idx, newPoint );
        output->SetPointData( idx, ItL.Value() );
        idx++;
        }
      ++It;
      ++ItL;
      }
    typedef itk::LabeledPointSetFileWriter<PointSetType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[4] );
    writer->SetInput( output );
    writer->Update();
    }
  catch(...)
    {
    typedef itk::ImageFileReader<RealImageType> ImageReaderType;
    typename ImageReaderType::Pointer imagereader = ImageReaderType::New();
    imagereader->SetFileName( argv[3] );
    imagereader->Update();

    typename PointSetType::PointsContainerIterator It
      = reader->GetOutput()->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator ItL
      = reader->GetOutput()->GetPointData()->Begin();

    typename PointSetType::Pointer output = PointSetType::New();
    output->Initialize();

    unsigned long idx = 0;
    while( It != reader->GetOutput()->GetPoints()->End() )
      {
      typename RealImageType::PointType point;
      point.CastFrom( It.Value() );

      itk::ContinuousIndex<double, ImageDimension> cidx;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        cidx[d] = point[d];
        }

      if( imagereader->GetOutput()->GetLargestPossibleRegion().IsInside( cidx ) )
        {
        imagereader->GetOutput()->TransformContinuousIndexToPhysicalPoint( cidx, point );
        output->SetPoint( idx, point );
        output->SetPointData( idx, ItL.Value() );
        idx++;
        }
      ++It;
      ++ItL;
      }
    typedef itk::LabeledPointSetFileWriter<PointSetType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[4] );
    writer->SetInput( output );
    writer->Update();
    }
  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << argv[0] << " ImageDimension inputPointSet deformationField outputPointSet" << std::endl;
    std::cout << "     Note: if deformation field is an image, we assume the input point set is in voxel coordinates and conversion is to physical coordinates of the image." << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     WarpPoints<2>( argc, argv );
     break;
   case 3:
     WarpPoints<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

