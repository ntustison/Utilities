#include "itkBSplineControlPointImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorLinearInterpolateImageFunction.h"

#include "itkImageRegionIterator.h"
#include "itkVector.h"

template <unsigned int ImageDimension>
int GenerateVectorFieldFromControlPointLattice( int argc, char *argv[] )
{
  typedef float PixelType;
  typedef float RealType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::Vector<PixelType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> VectorImageType;

  typedef itk::ImageFileReader<VectorImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  typename ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName( argv[4] );
  imageReader->Update();

  typename VectorImageType::RegionType::SizeType size
    = imageReader->GetOutput()->GetLargestPossibleRegion().GetSize();
  typename VectorImageType::PointType origin
    = imageReader->GetOutput()->GetOrigin();
  typename VectorImageType::SpacingType spacing =
    imageReader->GetOutput()->GetSpacing();

  RealType factor = 0.0;
  if ( argc > 6 )
    {
    factor = atof( argv[6] );
    }
  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    origin[d] -= factor * spacing[d];
    size[d] += static_cast<unsigned int>( 2*factor );
    }

  typedef itk::BSplineControlPointImageFilter
    <VectorImageType, VectorImageType> BSplineControlPointsFilterType;
  typename BSplineControlPointsFilterType::Pointer bspliner
    = BSplineControlPointsFilterType::New();

  bspliner->SetSplineOrder( atoi( argv[5] ) );
  bspliner->SetInput( reader->GetOutput() );
  bspliner->SetDirection( reader->GetOutput()->GetDirection() );
  bspliner->SetOrigin( origin );
  bspliner->SetSize( size );
  bspliner->SetSpacing( spacing );
  bspliner->Update();

  if( factor == 0.0 )
    {
    typedef itk::ImageFileWriter<VectorImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[3] );
    writer->SetInput( bspliner->GetOutput() );
    writer->Update();
    }
  else
    {
    typename VectorImageType::Pointer output = VectorImageType::New();
    output->SetOrigin( imageReader->GetOutput()->GetOrigin() );
    output->SetSpacing( imageReader->GetOutput()->GetSpacing() );
    output->SetRegions( imageReader->GetOutput()->GetLargestPossibleRegion() );
    output->Allocate();

    typedef itk::VectorLinearInterpolateImageFunction<VectorImageType, RealType>
      InterpolatorType;
    typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
    interpolator->SetInputImage( bspliner->GetOutput() );

    itk::ImageRegionIteratorWithIndex<VectorImageType> It( output,
      output->GetLargestPossibleRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      typename VectorImageType::PointType point;
      output->TransformIndexToPhysicalPoint( It.GetIndex(), point );
      typename InterpolatorType::PointType ipoint;
      ipoint.CastFrom( point );

      It.Set( interpolator->Evaluate( ipoint ) );
      }
    typedef itk::ImageFileWriter<VectorImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[3] );
    writer->SetInput( output );
    writer->Update();
    }


  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cout << "Usage: " << argv[0]
      << " imageDimension inputCPField outputField referenceImage"
      << " bsplineOrder [expansionFactor]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GenerateVectorFieldFromControlPointLattice<2>( argc, argv );
     break;
   case 3:
     GenerateVectorFieldFromControlPointLattice<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

