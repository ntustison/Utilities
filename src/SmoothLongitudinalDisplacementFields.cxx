#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImportImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkPointSet.h"

#include "Common.h"

template <unsigned int ImageDimension>
int SmoothLongitudinalDisplacementFields( int argc, char *argv[] )
{

  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;
  typedef itk::Image<VectorType, ImageDimension + 1> TimeVaryingDisplacementFieldType;


  std::vector<unsigned int> meshSize = ConvertVector<unsigned int>( std::string( argv[2] ) );
  if( meshSize.size() != ImageDimension + 1 )
    {
    std::cerr << "Mesh size needs to be specified as a vector of size ImageDimension + 1, "
              << "e.g. MxNxO,  where the last dimension is the number of mesh elements "
              << "in the temporal dimension." << std::endl;

    }

  unsigned int numberOfLevels = Convert<unsigned int>( std::string( argv[3] ) );

  // Read in the displacement fields

  std::vector<typename DisplacementFieldType::Pointer> inputFields;
  for( unsigned int n = 5; n < argc; n++ )
    {
    typedef itk::ImageFileReader<DisplacementFieldType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[2] );

    typename DisplacementFieldType::Pointer field = reader->GetOutput();
    field->Update();

    inputFields.push_back( field );
    }

  // accumulate the point set from all displacement fields

  typedef itk::PointSet<VectorType, ImageDimension + 1> PointSetType;
  typename PointSetType::Pointer pointSet = PointSetType::New();
  pointSet->Initialize();

  itk::SizeValueType numberOfPixels = 1;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    numberOfPixels *= inputFields[0]->GetBufferedRegion().GetSize()[d];
    }
  numberOfPixels *= inputFields.size();

  typename DisplacementFieldType::DirectionType identity;
  identity.SetIdentity();

  unsigned long count = 0;
  for( unsigned int n = 0; n < inputFields.size(); n++ )
    {
    typedef itk::ImportImageFilter<VectorType, ImageDimension> DisplacementFieldImporterType;
    typename DisplacementFieldImporterType::Pointer displacementFieldImporter = DisplacementFieldImporterType::New();
    displacementFieldImporter->SetImportPointer( inputFields[n]->GetBufferPointer(), numberOfPixels, false );
    displacementFieldImporter->SetRegion( inputFields[n]->GetBufferedRegion() );
    displacementFieldImporter->SetOrigin( inputFields[n]->GetOrigin() );
    displacementFieldImporter->SetSpacing( inputFields[n]->GetSpacing() );
    displacementFieldImporter->SetDirection( identity );
    displacementFieldImporter->Update();

    itk::ImageRegionIteratorWithIndex<DisplacementFieldType> It( displacementFieldImporter->GetOutput(),
      displacementFieldImporter->GetRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      typename DisplacementFieldType::PointType point;
      displacementFieldImporter->GetOutput()->TransformIndexToPhysicalPoint( It.GetIndex(), point );

      typename PointSetType::PointType spatioTemporalPoint;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        spatioTemporalPoint[d] = point[d];
        }
      spatioTemporalPoint[ImageDimension] = static_cast<float>( n );

      pointSet->SetPoint( count, spatioTemporalPoint );
      pointSet->SetPointData( count++, It.Get() );
      }
    }

  // set up the b-spline filter

  typename TimeVaryingDisplacementFieldType::PointType      timeVaryingDisplacementFieldOrigin;
  typename TimeVaryingDisplacementFieldType::SpacingType    timeVaryingDisplacementFieldSpacing;
  typename TimeVaryingDisplacementFieldType::SizeType       timeVaryingDisplacementFieldSize;
  typename TimeVaryingDisplacementFieldType::DirectionType  timeVaryingDisplacementFieldDirection;

  timeVaryingDisplacementFieldDirection.SetIdentity();

  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    timeVaryingDisplacementFieldOrigin[d] = inputFields[0]->GetOrigin()[d];
    timeVaryingDisplacementFieldSpacing[d] = inputFields[0]->GetSpacing()[d];
    timeVaryingDisplacementFieldSize[d] = inputFields[0]->GetBufferedRegion().GetSize()[d];
    }
  timeVaryingDisplacementFieldOrigin[ImageDimension] = 0.0;
  timeVaryingDisplacementFieldSpacing[ImageDimension] = 1.0;
  timeVaryingDisplacementFieldSize[ImageDimension] = inputFields.size();

  typedef itk::BSplineScatteredDataPointSetToImageFilter<PointSetType,
    TimeVaryingDisplacementFieldType> BSplineFilterType;

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
  bspliner->SetInput( pointSet );
  bspliner->SetOrigin( timeVaryingDisplacementFieldOrigin );
  bspliner->SetSpacing( timeVaryingDisplacementFieldSpacing );
  bspliner->SetSize( timeVaryingDisplacementFieldSize );
  bspliner->SetDirection( timeVaryingDisplacementFieldDirection );
  bspliner->SetNumberOfLevels( numberOfLevels );
  bspliner->SetSplineOrder( 3 );

  typename BSplineFilterType::ArrayType numberOfControlPoints;
  for( unsigned int d = 0; d <= ImageDimension; d++ )
    {
    numberOfControlPoints[d] = meshSize[d] + 3;
    }
  bspliner->SetNumberOfControlPoints( numberOfControlPoints );

  typename BSplineFilterType::ArrayType closeDimensions;
  closeDimensions.Fill( false );
  closeDimensions[ImageDimension] = true;
  bspliner->SetCloseDimension( closeDimensions );

  bspliner->Update();

  // We set up the size of the sampled b-spline filter so that we can simply
  // extract each slice


  typename TimeVaryingDisplacementFieldType::RegionType region;
  typename TimeVaryingDisplacementFieldType::RegionType::SizeType size =
    bspliner->GetOutput()->GetRequestedRegion().GetSize();
  typename TimeVaryingDisplacementFieldType::IndexType index =
    bspliner->GetOutput()->GetRequestedRegion().GetIndex();

  size[ImageDimension] = 0;

  for( unsigned int n = 0; n < inputFields.size(); n++ )
    {
    index[ImageDimension] = n;

    region.SetIndex( index );
    region.SetSize( size );

    typedef itk::ExtractImageFilter<TimeVaryingDisplacementFieldType, DisplacementFieldType> ExtracterType;
    typename ExtracterType::Pointer extracter = ExtracterType::New();
    extracter->SetInput( bspliner->GetOutput() );
    extracter->SetExtractionRegion( region );
    extracter->SetDirectionCollapseToIdentity();
    extracter->Update();

    extracter->GetOutput()->SetDirection( inputFields[n]->GetDirection() );

    std::ostringstream whichInput;
    whichInput << n;
    std::string outputFileName = std::string( argv[4] ) + whichInput.str() + std::string( ".nii.gz" );

    typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outputFileName.c_str() );
    writer->SetInput( extracter->GetOutput() );
    writer->Update();
    }

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cerr << argv[0] << " imageDimension meshSize numberOfLevels outputPrefix inputFiles" << std::endl;
    exit( 0 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     SmoothLongitudinalDisplacementFields<2>( argc, argv );
     break;
   case 3:
     SmoothLongitudinalDisplacementFields<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

