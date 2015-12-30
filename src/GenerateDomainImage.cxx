#include "itkBoundingBox.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkVector.h"
#include "itkVectorContainer.h"

#include "vnl/algo/vnl_svd.h"

#include "Common.h"

int GenerateDomainImage( unsigned int argc, char *argv[] )
{
  const unsigned int ImageDimension = 3;

  typedef float                                   RealType;
  typedef itk::Image<RealType, ImageDimension>    ImageType;
  typedef itk::Point<RealType, ImageDimension>    PointType;
  typedef itk::Vector<double, ImageDimension>   VectorType;

  typedef itk::VectorContainer<int, PointType> VectorContainerType;
  VectorContainerType::Pointer points = VectorContainerType::New();
  points->Initialize();

  ImageType::SpacingType spacing;

  std::vector<RealType> imageSpacing = ConvertVector<RealType>( std::string( argv[2] ) );
  if ( imageSpacing.size() == 1 )
    {
    spacing.Fill( imageSpacing[0] );
    }
  else if ( imageSpacing.size() == ImageDimension )
    {
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      spacing[d] = imageSpacing[d];
      }
    }
  else
    {
    std::cerr << "Invalid image spacing format." << std::endl;
    return EXIT_FAILURE;
    }


  // Find the point of intersection between the two LA planes and one of the SA planes

  VectorType zVector;
  zVector[0] = 0.0;
  zVector[1] = 0.0;
  zVector[2] = 1.0;

  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer readerLA1 = ReaderType::New();
  readerLA1->SetFileName( argv[argc-2] );
  readerLA1->Update();

  ImageType::DirectionType directionLA1 = readerLA1->GetOutput()->GetDirection();
  ImageType::PointType originLA1 = readerLA1->GetOutput()->GetOrigin();
  VectorType normalLA1 = directionLA1 * zVector;
  RealType kLA1 = normalLA1[0] * originLA1[0] + normalLA1[1] * originLA1[1] + normalLA1[2] * originLA1[2];

  ReaderType::Pointer readerLA2 = ReaderType::New();
  readerLA2->SetFileName( argv[argc-1] );
  readerLA2->Update();

  ImageType::DirectionType directionLA2 = readerLA2->GetOutput()->GetDirection();
  ImageType::PointType originLA2 = readerLA2->GetOutput()->GetOrigin();
  VectorType normalLA2 = directionLA2 * zVector;
  RealType kLA2 = normalLA2[0] * originLA2[0] + normalLA2[1] * originLA2[1] + normalLA2[2] * originLA2[2];

  ImageType::PointType intersectionPoint;
  intersectionPoint.Fill( 0.0 );

  unsigned int N = 0;
  for( unsigned int n = 3; n < argc-2; n++ )
    {
    ReaderType::Pointer readerSA = ReaderType::New();
    readerSA->SetFileName( argv[n] );
    readerSA->Update();

    ImageType::DirectionType directionSA = readerSA->GetOutput()->GetDirection();
    ImageType::PointType originSA = readerSA->GetOutput()->GetOrigin();
    VectorType normalSA = directionSA * zVector;
    RealType kSA = normalSA[0] * originSA[0] + normalSA[1] * originSA[1] + normalSA[2] * originSA[2];

    vnl_matrix<double> A( ImageDimension, ImageDimension );
    A.set_row( 0, normalSA.GetVnlVector() );
    A.set_row( 1, normalLA1.GetVnlVector() );
    A.set_row( 2, normalLA2.GetVnlVector() );

    vnl_vector<double> b( ImageDimension );
    b[0] = kSA;
    b[1] = kLA1;
    b[2] = kLA2;

    vnl_vector<double> x = vnl_svd<double>( A ).solve( b );

    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      intersectionPoint[d] = ( x[d] + static_cast<double>( N ) * intersectionPoint[d] ) / static_cast<double>( N + 1 );
      }

    if( ( intersectionPoint.GetVectorFromOrigin() ).GetNorm() > 1e8 )
      {
      std::cerr << "Planes are not parallel.  Make sure that you have specified a SA";
      std::cerr << "slice and two LA slices as the last in the image list." << std::endl;
      return EXIT_FAILURE;
      }
    N++;
    }

  ImageType::DirectionType direction;
  direction.SetIdentity();
  ImageType::DirectionType inverseDirection;
  inverseDirection.SetIdentity();

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[3] );

  ImageType::Pointer inputImage = reader->GetOutput();
  inputImage->Update();
  inputImage->DisconnectPipeline();

  direction = inputImage->GetDirection();
  inverseDirection = direction.GetInverse();

  intersectionPoint = inverseDirection * intersectionPoint;

  // Get the bounding box along the long axis

  N = 0;
  for( unsigned int n = 3; n < argc-2; n++ )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[n] );

    ImageType::Pointer inputImage = reader->GetOutput();
    inputImage->Update();
    inputImage->DisconnectPipeline();

    itk::ImageRegionConstIteratorWithIndex<ImageType> It( inputImage, inputImage->GetRequestedRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      ImageType::PointType imagePoint;
      inputImage->TransformIndexToPhysicalPoint( It.GetIndex(), imagePoint );

      imagePoint = inverseDirection * imagePoint;

      points->InsertElement( N++, imagePoint );
      }
    }

  typedef  itk::BoundingBox<unsigned int, ImageDimension, RealType, VectorContainerType> BoundingBoxType;
  BoundingBoxType::Pointer boundingBoxLA = BoundingBoxType::New();
  boundingBoxLA->SetPoints( points );
  boundingBoxLA->ComputeBoundingBox();
  BoundingBoxType::BoundsArrayType boundsLA = boundingBoxLA->GetBounds();

  // Assume heart is less than 120 mm in diameter

  RealType manhattanRadius = 60.0;

  N = 0;
  points->Initialize();
  for( unsigned int n = 3; n < argc; n++ )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[n] );

    ImageType::Pointer inputImage = reader->GetOutput();
    inputImage->Update();
    inputImage->DisconnectPipeline();

    itk::ImageRegionConstIteratorWithIndex<ImageType> It( inputImage, inputImage->GetRequestedRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      ImageType::PointType imagePoint;
      inputImage->TransformIndexToPhysicalPoint( It.GetIndex(), imagePoint );

      imagePoint = inverseDirection * imagePoint;

      if( ( imagePoint[0] < intersectionPoint[0] - manhattanRadius ) ||
          ( imagePoint[0] > intersectionPoint[0] + manhattanRadius ) ||
          ( imagePoint[1] < intersectionPoint[1] - manhattanRadius ) ||
          ( imagePoint[1] > intersectionPoint[1] + manhattanRadius )
        )
        {
        continue;
        }

      if( n >= argc - 2 )
        {
        if( imagePoint[2] < boundsLA[4] || imagePoint[2] > boundsLA[5] )
          {
          continue;
          }
        }

      points->InsertElement( N++, imagePoint );
      }
    }

  BoundingBoxType::Pointer boundingBox = BoundingBoxType::New();
  boundingBox->SetPoints( points );
  boundingBox->ComputeBoundingBox();
  BoundingBoxType::BoundsArrayType bounds = boundingBox->GetBounds();

  ImageType::PointType origin;

//   RealType minDistance = itk::NumericTraits<RealType>::max();
//   PointType minPoint;
//
//   const BoundingBoxType::PointsContainer* corners = boundingBox->GetCorners();
//   VectorContainerType::ConstIterator corner = corners->Begin();
//   while( corner != corners->End() )
//     {
//     RealType distance = originVolume.EuclideanDistanceTo( corner->Value() );
//     if( distance < minDistance )
//       {
//       minPoint.CastFrom( corner->Value() );
//       }
//     corner++;
//     }

  origin.CastFrom( boundingBox->GetMinimum() );
//   origin.CastFrom( minPoint );

  ImageType::PointType center;
  center.CastFrom( boundingBox->GetCenter() );
  ImageType::SizeType size;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    size[d] = static_cast<unsigned int>( ( bounds[2*d + 1] - bounds[2*d] ) / spacing[d] ) + 1;
    RealType epsilon = 0.5 * spacing[d] * static_cast<RealType>( size[d] - 1 ) - ( center[d] - origin[d] );
    origin[d] -= epsilon;
    }

  origin = direction * origin;

  ImageType::Pointer domainImage = ImageType::New();
  domainImage->SetOrigin( origin );
  domainImage->SetSpacing( spacing );
  domainImage->SetDirection( direction );
  domainImage->SetRegions( size );
  domainImage->Allocate();
  domainImage->FillBuffer( 0.0 );

  ImageType::IndexType intersectionIndex;
  intersectionPoint = direction * intersectionPoint;
  domainImage->TransformPhysicalPointToIndex( intersectionPoint, intersectionIndex );
  domainImage->SetPixel( intersectionIndex, 1.0 );

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[1] );
  writer->SetInput( domainImage );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cerr << argv[0] << " outputImage outputImageSpacing inputImage1* inputImage2 ... inputImageN" << std::endl;
    std::cerr << "*Note that we use inputImage1 to determine orientation.  Also, the last two images are LA." << std::endl;
    return EXIT_FAILURE;
    }

  return GenerateDomainImage( argc, argv );
}

