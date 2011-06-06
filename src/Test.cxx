#include "itkImage.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkBSplineDeformableTransform.h"
#include "itkBSplineDeformableTransformInitializer.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkTimeProbe.h"
#include "itkVector.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

#define NEWCLASS

template<unsigned int Dimension, unsigned int SplineOrder>
int Test(int argc, char* argv[])
{
  typedef itk::BSplineDeformableTransform<float, Dimension, SplineOrder> TransformType;
  typename TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  typedef typename TransformType::ImageType ImageType;

//  typedef itk::ImageFileReader<ImageType> ReaderType;
//  typename ReaderType::Pointer reader = ReaderType::New();
//  reader->SetFileName( argv[4] );
//  reader->Update();
//
//  transform->SetTransformDomainOrigin( reader->GetOutput()->GetOrigin() );
//  transform->SetTransformDomainSpacing( reader->GetOutput()->GetSpacing() );
//  transform->SetTransformDomainSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
//  transform->SetTransformDomainDirection( reader->GetOutput()->GetDirection() );
//
//  typename BSplineTransformType::SizeType meshSize;
//
  typedef itk::Vector<float, Dimension> VectorType;
  typedef itk::Image<VectorType, Dimension> VectorImageType;

  typedef itk::ImageFileReader<VectorImageType> VReaderType;
  typename VReaderType::Pointer vreader = VReaderType::New();
  vreader->SetFileName( argv[3] );
  vreader->Update();

  typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ImageType>
  SelectorType;
  typename SelectorType::Pointer selector = SelectorType::New();
  selector->SetInput( vreader->GetOutput() );

  typename TransformType::CoefficientImageArray images;
  for( unsigned int d = 0; d < Dimension; d++ )
    {
    selector->SetIndex( d );
    selector->Update();

    images[d] = selector->GetOutput();
    images[d]->DisconnectPipeline();
    }

  transform->SetCoefficientImage( images );

#ifdef NEWCLASS
  std::cout << "PhysDim: " << transform->GetTransformDomainPhysicalDimensions() << std::endl;
  std::cout << "Mesh   : " << transform->GetTransformDomainMeshSize() << std::endl;
  std::cout << "Direction: " << transform->GetTransformDomainDirection() << std::endl;
#endif

  typename VectorImageType::Pointer vectorImage = VectorImageType::New();
  typename VectorImageType::PointType origin;
  origin.Fill( 0.0 );
  typename VectorImageType::SpacingType spacing;
  spacing.Fill( 1.0 );
  typename VectorImageType::RegionType::SizeType size;
  size.Fill( 256 );
  typename VectorImageType::DirectionType direction;
  direction.SetIdentity();

  vectorImage->SetOrigin( origin );
  vectorImage->SetSpacing( spacing );
  vectorImage->SetRegions( size );
  vectorImage->SetDirection( direction );
  vectorImage->Allocate();

  itk::TimeProbe timer;

  itk::ImageRegionIteratorWithIndex<VectorImageType> It( vectorImage,
    vectorImage->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    typename VectorImageType::PointType point;
    vectorImage->TransformIndexToPhysicalPoint( It.GetIndex(), point );

    typename TransformType::InputPointType inputPoint;
    inputPoint.CastFrom( point );

    timer.Start();
    typename TransformType::OutputPointType outputPoint =
      transform->TransformPoint( inputPoint );
    timer.Stop();

    It.Set( outputPoint - inputPoint );
    }

  std::cout << "Mean time: " << timer.GetMeanTime() << std::endl;

  typedef itk::ImageFileWriter<VectorImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( vectorImage );
  writer->Update();


//  typedef float PixelType;
//
//  typedef itk::Image< PixelType, Dimension > FixedVolumeType;
//
//  typedef itk::ImageFileReader< FixedVolumeType > FixedVolumeReaderType;
//  typename FixedVolumeReaderType::Pointer fixedVolumeReader = FixedVolumeReaderType::New();
//  fixedVolumeReader->SetFileName( argv[3] );
//  fixedVolumeReader->Update();
//
//  typename FixedVolumeType::IndexType index;
//  index.Fill( 0 );
//  typename FixedVolumeType::RegionType region = fixedVolumeReader->GetOutput()->GetLargestPossibleRegion();
//  region.SetIndex( index );
//
//  fixedVolumeReader->GetOutput()->SetLargestPossibleRegion( region );
//
//  std::vector<typename FixedVolumeType::PointType> corners;
//
//  typename FixedVolumeType::PointType point;
//
//  // 0: 0,0,0
//  fixedVolumeReader->GetOutput()->TransformIndexToPhysicalPoint( index, point );
//  corners.push_back( point );
//  // 1: 0,0,1
//  index[0] = 0;
//  index[1] = 0;
//  index[2] = fixedVolumeReader->GetOutput()->GetLargestPossibleRegion().GetSize()[2] - 1;
//  fixedVolumeReader->GetOutput()->TransformIndexToPhysicalPoint( index, point );
//  corners.push_back( point );
//  // 2: 0,1,0
//  index[0] = 0;
//  index[1] = fixedVolumeReader->GetOutput()->GetLargestPossibleRegion().GetSize()[1] - 1;
//  index[2] = 0;
//  fixedVolumeReader->GetOutput()->TransformIndexToPhysicalPoint( index, point );
//  corners.push_back( point );
//  // 3: 0,1,1
//  index[0] = 0;
//  index[1] = fixedVolumeReader->GetOutput()->GetLargestPossibleRegion().GetSize()[1] - 1;
//  index[2] = fixedVolumeReader->GetOutput()->GetLargestPossibleRegion().GetSize()[2] - 1;
//  fixedVolumeReader->GetOutput()->TransformIndexToPhysicalPoint( index, point );
//  corners.push_back( point );
//  // 4: 1,0,0
//  index[0] = fixedVolumeReader->GetOutput()->GetLargestPossibleRegion().GetSize()[0] - 1;
//  index[1] = 0;
//  index[2] = 0;
//  fixedVolumeReader->GetOutput()->TransformIndexToPhysicalPoint( index, point );
//  corners.push_back( point );
//  // 5: 1,0,1
//  index[0] = fixedVolumeReader->GetOutput()->GetLargestPossibleRegion().GetSize()[0] - 1;
//  index[1] = 0;
//  index[2] = fixedVolumeReader->GetOutput()->GetLargestPossibleRegion().GetSize()[2] - 1;
//  fixedVolumeReader->GetOutput()->TransformIndexToPhysicalPoint( index, point );
//  corners.push_back( point );
//  // 6: 1,1,0
//  index[0] = fixedVolumeReader->GetOutput()->GetLargestPossibleRegion().GetSize()[0] - 1;
//  index[1] = fixedVolumeReader->GetOutput()->GetLargestPossibleRegion().GetSize()[1] - 1;
//  index[2] = 0;
//  fixedVolumeReader->GetOutput()->TransformIndexToPhysicalPoint( index, point );
//  corners.push_back( point );
//  // 7: 1,1,1
//  index[0] = fixedVolumeReader->GetOutput()->GetLargestPossibleRegion().GetSize()[0] - 1;
//  index[1] = fixedVolumeReader->GetOutput()->GetLargestPossibleRegion().GetSize()[1] - 1;
//  index[2] = fixedVolumeReader->GetOutput()->GetLargestPossibleRegion().GetSize()[2] - 1;
//  fixedVolumeReader->GetOutput()->TransformIndexToPhysicalPoint( index, point );
//  corners.push_back( point );
//
//  std::ofstream strIm( "imageOutline.txt" );
//  strIm << "0 0 0 0" << std::endl;
//  strIm << corners[0][0] << " " << corners[0][1] << " " << corners[0][2] << " 4" << std::endl;
//  strIm << corners[4][0] << " " << corners[4][1] << " " << corners[4][2] << " 4" << std::endl;
//  strIm << corners[6][0] << " " << corners[6][1] << " " << corners[6][2] << " 4" << std::endl;
//  strIm << corners[2][0] << " " << corners[2][1] << " " << corners[2][2] << " 4" << std::endl;
//  strIm << corners[0][0] << " " << corners[0][1] << " " << corners[0][2] << " 4" << std::endl;
//
//  strIm << corners[1][0] << " " << corners[1][1] << " " << corners[1][2] << " 4" << std::endl;
//  strIm << corners[5][0] << " " << corners[5][1] << " " << corners[5][2] << " 4" << std::endl;
//  strIm << corners[7][0] << " " << corners[7][1] << " " << corners[7][2] << " 4" << std::endl;
//  strIm << corners[3][0] << " " << corners[3][1] << " " << corners[3][2] << " 4" << std::endl;
//  strIm << corners[1][0] << " " << corners[1][1] << " " << corners[1][2] << " 4" << std::endl;
//
//  strIm << corners[0][0] << " " << corners[0][1] << " " << corners[0][2] << " 4" << std::endl;
//  strIm << corners[4][0] << " " << corners[4][1] << " " << corners[4][2] << " 4" << std::endl;
//  strIm << corners[5][0] << " " << corners[5][1] << " " << corners[5][2] << " 4" << std::endl;
//  strIm << corners[1][0] << " " << corners[1][1] << " " << corners[1][2] << " 4" << std::endl;
//  strIm << corners[0][0] << " " << corners[0][1] << " " << corners[0][2] << " 4" << std::endl;
//
//  strIm << corners[2][0] << " " << corners[2][1] << " " << corners[2][2] << " 4" << std::endl;
//  strIm << corners[3][0] << " " << corners[3][1] << " " << corners[3][2] << " 4" << std::endl;
//  strIm << corners[7][0] << " " << corners[7][1] << " " << corners[7][2] << " 4" << std::endl;
//  strIm << corners[6][0] << " " << corners[6][1] << " " << corners[6][2] << " 4" << std::endl;
//
//  strIm << corners[0][0] << " " << corners[0][1] << " " << corners[0][2] << " 1" << std::endl;
//  strIm << corners[4][0] << " " << corners[4][1] << " " << corners[4][2] << " 1" << std::endl;
//
//  strIm << corners[0][0] << " " << corners[0][1] << " " << corners[0][2] << " 2" << std::endl;
//  strIm << corners[2][0] << " " << corners[2][1] << " " << corners[2][2] << " 2" << std::endl;
//
//  strIm << corners[0][0] << " " << corners[0][1] << " " << corners[0][2] << " 3" << std::endl;
//  strIm << corners[1][0] << " " << corners[1][1] << " " << corners[1][2] << " 3" << std::endl;
//
//  strIm << "0 0 0 0" << std::endl;
//  strIm.close();
//
//  typedef itk::BSplineDeformableTransform<float, Dimension, SplineOrder> BSplineTransformType;
//  typename BSplineTransformType::Pointer initialBSplineTransform =
//    BSplineTransformType::New();
//  initialBSplineTransform->SetIdentity();
//
//  typedef typename BSplineTransformType::RegionType TransformRegionType;
//  typedef typename TransformRegionType::SizeType    TransformSizeType;
//
//  typedef itk::BSplineDeformableTransformInitializer
//    <BSplineTransformType, FixedVolumeType> InitializerType;
//  typename InitializerType::Pointer transformInitializer = InitializerType::New();
//
//  transformInitializer->SetTransform( initialBSplineTransform );
//
//  typename InitializerType::TransformDomainMeshSizeType meshSize;
//  meshSize[0] = 1;
//  meshSize[1] = 2;
//  meshSize[2] = 3;
//
//  transformInitializer->SetTransformDomainMeshSize( meshSize );
//  transformInitializer->SetImage( fixedVolumeReader->GetOutput() );
//
//  try
//    {
//    transformInitializer->InitializeTransform();
//    }
//  catch( itk::ExceptionObject & excp )
//    {
//    std::cerr << "Exception thrown " << std::endl;
//    std::cerr << excp << std::endl;
//    return EXIT_FAILURE;
//    }
//
//  typename BSplineTransformType::ParametersType temp = initialBSplineTransform->GetParameters();
//  initialBSplineTransform->SetParameters( temp );
//
//  typename BSplineTransformType::ImagePointer coefficientImage =
//    initialBSplineTransform->GetCoefficientImage()[0];
//
//  typedef itk::ImageRegionIterator<typename BSplineTransformType::ImageType> IteratorType;
//
//  IteratorType coefItr( coefficientImage, coefficientImage->GetLargestPossibleRegion() );
//  coefItr.GoToBegin( );
//
//  std::cout << "0 0 0 0" << std::endl;
//
//  while(!coefItr.IsAtEnd())
//    {
//    typename BSplineTransformType::ImageType::PointType currentPoint;
//    coefficientImage->TransformIndexToPhysicalPoint(coefItr.GetIndex(),currentPoint);
//
//    std::cout << currentPoint[0] << " " << currentPoint[1] << " " << currentPoint[2] << " 1" << std::endl;
//
//    ++coefItr;
//    }
//
//  std::cout << "0 0 0 0" << std::endl;
//
  return 0;
}

int main( unsigned int argc, char * argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] <<
      " imageDimension splineOrder inputCPImage [outputImage]" << std::endl;
    return 0;
    }


  switch( atoi( argv[1] ) )
    {
    case 2:
      {
      switch( atoi( argv[2] ) )
        {
        case 1:
          {
          Test<2,1>( argc, argv );
          break;
          }
        case 2:
          {
          Test<2,2>( argc, argv );
          break;
          }
        case 3:
          {
          Test<2,3>( argc, argv );
          break;
          }
        case 4:
          {
          Test<2,4>( argc, argv );
          break;
          }
        default:
          {
          std::cerr << "Unrecognized spline order." << std::endl;
          break;
          }
        }
      break;
      }
    case 3:
      {
      switch( atoi( argv[2] ) )
        {
        case 1:
          {
          Test<3,1>( argc, argv );
          break;
          }
        case 2:
          {
          Test<3,2>( argc, argv );
          break;
          }
        case 3:
          {
          Test<3,3>( argc, argv );
          break;
          }
        case 4:
          {
          Test<3,4>( argc, argv );
          break;
          }
        default:
          {
          std::cerr << "Unrecognized spline order." << std::endl;
          break;
          }
        }
      break;
      }
    default:
      {
      std::cerr << "Unrecognized dimension." << std::endl;
      break;
      }
    }
}
