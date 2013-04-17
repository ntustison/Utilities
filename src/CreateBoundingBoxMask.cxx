#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBoxSpatialObject.h"
#include "itkBresenhamLine.h"
#include "itkGroupSpatialObject.h"
#include "itkSpatialObjectToImageFilter.h"

#include "itkVector.h"

#include <string>
#include <vector>

template<class TValue>
TValue Convert( std::string optionString )
{
  TValue value;
  std::istringstream iss( optionString );
  iss >> value;
  return value;
}

template<class TValue>
std::vector<TValue> ConvertVector( std::string optionString )
{
  std::vector<TValue> values;
  std::string::size_type crosspos = optionString.find( 'x', 0 );

  if ( crosspos == std::string::npos )
   {
   values.push_back( Convert<TValue>( optionString ) );
   }
  else
   {
   std::string element = optionString.substr( 0, crosspos ) ;
   TValue value;
   std::istringstream iss( element );
   iss >> value;
   values.push_back( value );
   while ( crosspos != std::string::npos )
     {
     std::string::size_type crossposfrom = crosspos;
     crosspos = optionString.find( 'x', crossposfrom + 1 );
     if ( crosspos == std::string::npos )
       {
       element = optionString.substr( crossposfrom + 1, optionString.length() );
       }
     else
       {
       element = optionString.substr( crossposfrom + 1, crosspos ) ;
       }
     std::istringstream iss2( element );
     iss2 >> value;
     values.push_back( value );
     }
   }
  return values;
}

std::vector<float> CrossProduct( std::vector<float> vector1, std::vector<float> vector2 )
{
  std::vector<float> crossProduct;
  crossProduct.push_back( vector1[1] * vector2[2] - vector1[2] * vector2[1] );
  crossProduct.push_back( -vector1[0] * vector2[2] + vector1[2] * vector2[0] );
  crossProduct.push_back( vector1[0] * vector2[1] - vector1[1] * vector2[0] );

  return crossProduct;
}

float DotProduct( std::vector<float> vector1, std::vector<float> vector2 )
{
  return ( vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2] );
}

std::vector<float> Subtract( std::vector<float> minuend, std::vector<float> subtrahend )
{
  std::vector<float> difference( 3 );
  difference[0] = minuend[0] - subtrahend[0];
  difference[1] = minuend[1] - subtrahend[1];
  difference[2] = minuend[2] - subtrahend[2];

  return difference;
}

std::vector<float> Add( std::vector<float> addend1, std::vector<float> addend2 )
{
  std::vector<float> sum( 3 );
  sum[0] = addend1[0] + addend2[0];
  sum[1] = addend1[1] + addend2[1];
  sum[2] = addend1[2] + addend2[2];

  return sum;
}

bool IsInsideBox( std::vector<float> point,
                  std::vector<float> corner0,
                  std::vector<float> rowVector,
                  std::vector<float> columnVector,
                  std::vector<float> normalVector,
                  std::vector<float> dimensions )
{
  std::vector<float> corner1( 3 );
  corner1[0] = corner0[0] + dimensions[0] * rowVector[0];
  corner1[1] = corner0[1] + dimensions[0] * rowVector[1];
  corner1[2] = corner0[2] + dimensions[0] * rowVector[2];

  std::vector<float> corner3( 3 );
  corner3[0] = corner0[0] + dimensions[1] * columnVector[0];
  corner3[1] = corner0[1] + dimensions[1] * columnVector[1];
  corner3[2] = corner0[2] + dimensions[1] * columnVector[2];

  std::vector<float> corner4( 3 );
  corner4[0] = corner0[0] + dimensions[2] * normalVector[0];
  corner4[1] = corner0[1] + dimensions[2] * normalVector[1];
  corner4[2] = corner0[2] + dimensions[2] * normalVector[2];

  if( DotProduct( Subtract( point, corner0 ), normalVector ) < 0 || DotProduct( Subtract( point, corner4 ), normalVector ) > 0 )
    {
    return false;
    }
  if( DotProduct( Subtract( point, corner0 ), rowVector ) < 0 || DotProduct( Subtract( point, corner1 ), rowVector ) > 0 )
    {
    return false;
    }
  if( DotProduct( Subtract( point, corner0 ), columnVector ) < 0 || DotProduct( Subtract( point, corner3 ), columnVector ) > 0 )
    {
    return false;
    }

  return true;
}



int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cout << argv[0] << " inputImage outputImage position dimensions row column" << std::endl;
    exit( 1 );
    }

  const unsigned int ImageDimension = 3;

  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::GroupSpatialObject<ImageDimension> SceneType;
  typedef itk::SpatialObjectToImageFilter<SceneType, ImageType> SpatialObjectToImageFilterType;
  SceneType::Pointer scene = SceneType::New();

  typedef itk::BoxSpatialObject<ImageDimension> BoxType;
  typedef BoxType::TransformType       TransformType;
  typedef TransformType::MatrixType    MatrixType;
  typedef TransformType::OffsetType    OffsetType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  std::vector<float> positionVector = ConvertVector<float>( std::string( argv[3] ) );
  std::vector<float> dimensions = ConvertVector<float>( std::string( argv[4] ) );
  std::vector<float> rowVector = ConvertVector<float>( std::string( argv[5] ) );
  std::vector<float> columnVector = ConvertVector<float>( std::string( argv[6] ) );

  // set up coordinate system of box
  // rowVector x columnVector = normalVector

  std::vector<float> crossProduct = CrossProduct( rowVector, columnVector );

  ImageType::Pointer output = ImageType::New();
  output->CopyInformation( reader->GetOutput() );
  output->SetRegions( reader->GetOutput()->GetRequestedRegion() );
  output->Allocate();
  output->FillBuffer( 0.0 );

  itk::ImageRegionIteratorWithIndex<ImageType> It( output, output->GetRequestedRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    ImageType::PointType point;
    output->TransformIndexToPhysicalPoint( It.GetIndex(), point );

    std::vector<float> pointVector( 3 );
    pointVector[0] = point[0] + 0.5 * ( 0 * dimensions[0] * rowVector[0] + 0 * dimensions[1] * columnVector[0] + dimensions[2] * crossProduct[0] );
    pointVector[1] = point[1] + 0.5 * ( 0 * dimensions[0] * rowVector[1] + 0 * dimensions[1] * columnVector[1] + dimensions[2] * crossProduct[1] );
    pointVector[2] = point[2] + 0.5 * ( 0 * dimensions[0] * rowVector[2] + 0 * dimensions[1] * columnVector[2] + dimensions[2] * crossProduct[2] );

    if( IsInsideBox( pointVector, positionVector, rowVector, columnVector, crossProduct, dimensions ) )
      {
      It.Set( 1 );
      }
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( output );
  writer->Update();

  return 0;
}















//   MatrixType matrixBoxCoordinateSystem;
//
//   matrixBoxCoordinateSystem[0][0] = rowVector[0];
//   matrixBoxCoordinateSystem[1][0] = rowVector[1];
//   matrixBoxCoordinateSystem[2][0] = rowVector[2];
//
//   matrixBoxCoordinateSystem[0][1] = columnVector[0];
//   matrixBoxCoordinateSystem[1][1] = columnVector[1];
//   matrixBoxCoordinateSystem[2][1] = columnVector[2];
//
//   matrixBoxCoordinateSystem[0][2] = crossProduct[0];
//   matrixBoxCoordinateSystem[1][2] = crossProduct[1];
//   matrixBoxCoordinateSystem[2][2] = crossProduct[2];
//
//   // Set up the box spatial object
//
//   BoxType::Pointer box = BoxType::New();
//
//   BoxType::SizeType dimensions2;
//   dimensions2[0] = dimensions[0];
//   dimensions2[1] = dimensions[1];
//   dimensions2[2] = dimensions[2];
//
//   BoxType::TransformType::MatrixType matrixInPlane;
//   matrixInPlane.SetIdentity();
//
//   float theta = atof( argv[7] );
//
//   if( theta != 0 )
//     {
//     float u = crossProduct[0];
//     float v = crossProduct[1];
//     float w = crossProduct[2];
//
//     float u2 = vnl_math_sqr( u );
//     float v2 = vnl_math_sqr( v );
//     float w2 = vnl_math_sqr( w );
//
//     float cos_theta = vcl_cos( theta );
//     float sin_theta = vcl_sin( theta );
//
//     matrixInPlane[0][0] = cos_theta + u2 * ( 1.0 - cos_theta );
//     matrixInPlane[1][0] = u * v * ( 1.0 - cos_theta ) + w * sin_theta;
//     matrixInPlane[2][0] = u * w * ( 1.0 - cos_theta ) - v * sin_theta;
//
//     matrixInPlane[0][1] = u * v * ( 1.0 - cos_theta ) - w * sin_theta;
//     matrixInPlane[1][1] = cos_theta + v2 * ( 1.0 - cos_theta );
//     matrixInPlane[2][1] = v * w * ( 1.0 - cos_theta ) + u * sin_theta;
//
//     matrixInPlane[0][2] = u * w * ( 1.0 - cos_theta ) + v * sin_theta;
//     matrixInPlane[1][2] = v * w * ( 1.0 - cos_theta ) - u * sin_theta;
//     matrixInPlane[2][2] = cos_theta + w2 * ( 1.0 - cos_theta );
//     }
//
//   OffsetType offset;
//   for( unsigned int i = 0; i < ImageDimension; i++ )
//     {
//     offset[i] = 0.0;
//     for( unsigned int j = 0; j < ImageDimension; j++ )
//       {
//       offset[i] -= matrixBoxCoordinateSystem[i][j] * ( 0.5 * dimensions[j] );
//       }
//     }
//
//   TransformType::Pointer transform = TransformType::New();
//   transform->SetMatrix( matrixBoxCoordinateSystem );
//   transform->SetOffset( offset );
//
//   OffsetType offset2;
//   for( unsigned int i = 0; i < ImageDimension; i++ )
//     {
//     offset2[i] = positionVector[i];
//     for( unsigned int j = 0; j < ImageDimension; j++ )
//       {
//       offset2[i] += matrixInPlane[i][j] * ( 0.5 * dimensions[j] );
//       }
//     }
//
//   TransformType::Pointer transform2 = TransformType::New();
//   transform2->SetMatrix( matrixInPlane );
//   transform2->SetOffset( offset2 );
//
//   transform->Compose( transform2, false );
//
//   box->SetSize( dimensions2 );
//   box->SetObjectToParentTransform( transform );
//   box->ComputeObjectToWorldTransform();
//
//   std::cout << "Object to parent transform" << std::endl;
//   box->GetObjectToParentTransform()->Print( std::cout, 3 );
//
//   std::cout << "Object to world transform" << std::endl;
//   box->GetObjectToWorldTransform()->Print( std::cout, 3 );
//
//   scene->AddSpatialObject( box );
//
//   SpatialObjectToImageFilterType::Pointer filter = SpatialObjectToImageFilterType::New();
//   filter->SetInput( scene );
//   filter->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
//   filter->SetOrigin( reader->GetOutput()->GetOrigin() );
//   filter->SetSpacing( reader->GetOutput()->GetSpacing() );
//   filter->SetDirection( reader->GetOutput()->GetDirection() );
//   filter->SetInsideValue( 1 );
//   filter->SetOutsideValue( 0 );
//
//   ImageType::Pointer output = filter->GetOutput();
//   output->Update();
//   output->DisconnectPipeline();
//
//   if( argc > 8 )
//     {
//     typedef itk::BresenhamLine<ImageDimension> LinerType;
//     typedef LinerType::IndexType IndexType;
//     LinerType liner;
//
//     ImageType::IndexType sourceIndex;
//     ImageType::IndexType targetIndex;
//
//     ImageType::PointType positionPoint;
//
//     positionPoint[0] = positionVector[0];
//     positionPoint[1] = positionVector[1];
//     positionPoint[2] = positionVector[2];
//
//     output->TransformPhysicalPointToIndex( positionPoint, sourceIndex );
//
//     {  // do row vector -> 3
//     unsigned int label = 3;
//
//     ImageType::PointType point;
//     point[0] = positionVector[0] + dimensions[0] * rowVector[0];
//     point[1] = positionVector[1] + dimensions[0] * rowVector[1];
//     point[2] = positionVector[2] + dimensions[0] * rowVector[2];
//
//     output->TransformPhysicalPointToIndex( point, targetIndex );
//
//     std::cout << "row: " << targetIndex << std::endl;
//
//     LinerType::IndexArray indices = liner.BuildLine( sourceIndex, targetIndex );
//
//     LinerType::IndexArray::const_iterator it;
//     for( it = indices.begin(); it != indices.end(); it++ )
//       {
//       output->SetPixel( *it, label );
//       }
//     }
//
//     {  // do column vector -> 4
//     unsigned int label = 4;
//
//     ImageType::PointType point;
//     point[0] = positionVector[0] + dimensions[1] * columnVector[0];
//     point[1] = positionVector[1] + dimensions[1] * columnVector[1];
//     point[2] = positionVector[2] + dimensions[1] * columnVector[2];
//
//     output->TransformPhysicalPointToIndex( point, targetIndex );
//
//     std::cout << "col: " << targetIndex << std::endl;
//
//     LinerType::IndexArray indices = liner.BuildLine( sourceIndex, targetIndex );
//
//     LinerType::IndexArray::const_iterator it;
//     for( it = indices.begin(); it != indices.end(); it++ )
//       {
//       output->SetPixel( *it, label );
//       }
//     }
//
//     {  // do normal vector -> 5
//     unsigned int label = 5;
//
//     ImageType::PointType point;
//     point[0] = positionVector[0] + dimensions[2] * crossProduct[0];
//     point[1] = positionVector[1] + dimensions[2] * crossProduct[1];
//     point[2] = positionVector[2] + dimensions[2] * crossProduct[2];
//
//     output->TransformPhysicalPointToIndex( point, targetIndex );
//
//     std::cout << "normal: " << targetIndex << std::endl;
//
//     LinerType::IndexArray indices = liner.BuildLine( sourceIndex, targetIndex );
//
//     LinerType::IndexArray::const_iterator it;
//     for( it = indices.begin(); it != indices.end(); it++ )
//       {
//       output->SetPixel( *it, label );
//       }
//     }
//
//     output->TransformPhysicalPointToIndex( positionPoint, sourceIndex );
//     output->SetPixel( sourceIndex, 2 );
//
//     typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension> StructuringElementType;
//     StructuringElementType  element;
//     element.SetRadius( 2 );
//     element.CreateStructuringElement();
//
//     typedef itk::BinaryDilateImageFilter<ImageType, ImageType, StructuringElementType> MorphFilterType;
//     MorphFilterType::Pointer morphFilter = MorphFilterType::New();
//     morphFilter->SetKernel( element );
//     morphFilter->SetInput( output );
//     morphFilter->SetBackgroundValue( 0 );
//     morphFilter->SetForegroundValue( 2 );
//     morphFilter->Update();
//
//     typedef itk::ImageFileWriter<ImageType> WriterType;
//     WriterType::Pointer writer = WriterType::New();
//     writer->SetFileName( argv[2] );
//     writer->SetInput( morphFilter->GetOutput() );
//     writer->Update();
//     }
//   else
//     {
