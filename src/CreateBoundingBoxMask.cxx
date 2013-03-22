#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkBoxSpatialObject.h"
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

int main( int argc, char *argv[] )
{
  if ( argc < 7 )
    {
    std::cout << argv[0] << " inputImage outputImage center dimensions normal inplaneRotationAngle" << std::endl;
    exit( 1 );
    }

  const unsigned int ImageDimension = 3;

  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::GroupSpatialObject<ImageDimension> SceneType;
  typedef itk::SpatialObjectToImageFilter<SceneType, ImageType> SpatialObjectToImageFilterType;
  SceneType::Pointer scene = SceneType::New();

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  std::vector<float> center = ConvertVector<float>( std::string( argv[3] ) );
  std::vector<float> dimensions = ConvertVector<float>( std::string( argv[4] ) );
  std::vector<float> normal = ConvertVector<float>( std::string( argv[5] ) );

  // Find the angle between initial (assume initial = 1,0,0) and normal
  std::vector<float> initial;
  initial.push_back( 1.0 );
  for( unsigned int i = 1; i < ImageDimension; i++ )
    {
    initial.push_back( 0 );
    }

  // Normalize normal vector and initial vector and find angle between them
  float magnitudeNormal = 0.0;
  float magnitudeInitial = 0.0;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    magnitudeNormal += vnl_math_sqr( normal[i] );
    magnitudeInitial += vnl_math_sqr( initial[i] );
    }
  magnitudeNormal = vcl_sqrt( magnitudeNormal );
  magnitudeInitial = vcl_sqrt( magnitudeInitial );

  float dotProduct = 0.0;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    normal[i] /= magnitudeNormal;
    initial[i] /= magnitudeInitial;

    dotProduct += normal[i] * initial[i];
    }
  float theta = vcl_acos( dotProduct );

  std::vector<float> crossProduct;
  crossProduct.push_back( initial[1] * normal[2] - initial[2] * normal[1] );
  crossProduct.push_back( -initial[0] * normal[2] + initial[2] * normal[0] );
  crossProduct.push_back( initial[0] * normal[1] - initial[1] * normal[0] );

  float magnitudeCrossProduct = 0.0;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    magnitudeCrossProduct += vnl_math_sqr( crossProduct[i] );
    }

  magnitudeCrossProduct = vcl_sqrt( magnitudeCrossProduct );
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    crossProduct[i] /= magnitudeCrossProduct;
    }

  // Set up the box spatial object

  typedef itk::BoxSpatialObject<ImageDimension> BoxType;
  BoxType::Pointer box = BoxType::New();

  scene->AddSpatialObject( box );

  BoxType::SizeType dimensions2;
  dimensions2[0] = dimensions[0];
  dimensions2[1] = dimensions[1];
  dimensions2[2] = dimensions[2];

  BoxType::TransformType::MatrixType matrixImage = reader->GetOutput()->GetDirection();

  BoxType::TransformType::MatrixType matrixNormal;
  BoxType::TransformType::OffsetType offset;

  matrixNormal.SetIdentity();

  if( theta != 0.0 )
    {
    float u = crossProduct[0];
    float v = crossProduct[1];
    float w = crossProduct[2];

    float u2 = vnl_math_sqr( u );
    float v2 = vnl_math_sqr( v );
    float w2 = vnl_math_sqr( w );

    float cos_theta = vcl_cos( theta );
    float sin_theta = vcl_sin( theta );

    matrixNormal[0][0] = cos_theta + u2 * ( 1.0 - cos_theta );
    matrixNormal[1][0] = u * v * ( 1.0 - cos_theta ) - w * sin_theta;
    matrixNormal[2][0] = u * w * ( 1.0 - cos_theta ) + v * sin_theta;

    matrixNormal[0][1] = u * v * ( 1.0 - cos_theta ) + w * sin_theta;
    matrixNormal[1][1] = cos_theta + v2 * ( 1.0 - cos_theta );
    matrixNormal[2][1] = v * w * ( 1.0 - cos_theta ) - u * sin_theta;

    matrixNormal[0][2] = u * w * ( 1.0 - cos_theta ) - v * sin_theta;
    matrixNormal[1][2] = v * w * ( 1.0 - cos_theta ) + u * sin_theta;
    matrixNormal[2][2] = cos_theta + w2 * ( 1.0 - cos_theta );
    }

  BoxType::TransformType::MatrixType matrixInPlane;
  matrixInPlane.SetIdentity();

  theta = atof( argv[6] );

  if( theta != 0 )
    {
    itk::Vector<float, 3> normalVector;

    normalVector[0] = normal[0];
    normalVector[1] = normal[1];
    normalVector[2] = normal[2];

//     normalVector = matrixImage * normalVector;

    float u = normalVector[0];
    float v = normalVector[1];
    float w = normalVector[2];

    float u2 = vnl_math_sqr( u );
    float v2 = vnl_math_sqr( v );
    float w2 = vnl_math_sqr( w );

    float cos_theta = vcl_cos( theta );
    float sin_theta = vcl_sin( theta );

    matrixInPlane[0][0] = cos_theta + u2 * ( 1.0 - cos_theta );
    matrixInPlane[1][0] = u * v * ( 1.0 - cos_theta ) - w * sin_theta;
    matrixInPlane[2][0] = u * w * ( 1.0 - cos_theta ) + v * sin_theta;

    matrixInPlane[0][1] = u * v * ( 1.0 - cos_theta ) + w * sin_theta;
    matrixInPlane[1][1] = cos_theta + v2 * ( 1.0 - cos_theta );
    matrixInPlane[2][1] = v * w * ( 1.0 - cos_theta ) - u * sin_theta;

    matrixInPlane[0][2] = u * w * ( 1.0 - cos_theta ) - v * sin_theta;
    matrixInPlane[1][2] = v * w * ( 1.0 - cos_theta ) + u * sin_theta;
    matrixInPlane[2][2] = cos_theta + w2 * ( 1.0 - cos_theta );
    }

  BoxType::TransformType::MatrixType matrix = matrixInPlane * matrixNormal;

  std::cout << matrixNormal << std::endl;
  std::cout << matrixInPlane << std::endl;
  std::cout << matrixImage << std::endl;
  std::cout << matrix << std::endl;

//   matrixInPlane *= matrix;

  for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      dimensions2[i] = dimensions[i];
      offset[i] = center[i];
      for( unsigned int j = 0; j < ImageDimension; j++ )
        {
        offset[i] -= matrix[i][j] * ( 0.5 * dimensions[j] );
        }
      }

  box->SetSize( dimensions2 );
  box->GetObjectToParentTransform()->SetMatrix( matrix );
  box->GetObjectToParentTransform()->SetOffset( offset );
  box->ComputeObjectToWorldTransform();

  box->GetObjectToParentTransform()->Print( std::cout, 3 );

  box->GetParent()->GetObjectToWorldTransform()->Print( std::cout, 3 );

  SpatialObjectToImageFilterType::Pointer filter = SpatialObjectToImageFilterType::New();
  filter->SetInput( scene );
  filter->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  filter->SetOrigin( reader->GetOutput()->GetOrigin() );
  filter->SetSpacing( reader->GetOutput()->GetSpacing() );
  filter->SetDirection( reader->GetOutput()->GetDirection() );
  filter->SetInsideValue( 1 );
  filter->SetOutsideValue( 0 );
  filter->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return 0;
}

