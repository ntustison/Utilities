#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkBoxSpatialObject.h"
#include "itkGroupSpatialObject.h"
#include "itkSpatialObjectToImageFilter.h"

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

template <unsigned int ImageDimension>
int BB( unsigned int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::GroupSpatialObject<ImageDimension> SceneType;
  typedef itk::SpatialObjectToImageFilter<SceneType, ImageType>
    SpatialObjectToImageFilterType;
  typename SceneType::Pointer scene = SceneType::New();

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  std::vector<float> center =
    ConvertVector<float>( std::string( argv[4] ) );
  std::vector<float> dimensions =
    ConvertVector<float>( std::string( argv[5] ) );
  std::vector<float> normal =
    ConvertVector<float>( std::string( argv[6] ) );

  typedef itk::BoxSpatialObject<ImageDimension> BoxType;
  typename BoxType::Pointer box = BoxType::New();

  scene->AddSpatialObject( box );

  typename BoxType::SizeType dimensions2;
  typename BoxType::TransformType::MatrixType matrix;
  typename BoxType::TransformType::OffsetType offset;

  matrix.SetIdentity();

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    double theta = vcl_acos( normal[i] );
    if( ImageDimension == 2 )
      {
      matrix[0][0] = matrix[1][1] = vcl_cos( theta );
      matrix[0][1] = vcl_sin( theta );
      matrix[1][0] = -vcl_sin( theta );
      break;
      }
    else
      {
      typename BoxType::TransformType::MatrixType dmatrix;
      dmatrix.SetIdentity();
      if( i == 0 )
        {
        dmatrix[1][1] = dmatrix[2][2] = vcl_cos( theta );
        dmatrix[1][2] = -vcl_sin( theta );
        dmatrix[2][1] = vcl_sin( theta );
        }
      else if( i == 1 )
        {
        dmatrix[0][0] = dmatrix[2][2] = vcl_cos( theta );
        dmatrix[0][2] = vcl_sin( theta );
        dmatrix[2][0] = -vcl_sin( theta );
        }
      else
        {
        dmatrix[0][0] = dmatrix[1][1] = vcl_cos( theta );
        dmatrix[0][1] = -vcl_sin( theta );
        dmatrix[1][0] = vcl_sin( theta );
        }
      matrix *= dmatrix;
      }
    }

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

//  box->GetParent()->GetObjectToWorldTransform()->Print( std::cout, 3 );



  typename SpatialObjectToImageFilterType::Pointer filter =
    SpatialObjectToImageFilterType::New();
  filter->SetInput( scene );
  filter->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  filter->SetOrigin( reader->GetOutput()->GetOrigin() );
  filter->SetSpacing( reader->GetOutput()->GetSpacing() );
  filter->SetDirection( reader->GetOutput()->GetDirection() );
  filter->SetInsideValue( 1 );
  filter->SetOutsideValue( 0 );
  filter->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 7 )
    {
    std::cout << argv[0] << " imageDimension inputImage outputImage center dimensions normal" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     BB<2>( argc, argv );
     break;
   case 3:
     BB<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

