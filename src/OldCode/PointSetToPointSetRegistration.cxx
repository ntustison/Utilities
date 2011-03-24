#include "itkLandmarkBasedTransformInitializer.h"

#include "itkStatisticsImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkArray.h"
#include "itkVector.h"
#include "itkPoint.h"

#include <string>
#include "global.h"

#include "itkResampleImageFilter.h"

int main(int argc, char* argv[] )
{
  if ( argc != 5 )
    {
    std::cout << argv[0] << " fixed_image moving_image image_to_be_transformed outputimagePrefix" << std::endl; 
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer fixedImageReader = ReaderType::New();
  fixedImageReader->SetFileName( argv[1] );
  fixedImageReader->Update();
  ReaderType::Pointer movingImageReader = ReaderType::New();
  movingImageReader->SetFileName( argv[2] );
  movingImageReader->Update();

  typedef itk::StatisticsImageFilter<ImageType> StatisticsFilterType;
  StatisticsFilterType::Pointer statfilter = StatisticsFilterType::New();
  statfilter->SetInput( fixedImageReader->GetOutput() );
  statfilter->Update();

  // Set the transform type..
  typedef itk::VersorRigid3DTransform< double > TransformType;
  TransformType::Pointer transform = TransformType::New();
  typedef itk::LandmarkBasedTransformInitializer< TransformType, 
          ImageType, ImageType > TransformInitializerType;
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();

  // Set fixed and moving landmarks
  TransformInitializerType::LandmarkPointContainer fixedLandmarks;
  TransformInitializerType::LandmarkPointContainer movingLandmarks;
  TransformInitializerType::LandmarkPointType point;
  TransformInitializerType::LandmarkPointType tmp;
    
  itk::ImageRegionIteratorWithIndex<ImageType> Itf
    ( fixedImageReader->GetOutput(), fixedImageReader->GetOutput()->GetLargestPossibleRegion() );
  for ( unsigned int i = 1; i <= 6; i++ )
    {
    unsigned int n = 0;
    for ( unsigned int j = 0; j < ImageDimension; j++ )
      {
      point[j] = 0;
      }
    for ( Itf.GoToBegin(); !Itf.IsAtEnd(); ++Itf )
      {
      if ( Itf.Get() == i )
        {
        ImageType::PointType pt;
        fixedImageReader->GetOutput()->TransformIndexToPhysicalPoint( Itf.GetIndex(), pt );
        for ( unsigned int j = 0; j < ImageDimension; j++ )
          {
          point[j] += pt[j];
          }
        n++;
        }
      }
    for ( unsigned int j = 0; j < ImageDimension; j++ )
      {
      point[j] /= n;
      }
    std::cout << "Fixed: " << point << std::endl;
    fixedLandmarks.push_back( point );  
    }  

  itk::ImageRegionIteratorWithIndex<ImageType> Itm
    ( movingImageReader->GetOutput(), movingImageReader->GetOutput()->GetLargestPossibleRegion() );

  for ( unsigned int i = 1; i <= 6; i++ )
    {
    unsigned int n = 0;
    for ( unsigned int j = 0; j < ImageDimension; j++ )
      {
      point[j] = 0;
      }
    for ( Itm.GoToBegin(); !Itm.IsAtEnd(); ++Itm )
      {
      if ( Itm.Get() == i )
        {
        ImageType::PointType pt;
        movingImageReader->GetOutput()->TransformIndexToPhysicalPoint( Itm.GetIndex(), pt );
        for ( unsigned int j = 0; j < ImageDimension; j++ )
          {
          point[j] += pt[j];
          }
        n++;
        }
      }
    for ( unsigned int j = 0; j < ImageDimension; j++ )
      {
      point[j] /= n;
      }
    std::cout << "Moving: " << point << std::endl;
    movingLandmarks.push_back( point );  
    }  

  initializer->SetFixedLandmarks( fixedLandmarks );
  initializer->SetMovingLandmarks( movingLandmarks );
  initializer->SetTransform( transform );
  initializer->InitializeTransform();
  
  std::cout << "Solution: " << transform->GetParameters() << std::endl;
  
  // Resample image filter

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[3] );
  reader->Update();

  typedef itk::LinearInterpolateImageFunction<ImageType, double> LinearInterpolatorType;
  LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
  interpolator->SetInputImage( reader->GetOutput() );

  typedef itk::ResampleImageFilter<ImageType, ImageType, double> ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();
  resampler->SetTransform( transform );
  resampler->SetInterpolator( interpolator );
  resampler->SetInput( reader->GetOutput() );
  resampler->SetOutputSpacing( reader->GetOutput()->GetSpacing() );
  resampler->SetOutputOrigin( reader->GetOutput()->GetOrigin() );
  resampler->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  resampler->Update();

  std::string filename = std::string( argv[4] ) + std::string( ".hdr" );

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput( resampler->GetOutput() );
  writer->Update();
  
  // Write deformation field  
  
  typedef itk::Vector<double, ImageDimension> VectorType;
  typedef itk::Vector<double, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;

  DeformationFieldType::Pointer field = DeformationFieldType::New();
  field->SetOrigin( reader->GetOutput()->GetOrigin() );
  field->SetSpacing( reader->GetOutput()->GetSpacing() );
  field->SetRegions( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  field->Allocate();

  DeformationFieldType::Pointer inversefield = DeformationFieldType::New();
  inversefield->SetOrigin( reader->GetOutput()->GetOrigin() );
  inversefield->SetSpacing( reader->GetOutput()->GetSpacing() );
  inversefield->SetRegions( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  inversefield->Allocate();

  TransformType::Pointer inversetransform = TransformType::New();
  transform->GetInverse( inversetransform );
  
  itk::ImageRegionIteratorWithIndex<DeformationFieldType> 
    It( field, field->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    DeformationFieldType::PointType point;
    field->TransformIndexToPhysicalPoint( It.GetIndex(), point );
    TransformType::InputPointType P, Ps;
    VectorType V;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      P[i] = point[i];
      }  

    Ps = transform->TransformPoint( P );
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      V[i] = Ps[i] - P[i];
      }
    It.Set( V );   

    Ps = inversetransform->TransformPoint( P );
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      V[i] = Ps[i] - P[i];
      }
    inversefield->SetPixel( It.GetIndex(), V );
    } 
  
  typedef itk::Image<double, ImageDimension> RealImageType;
  typedef itk::VectorImageFileWriter<DeformationFieldType, RealImageType>
    VectorImageFileWriterType;

  VectorImageFileWriterType::Pointer fieldWriter = VectorImageFileWriterType::New();
  fieldWriter->SetUseAvantsNamingConvention( true );
  fieldWriter->SetInput( field );
  fieldWriter->SetFileName( filename.c_str() );
  fieldWriter->Update();

  std::string filename_inverse = std::string( argv[4] ) + std::string( "inverse.hdr" );

  VectorImageFileWriterType::Pointer inversefieldWriter = VectorImageFileWriterType::New();
  inversefieldWriter->SetUseAvantsNamingConvention( true );
  inversefieldWriter->SetInput( inversefield );
  inversefieldWriter->SetFileName( filename_inverse.c_str() );
  inversefieldWriter->Update();

}

