/*=========================================================================
  
  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: RigidVersor3DRegistration.cxx,v $
  Language:  C++      
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
  
=========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageToImageMetric.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkMutualInformationImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkVectorImageFileWriter.h"
#include "itkVectorImageFileReader.h"
#include "itkVersorRigid3DTransform.h"
#include "itkVersorRigid3DTransformOptimizer.h"

#include "global.h"

int main(int argc, char *argv[])        
{
  if ( argc < 7 )
    { 
      std::cout << "Usage: " << argv[0] << " fixedImage movingImage outputImage" 
                << " metric interpolator numberLevels " << std::endl;
      std::cout << "  metric: " << std::endl;
      std::cout << "    0 -> mean squares" << std::endl;
      std::cout << "    1 -> mattes mutual information (default)" << std::endl;
      std::cout << "    2 -> normalized correlation" << std::endl;  
      std::cout << "  interpolator: " << std::endl;
      std::cout << "    0 -> nearest neighbor" << std::endl;
      std::cout << "    1 -> linear (default)" << std::endl;
      exit( 1 ); 
    }
  
  typedef double RealType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer fixedReader = ReaderType::New();
  fixedReader->SetFileName( argv[1] );
  fixedReader->Update();
  ReaderType::Pointer movingReader = ReaderType::New();
  movingReader->SetFileName( argv[2] );
  movingReader->Update();
  
  typedef itk::VersorRigid3DTransformOptimizer OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();
  double steplength = 0.0;
  for ( unsigned int i = 0; i < ImageDimension; i++ ) 
    {
    steplength += ( fixedReader->GetOutput()->GetSpacing()[i]
      *movingReader->GetOutput()->GetSpacing()[i] );
    }
  optimizer->SetMaximumStepLength( vcl_sqrt( steplength ) ); 
  optimizer->SetMinimumStepLength( 0.1*optimizer->GetMaximumStepLength() );
  optimizer->SetNumberOfIterations( 100 );

  typedef itk::VersorRigid3DTransform<RealType> TransformType;
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  typedef itk::MultiResolutionPyramidImageFilter<ImageType, ImageType> PyramidType;
  PyramidType::Pointer fixedImagePyramid = PyramidType::New();
  PyramidType::Pointer movingImagePyramid = PyramidType::New();

  typedef itk::MultiResolutionImageRegistrationMethod<ImageType, ImageType> RegistrationType;
  RegistrationType::Pointer registration = RegistrationType::New();

  registration->SetNumberOfLevels( atoi( argv[5] ) );
  registration->SetFixedImagePyramid( fixedImagePyramid );
  registration->SetMovingImagePyramid( movingImagePyramid );
  registration->SetOptimizer( optimizer );
  registration->SetTransform( transform );
  registration->SetInitialTransformParameters( transform->GetParameters() );

  switch ( atoi( argv[3] ) )
    {
    case 0:
      {
      typedef itk::MeanSquaresImageToImageMetric<ImageType, ImageType> MetricType;
      MetricType::Pointer metric = MetricType::New();
      registration->SetMetric( metric );
      optimizer->MinimizeOn(); 
      break;
      }
    case 1: default:
      {
      typedef itk::MattesMutualInformationImageToImageMetric
        <ImageType, ImageType> MetricType;
      MetricType::Pointer metric = MetricType::New();
      metric->SetNumberOfHistogramBins( 32 );
      metric->SetNumberOfSpatialSamples( 5000 );
      registration->SetMetric( metric );
      optimizer->MinimizeOff(); 
      break;
      } 
    case 2:
      {
      typedef itk::NormalizedCorrelationImageToImageMetric
        <ImageType, ImageType> MetricType;
      MetricType::Pointer metric = MetricType::New();
      registration->SetMetric( metric );
      optimizer->MinimizeOn(); 
      break;
      } 
    }  

  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, RealType> InterpolatorType0;
  InterpolatorType0::Pointer interpolator0 = InterpolatorType0::New();
  typedef itk::LinearInterpolateImageFunction<ImageType, RealType> InterpolatorType1;
  InterpolatorType1::Pointer interpolator1 = InterpolatorType1::New();

  if ( atoi( argv[4] ) == 0 )
    {
    registration->SetInterpolator( interpolator0 );
    }
  else
    {  
    registration->SetInterpolator( interpolator1 );
    }

  registration->SetFixedImage( fixedReader->GetOutput() );
  registration->SetMovingImage( movingReader->GetOutput() );
  registration->SetFixedImageRegion(  fixedReader->GetOutput()->GetLargestPossibleRegion() );
  registration->DebugOn();
  registration->StartRegistration(); 
  transform->SetParameters( registration->GetLastTransformParameters() );  

  std::cout << "Transformation parameters: " 
            << registration->GetLastTransformParameters() << std::endl;
  std::cout << "  Rotation matrix: " << std::endl << transform->GetRotationMatrix() << std::endl;
  std::cout << "  Offset: " << std::endl << transform->GetOffset() << std::endl;


  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetTransform( transform );
  resampler->SetInput( movingReader->GetOutput() );
  resampler->SetSize( fixedReader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin( fixedReader->GetOutput()->GetOrigin() );
  resampler->SetOutputSpacing( fixedReader->GetOutput()->GetSpacing() ); 
  resampler->SetDefaultPixelValue( 0 );
  resampler->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( resampler->GetOutput() ); 
  writer->Write();    

/*  
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;
  DeformationFieldType::Pointer field = DeformationFieldType::New();
  field->SetRegions( reader1->GetOutput()->GetLargestPossibleRegion() );
  field->SetOrigin( reader1->GetOutput()->GetOrigin() );
  field->SetSpacing( reader1->GetOutput()->GetSpacing() );
  field->Allocate();  
  
  TransformType::InputPointType in;
  TransformType::OutputPointType out;
  VectorType vec;
  
  itk::ImageRegionIteratorWithIndex<DeformationFieldType> Id( 
    field, field->GetLargestPossibleRegion() );
  for ( Id.GoToBegin(); !Id.IsAtEnd(); ++Id )
    {
    field->TransformIndexToPhysicalPoint( Id.GetIndex(), in );
    out = transform->TransformPoint( in );
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      vec[i] = out[i] - in[i];
      }
    Id.Set( vec );  
    }

  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::VectorImageFileWriter<DeformationFieldType, RealImageType> FieldWriterType;
  FieldWriterType::Pointer fieldwriter = FieldWriterType::New();
  fieldwriter->SetInput( field );
  fieldwriter->SetFileName( argv[3] );
  fieldwriter->Update(); 
*/
  return 0;
 
}     

