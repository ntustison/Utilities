/*=========================================================================
  
  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: TranslationRegistration.cxx,v $
  Language:  C++      
  Date:      $Date: 2008/05/03 01:47:04 $
  Version:   $Revision: 1.1 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
  
=========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkTranslationTransform.h"
#include "itkResampleImageFilter.h"

#include "itkVectorImageFileWriter.h"
#include "itkVectorImageFileReader.h"

#include "global.h"

int main(int argc, char *argv[])        
{
  if ( argc != 6 )
    { 
      std::cout << "Usage: TranslationRegistration fixed_image moving_image output_image interpolator_type(0=NearestNeighbor, 1=Linear) number_of_levels " << std::endl; 
      exit( 1 ); 
    }
  
  typedef double RealType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[1] );
  reader1->Update();
  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  reader2->Update();
  
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();
  double steplength = 0.0;
  for ( unsigned int i = 0; i < ImageDimension; i++ ) 
    {
    steplength += reader1->GetOutput()->GetSpacing()[i]*reader2->GetOutput()->GetSpacing()[i];
    }
  optimizer->SetMaximumStepLength( sqrt( steplength ) ); 
  optimizer->SetMinimumStepLength( 0.01 );
  optimizer->SetNumberOfIterations( 10000 );
  optimizer->MinimizeOn(); 
  
  typedef itk::MattesMutualInformationImageToImageMetric<ImageType, ImageType > MetricType;
  MetricType::Pointer metric = MetricType::New();
  metric->SetNumberOfHistogramBins( 32 );
  metric->SetNumberOfSpatialSamples( 5000 );
  
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, RealType> InterpolatorType0;
  InterpolatorType0::Pointer interpolator0 = InterpolatorType0::New();
  typedef itk::LinearInterpolateImageFunction<ImageType, RealType> InterpolatorType1;
  InterpolatorType1::Pointer interpolator1 = InterpolatorType1::New();

  typedef itk::TranslationTransform<RealType, ImageDimension> TransformType;
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  typedef itk::MultiResolutionPyramidImageFilter<ImageType, ImageType >   PyramidType;
  PyramidType::Pointer fixedImagePyramid = PyramidType::New();
  PyramidType::Pointer movingImagePyramid = PyramidType::New();

  typedef itk::MultiResolutionImageRegistrationMethod<ImageType, ImageType> RegistrationType;
  RegistrationType::Pointer registration = RegistrationType::New();

  registration->SetNumberOfLevels( atoi( argv[5] ) );
  registration->SetFixedImagePyramid( fixedImagePyramid );
  registration->SetMovingImagePyramid( movingImagePyramid );
  registration->SetMetric( metric );
  registration->SetOptimizer( optimizer );
  registration->SetTransform( transform );
  registration->SetInitialTransformParameters( transform->GetParameters() );
  if ( atoi( argv[4] ) == 0 )
    {
    registration->SetInterpolator( interpolator0 );
    }
  else
    {  
    registration->SetInterpolator( interpolator1 );
    }
  registration->SetFixedImage( reader1->GetOutput() );
  registration->SetMovingImage( reader2->GetOutput() );
  registration->SetFixedImageRegion(  reader1->GetOutput()->GetLargestPossibleRegion() );
  registration->StartRegistration(); 
  transform->SetParameters( registration->GetLastTransformParameters() );  

  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetTransform( transform );
  resampler->SetInput( reader2->GetOutput() );
  resampler->SetSize( reader1->GetOutput()->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin(  reader1->GetOutput()->GetOrigin() );
  resampler->SetOutputSpacing( reader1->GetOutput()->GetSpacing() ); 
  resampler->SetDefaultPixelValue( 0 );
  resampler->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( resampler->GetOutput() ); 
  writer->Write();    
  
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

  return 0;
 
}     

