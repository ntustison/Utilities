/*=========================================================================
  
  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: TranslateImage.cxx,v $
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
#include "itkResampleImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "global.h"

int main(int argc, char *argv[])        
{
  if ( argc != ImageDimension + 3 )
    { 
      std::cout << "Usage: "<< argv[0] << " input_image output_image dx dy dz(if ImageDimension == 3) " << std::endl; 
      exit( 1 ); 
    }
  
  typedef double RealType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::TranslationTransform<RealType, ImageDimension> TransformType;
  TransformType::Pointer transform = TransformType::New();
  TransformType::ParametersType parameters( ImageDimension );
  parameters[0] = atof( argv[3] );
  parameters[1] = atof( argv[4] );
  if ( ImageDimension == 3 )
    {
    parameters[2] = atof( argv[5] );
    }
  transform->SetParameters( parameters);
  
  typedef itk::LinearInterpolateImageFunction<ImageType, RealType> LinearInterpolatorType;
  LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
  interpolator->SetInputImage( reader->GetOutput() );

  typedef itk::ResampleImageFilter<ImageType, ImageType, RealType> ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();
  ResamplerType::SpacingType spacing; 
  
  resampler->SetTransform( transform );
  resampler->SetInterpolator( interpolator );
  resampler->SetInput( reader->GetOutput() );
  resampler->SetOutputSpacing( reader->GetOutput()->GetSpacing() );
  resampler->SetOutputOrigin( reader->GetOutput()->GetOrigin() );
  resampler->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  resampler->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( resampler->GetOutput() );
  writer->Update();

 return 0;
}
