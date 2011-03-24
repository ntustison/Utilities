/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: RegisterWithMultiResDemons.cxx,v $
  Language:  C++
  Date:      $Date: 2008/05/03 01:47:04 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include "itkMultiResolutionPDEDeformableRegistration.h"
#include "itkIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkWarpImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkCommand.h"
#include "vnl/vnl_math.h"

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkVectorImageFileWriter.h"

#include "global.h"

#include <iostream>
#include <string>

namespace
{
  
// The following class is used to support callbacks
// on the filter in the pipeline that follows later
class ShowProgressPDEObject
{
public:
  ShowProgressPDEObject(itk::ProcessObject* o)
    {m_Process = o; m_Prefix="";}
  void ShowProgress()
    {
    std::cout <<  m_Prefix ;
    std::cout << "Progress " << m_Process->GetProgress() << std::endl;
    }
  void ShowIteration()
    {
    std::cout << "Level Completed" << std::endl;
    }
  itk::ProcessObject::Pointer m_Process;
  std::string m_Prefix;
};

template<typename TRegistration>
class PDERegistrationController
{
public:
  PDERegistrationController(TRegistration* o)
    {m_Process = o;}
  void ShowProgress()
    {
    if ( m_Process->GetCurrentLevel() == 3 )
      { 
      m_Process->GetRegistrationFilter()->StopRegistration();
      }
    }
  typename TRegistration::Pointer m_Process;
};

}

int main(int argc, char* argv[] )
{
  if ( argc <= 6 || argc != 5 + atoi( argv[4] ) )
    {
    std::cout << "Usage: RegisterWithMultiResDemons fixed_image moving_image output_image number_of_levels number_of_iterations_level0 ... number_of_iterations_levelN" << std::endl;
    exit( 0 );
    }

  typedef itk::Image<PixelType,ImageDimension> ImageType;  

  typedef float RealType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[1] );
  reader1->Update();
  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  reader2->Update();

  

  std::cout << "Run registration." << std::endl;

  typedef itk::MultiResolutionPDEDeformableRegistration<ImageType,
    ImageType, DeformationFieldType> RegistrationType;
  RegistrationType::Pointer registrator = RegistrationType::New();

  registrator->SetMovingImage( reader2->GetOutput() );
  registrator->SetFixedImage( reader1->GetOutput() );
 
  unsigned int numLevel = atoi( argv[4] );
  unsigned int numIterations[numLevel];
  for ( unsigned int i = 0; i < numLevel; i++ )
    {
    numIterations[i] = atoi( argv[i+5] );
    }
  
  registrator->SetNumberOfLevels( numLevel );
  registrator->SetNumberOfIterations( numIterations );

  registrator->Print(std::cout);

  typedef itk::SimpleMemberCommand<ShowProgressPDEObject> CommandType;

  ShowProgressPDEObject progressWatch(registrator);
  CommandType::Pointer command = CommandType::New();
  command->SetCallbackFunction(&progressWatch,
                               &ShowProgressPDEObject::ShowIteration);
  registrator->AddObserver(itk::IterationEvent(), command );

  PDERegistrationController<RegistrationType> controller(registrator);
  typedef itk::SimpleMemberCommand< PDERegistrationController<RegistrationType> > 
    ControllerType;
  ControllerType::Pointer controllerCommand = ControllerType::New();
  controllerCommand->SetCallbackFunction( 
    &controller, &PDERegistrationController<RegistrationType>::ShowProgress );
  registrator->AddObserver(itk::ProgressEvent(), controllerCommand );

  ShowProgressPDEObject innerWatch(registrator->GetRegistrationFilter() );
  innerWatch.m_Prefix = "    ";
  CommandType::Pointer innerCommand = CommandType::New();
  innerCommand->SetCallbackFunction(&innerWatch,
                               &ShowProgressPDEObject::ShowProgress);
  registrator->GetRegistrationFilter()->
    AddObserver(itk::ProgressEvent(), innerCommand);

  // make registration inplace
  registrator->GetRegistrationFilter()->InPlaceOn();
  registrator->Update();

 
  // -------------------------------------------------------
  std::cout << "Warp moving image" << std::endl;

  typedef itk::WarpImageFilter<ImageType, ImageType, DeformationFieldType> WarperType;
  WarperType::Pointer warper = WarperType::New();

  typedef WarperType::CoordRepType CoordRepType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType,CoordRepType>
    InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  warper->SetInput( reader2->GetOutput() );
  warper->SetDeformationField( registrator->GetOutput() );
  warper->SetInterpolator( interpolator );
  warper->SetOutputSpacing( reader1->GetOutput()->GetSpacing() );
  warper->SetOutputOrigin( reader1->GetOutput()->GetOrigin() );

  warper->Update();

 //-------------------------------------------------------

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( warper->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Write();
 
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::VectorImageFileWriter<DeformationFieldType, RealImageType> FieldWriterType;
  FieldWriterType::Pointer fieldwriter = FieldWriterType::New();
  fieldwriter->SetInput( registrator->GetOutput() );
  fieldwriter->SetFileName( argv[3] );
  fieldwriter->Update();
 
 
  return 0;
}
