/*=========================================================================
 
  Program:   Insight Segmentation & Registration Toolkit 
  Module:    $RCSfile: DeformableRegistration1.cxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:51 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/


#include "itkImageFileReader.h" 
#include "itkImageFileWriter.h" 
#include "itkRescaleIntensityImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkFEM.h"
#include "itkFEMRegistrationFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkDeformationFieldJacobianDeterminantFilter.h"
#include "itkVector.h"
#include "itkAddImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkSquaredDifferenceImageFilter.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMeanSquaresImageToImageMetric.h"


const unsigned int ImageDimension = 2;
 
typedef itk::Image<double, ImageDimension>                              RealImageType;
typedef itk::fem::FEMRegistrationFilter<RealImageType,RealImageType>    RegistrationType;
typedef RegistrationType::PDEDeformableMetricLoadType                    ImageLoadType;
template class itk::fem::ImageMetricLoadImplementation<ImageLoadType>;

typedef itk::fem::Element2DC0LinearQuadrilateralMembrane                ElementType;
typedef ElementType::LoadImplementationFunctionPointer                  LoadImpFP;
typedef ElementType::LoadType                                           ElementLoadType;
typedef itk::fem::VisitorDispatcher<ElementType, ElementLoadType, LoadImpFP>   
                                                           DispatcherType;
          
typedef itk::fem::Element2DBSplinePatchLinearElasticity                 BSplineElementType;
typedef BSplineElementType::LoadImplementationFunctionPointer           BSplineLoadImpFP;
typedef BSplineElementType::LoadType                                    BSplineElementLoadType;
typedef itk::fem::VisitorDispatcher<BSplineElementType, BSplineElementLoadType, BSplineLoadImpFP>   
                                                           BSplineDispatcherType;

int main(int argc, char *argv[])
{
  RegistrationType::Pointer registrationFilter = RegistrationType::New(); 

  ElementType::LoadImplementationFunctionPointer fp = 
    &itk::fem::ImageMetricLoadImplementation<ImageLoadType>::ImplementImageMetricLoad;
  DispatcherType::RegisterVisitor((ImageLoadType*)0, fp);

  BSplineElementType::LoadImplementationFunctionPointer bfp = 
    &itk::fem::ImageMetricLoadImplementation<ImageLoadType>::ImplementImageMetricLoad;
  BSplineDispatcherType::RegisterVisitor((ImageLoadType*)0, bfp);
 
  typedef itk::ImageFileReader<RealImageType>      FileSourceType;
  typedef RealImageType::PixelType PixType;

  FileSourceType::Pointer movingfilter = FileSourceType::New();
  movingfilter->SetFileName("BrianImages/r16slice.hdr");
  movingfilter->Update();
  FileSourceType::Pointer fixedfilter = FileSourceType::New();
  fixedfilter->SetFileName("BrianImages/r85slice.hdr");
  fixedfilter->Update();

  // Rescale the image intensities so that they fall between 0 and 255

  typedef itk::RescaleIntensityImageFilter<RealImageType,RealImageType> FilterType;
  FilterType::Pointer movingrescalefilter = FilterType::New();
  FilterType::Pointer fixedrescalefilter = FilterType::New();

  movingrescalefilter->SetInput(movingfilter->GetOutput());
  fixedrescalefilter->SetInput(fixedfilter->GetOutput());

  const double desiredMinimum =  0.0;
  const double desiredMaximum =  255.0;

  movingrescalefilter->SetOutputMinimum( desiredMinimum );
  movingrescalefilter->SetOutputMaximum( desiredMaximum );
  movingrescalefilter->UpdateLargestPossibleRegion();
  fixedrescalefilter->SetOutputMinimum( desiredMinimum );
  fixedrescalefilter->SetOutputMaximum( desiredMaximum );
  fixedrescalefilter->UpdateLargestPossibleRegion();
  
  // Histogram match the images
  typedef itk::HistogramMatchingImageFilter<RealImageType,RealImageType> HEFilterType;
  HEFilterType::Pointer IntensityEqualizeFilter = HEFilterType::New();

  IntensityEqualizeFilter->SetReferenceImage( fixedrescalefilter->GetOutput() );
  IntensityEqualizeFilter->SetInput( movingrescalefilter->GetOutput() );
  IntensityEqualizeFilter->SetNumberOfHistogramLevels( 100 );
  IntensityEqualizeFilter->SetNumberOfMatchPoints( 15 );
  IntensityEqualizeFilter->ThresholdAtMeanIntensityOn();

  IntensityEqualizeFilter->Update();

  registrationFilter->SetMovingImageInput(IntensityEqualizeFilter->GetOutput());
  registrationFilter->SetFixedImageInput(fixedrescalefilter->GetOutput());

  
//  Software Guide : BeginCodeSnippet
  // Create the material properties


  itk::fem::MaterialLinearElasticity::Pointer m;
  m = itk::fem::MaterialLinearElasticity::New();
  m->GN = 0;                  // Global number of the material
  m->A = 1.0;                 // Cross-sectional area
  m->h = 1.0;                 // Thickness
  m->I = 1.0;                 // Moment of inertia
  m->nu = atof(argv[1]);                 // Poisson's ratio -- DONT CHOOSE 1.0!!
  m->RhoC = 1.0;              // Density
  m->E = atof(argv[2]);


  if (argc < 10) 
  {
    BSplineElementType::Pointer eb1=BSplineElementType::New();
    eb1->m_mat=dynamic_cast<itk::fem::MaterialLinearElasticity*>( m );
    eb1->SetBSplineOrder(3);
    registrationFilter->SetElement(eb1);
  }
  else
  { 
    ElementType::Pointer e1 = ElementType::New();
    e1->m_mat=dynamic_cast<itk::fem::MaterialLinearElasticity*>( m );
    registrationFilter->SetElement(e1);
  }  

  registrationFilter->SetMaterial(m);
  registrationFilter->SetElasticity(m->E);

  registrationFilter->SetTimeStep(atof(argv[3]));
  registrationFilter->SetNumberOfLevels(atoi(argv[6]));
  registrationFilter->SetMeshResolution(atoi(argv[4]), 0); 
  registrationFilter->SetMaximumNumberOfIterations(atoi(argv[5]));
  registrationFilter->SetNumberOfIntegrationPoints(4);

/*
Linear triangle:  1 point
Quadratic triangle: 3 points
Linear quadrilateral: 4 points
Quadratic quadrilateral: 4 points
Linear brick: 8 points
Quadratic brick: 27 points
*/

//  typedef itk::MattesMutualInformationImageToImageMetric
//  typedef itk::MeanSquaresImageToImageMetric
//                      <RealImageType, RealImageType> ImageToImageMetricType;  
//  ImageToImageMetricType::Pointer ImageToImageMetric = ImageToImageMetricType::New(); 

//  ImageToImageMetric->UseAllPixelsOn();
//  ImageToImageMetric->DebugOn();

//  registrationFilter->SetImageToImageMetric(ImageToImageMetric);
  registrationFilter->SetMaximizeMetric(false);
  
  registrationFilter->DebugOn();
 
  try 
  {
    registrationFilter->Update();
  }
  catch (itk::ExceptionObject & exp)
  {
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << exp << std::endl; 
  }  
  
  // Write the warped image
  
  typedef itk::ImageFileWriter<RealImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName("results.hdr");
  writer->SetInput(registrationFilter->GetWarpedImage());
  writer->Update();    

  // Write the jacobian of the displacement field
  typedef RegistrationType::VectorType  VectorType;
  typedef RegistrationType::DeformationFieldType  DeformationFieldType;
  typedef double RealType;
  typedef itk::DeformationFieldJacobianDeterminantFilter
     <DeformationFieldType, RealType> JacobianFilterType;
  JacobianFilterType::Pointer JacobianFilter = JacobianFilterType::New();
  JacobianFilter->SetInput(registrationFilter->GetDeformationField());
  RealImageType::Pointer Ones = RealImageType::New();
  Ones->SetRegions(registrationFilter->GetWarpedImage()->GetLargestPossibleRegion());
  Ones->Allocate();
  Ones->FillBuffer(1.0);
  typedef itk::AddImageFilter<RealImageType, RealImageType, RealImageType> AddFilterType;
  AddFilterType::Pointer AddFilter = AddFilterType::New();
  AddFilter->SetInput1(Ones);
  AddFilter->SetInput2(JacobianFilter->GetOutput());
  AddFilter->Update();  
  /*
  typedef itk::ImageFileWriter<RegistrationType::RealImageType> RealWriterType;
  RealWriterType::Pointer realwriter = RealWriterType::New();
  realwriter->SetFileName("jacobian.hdr");
  realwriter->SetInput(AddFilter->GetOutput());
  realwriter->Update();
  */

  typedef itk::StatisticsImageFilter<RealImageType> StatisticsFilterType;
  StatisticsFilterType::Pointer statfilter = StatisticsFilterType::New();
  statfilter->SetInput(AddFilter->GetOutput());
  statfilter->Update();  
  std::cout << "   Minimum Jacobian = " << statfilter->GetMinimum() << std::endl;
  std::cout << "   Maximum Jacobian = " << statfilter->GetMaximum() << std::endl;
  std::cout << "   Average Jacobian = " << statfilter->GetMean() << std::endl;
  std::cout << "   std Jacobian     = " << statfilter->GetSigma() << std::endl;

  typedef itk::SquaredDifferenceImageFilter<RealImageType, RealImageType, RealImageType> SquaredDistanceFilterType;
  SquaredDistanceFilterType::Pointer sqdfilter = SquaredDistanceFilterType::New();
  sqdfilter->SetInput1(registrationFilter->GetOutput());
  sqdfilter->SetInput2(fixedrescalefilter->GetOutput());
  statfilter->SetInput(sqdfilter->GetOutput());
  statfilter->Update();  
  std::cout << "   Mean squares     = " << statfilter->GetMean() << std::endl;
  std::cout << "   std              = " << statfilter->GetSigma() << std::endl;


  // Write the components of the displacement field
  typedef RegistrationType::DeformationFieldType DeformationFieldType;
  typedef itk::VectorIndexSelectionCastImageFilter
          <DeformationFieldType, RealImageType>      IndexSelectCasterType;
  IndexSelectCasterType::Pointer fieldCaster = IndexSelectCasterType::New();   
  fieldCaster->SetInput( registrationFilter->GetDeformationField() );
  fieldCaster->SetIndex( 0 );
  writer->SetFileName("results_x.hdr");
  writer->SetInput(fieldCaster->GetOutput(0));
  writer->Update();
  fieldCaster->SetIndex( 1 );
  writer->SetFileName("results_y.hdr");
  writer->SetInput(fieldCaster->GetOutput(0));
  writer->Update();


  return 0;
}


