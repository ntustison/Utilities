/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkFEMRegistrationFilter.hxx,v $
Language:  C++

Date:      $Date: 2008/10/18 00:22:58 $
Version:   $Revision: 1.1.1.1 $

Copyright (c) Insight Software Consortium. All rights reser
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for detail.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef _itkFEMRegistrationFilter_hxx_
#define _itkFEMRegistrationFilter_hxx_

// disable debug warnings in MS compiler
#ifdef _MSC_VER
#pragma warning(disable: 4786)
#endif
 
#include "itkFEMRegistrationFilter.h"
#include "itkFEMElements.h"
#include "itkFEMLoadBC.h"
#include "itkFEMGenerateMesh.h"
#include "itkFEMGenerateBSplineMesh.h"

#include "itkBSplineKernelFunction.h"

#include "itkMeanSquareRegistrationFunction.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"
#include "itkVectorExpandImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"

#include "vnl/vnl_math.h"
#include "vnl/algo/vnl_determinant.h"
#include "vnl/algo/vnl_matrix_inverse.h"

namespace itk {
namespace fem {

template<class TMovingImage, class TFixedImage, class TWarpedImage>
FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::FEMRegistrationFilter( )
{
  this->SetNumberOfRequiredInputs(2);

  m_MeshResolution = ArrayContainerType::New();
  m_MeshResolution->Initialize();

  // Default values  
  this->SetDoLineSearchOnImageEnergy(false);
  this->SetLineSearchMaximumIterations(100);
  this->SetEnergyReductionFactor(0.5);
  this->SetAlpha(0.5);
  this->SetTimeStep(1.0);
  this->SetNumberOfLevels(1);
  this->SetElasticity(1.0e4); 
  this->SetRho(1.0e4);
  this->SetGamma(1.0); 
  this->SetMaximumNumberOfIterations(20);
  this->SetMeshResolution(64);
  this->SetNumberOfIntegrationPoints(4);
  this->SetMetricRegionWidth(1);  
  this->SetUseMassMatrix(true);

  // Default element and material
  typename MaterialLinearElasticity::Pointer 
     material = MaterialLinearElasticity::New();
  material->GN = 0;                  // Global number of the m_Material
  material->A = 1.0;                 // Cross-sectional area
  material->h = 1.0;                 // Thickness
  material->I = 1.0;                 // Moment of inertia
  material->nu = 0.0;                 // Poisson's ratio -- DONT CHOOSE 1.0!!
  material->RhoC = 1.0;              // Density
  material->E = this->GetElasticity();
  this->SetMaterial(material);

  m_UseImageToImageMetric = false;

  if (ImageDimension == 2)
  {
    typedef Element2DC0LinearQuadrilateralMembrane ElementType;
    typedef VisitorDispatcher<ElementType, 
       typename ElementType::LoadType, 
       typename ElementType::LoadImplementationFunctionPointer> DispatcherType;
    typename ElementType::LoadImplementationFunctionPointer fp = 
        &ImageMetricLoadImplementation<PDEDeformableMetricLoadType>::ImplementImageMetricLoad;
    DispatcherType::RegisterVisitor((PDEDeformableMetricLoadType*)0, fp);

    typename ElementType::Pointer element = ElementType::New();
    element->m_mat = dynamic_cast<MaterialLinearElasticity*>(m_Material);
    this->SetElement(element);
  }
  if (ImageDimension == 3)
  {
    typedef Element3DC0LinearHexahedronMembrane ElementType;
    typedef VisitorDispatcher<ElementType, 
       typename ElementType::LoadType, 
       typename ElementType::LoadImplementationFunctionPointer> DispatcherType;
    typename ElementType::LoadImplementationFunctionPointer fp = 
        &ImageMetricLoadImplementation<PDEDeformableMetricLoadType>::ImplementImageMetricLoad;
    DispatcherType::RegisterVisitor((PDEDeformableMetricLoadType*)0, fp);

    typename ElementType::Pointer element = ElementType::New();
    element->m_mat = dynamic_cast<MaterialLinearElasticity*>(m_Material);
    this->SetElement(element);
  }
 
  // Set up the default metric type (mean squares)
  typedef MeanSquareRegistrationFunction
        <FixedImageType, MovingImageType, DeformationFieldType> DefaultMetricType;  
  typename DefaultMetricType::Pointer DefaultMetric
          = DefaultMetricType::New();
  m_PDEDeformableMetric = DefaultMetric;
  m_PDEDeformableMetric->SetNormalizeGradient(false);
  this->SetMaximizeMetric(false);

  // Set up the default image interpolator
  typedef LinearInterpolateImageFunction
          <MovingImageType, RealType> DefaultImageInterpolatorType;
  typename DefaultImageInterpolatorType::Pointer Interpolator
          = DefaultImageInterpolatorType::New();
  m_ImageInterpolator = Interpolator;
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void 
FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::SetNumberOfLevels(unsigned int n)
{
  if (n > 0)
  {  
    RealType tmp_E, tmp_G, tmp_R;
    unsigned int tmp_MI, tmp_IP, tmp_MW;
    ArrayType tmp_MR;
    
    if (m_NumberOfLevels != 0)
    {
      tmp_E = m_Elasticity[0];
      tmp_G = m_Gamma[0];
      tmp_R = m_Rho[0];
   
      tmp_MI = m_MaximumNumberOfIterations[0];
      tmp_MR = m_MeshResolution->ElementAt(0);
      tmp_IP = m_NumberOfIntegrationPoints[0];
      tmp_MW = m_MetricRegionWidth[0];
    }
    
    m_Elasticity.set_size(n);
    m_Gamma.set_size(n);
    m_Rho.set_size(n);
    m_MaximumNumberOfIterations.set_size(n);
    m_NumberOfIntegrationPoints.set_size(n);
    m_MetricRegionWidth.set_size(n);

    if (m_NumberOfLevels != 0)
    {
      m_Elasticity.fill(tmp_E);
      m_Gamma.fill(tmp_G);
      m_Rho.fill(tmp_R);
      m_MaximumNumberOfIterations.fill(tmp_MI);
      m_NumberOfIntegrationPoints.fill(tmp_IP);
      m_MetricRegionWidth.fill(tmp_MW);
      for (unsigned int i = 0; i < ImageDimension; i++)
      {
        m_MeshResolution->InsertElement(i, tmp_MR);
      }
    }
    m_NumberOfLevels = n;
    this->Modified();
  }  
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void 
FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::GenerateData()
{
  if (this->GetNumberOfInputs() < 2)
  {
    itkExceptionMacro(<< "Images are not specified.");
  }

  this->RegisterImages();

  if (m_DeformationField)
  {
    typedef WarpImageFilter<MovingImageType, 
                            WarpedImageType, 
                            DeformationFieldType> WarperType;
    typename WarperType::Pointer warper = WarperType::New();

    warper->SetInput(this->GetInput(0));
    warper->SetDeformationField(m_DeformationField);
    warper->SetInterpolator(m_ImageInterpolator);
    warper->SetOutputSpacing(this->GetInput(1)->GetSpacing());
    warper->SetOutputOrigin(this->GetInput(1)->GetOrigin());
    warper->Update();
    m_WarpedImage = warper->GetOutput();  
    this->GraftNthOutput(0, m_WarpedImage);
  }
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void 
FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::RegisterImages()
{
  typedef RecursiveMultiResolutionPyramidImageFilter
                          <FixedImageType, FixedImageType>   FixedPyramidType;
  typedef RecursiveMultiResolutionPyramidImageFilter
                          <MovingImageType, MovingImageType> MovingPyramidType;
  
  typename MovingPyramidType::Pointer MovingPyramid;
  MovingPyramid = MovingPyramidType::New();
  MovingPyramid->SetInput(this->GetInput(0));
  MovingPyramid->SetNumberOfLevels(m_NumberOfLevels);

  typename FixedPyramidType::Pointer FixedPyramid;
  FixedPyramid = FixedPyramidType::New();
  FixedPyramid->SetInput(this->GetInput(1));
  FixedPyramid->SetNumberOfLevels(m_NumberOfLevels);

  typename FixedPyramidType::ScheduleType schedule(m_NumberOfLevels, ImageDimension);
  for (int i = m_NumberOfLevels-1; i >= 0; i--)
  {
    RealType factor = pow(2.0, static_cast<RealType>(m_NumberOfLevels-i-1));
    for (unsigned int j = 0; j < ImageDimension; j++)
    {
      if (static_cast<RealType>(this->GetInput(1)->GetRequestedRegion().GetSize()[j])/factor >= 32.0)
      { 
        schedule(i, j) = static_cast<unsigned int>(factor);
      }
      else
      {
        schedule(i, j) = (i < m_NumberOfLevels-1) ? schedule(i+1, j) : 1;
      }
    }
  }
  MovingPyramid->SetSchedule(schedule); 
  MovingPyramid->Update();
  FixedPyramid->SetSchedule(schedule); 
  FixedPyramid->Update();

  for (unsigned int i = 0; i < m_NumberOfLevels; i++)
  {
    m_PyramidLevelImageSizes.push_back(MovingPyramid->GetOutput(i)->GetLargestPossibleRegion().GetSize());
  } 
  
  VectorType V;  V.Fill(0.0);  
  m_DeformationField = DeformationFieldType::New();
  m_DeformationField->SetOrigin(this->GetInput(1)->GetOrigin());
  m_DeformationField->SetSpacing(this->GetInput(1)->GetSpacing());
  m_DeformationField->SetRegions(m_PyramidLevelImageSizes[0]);
  m_DeformationField->Allocate(); 
  m_DeformationField->FillBuffer(V);

  for (m_CurrentLevel = 0; m_CurrentLevel < m_NumberOfLevels; m_CurrentLevel++)
  {
    itkDebugMacro("Current Level = " << m_CurrentLevel << ", Image Size = " << m_PyramidLevelImageSizes[m_CurrentLevel]);

    if (m_MaximumNumberOfIterations[m_CurrentLevel] > 0)
    {
      RealType scaling[ImageDimension];
      for (unsigned int i = 0; i < ImageDimension; i++)
      {
        scaling[i] = static_cast<RealType>(MovingPyramid->GetSchedule()[m_CurrentLevel][i]); 
      }

      m_Solver.SetDeltatT(m_TimeStep);  
      m_Solver.SetRho(m_Rho[m_CurrentLevel]);     
      m_Solver.SetAlpha(m_Alpha);    

      this->GenerateFEMMesh(m_MeshResolution->ElementAt(m_CurrentLevel)); 
      this->ApplyLoads(scaling);
      this->ApplyImageLoads(MovingPyramid->GetOutput(m_CurrentLevel), FixedPyramid->GetOutput(m_CurrentLevel));

      LinearSystemWrapperItpack itpackWrapper; 
      itpackWrapper.SetTolerance(0.1);
      itpackWrapper.JacobianConjugateGradient(); 
      itpackWrapper.SetMaximumNumberIterations(2*m_Solver.GetNumberOfDegreesOfFreedom()); 
      itpackWrapper.SetMaximumNonZeroValuesInMatrix(m_Solver.GetMaximumNumberOfNonZeroElements());
      m_Solver.SetLinearSystemWrapper(&itpackWrapper); 

      if (m_UseMassMatrix) 
      {
        m_Solver.AssembleKandM(); 
      }
      else
      {
        m_Solver.InitializeForSolution();
        m_Solver.AssembleK();
      }  
      if (m_CurrentLevel > 0)  
      {
        this->PropogateSolutionToTheNextLevel();
      }
      this->IterativeSolve();

      if (m_UseBSplines && m_NumberOfLevels > 0 && m_CurrentLevel < m_NumberOfLevels-1)
      {
        this->GetBSplineNodalValuesForNextLevel();
      }
    }

    if (m_CurrentLevel < m_NumberOfLevels-1 && m_DeformationField) 
    {
      typedef VectorExpandImageFilter<DeformationFieldType, 
                                      DeformationFieldType> ExpanderType;
      typename ExpanderType::Pointer DeformationFieldExpander = ExpanderType::New(); 
      typename ExpanderType::ExpandFactorsType factors[ImageDimension];
      for (unsigned int i = 0; i < ImageDimension; i++)
      {
        factors[i] = static_cast<RealType>(FixedPyramid->GetOutput(m_CurrentLevel+1)->GetLargestPossibleRegion().GetSize()[i])
                     /static_cast<RealType>(FixedPyramid->GetOutput(m_CurrentLevel)->GetLargestPossibleRegion().GetSize()[i]);
      }
      DeformationFieldExpander->SetInput(m_DeformationField);
      DeformationFieldExpander->SetExpandFactors(factors);
      DeformationFieldExpander->Update();
      m_DeformationField = DeformationFieldExpander->GetOutput();
    }
  }
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void 
FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::IterativeSolve()
{
  unsigned int iters = 0;
  RealType MinimumEnergy = 10.e99;
  RealType LastE, deltE = 1.0;

  while (deltE > 0.0 && iters++ < m_MaximumNumberOfIterations[m_CurrentLevel])
  {
    itkDebugMacro(<< "Iteration = " << iters);

    if (m_UseImageToImageMetric) 
    {
      LastE = m_ImageToImageMetricLoad->GetCurrentEnergy();
      m_ImageToImageMetricLoad->SetCurrentEnergy(0.0);
      m_ImageToImageMetricLoad->InitializeMetric();
    }  
    else
    {
      LastE = m_PDEDeformableMetricLoad->GetCurrentEnergy();
      m_PDEDeformableMetricLoad->SetCurrentEnergy(0.0);
      m_PDEDeformableMetricLoad->GetMetric()->InitializeIteration();
    } 
    itkDebugMacro(<< "Current energy = " << LastE);

    (m_UseMassMatrix) ? m_Solver.AssembleFforTimeStep() : m_Solver.AssembleF();
    m_Solver.Solve();                  

    itkDebugMacro(<< "Linear system solved");

    if (m_UseImageToImageMetric)
    {
      deltE = vnl_math_max(LastE - m_ImageToImageMetricLoad->GetCurrentEnergy(), 
                         m_ImageToImageMetricLoad->GetCurrentEnergy() - LastE);
    }
    else
    {
      deltE = vnl_math_max(LastE - m_PDEDeformableMetricLoad->GetCurrentEnergy(), 
                         m_PDEDeformableMetricLoad->GetCurrentEnergy() - LastE);
    } 

    itkDebugMacro(<< "Energy change = " << deltE);

    if (m_DoLineSearchOnImageEnergy && deltE < 0.0)
    {
      LastE = this->GoldenSection(m_LineSearchMaximumIterations);
      deltE = (MinimumEnergy - LastE);
    }

    float curmaxsol = m_Solver.GetCurrentMaxSolution();
    if (curmaxsol == 0.0)
    {
      curmaxsol = 1.0;
    }
    m_Solver.AddToDisplacements(vnl_math_min(1.0, m_Gamma[m_CurrentLevel]/curmaxsol));
    MinimumEnergy = LastE; 

    this->InterpolateVectorField();
  }
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void 
FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::PropogateSolutionToTheNextLevel()
{

  if (m_UseBSplines)
  {
    this->DoubleResolutionOfFEMBSplineMesh();
  }
  else
  {
    // For the general case in which the displacement field interpolates
    // the nodal values, we simply interpolate the new nodal values based
    // on the expanded displacement field.

    typedef VectorLinearInterpolateImageFunction
          <DeformationFieldType, RealType>  DeformationFieldInterpolatorType; 
    typename DeformationFieldInterpolatorType::Pointer DeformationFieldInterpolator = 
           DeformationFieldInterpolatorType::New();

    DeformationFieldInterpolator->SetInputImage(m_DeformationField);

    Node::ArrayType* nodes = &(m_Solver.node);
    Node::ArrayType::iterator it; 
   
    for (it = nodes->begin(); it != nodes->end(); it++) 
    {
      Element::VectorType coord = (*it)->GetCoordinates();

      typename DeformationFieldInterpolatorType::ContinuousIndexType idx;
      typename DeformationFieldInterpolatorType::OutputType value;
      
      for (unsigned int i = 0; i < ImageDimension; i++) 
      {
        idx[i] = static_cast<RealType>(coord[i]);
      }
      value = DeformationFieldInterpolator->EvaluateAtContinuousIndex(idx);

      for (unsigned int i = 0; i < ImageDimension; i++)
      { 
        m_Solver.GetLinearSystemWrapper()->
          SetSolutionValue((*it)->GetDegreeOfFreedom(i), value[i], m_Solver.TotalSolutionIndex);    
        m_Solver.GetLinearSystemWrapper()->
          SetSolutionValue((*it)->GetDegreeOfFreedom(i), value[i], m_Solver.SolutionTMinus1Index);    
      }
    }
  } 
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void 
FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::GenerateFEMMesh(ArrayType resolution)
{
 if (!m_MeshFileName.empty())
  {    
    itkDebugMacro(<< "Reading mesh file. " );
    
    std::ifstream meshstream; 
    meshstream.open(m_MeshFileName.c_str());
    if (!meshstream)
    {
      itkExceptionMacro(<< "The mesh file, " << m_MeshFileName << ", was not found!" );
    }
    m_Solver.Read(meshstream); 
    m_Solver.GenerateGFN();

    itk::fem::MaterialLinearElasticity::Pointer m = dynamic_cast<MaterialLinearElasticity*>(m_Solver.mat.Find(0));   
    m_Material->E = this->GetElasticity(m_CurrentLevel);

    // now scale the mesh to the current scale
    Element::VectorType coord;  
    Node::ArrayType* nodes = &(m_Solver.node);
    Node::ArrayType::iterator node=nodes->begin();
    m_Element=( *((*node)->m_elements.begin()));
    for(node = nodes->begin(); node != nodes->end(); node++) 
    {
      coord = (*node)->GetCoordinates();    
      for (unsigned int i = 0; i < ImageDimension; i++)
      { 
        coord[i] /= pow(2.0, static_cast<RealType>(m_NumberOfLevels - m_CurrentLevel - 1));
      } 
      (*node)->SetCoordinates(coord);  
    }
    return;
  }

  vnl_vector<RealType> MeshOriginV(ImageDimension); 
  vnl_vector<RealType> MeshSizeV(ImageDimension); 
  vnl_vector<RealType> ImageSizeV(ImageDimension); 
  vnl_vector<RealType> ElementsPerDim(ImageDimension); 
  
  for (unsigned int i = 0; i < ImageDimension; i++)
  { 
    MeshSizeV[i] = static_cast<RealType>(m_PyramidLevelImageSizes[m_CurrentLevel][i]) - 1.0;
    MeshOriginV[i] = static_cast<RealType>(this->GetInput(1)->GetOrigin()[i]);
    ImageSizeV[i] = static_cast<RealType>(m_PyramidLevelImageSizes[m_CurrentLevel][i]);  
    ElementsPerDim[i] = resolution[i];
  }

  itkDebugMacro(<< "Generating FEM mesh. ");
  itkDebugMacro(<< "Number of mesh elements per dimension = " << ElementsPerDim);

  m_Material->E = this->GetElasticity(m_CurrentLevel);

  if (ImageDimension == 2)
  {
    (m_UseBSplines) ?
       Generate2DRectilinearBSplineMesh(m_Element, m_Solver, MeshOriginV, MeshSizeV, ElementsPerDim, 
                                 dynamic_cast<Element2DBSplinePatch*>(m_Element)->GetBSplineOrder())
     : Generate2DRectilinearMesh(m_Element, m_Solver, MeshOriginV, MeshSizeV, ElementsPerDim); 
  } 
  else if (ImageDimension == 3)  
  {
    (m_UseBSplines) ?
       Generate3DRectilinearBSplineMesh(m_Element, m_Solver, MeshOriginV, MeshSizeV, ElementsPerDim, 
                                 dynamic_cast<Element3DBSplinePatch*>(m_Element)->GetBSplineOrder())
     : Generate3DRectilinearMesh(m_Element, m_Solver, MeshOriginV, MeshSizeV, ElementsPerDim); 
  }
  else 
  {  
    throw FEMException(__FILE__, __LINE__, "GenerateFEMMesh ");
    return;
  }
  m_Solver.GenerateGFN();
  
  itkDebugMacro(<< "Initialize interpolation grid: (image size = " << ImageSizeV << ", MeshSize = " << MeshSizeV << ")" );
  m_Solver.InitializeInterpolationGrid(ImageSizeV, MeshOriginV, MeshSizeV);
}


template<class TMovingImage, class TFixedImage, class TWarpedImage>
void 
FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::ApplyImageLoads(MovingImageType*  MovingImage, FixedImageType* FixedImage)
{ 
  itkDebugMacro(<< "Applying the image loads.")

  ImageSizeType r;
  r.Fill(m_MetricRegionWidth[m_CurrentLevel]);

  if (m_UseImageToImageMetric)
  { 
    m_ImageToImageMetricLoad = ImageToImageMetricLoadType::New();
    m_ImageToImageMetricLoad->SetMovingImage(MovingImage); 
    m_ImageToImageMetricLoad->SetFixedImage(FixedImage);  
    m_ImageToImageMetricLoad->SetMetric(m_ImageToImageMetric);
    m_ImageToImageMetricLoad->SetDeformationField(m_DeformationField);
    m_ImageToImageMetricLoad->InitializeMetric();
    m_ImageToImageMetricLoad->SetGamma(m_Gamma[m_CurrentLevel]);
    m_ImageToImageMetricLoad->SetMaximizeMetric(m_MaximizeMetric);
    m_ImageToImageMetricLoad->SetMetricRadius(r);
    m_ImageToImageMetricLoad->SetNumberOfIntegrationPoints(m_NumberOfIntegrationPoints[m_CurrentLevel]);
    m_ImageToImageMetricLoad->GN = m_Solver.load.size()+1; 
    m_Solver.load.push_back(FEMP<Load>(&*m_ImageToImageMetricLoad));    
    m_ImageToImageMetricLoad = dynamic_cast<ImageToImageMetricLoadType*> 
      (&*m_Solver.load.Find(m_Solver.load.size()));  
  }
  else
  {
    m_PDEDeformableMetric->SetGradientStep(m_Gamma[m_CurrentLevel]);
    m_PDEDeformableMetricLoad = PDEDeformableMetricLoadType::New();
    m_PDEDeformableMetricLoad->SetMovingImage(MovingImage); 
    m_PDEDeformableMetricLoad->SetFixedImage(FixedImage);  
    m_PDEDeformableMetricLoad->SetDeformationField(m_DeformationField);
    m_PDEDeformableMetricLoad->SetMetric(m_PDEDeformableMetric);
    m_PDEDeformableMetricLoad->GetMetric()->InitializeIteration();
    m_PDEDeformableMetricLoad->SetGamma(m_Gamma[m_CurrentLevel]);
    m_PDEDeformableMetricLoad->SetMetricRadius(r);
    m_PDEDeformableMetricLoad->SetNumberOfIntegrationPoints(m_NumberOfIntegrationPoints[m_CurrentLevel]);
    m_PDEDeformableMetricLoad->GN = m_Solver.load.size()+1; 
    m_Solver.load.push_back(FEMP<Load>(&*m_PDEDeformableMetricLoad));    
    m_PDEDeformableMetricLoad = dynamic_cast<PDEDeformableMetricLoadType*> 
      (&*m_Solver.load.Find(m_Solver.load.size()));  
  }
  itkDebugMacro(<< "Finished applying the image loads.")
}   


template<class TMovingImage, class TFixedImage, class TWarpedImage>
void 
FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::ApplyLoads(RealType* scaling)
{
  itkDebugMacro(<< "Applying the non-image loads." );

  vnl_vector<RealType> pd; pd.set_size(ImageDimension);
  vnl_vector<RealType> pu; pu.set_size(ImageDimension);

  if (!m_LandmarkFileName.empty())
  {
    std::ifstream f(m_LandmarkFileName.c_str());
    if (!f)
    {
      itkExceptionMacro(<< "Cannot read landmark file " << m_LandmarkFileName );
    } 

    itkDebugMacro(<< "Loading landmarks from file." );

    Load::ArrayType::iterator it;
    LoadLandmark::Pointer l3;
    m_LandmarkArray.clear();

    try
    { 
      m_Solver.load.clear(); // NOTE: CLEARING ALL LOADS - LMS MUST BE APPLIED FIRST
      m_Solver.Read(f);
    }
    catch (itk::ExceptionObject &err)
    { 
      itkExceptionMacro(<< "Cannot load landmarks " << m_LandmarkFileName );
    }
    f.close();
        
    m_LandmarkArray.resize(m_Solver.load.size());
    unsigned int ct = 0;
    for(it = m_Solver.load.begin(); it !=m_Solver.load.end(); it++) 
    {
      if (l3 = dynamic_cast<LoadLandmark*>(&(*(*it)))) 
      {
        m_LandmarkArray[ct++] = dynamic_cast<LoadLandmark*>(l3->Clone());
      }
    }
    m_Solver.load.clear(); // NOTE: CLEARING ALL LOADS - LMS MUST BE APPLIED FIRST

    if ( !m_LandmarkArray.empty())
    {
      itkDebugMacro(<< "Scaling the landmarks." );
      
      for(unsigned int i = 0; i < m_LandmarkArray.size(); i++) 
      {
        m_LandmarkArray[i]->el[0]=NULL;
        if (scaling) 
        {
          m_LandmarkArray[i]->ScalePointAndForce(scaling, m_EnergyReductionFactor);
        }           
        pu = m_LandmarkArray[i]->GetSource();
        pd = m_LandmarkArray[i]->GetPoint();

        Element::ArrayType::const_iterator it;
        for (it = m_Solver.el.begin(); it != m_Solver.el.end(); it++) 
        {
          if ((*it)->GetLocalFromGlobalCoordinates(pu, pd))
          { 
            m_LandmarkArray[i]->SetPoint(pd);
            m_LandmarkArray[i]->el[0] = (&**it);
     break;
          }
        }
        m_LandmarkArray[i]->GN = i;            
        LoadLandmark::Pointer l5 = (dynamic_cast<LoadLandmark::Pointer>(m_LandmarkArray[i]->Clone()));
        m_Solver.load.push_back(FEMP<Load>(l5)); 
      } 
    }
  }

  itkDebugMacro(<< "Applying the boundary condition loads (i.e. pinning the image corner)." );

  LoadBC::Pointer l1;
  Node::ArrayType* nodes = &(m_Solver.node);

  unsigned int CornerCounter, EdgeCounter = 0;
  Element::VectorType coord; 
  Node::ArrayType::iterator node=nodes->begin();
  unsigned int nodect = 0;
  unsigned int nloads = 0;

  while(node++ != nodes->end() && EdgeCounter < ImageDimension) 
  {
    coord = (*node)->GetCoordinates();
    CornerCounter=0;
    for (unsigned int i = 0; i < ImageDimension; i++)
    { 
      if (coord[i] == this->GetInput(1)->GetOrigin()[i] || coord[i] == m_PyramidLevelImageSizes[m_CurrentLevel][i]-1) 
      {
        CornerCounter++;
      } 
    }
    if (CornerCounter == ImageDimension)  // the node is located at a true corner
    {
      unsigned int ndofpernode = (*((*node)->m_elements.begin()))->GetNumberOfDegreesOfFreedomPerNode();
      unsigned int numnodesperelt = (*((*node)->m_elements.begin()))->GetNumberOfNodes();
      
      typedef typename Node::SetOfElements NodeEltSetType;
      NodeEltSetType::iterator it;
      for (it = (*node)->m_elements.begin(); it!=(*node)->m_elements.end(); it++) 
      {
        for (unsigned int whichnode = 0; whichnode <= numnodesperelt-1; whichnode++)
        {
          coord = (*it)->GetNode(whichnode)->GetCoordinates();
          CornerCounter = 0;          
          for (unsigned int i = 0; i < ImageDimension; i++)
   {
            if (coord[i] == this->GetInput(1)->GetOrigin()[i] || coord[i] == m_PyramidLevelImageSizes[m_CurrentLevel][i]-1 ) 
     {  
              CornerCounter++;
     }  
   }        
          if (CornerCounter == ImageDimension - 1) // edge found
          {
            for (unsigned int j = 0; j < ndofpernode; j++)
            {            
              l1=LoadBC::New();
              // now we get the element from the node -- 
       // we assume we need fix the dof only once
              // even if more than one element shares it.
              l1->m_element = (*it);  
              l1->m_dof = whichnode*ndofpernode + j; 
              l1->m_value = vnl_vector<RealType>(1, 0.0);
              m_Solver.load.push_back(FEMP<Load>(&*l1));
            }
            EdgeCounter++;
          }
        }
      }
    }
    nodect++;
  }
  itkDebugMacro(<< "Finished applying the non-image loads." );
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void 
FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::InterpolateVectorField()
{
  vnl_vector<RealType> Pos(ImageDimension);  // solution at the point
  vnl_vector<RealType> Gpt(ImageDimension);  // global position given by local point

  DeformationFieldIteratorType It(m_DeformationField, m_DeformationField->GetLargestPossibleRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    for (unsigned int d = 0; d < ImageDimension; d++)
    { 
      Gpt[d] = static_cast<RealType>(It.GetIndex()[d]);
    }      
    Element::ConstPointer eltp = m_Solver.GetElementAtPoint(Gpt);   

    if (eltp)
    {
      eltp->GetLocalFromGlobalCoordinates(Gpt, Pos);
      typename Element::VectorType shapef(eltp->GetNumberOfNodes());
      shapef = eltp->ShapeFunctions(Pos);

      VectorType Sol; 
      for(unsigned int f = 0; f < ImageDimension; f++)
      {
 Sol[f] = 0.0;
 for(unsigned int n = 0; n < eltp->GetNumberOfNodes(); n++)
 {
          Sol[f] += shapef[n] * m_Solver.GetLS()->GetSolutionValue( 
            eltp->GetNode(n)->GetDegreeOfFreedom(f), m_Solver.TotalSolutionIndex);
 }
      }
      It.Set(Sol);
    }
  }
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void 
FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::GetBSplineNodalValuesForNextLevel()
{
  unsigned int order = ((ImageDimension == 2) 
                     ? dynamic_cast<Element2DBSplinePatch*>(m_Element)->GetBSplineOrder()
                     : dynamic_cast<Element3DBSplinePatch*>(m_Element)->GetBSplineOrder()) - 1;

  typename DeformationFieldType::RegionType::SizeType size;
  ArrayType array;

  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    size[i] = m_MeshResolution->ElementAt(m_CurrentLevel)[i] + order;
   
    array[i] = (m_PyramidLevelImageSizes[m_CurrentLevel][i] == 
                m_PyramidLevelImageSizes[m_CurrentLevel+1][i])
                ? m_MeshResolution->ElementAt(m_CurrentLevel)[i]
                : 2*m_MeshResolution->ElementAt(m_CurrentLevel)[i];
  }
  this->SetMeshResolution(array, m_CurrentLevel+1);

  VectorType V(0.0);
  m_BSplineNodalValues = DeformationFieldType::New();
  m_BSplineNodalValues->SetRegions(size);
  m_BSplineNodalValues->Allocate();
  m_BSplineNodalValues->FillBuffer(V);  

  DeformationFieldIteratorType It(m_BSplineNodalValues, m_BSplineNodalValues->GetBufferedRegion());

  Node::ArrayType* nodes = &(m_Solver.node);
  Node::ArrayType::iterator it; 
 
  for (it = nodes->begin(), It.GoToBegin(); it != nodes->end(); it++ ,++It) 
  {
    for (unsigned int i = 0; i < ImageDimension; i++)
    { 
      V[i] = m_Solver.GetLinearSystemWrapper()->
             GetSolutionValue((*it)->GetDegreeOfFreedom(i), m_Solver.TotalSolutionIndex);    
    } 
    It.Set(V);
  }
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void
FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::DoubleResolutionOfFEMBSplineMesh() 
{
  unsigned int order = ((ImageDimension == 2) 
                     ? dynamic_cast<Element2DBSplinePatch*>(m_Element)->GetBSplineOrder()
                     : dynamic_cast<Element3DBSplinePatch*>(m_Element)->GetBSplineOrder()) - 1;

  ArrayType NumberOfControlPoints = m_MeshResolution->ElementAt(m_CurrentLevel);

  ImageSizeType size;
  bool expandDimension[ImageDimension];
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    size[i] = m_MeshResolution->ElementAt(m_CurrentLevel)[i] + order;
    expandDimension[i] = (m_PyramidLevelImageSizes[m_CurrentLevel][i] == 
                          m_PyramidLevelImageSizes[m_CurrentLevel-1][i]) 
                       ? false : true;
  }

  DeformationFieldPointer RefinedLattice = DeformationFieldType::New();
  typename DeformationFieldType::PixelType V(0.0);
  RefinedLattice->SetRegions(size);
  RefinedLattice->Allocate();
  RefinedLattice->FillBuffer(V);  

  DeformationFieldIteratorType It(RefinedLattice, RefinedLattice->GetBufferedRegion());
  typename DeformationFieldType::IndexType idx, idx_Psi, tmp, tmp_Psi;
  typename DeformationFieldType::IndexType off, off_Psi;
  typename DeformationFieldType::RegionType::SizeType size_Psi;

  size.Fill(2);
  size_Psi.Fill(order+1);  

  typedef BSplineKernelFunction<3> KernelType;  
  typename KernelType::Pointer Kernel = KernelType::New();
  Kernel->SetSplineOrder(order);
  typename KernelType::MatrixType R, C, Coefficients;

  C = Kernel->GetShapeFunctionsInZeroToOneInterval();
  R = C;
  for (int i = 0; i < C.cols(); i++)
  {
    RealType c = pow(2.0, static_cast<int>(C.cols()-i-1));
    for (int j = 0; j < C.rows(); j++)
    {
      R(j, i) *= c;
    }
  }  
  R = R.transpose();  R.flipud();
  C = C.transpose();  C.flipud();
  Coefficients = vnl_matrix_inverse<RealType>(R)*C;

  It.GoToBegin();
  while (!It.IsAtEnd())
  {
    idx = It.GetIndex();
    for (unsigned int i = 0; i < ImageDimension; i++)
    {
      idx_Psi[i] = (expandDimension[i]) 
                 ? static_cast<unsigned int>(0.5*idx[i])
                 : static_cast<unsigned int>(idx[i]);    
    }
    
    for (unsigned int i = 0; i < pow(2, ImageDimension); i++)
    {
      VectorType Sum(0.0);
      off = this->IndexToSubscript(i, size);
      
      bool OutOfBoundary = false;
      for (unsigned int j = 0; j < ImageDimension; j++)
      {
        tmp[j] = idx[j] + off[j];
     if (tmp[j] >= NumberOfControlPoints[j])
     {
   OutOfBoundary = true;
          break;
     }
      } 
      if (OutOfBoundary)
      {
        continue;
      }      
      
      int N = static_cast<int>(pow(static_cast<int>(order+1), 
                                   static_cast<int>(ImageDimension)));
      for (unsigned int j = 0; j < N; j++)
      {
        off_Psi = this->IndexToSubscript(j, size_Psi);
        bool OutOfBoundary = false;
 for (unsigned int k = 0; k < ImageDimension; k++)
 {
   tmp_Psi[k] = idx_Psi[k] + off_Psi[k];
   if (tmp_Psi[k] >= m_MeshResolution->ElementAt(m_CurrentLevel)[k])
   {
     OutOfBoundary = true;
            break;
   }
 } 
 if (OutOfBoundary)
 {
          continue;
  }      

        RealType coeff = 1.0;
        for (unsigned int k = 0; k < ImageDimension; k++)
 {
          coeff *= Coefficients(off[k], off_Psi[k]);
 }  
        V = m_BSplineNodalValues->GetPixel(tmp_Psi);
        V *= coeff;
        Sum += V;
      }
      RefinedLattice->SetPixel(tmp, Sum);
    }  

    bool IsEvenIndex = false;
    while (!IsEvenIndex && !It.IsAtEnd())
    {      
      ++It;  
      idx = It.GetIndex();
      IsEvenIndex = true;
      for (unsigned int i = 0; i < ImageDimension; i++)
      {
        if (idx[i] % 2 != 0)
        {
          IsEvenIndex = false;
        }
      }
    }
  }
  Node::ArrayType* nodes = &(m_Solver.node);
  Node::ArrayType::iterator it; 
 
  for (it = nodes->begin(), It.GoToBegin(); it != nodes->end(); it++, ++It) 
  {
    for (unsigned int i = 0; i < ImageDimension; i++)
    { 
      m_Solver.GetLinearSystemWrapper()->
             SetSolutionValue((*it)->GetDegreeOfFreedom(i), It.Get()[i], m_Solver.TotalSolutionIndex);    
      m_Solver.GetLinearSystemWrapper()->
             SetSolutionValue((*it)->GetDegreeOfFreedom(i), It.Get()[i], m_Solver.SolutionTMinus1Index);  
    } 
  }
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::PrintSelf(std::ostream& os, Indent indent) const 
{ 
  Superclass::PrintSelf( os, indent );

  (!m_LandmarkFileName.empty()) ?
       std::cout << "   landmark file:    " << m_LandmarkFileName << std::endl
     : std::cout << "   No landmark file specified. " << std::endl;
  (!m_MeshFileName.empty()) ?
       std::cout << "   mesh file:    " << m_MeshFileName << std::endl
     : std::cout << "   No mesh file specified. " << std::endl;

  std::cout << "Variables: " << std::endl;
  std::cout << "   Alpha = " << m_Alpha << std::endl;
  std::cout << "   Energy reduction factor = " << m_EnergyReductionFactor << std::endl;
  std::cout << "   Time Step = " << m_TimeStep << std::endl;
  std::cout << "   Elasticity = " << m_Elasticity << std::endl;
  std::cout << "   Rho = " << m_Rho << std::endl;
  std::cout << "   Gamma = " << m_Gamma << std::endl;
  std::cout << "   Maximum iterations = " << m_MaximumNumberOfIterations << std::endl;
  std::cout << "   Number of integration points = " << m_NumberOfIntegrationPoints << std::endl;
  std::cout << "   Number of levels = " << m_NumberOfLevels << std::endl;
  std::cout << "   Do line search on image energy = " << m_DoLineSearchOnImageEnergy << std::endl;
  if (m_DoLineSearchOnImageEnergy)
  {
    std::cout << "     Line search maximum iterations = " << m_LineSearchMaximumIterations << std::endl; 
  }  
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
typename FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>::RealType 
FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::EvaluateResidual(RealType t)
{
  RealType SimE = (m_UseImageToImageMetric)
                ? m_ImageToImageMetricLoad->EvaluateMetricGivenSolution(&(m_Solver.el), t)
                : m_PDEDeformableMetricLoad->EvaluateMetricGivenSolution(&(m_Solver.el), t);
  return vnl_math_abs(static_cast<RealType>(SimE));
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void 
FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::FindBracketingTriplet(RealType* a, RealType* b, RealType* c)
{
  const RealType Gold = 1.618034;
  const RealType Glimit = 100.0;
  const RealType Tiny = 1.e-20;
  RealType ax = 0.0;
  RealType bx = 1.0;
  RealType fa = EvaluateResidual(ax);
  RealType fb = EvaluateResidual(bx);

  RealType dum;
  if ( fb > fa )
  {
    dum = ax; ax = bx; bx = dum;
    dum = fb; fb = fa; fa = dum;
  }

  RealType cx = bx+Gold*(bx-ax);  // first guess for c - the 3rd pt needed to bracket the min
  RealType fc = vnl_math_abs(EvaluateResidual(cx));

  RealType ulim,u,r,q,fu;
  while (fb > fc)
  {
    r = (bx-ax)*(fb-fc);
    q = (bx-cx)*(fb-fa);
    RealType denom = (2.0*m_Solver.GSSign(m_Solver.GSMax(vnl_math_abs(q-r),Tiny),q-r));
    u = (bx)-((bx-cx)*q-(bx-ax)*r)/denom;
    ulim = bx + Glimit*(cx-bx);
    if ((bx-u)*(u-cx) > 0.0)
    {
      fu = EvaluateResidual(u);
      if (fu < fc)
      {
        ax = bx;
        bx = u;
        *a = ax; 
 *b = bx; 
 *c = cx;
        return;
      }
      else if (fu > fb)
      {
        cx = u;
        *a = ax; 
 *b = bx; 
 *c = cx;
        return;
      }

      u = cx+Gold*(cx-bx);
      fu = EvaluateResidual(u);

    }
    else if ( (cx-u)*(u-ulim) > 0.0)
    {
      fu = EvaluateResidual(u);
      if (fu < fc)
      {
        bx = cx; 
 cx = u; 
 u = cx+Gold*(cx-bx);
        fb = fc; 
 fc = fu; 
 fu = EvaluateResidual(u);
      }

    }
    else if ( (u-ulim)*(ulim-cx) >=  0.0)
    {
      u = ulim;
      fu = EvaluateResidual(u);
    }
    else
    {
      u = cx+Gold*(cx-bx);
      fu = EvaluateResidual(u);
    }

    ax = bx; 
    bx = cx; 
    cx = u;
    fa = fb; 
    fb = fc; 
    fc = fu;

  }

  if ( vnl_math_abs(ax) > 1.e3  || vnl_math_abs(bx) > 1.e3 || vnl_math_abs(cx) > 1.e3)
  {
    ax = -2.0;  
    bx = 1.0;  
    cx = 2.0;
  }               

  *a = ax; 
  *b = bx; 
  *c = cx;
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
typename FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>::RealType 
FEMRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::GoldenSection(unsigned int MaxIters)
{
  // We should now have a, b and c, as well as f(a), f(b), f(c), 
  // where b gives the minimum energy position;
  RealType ax, bx, cx;
  FindBracketingTriplet(&ax, &bx, &cx);

  const RealType R = 0.6180339;
  const RealType C = (1.0 - R);
  const RealType tol = 1.0;

  RealType x0 = ax;
  RealType x1;
  RealType x2;
  RealType x3 = cx;
  if (vnl_math_abs(cx-bx) > vnl_math_abs(bx-ax))
  {
    x1 = bx;
    x2 = bx+C*(cx-bx);
  }
  else
  {
    x2 = bx;
    x1 = bx-C*(bx-ax);
  }

  RealType f1 = EvaluateResidual(x1);
  RealType f2 = EvaluateResidual(x2);
  unsigned int iters=0;
  while (vnl_math_abs(x3-x0) > tol*(vnl_math_abs(x1)+vnl_math_abs(x2)) && iters++ < MaxIters)
  {
    if (f2 < f1)
    {
      x0 = x1; x1 = x2; x2=R*x1+C*x3;
      f1 = f2; f2 = EvaluateResidual(x2);
    } 
    else
    {
      x3 = x2; 
      x2 = x1; 
      x1 = R*x2+C*x0;
      f2 = f1; 
      f1 = EvaluateResidual(x1);
    }
  }
  RealType xmin,fmin;
  if (f1<f2)
  {
    xmin=x1;
    fmin=f1;
  }
  else
  {
    xmin=x2;
    fmin=f2;
   }

  m_Solver.SetEnergyToMin(xmin);
  return vnl_math_abs(fmin);
}


}} // end namespace itk::fem

#endif
