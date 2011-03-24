/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMRegistrationFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:57 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef _itkFEMRegistrationFilter_h_
#define _itkFEMRegistrationFilter_h_

#include "itkFEMLinearSystemWrapperItpack.h"
#include "itkFEMSolverCrankNicolson.h"
#include "itkFEMMaterialLinearElasticity.h"
#include "itkFEMImageMetricLoad.h"
#include "itkFEMFiniteDifferenceFunctionLoad.h"
#include "itkFEMLoadLandmark.h"

#include "itkVector.h"
#include "itkVectorContainer.h"
#include "itkFixedArray.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkWarpImageFilter.h"
 
#include "vnl/vnl_vector.h"

#include <iostream>
#include <string>

/** Set/Get built-in type to handle elements of arrays; */
#define itkSetElementMacro(name,type) \
  void Set##name (const type _arg, const unsigned int _i) \
  { \
    itkDebugMacro("setting element " << _i << " of " #name " to " << _arg); \
    if (this->m_##name[_i] != _arg) \
      { \
      this->m_##name[_i] = _arg; \
      this->Modified(); \
      } \
  } 
#define itkSetAllElementsMacro(name,type) \
  void Set##name (const type _arg) \
  { \
    itkDebugMacro("setting the entire vector to " << _arg); \
    this->m_##name.fill(_arg); \
    this->Modified(); \
  } 
#define itkGetElementConstMacro(name,type) \
  virtual type Get##name (const unsigned int _i = 0) const \
  { \
    itkDebugMacro("returning element " << _i << " of " << #name " of " << this->m_##name ); \
    return this->m_##name[_i]; \
  }

namespace itk {
namespace fem {

/** \class FEMRegistrationFilter 
    \brief FEM Image registration filter.

     The image registration problem is modeled here with the finite element method.
     Image registration is, in general, an ill-posed problem.  Thus, we use an optimization
     scheme where the optimization criterion is given by a regularized variational energy.
     The variational energy arises from modeling the image as a physical body on which 
     external forces act.  The body is allowed to deform in response to the applied force.  
     The resistance of the physical body to deformation, determined by the physics 
     associated with the body, serves to regularize the solution.  The forces applied 
     to the body are, generally, highly non-linear and so the body is allowed to deform 
     slowly and incrementally.  The direction it deforms follows the gradient of the 
     potential energy (the force) we define.  The potential energies we may choose from 
     are given by the itk image-to-image metrics or those derived from the PDE deformable
     registration functions.  

     \par   
     The forces driving the problem may also be given by user-supplied landmarks.  
     The corners of the image, in this example, are always pinned.  This example is 
     designed for 2D or 3D images.  A rectilinear mesh is generated automatically 
     given the correct element type (Quadrilateral or Hexahedral).

     \par 
     Our specific Solver for this example uses trapezoidal time stepping.  This is 
     a method for solving a second-order PDE in time.  The solution is penalized 
     by the zeroth (mass matrix) and first derivatives (stiffness matrix) of the 
     shape functions.  There is an option to perform a line search on the energy 
     after each iteration.  Optimal parameter settings require experimentation.
     The following approach tends to work well:  
        -> Choose the relative size of density  to elasticity (e.g. Rho / E ~= 1.)
           such that the image deforms locally and slowly.  This also affects the 
    stability of the solution.
        -> Choose the time step to control the size of the deformation at each step.
        -> Choose enough iterations to allow the solution to converge (this may be 
    automated).

     \par
     Reading images is up to the user.  Either set the images using 
     SetMoving/FixedImage or see the ReadImages function.  
     
     \note This code works for only 2 or 3 dimensions b/c we do not have > 3D elements.

  \note TODO :  Keep the full field around (if using re-gridding).
                Introduce compensation for kinematic non-linearity in time (if using Eulerian frame).
*/

template<class TMovingImage, class TFixedImage, class TWarpedImage = TFixedImage> 
class ITK_EXPORT FEMRegistrationFilter 
: public ImageToImageFilter<TFixedImage, TWarpedImage>
{
public:
  typedef FEMRegistrationFilter                              Self;
  typedef ImageToImageFilter<TMovingImage, TWarpedImage>     Superclass;
  typedef SmartPointer<Self>                                 Pointer;
  typedef SmartPointer<const Self>                           ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  
  
  /** Run-time type information (and related methods) */
  itkTypeMacro(FEMRegistrationFilter, ImageToImageFilter);
  
  typedef TMovingImage                                       MovingImageType;
  typedef TFixedImage                                        FixedImageType;
  typedef TWarpedImage                                       WarpedImageType;
  typedef typename FixedImageType::PixelType                 PixelType;
  typedef typename FixedImageType::SizeType                  ImageSizeType;
  
  /** Dimensionality of input and output data is assumed to be the same. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      FixedImageType::ImageDimension);

  typedef double                                             RealType;
  typedef Image<RealType, 
          itkGetStaticConstMacro(ImageDimension)>            RealImageType;    
  typedef FixedArray<unsigned int, 
          itkGetStaticConstMacro(ImageDimension)>            ArrayType;
  typedef VectorContainer<unsigned, ArrayType>               ArrayContainerType;     
   
  /** Typedefs for the deformation field */   
  typedef Vector<RealType, 
          itkGetStaticConstMacro(ImageDimension)>            VectorType;
  typedef Image<VectorType, 
          itkGetStaticConstMacro(ImageDimension)>            DeformationFieldType;
  typedef typename DeformationFieldType::Pointer             DeformationFieldPointer;
  typedef ImageRegionIteratorWithIndex<DeformationFieldType> DeformationFieldIteratorType; 

  typedef LinearSystemWrapperItpack                          LinearSystemSolverType;
  typedef SolverCrankNicolson                                SolverType;
  typedef MaterialLinearElasticity                           MaterialType;
  typedef std::vector<typename LoadLandmark::Pointer>        LandmarkArrayType;

  /** Typedef support for the interpolation function */
  typedef InterpolateImageFunction<MovingImageType, 
                                   RealType>                 ImageInterpolatorType;
  typedef typename ImageInterpolatorType::Pointer            ImageInterpolatorTypePointer; 
       
  /** Typedefs for image metrics (PDEDeformableRegistrationFunction) */       
  typedef PDEDeformableRegistrationFunction
     <FixedImageType, MovingImageType, DeformationFieldType> PDEDeformableMetricType; 
  typedef typename PDEDeformableMetricType::Pointer          PDEDeformableMetricTypePointer;
  typedef FiniteDifferenceFunctionLoad<MovingImageType, 
                                       FixedImageType>       PDEDeformableMetricLoadType;
  typedef typename PDEDeformableMetricLoadType::Pointer      PDEDeformableMetricLoadTypePointer;

  /** Typedefs for image metrics (ImageToImageMetric) */       
  typedef ImageToImageMetric<FixedImageType, 
                             MovingImageType>                ImageToImageMetricType; 
  typedef typename ImageToImageMetricType::Pointer           ImageToImageMetricTypePointer;
  typedef ImageMetricLoad<MovingImageType, FixedImageType>   ImageToImageMetricLoadType;
  typedef typename ImageToImageMetricLoadType::Pointer       ImageToImageMetricLoadTypePointer;

 
  /** Main functions */
  
  void RunRegistration()
    { this->Update(); }

  void IterativeSolve();  
  void RegisterImages();

  /** Set/Get functions (variables are explained below). */
  itkSetInputMacro(MovingImage, MovingImageType, 0); 
  itkSetInputMacro(FixedImage, FixedImageType, 1); 
  itkGetConstMacro(WarpedImage, typename WarpedImageType::Pointer);  
  itkGetConstMacro(DeformationField, typename DeformationFieldType::Pointer);

  itkSetStringMacro(LandmarkFileName); 
  itkGetStringMacro(LandmarkFileName);
  itkSetStringMacro(MeshFileName); 
  itkGetStringMacro(MeshFileName);

  itkSetMacro(Alpha, RealType);
  itkGetConstMacro(Alpha, RealType);
  itkSetMacro(TimeStep, RealType);
  itkGetConstMacro(TimeStep, RealType);
  itkSetMacro(EnergyReductionFactor, RealType);
  itkGetConstMacro(EnergyReductionFactor, RealType);
  itkSetMacro(Material, typename MaterialType::Pointer);
  itkGetConstMacro(Material, typename MaterialType::Pointer);
  itkSetMacro(PDEDeformableMetricLoad, PDEDeformableMetricLoadTypePointer);
  itkGetConstMacro(PDEDeformableMetricLoad, PDEDeformableMetricLoadTypePointer);
  itkSetMacro(ImageToImageMetricLoad, ImageToImageMetricLoadTypePointer);
  itkGetConstMacro(ImageToImageMetricLoad, ImageToImageMetricLoadTypePointer);
  void SetElement(typename Element::Pointer element)
    {
      m_Element = element; 
      m_UseBSplines = 
        (dynamic_cast<Element2DBSplinePatch*>(m_Element) != NULL ||
         dynamic_cast<Element3DBSplinePatch*>(m_Element) != NULL) ?
         true : false;
    }
  itkGetConstMacro(Element, typename Element::Pointer);

  itkSetElementMacro(Elasticity, RealType);
  itkSetAllElementsMacro(Elasticity, RealType);
  itkGetElementConstMacro(Elasticity, RealType);
  itkSetElementMacro(Rho, RealType);
  itkSetAllElementsMacro(Rho, RealType);
  itkGetElementConstMacro(Rho, RealType);
  itkSetElementMacro(Gamma, RealType);
  itkSetAllElementsMacro(Gamma, RealType);
  itkGetElementConstMacro(Gamma, RealType);
  itkSetElementMacro(MaximumNumberOfIterations, unsigned int);
  itkSetAllElementsMacro(MaximumNumberOfIterations, unsigned int);
  itkGetElementConstMacro(MaximumNumberOfIterations, unsigned int);
  itkSetElementMacro(NumberOfIntegrationPoints, unsigned int);
  itkSetAllElementsMacro(NumberOfIntegrationPoints, unsigned int);
  itkGetElementConstMacro(NumberOfIntegrationPoints, unsigned int);
  itkSetElementMacro(MetricRegionWidth, unsigned int);
  itkSetAllElementsMacro(MetricRegionWidth, unsigned int);
  itkGetElementConstMacro(MetricRegionWidth, unsigned int);

  void SetMeshResolution(ArrayType array, unsigned int i)
    {
      m_MeshResolution->InsertElement(i, array);  
    } 
  void SetMeshResolution(unsigned int N, unsigned int i)
    {
      ArrayType array;
      array.Fill(N);
      m_MeshResolution->InsertElement(i, array);  
    } 
  void SetMeshResolution(unsigned int N)
    {
      ArrayType array;
      array.Fill(N);
      for (unsigned int i = 0; i < itkGetStaticConstMacro(ImageDimension); i++)
      {   
        m_MeshResolution->InsertElement(i, array);  
      }
    } 

  void SetNumberOfLevels(unsigned int);
  itkGetConstMacro(NumberOfLevels, unsigned int);

  itkSetMacro(DoLineSearchOnImageEnergy, bool);
  itkGetConstMacro(DoLineSearchOnImageEnergy, bool);
  itkSetMacro(UseMassMatrix, bool);
  itkGetConstMacro(UseMassMatrix, bool);
  itkSetMacro(LineSearchMaximumIterations, unsigned int);
  itkGetConstMacro(LineSearchMaximumIterations, unsigned int);
  itkSetMacro(MaximizeMetric, bool);
  itkGetConstMacro(MaximizeMetric, bool);

  /** Set/Get the image interpolator. */ 
  itkSetObjectMacro(ImageInterpolator, ImageInterpolatorType); 
  itkGetObjectMacro(ImageInterpolator, ImageInterpolatorType);

  /** Set/Get the image metric.  (PDEDeformableRegistrationFunction) */ 
  void SetPDEDeformableMetric(PDEDeformableMetricType* M)
    {
      m_PDEDeformableMetric = M;
      m_UseImageToImageMetric = false;
    } 
  itkGetObjectMacro(PDEDeformableMetric, PDEDeformableMetricType);

  /** Set/Get the image metric. (ImageToImageMetric) */ 
  void SetImageToImageMetric(ImageToImageMetricType* M)
    {
      m_ImageToImageMetric = M; 
      m_UseImageToImageMetric = true;
    } 
  itkGetObjectMacro(ImageToImageMetric, ImageToImageMetricType);

protected :
  /** de/constructor */
  FEMRegistrationFilter( ); 
  ~FEMRegistrationFilter() {} 

  void PrintSelf(std::ostream& os, Indent indent) const;
  
  void GenerateData();

  /** Easy access to the FEMObjectFactory (short and not templated) */
  class FEMOF 
  : public FEMObjectFactory<FEMLightObject>
  {
    protected:
      FEMOF();
      ~FEMOF();
  };

private : 
  
  FEMRegistrationFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
 
  void Initialize();
 
  /**
   * If an external mesh is not specified, the FEM filter can generate its 
   * own rectilinear mesh for 2 or 3 dimensions. The mesh can be generated 
   * for linear elements (quadrilaterals or hexahedra) or B-spline patches.  
   * To generate 10 elements per dimension in the 1st resolution, 
   * call the function this->SetMeshResolution(10, 0).
   */
  void GenerateFEMMesh(ArrayType);

  /** This is used for changing between mesh resolutions. */
  void SampleVectorFieldAtNodes();

  itkSetMacro(DeformationField, typename DeformationFieldType::Pointer);

  void ApplyLoads(RealType* spacing = NULL); 
  void ApplyImageLoads(MovingImageType*, FixedImageType*); 
  
  /** Evaluates the image similarity energy by calling the image metric */
  RealType EvaluateEnergy();
 
  /** Interpolates the vector field over the domain.  
    * Our convention is to always keep the vector field
    * at the scale of the original images.
    */
  void InterpolateVectorField(); 

  void PropogateSolutionToTheNextLevel();
  void DoubleResolutionOfFEMBSplineMesh();
  void GetBSplineNodalValuesForNextLevel();

  inline typename RealImageType::IndexType 
  IndexToSubscript(unsigned int index, typename RealImageType::SizeType size)
    {
      typename RealImageType::IndexType k;     k[0] = 1;    
      for (unsigned int i = 1; i < ImageDimension; i++)
      {
        k[i] = size[i]*k[i-1];
      }  
      typename RealImageType::IndexType sub;
      for (unsigned int i = 0; i < ImageDimension; i++)
      {
        sub[ImageDimension-i-1] = static_cast<unsigned int>(index/k[ImageDimension-i-1]);
        index %= k[ImageDimension-i-1];
      }
      return sub;
    }
  
  RealType EvaluateResidual(RealType);
  void FindBracketingTriplet(RealType*, RealType*, RealType*);
  RealType GoldenSection(unsigned int MaxIters = 25);

private :

  /** Filter Outputs */  
  typename WarpedImageType::Pointer      m_WarpedImage;  

  DeformationFieldPointer                m_DeformationField;
  DeformationFieldPointer                m_BSplineNodalValues;
 
  /** File names for input/ouput */
  std::string                            m_LandmarkFileName;
  std::string                            m_DisplacementsFileName;
  std::string                            m_MeshFileName;
  
  /** Set/get variables */
  unsigned int                           m_NumberOfLevels;                        // Number of resolution levels
  unsigned int                           m_LineSearchMaximumIterations;           // 
  bool                                   m_DoLineSearchOnImageEnergy;             //
  bool                                   m_UseMassMatrix;            
  bool                                   m_UseBSplines; 

  RealType                               m_Alpha;                                 // Solver parameter 
  RealType                               m_EnergyReductionFactor;     // Convergence threshold 
  RealType                               m_TimeStep;                              // Time step
  typename Element::Pointer              m_Element;                               //
  typename MaterialType::Pointer         m_Material;                              // 

  vnl_vector<RealType>                   m_Elasticity;                            // elasticity 
  vnl_vector<RealType>                   m_Rho;                                   // mass matrix weight
  vnl_vector<RealType>                   m_Gamma;                                 // image similarity weight
  vnl_vector<unsigned int>               m_MaximumNumberOfIterations;             // 
  vnl_vector<unsigned int>               m_NumberOfIntegrationPoints;             // resolution of integration
  vnl_vector<unsigned int>               m_MetricRegionWidth;                     //  
  typename ArrayContainerType::Pointer   m_MeshResolution;

  unsigned int                           m_CurrentLevel;
  std::vector<ImageSizeType>             m_PyramidLevelImageSizes;
  SolverType                             m_Solver;                                // Defines the solver (currently Crank-Nicolson)  
  LandmarkArrayType                      m_LandmarkArray;                         // Contains pointers to the landmark load pointers 

  bool                                   m_MaximizeMetric;
  bool                                   m_UseImageToImageMetric;  
  ImageToImageMetricTypePointer          m_ImageToImageMetric;
  ImageToImageMetricLoadTypePointer      m_ImageToImageMetricLoad;                // Defines the load to use
  ImageInterpolatorTypePointer           m_ImageInterpolator;                    
  PDEDeformableMetricTypePointer         m_PDEDeformableMetric;
  PDEDeformableMetricLoadTypePointer     m_PDEDeformableMetricLoad;               // Defines the load to use

};

}} // end namespace fem

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFEMRegistrationFilter.txx"
#endif

#endif

