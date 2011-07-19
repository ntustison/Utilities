/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBoykovAlphaExpansionMRFImageFilter.h,v $
  Language:  C++
  Date:      
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBoykovAlphaExpansionMRFImageFilter_h
#define __itkBoykovAlphaExpansionMRFImageFilter_h

#include "itkGraph.h"
#include "itkImage.h"
#include "itkMRFImageFilter.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkBoykovGraphTraits.h"
#include "itkBoykovImageToGraphFunctor.h"

#include <vector>

namespace itk {

/** \class itkBoykovAlphaExpansionMRFImageFilter
 * \brief Classes an input image based on graph-based energy minimization.
 * 
 * \par
 * This class implements the alpha-expansion algorithm of Boykov et al. 
 * given in the reference below.  Applications include vision problems
 * for which there are large numbers of labels, e.g. stereo, motion, image
 * restoration, segmentation, and scene reconstruction.  
 *
 * \par
 * This method seeks to minimize the energy function derived as a
 * posterior probability in a Markov Random Field (MRF).  While 
 * finding the global minimum is NP-hard, the alpha-expansion algorithm
 * is guaranteed to find a local minimum that is within a certain
 * factor of the global minimum.  
 * 
 * \par INPUTS
 * Suppose that a particular image has N labels.  The input consists
 * partly of N likelihood images of the same dimension and size 
 * as the input image.  These can be derived from non-parametric 
 * techniques such as Parzen windowing, etc.  Optional: N lists of 
 * pixel indices can be specified to "hard-constrain" those pixels 
 * to be of a specific labeling.
 *
 * \par REFERENCE
 * Y. Boykov, O. Veksler, and R. Zabih, "Fast Approximate Energy 
 * Minimization via Graph Cuts," IEEE-PAMI, 23(11):1222-1239, 2001.
 *
 **/
  
template<typename TInputImage, 
         typename TGraphTraits,
         typename TClassifiedImage = 
           Image<int, GetImageDimension<TInputImage>::ImageDimension> >
class ITK_EXPORT BoykovAlphaExpansionMRFImageFilter 
: public MRFImageFilter<TInputImage, TClassifiedImage>
{
public:
  /** Extract dimension from input and output image. */

  /** Template classes typedef */
  typedef TInputImage                                       InputImageType;
  typedef TGraphTraits                                      GraphTraitsType;  
  typedef TClassifiedImage                                  OutputImageType;  

  /** Standard class typedefs */
  typedef BoykovAlphaExpansionMRFImageFilter                Self;
  typedef MRFImageFilter<InputImageType, 
                         OutputImageType>                   Superclass;
  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Image dimension */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage::ImageDimension );

  /** Image typedefs */
  typedef typename InputImageType::PixelType                InputPixelType;  
  typedef typename InputImageType::IndexType                IndexType;  
  typedef typename OutputImageType::PixelType               OutputPixelType;  

  /** Graph typedefs */
  typedef typename GraphTraitsType::NodePointerType         NodePointerType;
  typedef typename GraphTraitsType::EdgePointerType         EdgePointerType;
  typedef typename GraphTraitsType::NodeIdentifierType      NodeIdentifierType;
  typedef typename GraphTraitsType::EdgeIdentifierType      EdgeIdentifierType;
  typedef typename GraphTraitsType::NodeWeightType          NodeWeightType;
  typedef typename GraphTraitsType::EdgeWeightType          EdgeWeightType;
  typedef Graph<GraphTraitsType>                            GraphType;  
  typedef typename GraphType::NodeIterator                  NodeIteratorType;
  typedef typename GraphType::EdgeIterator                  EdgeIteratorType;

  /** Other related image typedefs */
  typedef ShapedNeighborhoodIterator<OutputImageType>  NeighborhoodIteratorType;

  /** Other typedefs */
  typedef double                                   RealType;
  typedef std::vector<IndexType>                   IndexContainerType;
  typedef std::vector<IndexContainerType>          IndexContainerContainerType;
  typedef Image<RealType, ImageDimension>          ProbabilityImageType;

 
  /** Image To Graph Functor Type */
  typedef BoykovImageToGraphFunctor
    <InputImageType, GraphType>               BoykovImageToGraphFunctorType;
  typedef typename 
    BoykovImageToGraphFunctorType::Pointer    BoykovImageToGraphFunctorPointer;

  /** Set/Get Functions */
  
  /** Set/Get BoykovImageToGraphFunctor */
  itkSetObjectMacro( ImageToGraphFunctor, BoykovImageToGraphFunctorType )
  itkGetObjectMacro( ImageToGraphFunctor, BoykovImageToGraphFunctorType )

  itkGetMacro(RandomizeInitialLabeling, bool);
  itkSetMacro(RandomizeInitialLabeling, bool);
  

  void SetLikelihoodImage( unsigned int i, ProbabilityImageType *image )
    { this->SetNthInput( i, const_cast<ProbabilityImageType *>( image ) ); }
  void SetIndexContainer( unsigned int i, IndexContainerType indices )
    { this->m_Indices[i] = indices; }
    
  ProbabilityImageType* GetPosteriorProbabilityImage( unsigned int i )
    { return dynamic_cast<ProbabilityImageType *>(
        this->ProcessObject::GetOutput( i ) ); }

protected:
  BoykovAlphaExpansionMRFImageFilter();
  virtual ~BoykovAlphaExpansionMRFImageFilter() {}
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

private:
  BoykovAlphaExpansionMRFImageFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  virtual void ApplyMRFImageFilter();  

  /** Private functions for applying the MRF Image filter */

  void AlphaExpansion();
  //void AlphaBetaSwap();

  void Initialize();
  RealType CalculateCurrentEnergy();
  void FindMinimumEnergyBinaryLabeling();
  void FindMinimumEnergyLabeling( unsigned int );
  void Relabel();
  
   /** 
    * The posterior probability at a certain pixel 'p' with intensity
    * Ip, i.e. Pr(label i|Ip), is proportional to exp[-U(f_i)] where 
    * U(f_i) is the energy associated with label i at pixel p.  We simply
    * generate the energy images and leave it to the user to formulate
    * a desired posterior probability. 
    */
  void GeneratePosteriorEnergyImages();
  
  EdgeWeightType CalculateSmoothnessPenaltyTerm( 
    IndexType, IndexType, unsigned int, unsigned int );  
  void AddUnaryTerm(
    GraphType *, NodePointerType, NodeWeightType, NodeWeightType );  
  void AddBinaryTerm( GraphType *, EdgePointerType, 
    EdgeWeightType, EdgeWeightType, EdgeWeightType, EdgeWeightType );  

  /** private data members */

  typename OutputImageType::Pointer   m_LabelImage;  
  bool                                m_RandomizeInitialLabeling;
  BoykovImageToGraphFunctorPointer    m_ImageToGraphFunctor;
  IndexContainerContainerType         m_Indices;
  RealType                            m_CurrentEnergy;

  /** 
   * private data members in base MRF class. 
   * Should be 'protected' in MRF base class? 
   */
  unsigned int m_NumberOfIterations;
  unsigned int m_MaximumNumberOfIterations;

} ; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBoykovAlphaExpansionMRFImageFilter.hxx"
#endif

#endif


