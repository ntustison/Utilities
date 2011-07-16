/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBoykovAlphaExpansionMRFImageFilter.hxx,v $
  Language:  C++
  Date:      $Date: 2008/11/11 03:09:13 $
  Version:   $ $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkBoykovAlphaExpansionMRFImageFilter_hxx
#define _itkBoykovAlphaExpansionMRFImageFilter_hxx

#include "itkBoykovAlphaExpansionMRFImageFilter.h"
#include "itkBoykovMinCutGraphFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageToGraphFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

#include "vnl/vnl_math.h"

namespace itk
{

template<typename TInputImage, typename TGraphTraits, typename TClassifiedImage>
BoykovAlphaExpansionMRFImageFilter<TInputImage, TGraphTraits, TClassifiedImage>
::BoykovAlphaExpansionMRFImageFilter()
{
  this->SetNumberOfClasses( 2 ); 
  this->m_RandomizeInitialLabeling = false;
  this->m_MaximumNumberOfIterations = 10;
  this->m_Indices.clear();
}

template<typename TInputImage, typename TGraphTraits, typename TClassifiedImage>
void
BoykovAlphaExpansionMRFImageFilter<TInputImage, TGraphTraits, TClassifiedImage>
::GenerateData()
{
  this->ApplyMRFImageFilter();
  this->GeneratePosteriorEnergyImages();
} 

template<typename TInputImage, typename TGraphTraits, typename TClassifiedImage>
void
BoykovAlphaExpansionMRFImageFilter<TInputImage, TGraphTraits, TClassifiedImage>
::ApplyMRFImageFilter()
{
  this->AlphaExpansion();
//  this->AlphaBetaSwap();
}

template<typename TInputImage, typename TGraphTraits, typename TClassifiedImage>
void
BoykovAlphaExpansionMRFImageFilter<TInputImage, TGraphTraits, TClassifiedImage>
::AlphaExpansion()
{
  this->Initialize();
  
  /** special case for binary labeling */
  if( this->GetNumberOfClasses() == 2 )
    {
    this->FindMinimumEnergyBinaryLabeling();
    this->Relabel();
    return;
    }
  this->m_CurrentEnergy = this->CalculateCurrentEnergy();    
  RealType E = 2.0*m_CurrentEnergy;

  itkDebugMacro( "Initial Energy = " << this->m_CurrentEnergy );

  this->m_NumberOfIterations = 0;
  while( this->m_CurrentEnergy < E && 
    this->m_NumberOfIterations++ < this->m_MaximumNumberOfIterations )
    {
    itkDebugMacro( "Iteration Number = " << m_NumberOfIterations );
    for( unsigned int i = 1; i <= this->GetNumberOfClasses(); i++ )
      {  
      this->FindMinimumEnergyLabeling( i );
      E = this->CalculateCurrentEnergy();
      if( E < this->m_CurrentEnergy ) 
        {
        this->m_CurrentEnergy = E;
        this->Relabel();
        itkDebugMacro( "New Energy = " << m_CurrentEnergy );
      }
    }
  }  
} 

template<typename TInputImage, typename TGraphTraits, typename TClassifiedImage>
void
BoykovAlphaExpansionMRFImageFilter<TInputImage, TGraphTraits, TClassifiedImage>
::Initialize()
{
  typename InputImageType::Pointer input0 
    = const_cast<InputImageType*>( static_cast<const InputImageType*>( 
    this->ProcessObject::GetInput( 0 ) ) );

  /** Create the functor which will determine how the 
   *  graph is constructed from the input image. 
   */
  this->m_ImageToGraphFunctor->SetInput( input0 );
  
  /** Create the output image */ 
  typename OutputImageType::Pointer output = this->GetOutput();  
  output->SetRegions( input0->GetBufferedRegion() ); 
  output->Allocate(); 
  
  /** Set the initial labeling.  */
  typedef Statistics::MersenneTwisterRandomVariateGenerator GeneratorType; 
  typename GeneratorType::Pointer generator = GeneratorType::New();

  ImageRegionIteratorWithIndex<OutputImageType> 
     it( output, output->GetRequestedRegion() );         
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    typename OutputImageType::IndexType idx = it.GetIndex();
    unsigned int label = 0;     
    if( this->m_ImageToGraphFunctor->IsPixelANode( idx ) )
      {
      if( this->m_RandomizeInitialLabeling )
        {
        label = generator->GetIntegerVariate( 
          this->GetNumberOfClasses()-1 ) + 1; 
        }
      else
        {    
        label = 1;
        }
      }
    it.Set( label );
    }  

  /** Duplicate the initially labeled output image
   *  and assign its duplication to the m_LabelImage.
   *  This temporary OutputImageType is used to 
   *  determine whether the energy found in the
   *  current iteration of the alpha-expansion 
   *  algorithm is the minimum.
   */
  typedef ImageDuplicator<OutputImageType> ImageDuplicatorType;
  typename ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New();  
  duplicator->SetInputImage( output );
  duplicator->Update();
  this->m_LabelImage = duplicator->GetOutput();  
}

template<typename TInputImage, typename TGraphTraits, typename TClassifiedImage>
typename BoykovAlphaExpansionMRFImageFilter
  <TInputImage, TGraphTraits, TClassifiedImage>::RealType
BoykovAlphaExpansionMRFImageFilter<TInputImage, TGraphTraits, TClassifiedImage>
::CalculateCurrentEnergy()
{
  typename InputImageType::Pointer input0 = const_cast<InputImageType*>(
    static_cast<const InputImageType*>( this->ProcessObject::GetInput( 0 ) ) );
  this->m_ImageToGraphFunctor->SetInput( input0 );

  NeighborhoodIteratorType It( this->m_ImageToGraphFunctor->GetRadius(), 
     this->m_LabelImage, this->m_LabelImage->GetBufferedRegion() );  
  It.ClearActiveList();      
  
  typename BoykovImageToGraphFunctorType::IndexListType::const_iterator it;
  for( it = this->m_ImageToGraphFunctor->GetActiveIndexList().begin(); 
       it != this->m_ImageToGraphFunctor->GetActiveIndexList().end(); ++it )
    {
    It.ActivateOffset( It.GetOffset( *it ) );        
    }  

  RealType energy = 0.0;
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    IndexType idx1 = It.GetIndex();

    if( this->m_ImageToGraphFunctor->IsPixelANode( idx1 ) )
      {
      OutputPixelType label1 = It.GetCenterPixel();
    
      typename ProbabilityImageType::Pointer input;
      input = const_cast<ProbabilityImageType*>( 
        static_cast<const ProbabilityImageType*>( 
        this->ProcessObject::GetInput( label1 ) ) ); 
      this->m_ImageToGraphFunctor->SetSinkLikelihoodImage( input );
    
      /** data term */
      energy += m_ImageToGraphFunctor->GetSinkDataTerm( idx1 );

      typename ShapedNeighborhoodIterator<OutputImageType>::Iterator it;
      for( it = It.Begin(); !it.IsAtEnd(); it++ )
        {
        unsigned int j = it.GetNeighborhoodIndex(); 
        bool IsInBounds;
        OutputPixelType label2 = It.GetPixel( j, IsInBounds );
        if( IsInBounds )
          {
          IndexType idx2 = It.GetIndex( j );
          
          /** smoothness term */
          energy += this->CalculateSmoothnessPenaltyTerm(
            idx1, idx2, label1, label2);
          }
        }
      }
    }   
    
  return energy;
}

template<typename TInputImage, typename TGraphTraits, typename TClassifiedImage>
void
BoykovAlphaExpansionMRFImageFilter<TInputImage, TGraphTraits, TClassifiedImage>
::FindMinimumEnergyBinaryLabeling()
{

  /** Prepare the image-to-graph functor */
  typename ProbabilityImageType::Pointer source = 
    const_cast<ProbabilityImageType*>( static_cast<const ProbabilityImageType*>(
       this->ProcessObject::GetInput( 1 ) ) ); 
  typename ProbabilityImageType::Pointer sink = 
    const_cast<ProbabilityImageType*>( static_cast<const ProbabilityImageType*>(
       this->ProcessObject::GetInput( 2 ) ) ); 
  this->m_ImageToGraphFunctor->SetSourceLikelihoodImage( source );  
  this->m_ImageToGraphFunctor->SetSinkLikelihoodImage( sink ); 

  if( this->m_Indices.size() > 1 )
    {
    this->m_ImageToGraphFunctor->SetSourceIndexContainer( this->m_Indices[1] );
    }
  if( this->m_Indices.size() > 2 )
    {
    this->m_ImageToGraphFunctor->SetSinkIndexContainer( this->m_Indices[2] );
    }   
  
  /** Create the graph from the labeled image */

  typedef ImageToGraphFilter<InputImageType, GraphType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( const_cast<InputImageType*>( 
    static_cast<const InputImageType*>( 
    this->ProcessObject::GetInput( 0 ) ) ) );    
  filter->SetImageToGraphFunctor( this->m_ImageToGraphFunctor ); 
  filter->Update();    
  
  /** Label the graph nodes as 'sink' or 'source' (alpha or not alpha) */
  
  typedef BoykovMinCutGraphFilter<GraphType> MinCutFilterType;
  typename MinCutFilterType::Pointer mincut = MinCutFilterType::New();  
  mincut->SetInput( filter->GetOutput() );
  mincut->Update();  
  typename GraphType::Pointer graph = mincut->GetOutput();

  /** Update the m_LabelImage with the new labeling */
  
  NodePointerType node;
  NodeIteratorType It( graph );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    node = It.GetPointer();
    IndexType idx = node->ImageIndex;       
    if( this->m_LabelImage->GetBufferedRegion().IsInside( idx ) )
      {
      if( !node->IsSink && node->Parent != NULL )
        {        
        this->m_LabelImage->SetPixel( idx, 1 ); 
        }  
      else
        {        
        this->m_LabelImage->SetPixel( idx, 2 ); 
        }  
      }  
    }
}

template<typename TInputImage, typename TGraphTraits, typename TClassifiedImage>
void
BoykovAlphaExpansionMRFImageFilter<TInputImage, TGraphTraits, TClassifiedImage>
::FindMinimumEnergyLabeling( unsigned int alpha )
{
  /** Create the alpha / not-alpha graph. */
  
  typename ProbabilityImageType::Pointer sink = ProbabilityImageType::New();
  sink->SetRegions( this->m_LabelImage->GetBufferedRegion() );
  sink->Allocate();

  typedef CastImageFilter<OutputImageType, InputImageType> CasterType;
  typename CasterType::Pointer caster = CasterType::New();
  caster->SetInput( this->m_LabelImage );
  caster->Update();

  /** Create the not-alpha sink image */

  ImageRegionIterator<InputImageType> 
     It( caster->GetOutput(), caster->GetOutput()->GetRequestedRegion() );  
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    unsigned int label = this->m_LabelImage->GetPixel( It.GetIndex() );
    if( label == alpha || label == 0 )
      {
      It.Set( this->m_ImageToGraphFunctor->GetBackgroundValue() );
      }
    else 
      {
      typename ProbabilityImageType::Pointer input 
        = const_cast<ProbabilityImageType*>(
        static_cast<const ProbabilityImageType*>(
        this->ProcessObject::GetInput( label ) ) );
      sink->SetPixel( It.GetIndex(), input->GetPixel( It.GetIndex() ) );
      }    
    }  

  /** Create the not-alpha Index vector */
  IndexContainerType sink_index;  
  typename IndexContainerType::const_iterator it;

  for( unsigned int i = 1; i <= this->GetNumberOfClasses(); i++ )
    {
    if( i != alpha && this->m_Indices.size() > i )
      {
      for( it = this->m_Indices[i].begin(); 
        it != this->m_Indices[i].end(); ++it ) 
        {
        sink_index.push_back( *it ); 
        }
      }
    }

  /** Prepare the image-to-graph functor */
  typename ProbabilityImageType::Pointer source = 
    const_cast<ProbabilityImageType*>( static_cast<const ProbabilityImageType*>(
       this->ProcessObject::GetInput( alpha ) ) ); 
  this->m_ImageToGraphFunctor->SetSourceLikelihoodImage( source );  
  this->m_ImageToGraphFunctor->SetSinkLikelihoodImage( sink ); 
  bool tmp_ExcludeBackground = m_ImageToGraphFunctor->GetExcludeBackground();
  this->m_ImageToGraphFunctor->SetExcludeBackground( true );
  if( this->m_Indices.size() > alpha )
    {
    this->m_ImageToGraphFunctor->SetSourceIndexContainer( 
      this->m_Indices[alpha] );
    }  
  if( sink_index.size() > 0 )
    {
    this->m_ImageToGraphFunctor->SetSinkIndexContainer( sink_index );   
    }

  /** Create the graph from the labeled image */
  typedef ImageToGraphFilter<InputImageType, GraphType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( caster->GetOutput() );    
  filter->SetImageToGraphFunctor( this->m_ImageToGraphFunctor ); 
  filter->Update();    
  typename GraphType::Pointer graph = filter->GetOutput();
  
  this->m_ImageToGraphFunctor->SetExcludeBackground( tmp_ExcludeBackground );
  this->m_ImageToGraphFunctor->SetInput( this->GetInput() );

  /** Reassign the edge and node weights */

  typedef Image<NodePointerType, ImageDimension> NodeImageType;
  typename NodeImageType::Pointer nodes = NodeImageType::New();  
  nodes->SetRegions( this->m_LabelImage->GetBufferedRegion() );
  nodes->Allocate();
  nodes->FillBuffer( NULL );

  NodePointerType node;
  NodeIteratorType NIt( graph );
  for( NIt.GoToBegin(); !NIt.IsAtEnd(); ++NIt )  
    {    
    node = NIt.GetPointer();
    nodes->SetPixel( node->ImageIndex, node );
    }  

  NeighborhoodIteratorType nit( this->m_ImageToGraphFunctor->GetRadius(), 
    this->m_LabelImage, this->m_LabelImage->GetBufferedRegion() );  
  nit.ClearActiveList();      
  typename BoykovImageToGraphFunctorType::IndexListType::const_iterator iter;
  for( iter = this->m_ImageToGraphFunctor->GetActiveIndexList().begin(); 
       iter != this->m_ImageToGraphFunctor->GetActiveIndexList().end(); ++iter )
    {
    nit.ActivateOffset( nit.GetOffset( *iter ) );        
    }  

  for( nit.GoToBegin(); !nit.IsAtEnd(); ++nit )
    {
    IndexType idx = nit.GetIndex();    
    node = nodes->GetPixel( idx );   
    if( node == NULL )
      {
      continue; 
      }  
    typename ShapedNeighborhoodIterator<OutputImageType>::Iterator it;
    for( it = nit.Begin(); !it.IsAtEnd(); it++ )
      {
      unsigned int d = it.GetNeighborhoodIndex(); 
      bool IsInBounds;
      nit.GetPixel( d, IsInBounds );
      IndexType nidx = nit.GetIndex( d );
      if( idx != nidx && IsInBounds )
        {
        if( nit.GetPixel( d ) != static_cast<int>( alpha ) 
            && nodes->GetPixel( nidx ) != NULL )
          {   
          EdgeWeightType A = this->CalculateSmoothnessPenaltyTerm(
            idx, nidx, alpha, alpha );
          EdgeWeightType B = this->CalculateSmoothnessPenaltyTerm(
            idx, nidx, alpha, nit.GetPixel( d ) );
          EdgeWeightType C = this->CalculateSmoothnessPenaltyTerm(
            idx, nidx, nit.GetCenterPixel(), alpha );
          EdgeWeightType D = this->CalculateSmoothnessPenaltyTerm(
            idx, nidx, nit.GetCenterPixel(), nit.GetPixel( d ) );
          this->AddBinaryTerm( graph, graph->GetEdgePointer(
                  nodes->GetPixel( idx ), nodes->GetPixel( nidx ) ), 
                  A, B, C, D );
          }
        else
          {
          EdgeWeightType A = this->CalculateSmoothnessPenaltyTerm(
            idx, nidx, alpha, nit.GetPixel( d ) );
          EdgeWeightType B = this->CalculateSmoothnessPenaltyTerm( 
            idx, nidx, nit.GetCenterPixel(), alpha );
          this->AddUnaryTerm( graph, node, static_cast<NodeWeightType>( A ), 
            static_cast<NodeWeightType>( B ) );
          }
        }    
      } 
    }
  
  /** Label the graph nodes as 'sink' or 'source' (alpha or not alpha) */
  typedef BoykovMinCutGraphFilter<GraphType> MinCutFilterType;
  typename MinCutFilterType::Pointer mincut = MinCutFilterType::New();  
  mincut->SetInput( filter->GetOutput() );
  mincut->Update();  
  graph = mincut->GetOutput();

  /** Update the m_LabelImage with the new labeling */
  for( NIt.GoToBegin(); !NIt.IsAtEnd(); ++NIt )
    {
    node = NIt.GetPointer();
    IndexType idx = node->ImageIndex;       
    if( this->m_LabelImage->GetBufferedRegion().IsInside( idx ) )
      {
      if( !node->IsSink && node->Parent != NULL )
        {        
        this->m_LabelImage->SetPixel(idx, alpha); 
        }  
      }  
    }
}

template<typename TInputImage, typename TGraphTraits, typename TClassifiedImage>
void 
BoykovAlphaExpansionMRFImageFilter<TInputImage, TGraphTraits, TClassifiedImage>
::Relabel()
{
  typedef ImageDuplicator<OutputImageType> ImageDuplicatorType;
  typename ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New();  
  duplicator->SetInputImage( this->m_LabelImage );
  duplicator->Update();
  this->SetNthOutput( 0, duplicator->GetOutput() );  
}

template<typename TInputImage, typename TGraphTraits, typename TClassifiedImage>
void 
BoykovAlphaExpansionMRFImageFilter<TInputImage, TGraphTraits, TClassifiedImage>
::GeneratePosteriorEnergyImages()
{
  for( unsigned int i = 1; i <= this->GetNumberOfClasses(); i++ )
    {
    typename ProbabilityImageType::Pointer output = ProbabilityImageType::New();
    output->SetRegions( this->m_LabelImage->GetBufferedRegion() );
    output->Allocate();
    output->FillBuffer( 0.0 );   
    this->SetNthOutput( i, output );  
    }  
  
  typename InputImageType::Pointer input0 = const_cast<InputImageType*>(
    static_cast<const InputImageType*>( this->ProcessObject::GetInput( 0 ) ) );
  this->m_ImageToGraphFunctor->SetInput( input0 );

  NeighborhoodIteratorType It( this->m_ImageToGraphFunctor->GetRadius(), 
     this->m_LabelImage, this->m_LabelImage->GetBufferedRegion() );  
  It.ClearActiveList();      
  typename BoykovImageToGraphFunctorType::IndexListType::const_iterator it;
  for( it = this->m_ImageToGraphFunctor->GetActiveIndexList().begin(); 
    it != this->m_ImageToGraphFunctor->GetActiveIndexList().end(); ++it )
    {
    It.ActivateOffset( It.GetOffset( *it ) );        
    }  

  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( this->m_ImageToGraphFunctor->IsPixelANode( It.GetIndex() ) )
      {
      typename ProbabilityImageType::Pointer input;
      input = const_cast<ProbabilityImageType*>(
        static_cast<const ProbabilityImageType*>(
        this->ProcessObject::GetInput( It.GetCenterPixel() ) ) ); 
      this->m_ImageToGraphFunctor->SetSinkLikelihoodImage( input );

      typename ProbabilityImageType::Pointer output;
      output = const_cast<ProbabilityImageType*>(
        static_cast<const ProbabilityImageType*>(
        this->ProcessObject::GetOutput( It.GetCenterPixel() ) ) ); 
      output->SetPixel( It.GetIndex(), this->m_CurrentEnergy );
      
      /** Calculate current neighborhood energy */

      /** data term */
      RealType Energy_N 
        = this->m_ImageToGraphFunctor->GetSinkDataTerm( It.GetIndex() );

      typename ShapedNeighborhoodIterator<OutputImageType>::Iterator it;
      for( it = It.Begin(); !it.IsAtEnd(); it++ )
        {
        unsigned int j = it.GetNeighborhoodIndex(); 
        bool IsInBounds;
        It.GetPixel( j, IsInBounds );
        if( IsInBounds )  
          {
          /** smoothness term */
          Energy_N += this->CalculateSmoothnessPenaltyTerm(
            It.GetIndex(), It.GetIndex( j ), 
            this->m_LabelImage->GetPixel( It.GetIndex() ), 
            this->m_LabelImage->GetPixel( It.GetIndex( j ) ) );
          } 
        }

      /** Calculate neighborhood energy for each label */

      for( unsigned int i = 1; i <= this->GetNumberOfClasses(); i++)
        {
        if( i == static_cast<unsigned int>( It.GetCenterPixel() ) )
          {
          continue;
          } 
        typename ProbabilityImageType::Pointer input;
        input = const_cast<ProbabilityImageType*>(
          static_cast<const ProbabilityImageType*>(
          this->ProcessObject::GetInput( i ) ) ); 
        this->m_ImageToGraphFunctor->SetSinkLikelihoodImage( input );

        /** data term */
        RealType Energy_Ni 
          = this->m_ImageToGraphFunctor->GetSinkDataTerm( It.GetIndex() );

        typename ShapedNeighborhoodIterator<OutputImageType>::Iterator it;
        for( it = It.Begin(); !it.IsAtEnd(); it++ )
          {
          unsigned int j = it.GetNeighborhoodIndex(); 
          bool IsInBounds;
          It.GetPixel( j, IsInBounds );
          if( IsInBounds )  
            {
            /** smoothness term */
            Energy_Ni += this->CalculateSmoothnessPenaltyTerm( It.GetIndex(), 
              It.GetIndex( j ), this->m_LabelImage->GetPixel( It.GetIndex() ), 
              this->m_LabelImage->GetPixel( It.GetIndex( j ) ) );
            }
          }

        typename ProbabilityImageType::Pointer output;
        output = const_cast<ProbabilityImageType*>(
          static_cast<const ProbabilityImageType*>( 
          this->ProcessObject::GetOutput( i ) ) ); 
        output->SetPixel( It.GetIndex(), 
          this->m_CurrentEnergy - Energy_N + Energy_Ni ); 
        }
      }             
    }
}

template<typename TInputImage, typename TGraphTraits, typename TClassifiedImage>
typename BoykovAlphaExpansionMRFImageFilter
  <TInputImage, TGraphTraits, TClassifiedImage>::EdgeWeightType
BoykovAlphaExpansionMRFImageFilter<TInputImage, TGraphTraits, TClassifiedImage>
::CalculateSmoothnessPenaltyTerm(
  IndexType idx1, IndexType idx2, unsigned int label1, unsigned int label2 )
{     
  if( label1 == label2 )
    {
    return static_cast<EdgeWeightType>( 0 );
    }
  return this->m_ImageToGraphFunctor->GetSmoothnessTerm( idx1, idx2 ); 
}

template<typename TInputImage, typename TGraphTraits, typename TClassifiedImage>
void 
BoykovAlphaExpansionMRFImageFilter<TInputImage, TGraphTraits, TClassifiedImage>
::AddUnaryTerm( GraphType *graph, NodePointerType node, 
  NodeWeightType A, NodeWeightType B )
{
  graph->AddNodeWeight( node, B-A );
}

template<typename TInputImage, typename TGraphTraits, typename TClassifiedImage>
void 
BoykovAlphaExpansionMRFImageFilter<TInputImage, TGraphTraits, TClassifiedImage>
::AddBinaryTerm( GraphType *graph, EdgePointerType edge, EdgeWeightType A, 
  EdgeWeightType B, EdgeWeightType C, EdgeWeightType D )
{
  this->AddUnaryTerm( graph, graph->GetNodePointer( edge->SourceIdentifier ), 
    static_cast<NodeWeightType>( A ), static_cast<NodeWeightType>( D ) );
  B -= A;
  C -= D;

  if( B + C < 0 )
    {
    itkExceptionMacro( "Energy function is not regular." );
    }

  if( B < 0 )
    {
    this->AddUnaryTerm( graph, graph->GetNodePointer( edge->SourceIdentifier ), 
      static_cast<NodeWeightType>( B ), static_cast<NodeWeightType>( 0 ) );
    this->AddUnaryTerm( graph, 
      graph->GetNodePointer( edge->TargetIdentifier ), 
      static_cast<NodeWeightType>( -B ), static_cast<NodeWeightType>( 0 ) );
    graph->SetEdgeWeight( edge, static_cast<EdgeWeightType>( 0 ) );
    graph->SetEdgeWeight( graph->GetReverseEdgePointer( edge ), B+C );
    }
  else if( C < 0 )
    {
    this->AddUnaryTerm( graph, graph->GetNodePointer( edge->SourceIdentifier ), 
      static_cast<NodeWeightType>( -C ), static_cast<NodeWeightType>( 0 ) );
    this->AddUnaryTerm(graph, graph->GetNodePointer( edge->TargetIdentifier ), 
      static_cast<NodeWeightType>( C ), static_cast<NodeWeightType>( 0 ) );
    graph->SetEdgeWeight( edge, B+C );
    graph->SetEdgeWeight( graph->GetReverseEdgePointer( edge ), 
      static_cast<EdgeWeightType>( 0 ) );
  }
  else
  {
    graph->SetEdgeWeight( edge, B );
    graph->SetEdgeWeight( graph->GetReverseEdgePointer( edge ), C );
  }
}


template<typename TInputImage, typename TGraphTraits, typename TClassifiedImage>
void 
BoykovAlphaExpansionMRFImageFilter<TInputImage, TGraphTraits, TClassifiedImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

} // end namespace itk

#endif
