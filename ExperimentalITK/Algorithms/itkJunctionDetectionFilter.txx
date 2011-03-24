/*=========================================================================
 * Automatic Junction Detection for Tubular Structures
 * Author:      Guanglei Xiong (guangleixiong at gmail.com)
 * Date:        March 23, 2009
 * Reference:   http://www.insight-journal.org/browse/publication/324
=========================================================================*/

#ifndef __itkJunctionDetectionFilter_txx
#define __itkJunctionDetectionFilter_txx
#include "itkJunctionDetectionFilter.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNeighborhoodIterator.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkListSample.h"
#include "itkSampleMeanShiftClusteringFilter.h"
#include "itkProgressReporter.h"
#include <queue>
#include <set>

namespace itk
{

template <typename T>
struct SetComp
{
  bool operator() (const T& lhs, const T& rhs) const
  {
    if( lhs[2] < rhs[2] )  return true;
    else if( lhs[2] == rhs[2] )
    {
      if( lhs[1] < rhs[1] )  return true;
      else if( lhs[1] == rhs[1] && lhs[0] < rhs[0] )  return true;
      else     return false;
    }
    else
    {
      return false;
    }
  }
};

template <class TInputImage>
JunctionDetectionFilter<TInputImage>
::JunctionDetectionFilter()
{
    m_BackgroundValue = NumericTraits<InputPixelType>::Zero;
    m_InnerRadius = 2.0;
    m_OuterRadius = 3.0;
    m_MinNumberOfPixel = 16;
}

template <class TInputImage>
JunctionDetectionFilter<TInputImage>
::~JunctionDetectionFilter()
{
}

template <class TInputImage>
void JunctionDetectionFilter<TInputImage>
::GenerateData()
{
  typedef itk::Image< float, ImageDimension>               FloatImageType;
  typedef itk::Image< IndexType, ImageDimension >          IndexImageType;
  typedef itk::ImageRegionConstIteratorWithIndex< FloatImageType >  FloatConstIteratorType;
  typedef itk::ImageRegionIterator< IndexImageType >                IndexIteratorType;
  typedef itk::ImageRegionIteratorWithIndex< FloatImageType >       FloatIteratorType;
  typedef itk::NeighborhoodIterator< IndexImageType >               IndexNeighborhoodIteratorType;
  typedef itk::NeighborhoodIterator< FloatImageType >               FloatNeighborhoodIteratorType;
  typedef itk::ImageRegionIteratorWithIndex< OutputImageType >      OutputIteratorType;
  typedef itk::NeighborhoodIterator< OutputImageType >              OutputNeighborIteratorType;

  typedef std::queue< IndexType >                     DistanceQueueType;
  typedef std::set< IndexType, SetComp<IndexType> >   DistanceSetType;
  typedef itk::Vector< typename IndexType::IndexValueType, ImageDimension >   JCIndexType;
  typedef itk::SignedMaurerDistanceMapImageFilter< InputImageType, FloatImageType >  DistanceFilterType;
  typedef itk::Statistics::ListSample< JCIndexType >    JCSampleType;
  typedef itk::Statistics::SampleMeanShiftClusteringFilter< JCSampleType > JCClusteringType;

  typename Superclass::InputImageConstPointer  inputPtr = this->GetInput();
  typename Superclass::OutputImagePointer outputPtr = this->GetOutput(0);

  // Compute the distance map of the input
  typename DistanceFilterType::Pointer wallDister = DistanceFilterType::New();
  wallDister->UseImageSpacingOn();
  wallDister->SquaredDistanceOff();
  wallDister->InsideIsPositiveOn();
  wallDister->SetBackgroundValue(m_BackgroundValue);
  wallDister->SetInput( inputPtr );
  wallDister->Update();

  // Define images to be used
  typename FloatImageType::Pointer wallDistImage = wallDister->GetOutput();
  typename FloatImageType::Pointer distImage = FloatImageType::New();
  distImage->SetRegions( wallDistImage->GetRequestedRegion() );
  distImage->CopyInformation( wallDistImage );
  distImage->Allocate();
  typename IndexImageType::Pointer labelImage = IndexImageType::New();
  labelImage->SetRegions( wallDistImage->GetRequestedRegion() );
  labelImage->CopyInformation( wallDistImage );
  labelImage->Allocate();
  outputPtr->SetRequestedRegion(inputPtr->GetRequestedRegion());
  outputPtr->SetBufferedRegion(inputPtr->GetBufferedRegion());
  outputPtr->SetLargestPossibleRegion(inputPtr->GetLargestPossibleRegion());
  outputPtr->Allocate();

  // Define some iterators
  FloatConstIteratorType wallDistIt(wallDistImage, wallDistImage->GetRequestedRegion());
  IndexIteratorType labelIt(labelImage, labelImage->GetRequestedRegion());
  FloatIteratorType distIt(distImage, distImage->GetRequestedRegion());
  typename IndexNeighborhoodIteratorType::RadiusType radius;
  radius.Fill(1);
  FloatNeighborhoodIteratorType wallneighborIt(radius, wallDistImage, wallDistImage->GetRequestedRegion());
  IndexNeighborhoodIteratorType labelneighborIt(radius, labelImage, labelImage->GetRequestedRegion());
  FloatNeighborhoodIteratorType distneighborIt(radius, distImage, distImage->GetRequestedRegion());
  OutputIteratorType outputIt(outputPtr, outputPtr->GetRequestedRegion());
  OutputNeighborIteratorType outputneighborIt(radius, outputPtr, outputPtr->GetRequestedRegion());

  // Clear the label and count pixels to process
  IndexType clearIndex;
  clearIndex.Fill(-1);
  unsigned long pixels = 0;
  for( labelIt.GoToBegin(), wallDistIt.GoToBegin(); !labelIt.IsAtEnd(); ++labelIt, ++wallDistIt )
  {
    labelIt.Set(clearIndex);
    if( wallDistIt.Get() > 0.0f)
      pixels++;
  }

  // Get image infos
  SpacingType voxelsize = distImage->GetSpacing();
  
  // Loop through all the pixels which has wallDist>0
  ProgressReporter progress(this, 0, pixels);
  typename JCSampleType::Pointer jcSample = JCSampleType::New();
  for( wallDistIt.GoToBegin(), distIt.GoToBegin(), labelIt.GoToBegin(); !wallDistIt.IsAtEnd(); ++wallDistIt, ++distIt, ++labelIt )
  {
    float wallDist = wallDistIt.Get();
    if( wallDist > 0.0f )    //Only consider pixels inside the object
    {
      IndexType currentLabel = distIt.GetIndex();
      distIt.Set(0.0f);
      labelIt.Set(currentLabel);
      DistanceQueueType distQueue;
      DistanceSetType distSet;
      distQueue.push(currentLabel);
      while( !distQueue.empty() )
      {
        IndexType centerIndex = distQueue.front();
        distQueue.pop();
        wallneighborIt.SetLocation(centerIndex);
        labelneighborIt.SetLocation(centerIndex);
        distneighborIt.SetLocation(centerIndex);
        float centerDist = distneighborIt.GetCenterPixel();
        for(unsigned int i=0; i<distneighborIt.Size(); i++)
        {
          if(i != distneighborIt.Size()/2)
          {
            OffsetType offset = distneighborIt.GetOffset(i);
            
            bool inbounds;
            float offwall = wallneighborIt.GetPixel(i, inbounds);
            if( offwall<=0.0f || !inbounds ) continue;

            IndexType offlabel = labelneighborIt.GetPixel(i);
            float offdist = distneighborIt.GetPixel(i);
            float newdist = 0.0f;
            for(unsigned int j=0; j<ImageDimension; j++) 
              newdist += offset[j]*offset[j]*voxelsize[j]*voxelsize[j];
            newdist = centerDist + sqrt(newdist);
            if( offlabel == currentLabel && offdist <= newdist )
            {
              continue;
            }
            else if( newdist < m_OuterRadius*wallDist )
            {
              if( offlabel != currentLabel )
              {
                labelneighborIt.SetPixel(i, currentLabel);
              }
              distneighborIt.SetPixel(i, newdist);
              IndexType offIndex = distneighborIt.GetIndex(i);
              distQueue.push(offIndex);
              if( newdist > m_InnerRadius*wallDist )
              {
                distSet.insert(offIndex);
              }
              else
              {
                typename DistanceSetType::iterator dsIt = distSet.find(offIndex);
                if( dsIt != distSet.end() )
                {
                  distSet.erase(dsIt);
                }
              }
            }

          }
        }
      }
      // The minimal number of pixels inside of the hollow sphere
      if( distSet.size() > m_MinNumberOfPixel )
      {
        size_t origdistsize = distSet.size();
        DistanceQueueType distQueue;
        int compCount = 0;
        while( compCount < 3 && distSet.size() > 0 )
        {
          size_t olddistsize = distSet.size();
          typename DistanceSetType::iterator dsIt = distSet.begin();
          distQueue.push(*dsIt);
          distSet.erase(dsIt);
          while( !distQueue.empty() )
          {
            IndexType oneIndex = distQueue.front();
            distQueue.pop();
            for(dsIt = distSet.begin(); dsIt!=distSet.end(); ++dsIt)
            {
              IndexType anotherIndex = *dsIt;
              OffsetType offset = anotherIndex - oneIndex;
              typename OffsetType::OffsetValueType offsetSumSquared = 0;
              for(int j=0; j<ImageDimension; j++)
                offsetSumSquared += offset[j]*offset[j];
              if( offsetSumSquared <= 3 )
              {
                distQueue.push(anotherIndex);
                distSet.erase(dsIt);
              }
            }
          }
          size_t newdistsize = distSet.size();
          if( olddistsize - newdistsize > origdistsize/m_MinNumberOfPixel)
          {
            compCount++;
          }
        }
        if( compCount >= 3 )
        {
          JCIndexType jcIndex;
          for(int j=0; j<ImageDimension; j++)
            jcIndex[j] = currentLabel[j];
          jcSample->PushBack(jcIndex);
          //std::cout << "Junction: " << currentLabel << std::endl;
        }
      }
      progress.CompletedPixel();
    }
  }

  typename JCClusteringType::Pointer clustering = JCClusteringType::New();
  clustering->SetInputSample(jcSample);
  clustering->SetThreshold(8);
  clustering->SetMinimumClusterSize(1);
  clustering->Update(); 

  m_JCLabelMap.clear();
  typename JCClusteringType::ClusterLabelsType& jcLabels = clustering->GetOutput();
  typename JCSampleType::Iterator jcsIt;
  typename JCClusteringType::ClusterLabelsType::iterator jclIt;
  for( jcsIt = jcSample->Begin(), jclIt = jcLabels.begin(); jcsIt!=jcSample->End(); ++jcsIt, ++jclIt )
  {
    typename JCSampleType::MeasurementVectorType jcVector = jcsIt.GetMeasurementVector();
    IndexType jcIndex;
    for(int j=0; j<ImageDimension; j++)
      jcIndex[j] = jcVector[j];
    long jcLabel = *jclIt;
    //std::cout << jcIndex << " " << jcLabel << std::endl;
    if( m_JCLabelMap.find(jcLabel) != m_JCLabelMap.end() )
    {
      JCLabelPairType jclPair = m_JCLabelMap[jcLabel];
      wallDistIt.SetIndex(jcIndex);
      if( wallDistIt.Get() > jclPair.second )
      {
        m_JCLabelMap[jcLabel] = JCLabelPairType(jcIndex, wallDistIt.Get());
      }
    }
    else
    {
      wallDistIt.SetIndex(jcIndex);
      m_JCLabelMap[jcLabel] = JCLabelPairType(jcIndex, wallDistIt.Get());
    }
  }

  // Clear the label
  for( outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt )
  {
    outputIt.Set(0);
  }
  
  // Render spheres on the output image to represent found junctions
  //std::cout << "Junctions (" << m_JCLabelMap.size() << "):" << std::endl;
  for(typename JCLabelMapType::iterator jclmIt=m_JCLabelMap.begin(); jclmIt!=m_JCLabelMap.end(); ++jclmIt)
  {
    long jcLabel = (*jclmIt).first;
    JCLabelPairType jclPair = (*jclmIt).second;
    IndexType jcIndex = jclPair.first;
    float jcRadius = jclPair.second;
    //std::cout << jcLabel << " " << jcIndex << " " << jcRadius << std::endl;
    outputIt.SetIndex(jcIndex);
    outputIt.Set(jcLabel);
    DistanceQueueType distQueue;
    distQueue.push(jcIndex);
    while( !distQueue.empty() )
    {
      IndexType centerIndex = distQueue.front();
      distQueue.pop();
      outputneighborIt.SetLocation(centerIndex);
      for(unsigned int i=0; i<outputneighborIt.Size(); i++)
      {
        if( i != outputneighborIt.Size()/2 )
        {
          bool inbounds;
          long offLabel = outputneighborIt.GetPixel(i, inbounds);
          if( offLabel == 0 && inbounds )
          {
            IndexType offIndex = outputneighborIt.GetIndex(i);
            OffsetType offset = offIndex - jcIndex;
            float offdist = 0.0f;
            for(int j=0; j<ImageDimension; j++)
            {
              offdist += offset[j]*offset[j]*voxelsize[j]*voxelsize[j];
            }
            offdist = sqrt(offdist);
            if( offdist <= jcRadius )
            {
              outputneighborIt.SetPixel(i, jcLabel);
              distQueue.push(offIndex);
            }
          }
        }
      }
    }
  }

}

template <class TInputImage>
void JunctionDetectionFilter<TInputImage>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  
  os << indent << "Inner radius coefficient: " << m_InnerRadius << std::endl;
  os << indent << "Outer radius coefficient: " << m_OuterRadius << std::endl;
  os << indent << "Miminal number of pixels: " << m_MinNumberOfPixel << std::endl;
}

} // end of namespace itk

#endif
