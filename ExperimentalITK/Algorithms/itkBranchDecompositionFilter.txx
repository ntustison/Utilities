#ifndef __itkBranchDecompositionFilter_txx
#define __itkBranchDecompositionFilter_txx
#include "itkBranchDecompositionFilter.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNeighborhoodIterator.h"
#include "itkListSample.h"
#include "itkProgressReporter.h"
#include <queue>

namespace itk
{

template <class TInputImage>
BranchDecompositionFilter<TInputImage>
::BranchDecompositionFilter()
{
    m_JCLabelMap = NULL;
    m_InnerRadius = 2.0;
    m_OuterRadius = 3.0;
    m_MinNumberOfPixel = 16;
}

template <class TInputImage>
BranchDecompositionFilter<TInputImage>
::~BranchDecompositionFilter()
{
}

template <class TInputImage>
void BranchDecompositionFilter<TInputImage>
::GenerateData()
{
  typedef itk::Image< float, ImageDimension>               FloatImageType;
  typedef itk::ImageRegionConstIteratorWithIndex< FloatImageType >  FloatConstIteratorType;
  typedef itk::ImageRegionIteratorWithIndex< FloatImageType >       FloatIteratorType;
  typedef itk::NeighborhoodIterator< FloatImageType >               FloatNeighborhoodIteratorType;
  typedef itk::ImageRegionConstIteratorWithIndex< InputImageType >  InputIteratorType;
  typedef itk::ImageRegionIteratorWithIndex< OutputImageType >      OutputIteratorType;
  typedef itk::NeighborhoodIterator< OutputImageType >              OutputNeighborIteratorType;

  typedef std::queue< IndexType >                     QueueType;

  typename Superclass::InputImageConstPointer  inputPtr = this->GetInput();
  typename Superclass::OutputImagePointer outputPtr = this->GetOutput(0);

  // Test if the junction label map is provided
  if( m_JCLabelMap == NULL || m_JCLabelMap->size() == 0 )
  {
    itkExceptionMacro(<< "An nonempty junction label map must be provided");
  }

  // Define images to be used
  typename FloatImageType::Pointer distImage = FloatImageType::New();
  distImage->SetRegions( inputPtr->GetRequestedRegion() );
  distImage->CopyInformation( inputPtr );
  distImage->Allocate();
  outputPtr->SetRequestedRegion(inputPtr->GetRequestedRegion());
  outputPtr->SetBufferedRegion(inputPtr->GetBufferedRegion());
  outputPtr->SetLargestPossibleRegion(inputPtr->GetLargestPossibleRegion());
  outputPtr->Allocate();

  // Get image infos
  SpacingType voxelsize = distImage->GetSpacing();

  //Define some iterators
  FloatIteratorType distIt(distImage, distImage->GetRequestedRegion());
  typename FloatNeighborhoodIteratorType::RadiusType radius;
  radius.Fill(1);
  FloatNeighborhoodIteratorType distneighborIt(radius, distImage, distImage->GetRequestedRegion());
  InputIteratorType inputIt(inputPtr, inputPtr->GetRequestedRegion());
  OutputIteratorType outputIt(outputPtr, outputPtr->GetRequestedRegion());
  OutputNeighborIteratorType outputneighborIt(radius, outputPtr, outputPtr->GetRequestedRegion());

  // Clear labels in the output image and initialize distances to be maximum in the distance image
  unsigned long pixels = 0;
  for( outputIt.GoToBegin(), inputIt.GoToBegin(), distIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt, ++inputIt, ++distIt )
  {
    outputIt.Set(0);
    distIt.Set(itk::NumericTraits<float>::max());
    if( inputIt.Get() != itk::NumericTraits<InputPixelType>::Zero )
      pixels++;
  }

  //Progress reporter
  ProgressReporter progress(this, 0, pixels);

  // Label pixels at the junctions
  QueueType distQueue;
  for(typename JCLabelMapType::const_iterator jclmIt=m_JCLabelMap->begin(); jclmIt!=m_JCLabelMap->end(); ++jclmIt)
  {
    long jcLabel = (*jclmIt).first;
    JCLabelPairType jclPair = (*jclmIt).second;
    IndexType jcIndex = jclPair.first;
    inputIt.SetIndex(jcIndex);
    if( inputIt.Get() != itk::NumericTraits<InputPixelType>::Zero )
    {
      distIt.SetIndex(jcIndex);
      distIt.Set(0.0f);
      outputIt.SetIndex(jcIndex);
      outputIt.Set(jcLabel);
      distQueue.push(jcIndex);
      progress.CompletedPixel();
    }
    else
    {
      itkDebugMacro(<< "Junction " << jcLabel << " is not within the foreground of the input image!");
    }
  }

  // Fill pixels within m_OuterRadius away from the junctions with corresponding junction labels
  while( !distQueue.empty() )
  {
    IndexType centerIndex = distQueue.front();
    distQueue.pop();
    outputneighborIt.SetLocation(centerIndex);

    long centerLabel = outputneighborIt.GetCenterPixel();
    JCLabelPairType centerlPair = (*m_JCLabelMap)[centerLabel];
    float centerRadius = centerlPair.second;
    distIt.SetIndex(centerIndex);
    float centerDist = distIt.Get();
    for(unsigned int i=0; i<outputneighborIt.Size(); i++)
    {
      if( i == outputneighborIt.Size()/2 )
        continue;
      IndexType offIndex = outputneighborIt.GetIndex(i);
      if( !inputIt.GetRegion().IsInside(offIndex) )
      {
        continue;
      }
      inputIt.SetIndex(offIndex);
      if( inputIt.Get() == itk::NumericTraits<InputPixelType>::Zero )
      {
        continue;
      }
      else
      {
        bool inbounds;
        long offLabel = outputneighborIt.GetPixel(i, inbounds);
        if( inbounds )
        {
          OffsetType offset = offIndex - centerIndex;
          float offdist = 0.0f;
          for(int j=0; j<ImageDimension; j++)
          {
            offdist += offset[j]*offset[j]*voxelsize[j]*voxelsize[j];
          }
          offdist = centerDist + sqrt(offdist);
          distIt.SetIndex(offIndex);
          float offDist = distIt.Get();

          if( offdist <= m_OuterRadius * centerRadius && offdist < offDist )
          {
            outputneighborIt.SetPixel(i, centerLabel);
            distIt.Set(offdist);
            distQueue.push(offIndex);
            if( offLabel == 0 )
              progress.CompletedPixel();
          }
        }
      }
    }
  }

  // Get the maximal junction label and Set the current available label
  typename JCLabelMapType::reverse_iterator jclmrIt = m_JCLabelMap->rbegin();
  long maxjcLabel = jclmrIt->first;
  long availLabel = maxjcLabel + 1;

  // Define some maps related to branches
  std::map< long, std::set<long> >  branch_junctionset_map; //branch to its connected junctions
  std::map< long, long >            branch_junction_map;    //branch to its generating junction
  std::map< long, unsigned long >   branch_pixels_map;      //branch to its number of pixels
  std::map< long, std::set<long> >  branch_connect_map;     //branch to its connecting branches
  std::map< long, long >            branch_replace_map;     //branch to its replacement branch

  // Assign new labels to pixels between m_InnerRadius and m_OuterRadius away from the junctions
  QueueType rgQueue;
  for( outputIt.GoToBegin(), distIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt, ++distIt )
  {
    IndexType jcIndex = outputIt.GetIndex();
    inputIt.SetIndex(jcIndex);
    if( inputIt.Get() != itk::NumericTraits<InputPixelType>::Zero )
    {
      long jcLabel = outputIt.Get();
      if( jcLabel > 0 && jcLabel <= maxjcLabel )
      {
        JCLabelPairType jclPair = (*m_JCLabelMap)[jcLabel];
        float jcRadius = jclPair.second;
        distIt.SetIndex(jcIndex);
        float jcDist = distIt.Get();
        float innerRadius = m_InnerRadius * jcRadius;
        float outerRadius = m_OuterRadius * jcRadius;
        if( jcDist > innerRadius && jcDist <= outerRadius )
        {
          long centerLabel = availLabel;
          outputIt.Set(centerLabel);
          branch_junction_map[centerLabel] = jcLabel;
          branch_pixels_map[centerLabel] = 1;
          branch_connect_map[centerLabel] = std::set<long>();
          branch_replace_map[centerLabel] = 0;
          distQueue.push(jcIndex);
          rgQueue.push(jcIndex);
          while( !distQueue.empty() )
          {
            IndexType centerIndex = distQueue.front();
            distQueue.pop();
            outputneighborIt.SetLocation(centerIndex);
            for(unsigned int i=0; i<outputneighborIt.Size(); i++)
            {
              if( i == outputneighborIt.Size()/2 )
                continue;
              IndexType offIndex = outputneighborIt.GetIndex(i);
              if( !inputIt.GetRegion().IsInside(offIndex) )
              {
                continue;
              }
              inputIt.SetIndex(offIndex);
              if( inputIt.Get() == itk::NumericTraits<InputPixelType>::Zero )
              {
                continue;
              }
              else
              {
                bool inbounds;
                long offLabel = outputneighborIt.GetPixel(i, inbounds);
                if( inbounds && offLabel == jcLabel )
                {
                  distIt.SetIndex(offIndex);
                  float offDist = distIt.Get();
                  if( offDist > innerRadius && offDist <= outerRadius )
                  {
                    outputneighborIt.SetPixel(i, centerLabel);
                    branch_pixels_map[centerLabel] = branch_pixels_map[centerLabel] + 1;
                    distQueue.push(offIndex);
                    rgQueue.push(offIndex);
                  }
                }
              }
            }
          }
          availLabel++;
        }
      }
    }
  }
  
  // Region growing with new labels
  while( !rgQueue.empty() ) 
  {
    IndexType centerIndex = rgQueue.front();
    rgQueue.pop();
    outputneighborIt.SetLocation(centerIndex);
    long centerLabel = outputneighborIt.GetCenterPixel();
    for(unsigned int i=0; i<outputneighborIt.Size(); i++)
    {
      if( i == outputneighborIt.Size()/2 )
        continue;
      IndexType offIndex = outputneighborIt.GetIndex(i);
      if( !inputIt.GetRegion().IsInside(offIndex) )
      {
        continue;
      }
      inputIt.SetIndex(offIndex);
      if( inputIt.Get() == itk::NumericTraits<InputPixelType>::Zero )
      {
        continue;
      }
      else
      {
        bool inbounds;
        long offLabel = outputneighborIt.GetPixel(i, inbounds);
        if( inbounds )
        {
          if( offLabel == 0 )
          {
            outputneighborIt.SetPixel(i, centerLabel);
            branch_pixels_map[centerLabel] = branch_pixels_map[centerLabel] + 1;
            rgQueue.push(offIndex);
            progress.CompletedPixel();
          }
          else if( offLabel > maxjcLabel && offLabel != centerLabel )
          {
            long minLabel = offLabel < centerLabel ? offLabel:centerLabel;
            long maxLabel = offLabel > centerLabel ? offLabel:centerLabel;
            std::set<long>& connectSet = branch_connect_map[maxLabel];
            connectSet.insert(minLabel);
          }
        }
      }
    }
  }

  // Generate branch_replace_map
  for(std::map< long, long >::reverse_iterator briter = branch_replace_map.rbegin(); briter != branch_replace_map.rend(); ++briter)
  {
    long brLabel = briter->first;
    long reLabel = briter->second;
    if( reLabel == 0 )
    {
      std::set<long> connectSet = branch_connect_map[brLabel];
      std::set<long> combineSet;
      combineSet.insert(brLabel);
      BRJCSetType jcSet; jcSet.insert(branch_junction_map[brLabel]);
      unsigned long brPixels = branch_pixels_map[brLabel];
      while( !connectSet.empty() )
      {
        long connLabel = *connectSet.rbegin();
        connectSet.erase( connLabel );
        connectSet.insert( branch_connect_map[connLabel].begin(), branch_connect_map[connLabel].end() );
        combineSet.insert( connLabel );
        jcSet.insert(branch_junction_map[connLabel]);
        brPixels = brPixels + branch_pixels_map[connLabel];
      }
      long rpLabel = *combineSet.begin();
      std::set<long>::iterator cbiter;
      if( brPixels >= m_MinNumberOfPixel )
      {
        for(std::set<long>::iterator cbiter = combineSet.begin(); cbiter != combineSet.end(); ++cbiter)
        {
          long mbLabel = *cbiter;
          branch_replace_map[mbLabel] = rpLabel;
          if( mbLabel != rpLabel )
          {
            branch_pixels_map.erase(mbLabel);
          }
          else
          {
            branch_junctionset_map[mbLabel] = jcSet;
            branch_pixels_map[mbLabel] = brPixels;
          }
        }
      }
      else
      {
        for(std::set<long>::iterator cbiter = combineSet.begin(); cbiter != combineSet.end(); ++cbiter)
        {
          long mbLabel = *cbiter;
          branch_replace_map[mbLabel] = rpLabel;
          branch_pixels_map.erase(mbLabel);
        }
      }
    }
  }

  // Resort the branch labels
  m_BRJCConnectMap.clear();
  std::map< long, long > branch_resort_map;
  std::map< long, std::set<long> >::iterator brjcsiter;
  long newLabel = maxjcLabel + 1;
  for(brjcsiter = branch_junctionset_map.begin(); brjcsiter != branch_junctionset_map.end(); ++brjcsiter, ++newLabel)
  {
    long oldLabel = brjcsiter->first;
    m_BRJCConnectMap[newLabel] = brjcsiter->second;
    branch_resort_map[oldLabel] = newLabel;
  }
  
  // Replace branch labels with their resorted representative labels
  for( outputIt.GoToBegin(), distIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt, ++distIt )
  {
    IndexType jcIndex = outputIt.GetIndex();
    inputIt.SetIndex(jcIndex);
    if( inputIt.Get() != itk::NumericTraits<InputPixelType>::Zero )
    {
      long brjcLabel = outputIt.Get();
      if( brjcLabel > maxjcLabel )
      {
        outputIt.Set( branch_resort_map[ branch_replace_map[brjcLabel] ] );
      }
    }
  }
  

}

template <class TInputImage>
void BranchDecompositionFilter<TInputImage>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  
  os << indent << "Inner radius coefficient: " << m_InnerRadius << std::endl;
  os << indent << "Outer radius coefficient: " << m_OuterRadius << std::endl;
  os << indent << "Miminal number of pixels: " << m_MinNumberOfPixel << std::endl;
}

} // end of namespace itk

#endif
