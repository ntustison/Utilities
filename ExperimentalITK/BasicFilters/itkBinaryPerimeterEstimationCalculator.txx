/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBinaryPerimeterEstimationCalculator.txx,v $
  Language:  C++
  Date:      $Date: 2004/12/21 22:47:30 $
  Version:   $Revision: 1.12 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

    This software is distributed WITHOUT ANY WARRANTY; without even 
    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
    PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBinaryPerimeterEstimationCalculator_txx
#define __itkBinaryPerimeterEstimationCalculator_txx

#include "itkBinaryPerimeterEstimationCalculator.h"
#include "itkProgressReporter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkSize.h"
#include "itkConnectedComponentAlgorithm.h"

namespace itk {

template <class TInputImage>
BinaryPerimeterEstimationCalculator<TInputImage>
::BinaryPerimeterEstimationCalculator()
{
  m_FullyConnected = false;
  m_ForegroundValue = NumericTraits<InputImagePixelType>::max();
  m_Perimeter = 0;
}


template<class TInputImage>
void
BinaryPerimeterEstimationCalculator<TInputImage>
::Compute()
{  
  m_Perimeter = 0;

  // reduce the region to avoid reading outside
  RegionType region = this->GetImage()->GetBufferedRegion();
  SizeType size = region.GetSize();
  for( int i=0; i<ImageDimension; i++ )
    {
    size[i]--;
    }
  region.SetSize( size );

  // the radius which will be used for all the shaped iterators
  Size< ImageDimension > radius;
  radius.Fill(1);

  // set up the iterator
  typedef ConstShapedNeighborhoodIterator< InputImageType > IteratorType;
  typename IteratorType::ConstIterator nIt;
  IteratorType iIt(radius, this->GetImage(), region);
  // we want to search the neighbors with offset >= 0
  // 2D -> 4 neighbors
  // 3D -> 8 neighbors
  typename IteratorType::OffsetType offset;
  unsigned int centerIndex = iIt.GetCenterNeighborhoodIndex();
  // store the offsets to reuse them to evaluate the contributions of the
  // configurations
  typename std::vector< IndexType > indexes;
  IndexType idx0;
  idx0.Fill(0);
  for ( unsigned int d = centerIndex; d < 2 * centerIndex + 1; d++ )
    {
    offset = iIt.GetOffset(d);
    bool deactivate = false;
    for ( unsigned int j = 0; j < ImageDimension && !deactivate; j++ )
      {
      if ( offset[j] < 0 )
        {
        deactivate = true;
        }
      }
    if ( deactivate )
      {
      iIt.DeactivateOffset(offset);
      }
    else
      {
      iIt.ActivateOffset(offset);
      indexes.push_back(idx0 + offset);
      }
    }
  
  // to store the configurations count
  typedef typename std::map< unsigned long, unsigned long > MapType;
  MapType confCount;

  for( iIt.GoToBegin(); !iIt.IsAtEnd(); ++iIt )
    {
    unsigned long conf = 0;
    int i=0;
    for ( nIt= iIt.Begin();
      nIt != iIt.End();
      nIt++, i++ )
      {
      if( nIt.Get() == m_ForegroundValue )
        {
        conf += 1 << i;
        }
      }
    confCount[ conf ]++;
    // progress.CompletedPixel();
    }

  // compute the participation to the perimeter for all the configurations
  double physicalSize = 1.0;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    physicalSize *= this->GetImage()->GetSpacing()[i];
    }
  typedef typename std::map< unsigned long, double > ContributionMapType;
  ContributionMapType contributions;
  const unsigned int  numberOfNeighbors      =
    static_cast< unsigned int >( vcl_pow( 2.0, static_cast< double >( ImageDimension ) ) );
  const unsigned int numberOfConfigurations =
    static_cast< unsigned int >( vcl_pow( 2.0, static_cast< double >( numberOfNeighbors ) ) );
  // create an image to store the neighbors
  typedef typename itk::Image< bool, ImageDimension > ImageType;
  typename ImageType::Pointer neighborsImage = ImageType::New();
  // typename ImageType::SizeType size;
  size.Fill(2);
  neighborsImage->SetRegions(size);
  neighborsImage->Allocate();
  for ( unsigned int i = 0; i < numberOfConfigurations; i++ )
    {
    neighborsImage->FillBuffer(false);
    for ( unsigned int j = 0; j < numberOfNeighbors; j++ )
      {
      if ( i & 1 << j )
        {
        neighborsImage->SetPixel(indexes[j], true);
        }
      }
    // the image is created - we can now compute the contributions of the pixels
    // for that configuration
    contributions[i] = 0;
    for ( unsigned int j = 0; j < numberOfNeighbors; j++ )
      {
      IndexType currentIdx = indexes[j];
      if ( neighborsImage->GetPixel(currentIdx) )
        {
        for ( unsigned int k = 0; k < ImageDimension; k++ )
          {
          IndexType idx = currentIdx;
          idx[k] = vcl_abs(idx[k] - 1);
          if ( !neighborsImage->GetPixel(idx) )
            {
            contributions[i] += physicalSize / this->GetImage()->GetSpacing()[k] / 2.0;
            }
          }
        }
      }
    contributions[i] /= ImageDimension;
    }

   for( typename MapType::const_iterator it=confCount.begin();
     it!=confCount.end();
     it++ )
     {
     m_Perimeter += contributions[it->first] * it->second;
     }

}



template<class TInputImage>
void
BinaryPerimeterEstimationCalculator<TInputImage>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  
  os << indent << "FullyConnected: "  << m_FullyConnected << std::endl;
  os << indent << "ForegroundValue: "  << static_cast<typename NumericTraits<InputImagePixelType>::PrintType>(m_ForegroundValue) << std::endl;
  os << indent << "Perimeter: " << static_cast<typename NumericTraits< double >::PrintType>(m_Perimeter) << std::endl;
}
  
}// end namespace itk
#endif
