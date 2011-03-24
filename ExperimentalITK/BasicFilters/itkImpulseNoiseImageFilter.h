/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImpulseNoiseImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:52 $
  Version:   $Revision: 1.1.1.1 $
  Author:    Gavin Baker <gavinb@cs.mu.oz.au>

  Copyright (c) 2004 Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkImpulseNoiseImageFilter_h
#define __itkImpulseNoiseImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkUnaryFunctorImageFilter.h"

#include <stdlib.h>

namespace itk
{

/** \class ImpulseNoiseFunctor
 * \brief Pixel functor that adds impulse noise as Salt & Pepper
 */
template < class InputPixelType >
class ImpulseNoiseFunctor
{
public:

  ImpulseNoiseFunctor()
  {
    m_Probability = 0.01f;
    this->m_Minimum = itk::NumericTraits< InputPixelType >::min();
    this->m_Maximum = itk::NumericTraits< InputPixelType >::max();
  }

  float GetProbability() const
  {
    return m_Probability;
  }

  void SetProbability( float prob )
  {
    m_Probability = prob;
  }

  InputPixelType GetMinimum() const
  {
    return m_Minimum;
  }

  void SetMinimum( InputPixelType value )
  {
    m_Minimum = value;
  }

  InputPixelType GetMaximum() const
  {
    return m_Maximum;
  }

  void SetMaximum( InputPixelType value )
  {
    m_Maximum = value;
  }

  float GetNextRandom()
  {
    // Get a random number from [0,1]
    return (float)rand()/(RAND_MAX+1.0);
  }

  InputPixelType operator()( InputPixelType input )
  {

    // There is a 50% chance of salt vs pepper
    float e = 0.0f;
    float var = GetNextRandom();
    float salt_pepper = GetNextRandom() > 0.5 ? this->m_Minimum : this->m_Maximum;

    if ( var < m_Probability )
      {
      e = 1.0f;
      }

    float output = (1-e) * static_cast<float>( input ) + e*salt_pepper;

    // Clamp the output value in valid range
    output = ( output < this->m_Minimum ? this->m_Minimum : output );
    output = ( output > this->m_Maximum ? this->m_Maximum : output );

    return static_cast< InputPixelType > ( output );
  }

private:

  float m_Probability;
  InputPixelType m_Minimum;
  InputPixelType m_Maximum;
};

/** \class ImpulseNoiseImageFilter
 * \brief Adds impulse noise to the input image
 *
 * Adds impulse noise to the input image, according to a probabilistic
 * model.  The noise is 50% salt and 50% pepper.
 *
 * \author Gavin Baker <gavinb at cs_mu_oz_au>
 *
 * \ingroup ImageToImageFilters
 */
template <class TInputImage >
class ITK_EXPORT ImpulseNoiseImageFilter :
    public ImageToImageFilter< TInputImage, TInputImage >
{
public:
  /** Standard class typedefs. */
  typedef ImpulseNoiseImageFilter Self;
  typedef ImageToImageFilter<TInputImage, TInputImage>  Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImpulseNoiseImageFilter, ImageToImageFilter);

  /** Superclass typedefs. */
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  typedef typename Superclass::OutputImagePointer    OutputImagePointer;

  /** Some convenient typedefs. */
  typedef TInputImage InputImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::RegionType     InputRegionType; 
  typedef typename InputImageType::PixelType      InputPixelType; 

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  virtual void GenerateData();

  // Accessor & Mutator methods

  /**
   *    Returns the probability density for the noise.
   *    The default is 0.01.
   */
  float GetProbability() const
  {
      return m_NoiseFilter->GetFunctor().GetProbability();
  }

  /**
   *    Specifies the probability [0,1] that a given voxel will be affected
   *    by noise.
   */
  void SetProbability( float prob )
  {
      return m_NoiseFilter->GetFunctor().SetProbability( prob );
  }

  /**
   *    Return the minimum value for the noise.
   */
  InputPixelType GetMinimum() const
  {
      return m_NoiseFilter->GetFunctor().GetMinimum();
  }

  /**
   *    Set the minimum value for the noise.
   */
  void SetMinimum( InputPixelType value )
  {
      return m_NoiseFilter->GetFunctor().SetMinimum( value );
  }


  /**
   *    Returns the maximum value for the noise.
   */
  InputPixelType GetMaximum() const
  {
      return m_NoiseFilter->GetFunctor().GetMaximum();
  }

  /**
   *    Set the maximum value for the noise.
   */
  void SetMaximum( InputPixelType value )
  {
      return m_NoiseFilter->GetFunctor().SetMaximum( value );
  }



  /**
   *    Specifies the seed for the normal variate generator.  The same seed
   *    will produce the same pseduo-random sequence, which can be used to
   *    reproduce results.  For a higher dose of entropy, initialise with
   *    the current system time (in ms).
   */
  void SetSeed( unsigned long seed )
  {
    m_NoiseFilter->GetFunctor().SetSeed( seed );
  }

protected:

  ImpulseNoiseImageFilter();

  virtual void PrintSelf(std::ostream& os, Indent indent) const;

private:

  ImpulseNoiseImageFilter(const Self&);  // intentionally not implemented
  void operator=(const Self&);      // intentionally not implemented

public:

  typedef UnaryFunctorImageFilter<
    InputImageType,
    InputImageType,
    ImpulseNoiseFunctor< typename InputImageType::PixelType > >
  NoiseFilterType;

private:

  typename NoiseFilterType::Pointer m_NoiseFilter;
};

} /* end namespace itk */

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImpulseNoiseImageFilter.txx"
#endif

#endif /* __itkImpulseNoiseImageFilter_h */
