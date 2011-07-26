/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAdditiveGaussianNoiseImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:49 $
  Version:   $Revision: 1.1.1.1 $
  Author:    Gavin Baker <gavinb@cs.mu.oz.au>

  Copyright (c) 2004 Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAdditiveGaussianNoiseImageFilter_h
#define __itkAdditiveGaussianNoiseImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkUnaryFunctorImageFilter.h"

namespace itk
{

/** \class NormalVariateMultiplierFunctor
 * \brief Pixel functor that adds Gaussian noise
 */
template < class InputPixelType >
class NormalVariateNoiseFunctor
{
public:

  NormalVariateNoiseFunctor()
  {
    m_Mean = 0.0f;
    m_StandardDeviation = 1.0f;

    m_Generator = Statistics::MersenneTwisterRandomVariateGenerator::New();
    m_Generator->SetSeed();
  }

  float GetMean() const
  {
    return m_Mean;
  }

  void SetMean( float mean )
  {
    m_Mean = mean;
  }

  float GetStandardDeviation() const
  {
    return m_StandardDeviation;
  }

  void SetStandardDeviation( float stddev )
  {
    m_StandardDeviation = stddev;
  }

  void SetSeed( unsigned long seed )
  {
    m_Generator->Initialize( seed );
  }

  InputPixelType operator()( InputPixelType input )
  {
    static const float min = itk::NumericTraits< InputPixelType >::min();
    static const float max = itk::NumericTraits< InputPixelType >::max();

    float output = static_cast<float>( input ) +
      m_Generator->GetNormalVariate( this->m_Mean,
      vnl_math_sqr( this->m_StandardDeviation ) );

    // Clamp the output value in valid range
    output = ( output < min ? min : output );
    output = ( output > max ? max : output );
    return static_cast< InputPixelType > ( output );
  }

private:

  float m_Mean;
  float m_StandardDeviation;
  typename Statistics
     ::MersenneTwisterRandomVariateGenerator::Pointer m_Generator;
};


/** \class AdditiveGaussianNoiseImageFilter
 * \brief Adds Gaussian noise to the input image
 *
 * Adds noise to the input image according to a Gaussian normal variate
 * distribution.  The user supplies the mean \f$\bar{x}\f$ and standard
 * deviation \f$\sigma\f$, such that the output is given by:
 *
 * \f[
 *     v_{out} = v_{in} + \bar{x} + \sigma * G(d)
 * \f]
 *
 * where G() is the Gaussian generator and d is the seed.  A particular seed
 * can be specified in order to perform repeatable tests.
 *
 * \author Gavin Baker <gavinb at cs_mu_oz_au>
 *
 * \ingroup ImageToImageFilters
 */
template <class TInputImage >
class ITK_EXPORT AdditiveGaussianNoiseImageFilter :
    public ImageToImageFilter< TInputImage, TInputImage >
{
public:
  /** Standard class typedefs. */
  typedef AdditiveGaussianNoiseImageFilter Self;
  typedef ImageToImageFilter<TInputImage, TInputImage>  Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(AdditiveGaussianNoiseImageFilter, ImageToImageFilter);

  /** Superclass typedefs. */
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  typedef typename Superclass::OutputImagePointer    OutputImagePointer;

  /** Some convenient typedefs. */
  typedef TInputImage InputImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::RegionType     InputImageRegionType;
  typedef typename InputImageType::PixelType      InputImagePixelType;

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  // virtual void GenerateOutputInformation();

  virtual void GenerateData();

  // Accessor & Mutator methods

  /**
   *    Specifies the average noise added to the image per pixel.
   *    The default is 0.
   */
  void SetMean( float mean )
  {
    m_NoiseFilter->GetFunctor().SetMean( mean );
  }

  /**
   *    Returns the average noise added to the image per pixel.
   *    The default is 0.
   */
  float GetMean() const
  {
    return m_NoiseFilter->GetFunctor().GetMean();
  }

  /**
   *    Specifies the standard deviation of the noise added to the image.
   *    The default is 1.
   */
  void SetStandardDeviation( float stddev )
  {
    m_NoiseFilter->GetFunctor().SetStandardDeviation( stddev );
  }

  /**
   *    Returns the standard deviation of the noise added to the image.
   *    The default is 1.
   */
  float GetStandardDeviation() const
  {
    return m_NoiseFilter->GetFunctor().GetStandardDeviation();
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

  AdditiveGaussianNoiseImageFilter();

  virtual void PrintSelf(std::ostream& os, Indent indent) const;

private:

  AdditiveGaussianNoiseImageFilter(const Self&);  // intentionally not implemented
  void operator=(const Self&);      // intentionally not implemented

public:

  typedef UnaryFunctorImageFilter<
    InputImageType,
    InputImageType,
    NormalVariateNoiseFunctor< typename InputImageType::PixelType > >
  NoiseFilterType;

private:

  typename NoiseFilterType::Pointer m_NoiseFilter;
};

} /* end namespace itk */

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAdditiveGaussianNoiseImageFilter.hxx"
#endif

#endif /* __itkAdditiveGaussianNoiseImageFilter_h */
