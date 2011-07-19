#ifndef __itkTopologyPreservingDigitalSurfaceEvolutionImageFilter_hxx
#define __itkTopologyPreservingDigitalSurfaceEvolutionImageFilter_hxx

#include "itkTopologyPreservingDigitalSurfaceEvolutionImageFilter.h"

#include "itkAddImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkMapContainer.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkTimeProbe.h"
#include "itkVariableSizeMatrix.h"
#include "itkWellComposedImageFilter.h"

#include "itkImageFileWriter.h"

namespace itk
{

template<class TImage>
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::TopologyPreservingDigitalSurfaceEvolutionImageFilter()
{

  this->SetNumberOfRequiredInputs( 1 );
  this->m_NumberOfIterations = 100;
  this->m_ThresholdValue = 0.5;
  this->m_BackgroundValue = NumericTraits<PixelType>::Zero;
  this->m_ForegroundValue = NumericTraits<PixelType>::One;
  this->m_GrowOnly = false;
  this->m_GlueObjects = false;
  this->m_FillCavities = false;
  this->m_GluingStrategy = 0;
  this->m_TargetImage = NULL;

  if( ImageDimension == 2 )
    {
    this->InitializeIndices2D();
    }
  else if( ImageDimension == 3 )
    {
    this->InitializeIndices3D();
    }
  else
    {
    itkExceptionMacro( "Image dimension must be equal to 2 or 3." );
    }
}

template<class TImage>
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::~TopologyPreservingDigitalSurfaceEvolutionImageFilter()
{
}

template<class TImage>
void
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::GenerateData()
{
  this->AllocateOutputs();

  if( !this->m_TargetImage )
    {
    itkExceptionMacro( "TargetImage not specified." );
    }

  typedef BinaryThresholdImageFilter<ImageType, ImageType> ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( this->GetInput() );
  thresholder->SetOutsideValue( this->m_ForegroundValue );
  thresholder->SetInsideValue( this->m_BackgroundValue );
  thresholder->SetLowerThreshold( this->m_BackgroundValue );
  thresholder->SetUpperThreshold( this->m_BackgroundValue );
  thresholder->Update();

  this->m_LabelSurfaceImage = thresholder->GetOutput();

  RealType totalDifference = 0.0;

  ImageRegionIterator<ImageType> ItL( this->m_LabelSurfaceImage,
    this->m_LabelSurfaceImage->GetRequestedRegion() );
  ImageRegionIterator<RealImageType> ItT( this->m_TargetImage,
    this->m_TargetImage->GetRequestedRegion() );
  for( ItL.GoToBegin(), ItT.GoToBegin(); !ItL.IsAtEnd(); ++ItL, ++ItT )
    {
    if( ItT.Get() >= this->m_ThresholdValue
      && ItL.Get() == this->m_BackgroundValue )
      {
      totalDifference += 1.0;
      }
    else if( !this->m_GrowOnly &&
      ( 1.0 - ItT.Get() >= this->m_ThresholdValue
      && ItL.Get() == this->m_ForegroundValue ) )
      {
      totalDifference += 1.0;
      }
    }

  this->m_SurfaceLabel = this->m_ForegroundValue + 1;

  unsigned int iterations = 0;
  bool changeDetected = true;
  bool inverseChangeDetected = true;
  if( this->m_GrowOnly )
    {
    inverseChangeDetected = false;
    }

  RealType totalNumberOfChanges = 0.0;
  std::cout << "/" << std::flush;

  TimeProbe timer;
  timer.Start();

  this->SetProgress( 0.0 );
  RealType oldProgress = this->GetProgress() * 100.0;
  unsigned int progressCount = 0;

  this->m_IsInversionStep = false;
  while( iterations++ < this->m_NumberOfIterations && ( changeDetected ||
    inverseChangeDetected ) )
    {
    if( this->m_IsInversionStep )
      {
      inverseChangeDetected = false;
      }
    else
      {
      changeDetected = false;
      }

    this->CreateLabelSurfaceImage();

    ImageRegionIterator<RealImageType> ItT( this->m_TargetImage,
      this->m_TargetImage->GetRequestedRegion() );
    ImageRegionIterator<ImageType> ItL( this->m_LabelSurfaceImage,
      this->m_LabelSurfaceImage->GetRequestedRegion() );

    ItT.GoToBegin();
    ItL.GoToBegin();
    while( !ItL.IsAtEnd() )
      {
      bool meetsThresholdCriteria = false;
      if( this->m_IsInversionStep )
        {
        meetsThresholdCriteria =
          ( 1.0 - ItT.Get() >= this->m_ThresholdValue );
        }
      else
        {
        meetsThresholdCriteria = ( ItT.Get() >= this->m_ThresholdValue );
        }
      if( ItL.Get() == this->m_SurfaceLabel && meetsThresholdCriteria &&
        this->IsChangeSafe( ItL.GetIndex() ) )
        {
        if( this->m_IsInversionStep )
          {
          ItL.Set( this->m_BackgroundValue );
          inverseChangeDetected = true;
          }
        else
          {
          ItL.Set( this->m_ForegroundValue );
          changeDetected = true;
          }

        totalNumberOfChanges += 1.0;
        this->SetProgress( totalNumberOfChanges / totalDifference );
        RealType newProgress = this->GetProgress() * 100.0;

        if( newProgress - oldProgress >= 1.0 )
          {
          oldProgress = newProgress;
          std::cout << "*" << std::flush;
          progressCount++;
          if( progressCount % 10 == 0 )
            {
            std::cout << progressCount << std::flush;
            }
          }
        }
      if( ItL.Get() == this->m_SurfaceLabel )
        {
        if( this->m_IsInversionStep )
          {
          ItL.Set( this->m_ForegroundValue );
          }
        else
          {
          ItL.Set( this->m_BackgroundValue );
          }
        }
      ++ItT;
      ++ItL;
      }
    if( !this->m_GrowOnly && ( this->m_IsInversionStep && changeDetected ) ||
      ( !this->m_IsInversionStep && inverseChangeDetected ) )
      {
      this->m_IsInversionStep = !this->m_IsInversionStep;
      }
    }
  timer.Stop();
  std::cout << "/ -> " << this->GetProgress() << " ("
    << timer.GetMeanTime() << " seconds)" << std::endl;

  typedef ConnectedComponentImageFilter
    <ImageType, ImageType> ConnecterType;
  typename ConnecterType::Pointer connecter = ConnecterType::New();
  connecter->SetInput( this->m_LabelSurfaceImage );
  connecter->Update();

  typedef MapContainer<PixelType, PixelType> MapType;
  typename MapType::Pointer map = MapType::New();
  map->Initialize();

  ImageRegionConstIterator<ImageType> ItI( this->GetInput(),
    this->GetInput()->GetRequestedRegion() );
  ImageRegionIterator<ImageType> ItC( connecter->GetOutput(),
    connecter->GetOutput()->GetRequestedRegion() );
  for( ItI.GoToBegin(), ItC.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItC )
    {
    PixelType valI = ItI.Get();
    PixelType valC = ItC.Get();
    if( valI != this->m_BackgroundValue && valC != this->m_BackgroundValue
      && !map->IndexExists( valC ) )
      {
      map->SetElement( valC, valI );
      }
    }

  ImageRegionIterator<ImageType> ItO( this->GetOutput(),
    this->GetOutput()->GetRequestedRegion() );
  for( ItO.GoToBegin(), ItC.GoToBegin(); !ItO.IsAtEnd(); ++ItO, ++ItC )
    {
    if( ItC.Get() != this->m_BackgroundValue )
      {
      ItO.Set( map->GetElement( ItC.Get() ) );
      }
    else
      {
      ItO.Set( this->m_BackgroundValue );
      }
    }

  if( this->m_GlueObjects && map->Size() > 1 )
    {
    this->FindLabelOrderingForObjectGluing();
    this->GlueObjects();

    if( this->m_FillCavities )
      {
      this->FillCavities();
      }
    }
}

template<class TImage>
void
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::CreateLabelSurfaceImage()
{

  if( this->m_IsInversionStep )
    {
    typedef BinaryContourImageFilter<ImageType, ImageType> ContourFilterType;
    typename ContourFilterType::Pointer contourFilter
      = ContourFilterType::New();
    contourFilter->SetInput( this->m_LabelSurfaceImage );
    contourFilter->SetFullyConnected( false );
    contourFilter->SetBackgroundValue( this->m_BackgroundValue );
    contourFilter->SetForegroundValue( this->m_ForegroundValue );
    contourFilter->Update();

    typedef AddImageFilter<ImageType, ImageType, ImageType> AdderType;
    typename AdderType::Pointer adder = AdderType::New();
    adder->SetInput1( contourFilter->GetOutput() );
    adder->SetInput2( this->m_LabelSurfaceImage );
    adder->Update();

    this->m_LabelSurfaceImage = adder->GetOutput();
    }
  else
    {
    typedef BinaryThresholdImageFilter<ImageType, ImageType> ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput( this->m_LabelSurfaceImage );
    thresholder->SetLowerThreshold( this->m_ForegroundValue );
    thresholder->SetUpperThreshold( this->m_ForegroundValue );
    thresholder->SetInsideValue( this->m_BackgroundValue );
    thresholder->SetOutsideValue( this->m_ForegroundValue );
    thresholder->Update();

    typedef BinaryContourImageFilter<ImageType, ImageType> ContourFilterType;
    typename ContourFilterType::Pointer contourFilter
      = ContourFilterType::New();
    contourFilter->SetInput( thresholder->GetOutput() );
    contourFilter->SetFullyConnected( false );
    contourFilter->SetBackgroundValue( this->m_BackgroundValue );
    contourFilter->SetForegroundValue( this->m_ForegroundValue );
    contourFilter->Update();

    typename ThresholderType::Pointer thresholder2 = ThresholderType::New();
    thresholder2->SetInput( contourFilter->GetOutput() );
    thresholder2->SetLowerThreshold( this->m_ForegroundValue );
    thresholder2->SetUpperThreshold( this->m_ForegroundValue );
    thresholder2->SetInsideValue( this->m_SurfaceLabel );
    thresholder2->SetOutsideValue( this->m_BackgroundValue );
    thresholder2->Update();

    typedef AddImageFilter<ImageType, ImageType, ImageType> AdderType;
    typename AdderType::Pointer adder = AdderType::New();
    adder->SetInput1( thresholder2->GetOutput() );
    adder->SetInput2( this->m_LabelSurfaceImage );
    adder->Update();

    this->m_LabelSurfaceImage = adder->GetOutput();
    }
}

template<class TImage>
bool
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::IsChangeSafe( IndexType idx )
{
  bool wellComposed = false;
  if( ImageDimension == 2 )
    {
    wellComposed = this->IsChangeWellComposed2D( idx );
    }
  else
    {
    wellComposed = this->IsChangeWellComposed3D( idx );
    }

  if( wellComposed && !this->IsCriticalTopologicalConfiguration( idx ) )
    {
    return true;
    }

  return false;
}


/*
 * 2-D
 */
template<class TImage>
bool
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::IsChangeWellComposed2D( IndexType idx )
{
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->m_LabelSurfaceImage,
    this->m_LabelSurfaceImage->GetRequestedRegion() );
  It.SetLocation( idx );

  Array<short> neighborhoodPixels( 9 );

  PixelType checkValue = this->m_ForegroundValue;
  if( this->m_IsInversionStep )
    {
    checkValue = this->m_BackgroundValue; 
    }

  // Check for critical configurations: 4 90-degree rotations

  for ( unsigned int i = 0; i < 4; i++ )
    {
    for ( unsigned int j = 0; j < 9; j++ )
      {
      neighborhoodPixels[j] =
        ( It.GetPixel( this->m_RotationIndices[i][j] ) != checkValue );
      if( this->m_RotationIndices[i][j] == 4 )
        {
        neighborhoodPixels[j] = !neighborhoodPixels[j];
        }
      }

    if( this->IsCriticalC1Configuration2D( neighborhoodPixels )
      || this->IsCriticalC2Configuration2D( neighborhoodPixels )
      || this->IsCriticalC3Configuration2D( neighborhoodPixels )
      || this->IsCriticalC4Configuration2D( neighborhoodPixels ) )
      {
      return false;
      }
    }

  // Check for critical configurations: 2 reflections
  //  Note that the reflections for the C1 and C2 cases
  //  are covered by the rotation cases above (except
  //  in the case of FullInvariance == false.

  for ( unsigned int i = 0; i < 2; i++ )
    {
    for ( unsigned int j = 0; j < 9; j++ )
      {
      neighborhoodPixels[j] =
        ( It.GetPixel( this->m_ReflectionIndices[i][j] ) != checkValue );
      if( this->m_ReflectionIndices[i][j] == 4 )
        {
        neighborhoodPixels[j] = !neighborhoodPixels[j];
        }
      }
//    if( !this->m_FullInvariance
//      && ( this->IsCriticalC1Configuration2D( neighborhoodPixels )
//        || this->IsCriticalC2Configuration2D( neighborhoodPixels ) ) )
//      {
//      return false;
//      }
    if( this->IsCriticalC3Configuration2D( neighborhoodPixels )
      || this->IsCriticalC4Configuration2D( neighborhoodPixels ) )
      {
      return false;
      }
    }
  return true;
}

template<class TImage>
bool
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::IsCriticalTopologicalConfiguration( IndexType idx )
{
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->m_LabelSurfaceImage,
    this->m_LabelSurfaceImage->GetRequestedRegion() );
  It.SetLocation( idx );

  PixelType checkValue = this->m_ForegroundValue;
  if( this->m_IsInversionStep )
    {
    checkValue = this->m_BackgroundValue; 
    }

  unsigned int numberOfCriticalC3Configurations = 0;
  unsigned int numberOfFaces = 0;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    if( It.GetNext( d ) == checkValue )
      {
      numberOfFaces++;
      }
    if( It.GetPrevious( d ) == checkValue )
      {
      numberOfFaces++;
      }
    if( It.GetNext( d ) == checkValue && It.GetPrevious( d ) == checkValue )
      {
      numberOfCriticalC3Configurations++;
      }
    }

  if( numberOfCriticalC3Configurations > 0 && numberOfFaces % 2 == 0
      && numberOfCriticalC3Configurations * 2 == numberOfFaces )
    {
    return true;
    }
  return false;
}

/*
 * 2-D
 */
template<class TImage>
bool
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::IsCriticalC1Configuration2D( Array<short> neighborhood )
{
  return ( !neighborhood[0] &&  neighborhood[1] &&
            neighborhood[3] && !neighborhood[4] &&
           !neighborhood[8] );
}

/*
 * 2-D
 */
template<class TImage>
bool
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::IsCriticalC2Configuration2D( Array<short> neighborhood )
{
  return ( !neighborhood[0] &&  neighborhood[1] &&
            neighborhood[3] && !neighborhood[4] &&
            neighborhood[8] &&
           ( neighborhood[5] || neighborhood[7] ) );
}

/*
 * 2-D
 */
template<class TImage>
bool
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::IsCriticalC3Configuration2D( Array<short> neighborhood )
{
  return ( !neighborhood[0] &&  neighborhood[1] &&
            neighborhood[3] && !neighborhood[4] &&
           !neighborhood[5] &&  neighborhood[6] &&
           !neighborhood[7] &&  neighborhood[8] );
}

/*
 * 2-D
 */
template<class TImage>
bool
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::IsCriticalC4Configuration2D( Array<short> neighborhood )
{
  return ( !neighborhood[0] &&  neighborhood[1] &&
            neighborhood[3] && !neighborhood[4] &&
           !neighborhood[5] && !neighborhood[6] &&
           !neighborhood[7] &&  neighborhood[8] );
}

/*
 * 2-D
 */
template<class TImage>
bool
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::IsSpecialCaseOfC4Configuration2D( PixelType label, IndexType idx,
                                    IndexType idx6, IndexType idx7 )
{
  IndexType idxa;
  IndexType idxb;
  for ( unsigned int j = 0; j < 2; j++ )
    {
    idxa[j] = idx7[j] + ( idx7[j] - idx[j] );
    idxb[j] = idx6[j] + ( idx7[j] - idx[j] );
    }
  return ( this->GetOutput()->GetRequestedRegion().IsInside( idxa ) &&
           this->GetOutput()->GetRequestedRegion().IsInside( idxb ) &&
           this->m_LabelSurfaceImage->GetPixel( idxa ) <  label &&
           this->m_LabelSurfaceImage->GetPixel( idxb ) >= label );
}

/*
 * 2-D
 */
template<class TImage>
void
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::InitializeIndices2D()
{
  this->m_RotationIndices[0].SetSize( 9 );
  this->m_RotationIndices[1].SetSize( 9 );
  this->m_RotationIndices[2].SetSize( 9 );
  this->m_RotationIndices[3].SetSize( 9 );

  this->m_RotationIndices[0][0] = 0;
  this->m_RotationIndices[0][1] = 1;
  this->m_RotationIndices[0][2] = 2;
  this->m_RotationIndices[0][3] = 3;
  this->m_RotationIndices[0][4] = 4;
  this->m_RotationIndices[0][5] = 5;
  this->m_RotationIndices[0][6] = 6;
  this->m_RotationIndices[0][7] = 7;
  this->m_RotationIndices[0][8] = 8;

  this->m_RotationIndices[1][0] = 2;
  this->m_RotationIndices[1][1] = 5;
  this->m_RotationIndices[1][2] = 8;
  this->m_RotationIndices[1][3] = 1;
  this->m_RotationIndices[1][4] = 4;
  this->m_RotationIndices[1][5] = 7;
  this->m_RotationIndices[1][6] = 0;
  this->m_RotationIndices[1][7] = 3;
  this->m_RotationIndices[1][8] = 6;

  this->m_RotationIndices[2][0] = 8;
  this->m_RotationIndices[2][1] = 7;
  this->m_RotationIndices[2][2] = 6;
  this->m_RotationIndices[2][3] = 5;
  this->m_RotationIndices[2][4] = 4;
  this->m_RotationIndices[2][5] = 3;
  this->m_RotationIndices[2][6] = 2;
  this->m_RotationIndices[2][7] = 1;
  this->m_RotationIndices[2][8] = 0;

  this->m_RotationIndices[3][0] = 6;
  this->m_RotationIndices[3][1] = 3;
  this->m_RotationIndices[3][2] = 0;
  this->m_RotationIndices[3][3] = 7;
  this->m_RotationIndices[3][4] = 4;
  this->m_RotationIndices[3][5] = 1;
  this->m_RotationIndices[3][6] = 8;
  this->m_RotationIndices[3][7] = 5;
  this->m_RotationIndices[3][8] = 2;

  this->m_ReflectionIndices[0].SetSize( 9 );
  this->m_ReflectionIndices[1].SetSize( 9 );

  this->m_ReflectionIndices[0][0] = 6;
  this->m_ReflectionIndices[0][1] = 7;
  this->m_ReflectionIndices[0][2] = 8;
  this->m_ReflectionIndices[0][3] = 3;
  this->m_ReflectionIndices[0][4] = 4;
  this->m_ReflectionIndices[0][5] = 5;
  this->m_ReflectionIndices[0][6] = 0;
  this->m_ReflectionIndices[0][7] = 1;
  this->m_ReflectionIndices[0][8] = 2;

  this->m_ReflectionIndices[1][0] = 2;
  this->m_ReflectionIndices[1][1] = 1;
  this->m_ReflectionIndices[1][2] = 0;
  this->m_ReflectionIndices[1][3] = 5;
  this->m_ReflectionIndices[1][4] = 4;
  this->m_ReflectionIndices[1][5] = 3;
  this->m_ReflectionIndices[1][6] = 8;
  this->m_ReflectionIndices[1][7] = 7;
  this->m_ReflectionIndices[1][8] = 6;
}

/*
 * 3-D
 */
template<class TImage>
bool
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::IsChangeWellComposed3D( IndexType idx )
{
  Array<short> neighborhoodPixels( 8 );

  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->m_LabelSurfaceImage,
    this->m_LabelSurfaceImage->GetRequestedRegion() );
  It.SetLocation( idx );

  PixelType checkValue = this->m_ForegroundValue;
  if( this->m_IsInversionStep )
    {
    checkValue = this->m_BackgroundValue; 
    }

  // Check for C1 critical configurations
  for ( unsigned int i = 0; i < 12; i++ )
    {
    for ( unsigned int j = 0; j < 4; j++ )
      {
      neighborhoodPixels[j]
        = ( It.GetPixel( this->m_C1Indices[i][j] ) == checkValue );
      if( this->m_C1Indices[i][j] == 13 )
        {
        neighborhoodPixels[j] = !neighborhoodPixels[j];
        }
      }
    if( this->IsCriticalC1Configuration3D( neighborhoodPixels ) )
      {
      return false;
      }
    }

  // Check for C2 critical configurations
  for ( unsigned int i = 0; i < 8; i++ )
    {
    for ( unsigned int j = 0; j < 8; j++ )
      {
      neighborhoodPixels[j]
        = ( It.GetPixel( this->m_C2Indices[i][j] ) == checkValue );
      if( this->m_C2Indices[i][j] == 13 )
        {
        neighborhoodPixels[j] = !neighborhoodPixels[j];
        }
      }
    if( this->IsCriticalC2Configuration3D( neighborhoodPixels ) )
      {
      return false;
      }
    }

  return true;
}


/*
 * 3-D
 */
template<class TImage>
bool
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::IsCriticalC1Configuration3D( Array<short> neighborhood )
{
  return ( (  neighborhood[0] &&  neighborhood[1] &&
             !neighborhood[2] && !neighborhood[3] ) ||
           ( !neighborhood[0] && !neighborhood[1] &&
              neighborhood[2] &&  neighborhood[3] ) );
}

/*
 * 3-D
 */
template<class TImage>
unsigned int
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::IsCriticalC2Configuration3D( Array<short> neighborhood )
{
  // Check if Type 1 or Type 2
  for ( unsigned int i = 0; i < 4; i++ )
    {
    bool isC2 = false;
    if( neighborhood[2*i] == neighborhood[2*i+1] )
      {
      isC2 = true;
      for ( unsigned int j = 0; j < 8; j++ )
        {
        if( neighborhood[j] == neighborhood[2*i] &&
               j != 2*i && j != 2*i+1 )
          {
          isC2 = false;
          }
        }
      }
    if( isC2 )
      {
      if( neighborhood[2*i] )
        {
        return 1;
        }
      else
        {
        return 2;
        }
      }
    }

  return 0;
}

/*
 * 3-D
 */
template<class TImage>
void
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::InitializeIndices3D()
{
  for ( unsigned int i = 0; i <  12; i++ )
    {
    this->m_C1Indices[i].SetSize( 4 );
    }
  for ( unsigned int i = 0; i <  8; i++ )
    {
    this->m_C2Indices[i].SetSize( 8 );
    }

  this->m_C1Indices[0][0] = 1;
  this->m_C1Indices[0][1] = 13;
  this->m_C1Indices[0][2] = 4;
  this->m_C1Indices[0][3] = 10;

  this->m_C1Indices[1][0] = 9;
  this->m_C1Indices[1][1] = 13;
  this->m_C1Indices[1][2] = 10;
  this->m_C1Indices[1][3] = 12;

  this->m_C1Indices[2][0] = 3;
  this->m_C1Indices[2][1] = 13;
  this->m_C1Indices[2][2] = 4;
  this->m_C1Indices[2][3] = 12;

  this->m_C1Indices[3][0] = 4;
  this->m_C1Indices[3][1] = 14;
  this->m_C1Indices[3][2] = 5;
  this->m_C1Indices[3][3] = 13;

  this->m_C1Indices[4][0] = 12;
  this->m_C1Indices[4][1] = 22;
  this->m_C1Indices[4][2] = 13;
  this->m_C1Indices[4][3] = 21;

  this->m_C1Indices[5][0] = 13;
  this->m_C1Indices[5][1] = 23;
  this->m_C1Indices[5][2] = 14;
  this->m_C1Indices[5][3] = 22;

  this->m_C1Indices[6][0] = 4;
  this->m_C1Indices[6][1] = 16;
  this->m_C1Indices[6][2] = 7;
  this->m_C1Indices[6][3] = 13;

  this->m_C1Indices[7][0] = 13;
  this->m_C1Indices[7][1] = 25;
  this->m_C1Indices[7][2] = 16;
  this->m_C1Indices[7][3] = 22;

  this->m_C1Indices[8][0] = 10;
  this->m_C1Indices[8][1] = 22;
  this->m_C1Indices[8][2] = 13;
  this->m_C1Indices[8][3] = 19;

  this->m_C1Indices[9][0] = 12;
  this->m_C1Indices[9][1] = 16;
  this->m_C1Indices[9][2] = 13;
  this->m_C1Indices[9][3] = 15;

  this->m_C1Indices[10][0] = 13;
  this->m_C1Indices[10][1] = 17;
  this->m_C1Indices[10][2] = 14;
  this->m_C1Indices[10][3] = 16;

  this->m_C1Indices[11][0] = 10;
  this->m_C1Indices[11][1] = 14;
  this->m_C1Indices[11][2] = 11;
  this->m_C1Indices[11][3] = 13;

  this->m_C2Indices[0][0] = 0;
  this->m_C2Indices[0][1] = 13;
  this->m_C2Indices[0][2] = 1;
  this->m_C2Indices[0][3] = 12;
  this->m_C2Indices[0][4] = 3;
  this->m_C2Indices[0][5] = 10;
  this->m_C2Indices[0][6] = 4;
  this->m_C2Indices[0][7] = 9;

  this->m_C2Indices[4][0] = 9;
  this->m_C2Indices[4][1] = 22;
  this->m_C2Indices[4][2] = 10;
  this->m_C2Indices[4][3] = 21;
  this->m_C2Indices[4][4] = 12;
  this->m_C2Indices[4][5] = 19;
  this->m_C2Indices[4][6] = 13;
  this->m_C2Indices[4][7] = 18;

  for ( unsigned int i = 1; i < 4; i++ )
    {
    int addend;
    if( i == 2 )
      {
      addend = 2;
      }
    else
      {
      addend = 1;
      }
    for ( unsigned int j = 0; j < 8; j++ )
      {
      this->m_C2Indices[i  ][j] = this->m_C2Indices[i-1][j] + addend;
      this->m_C2Indices[i+4][j] = this->m_C2Indices[i+3][j] + addend;
      }
    }
}

template <class TImage>
void
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::FindLabelOrderingForObjectGluing()
{
  LabelArrayType labels;

  typedef NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );

  NeighborhoodIteratorType ItN( radius, this->GetOutput(),
    this->GetOutput()->GetLargestPossibleRegion() );

  for( ItN.GoToBegin(); !ItN.IsAtEnd(); ++ItN )
    {
    PixelType centerPixel = ItN.GetCenterPixel();
    if( centerPixel != this->m_BackgroundValue &&
      find( labels.begin(), labels.end(), centerPixel )
      == labels.end() )
      {
      labels.push_back( centerPixel );
      }
    }

  if( labels.size() <= 2 )
    {
    this->m_OrderedObjectLabels = labels;
    return;
    }

  /**
   * Iterate through the image and count the neighbors.  This won't
   * give an exact count but we're only after a relative approximation.
   */

  VariableSizeMatrix<unsigned long> labelNeighborCount;
  labelNeighborCount.SetSize( labels.size(), labels.size() );
  labelNeighborCount.Fill( 0 );

  ImageRegionIterator<RealImageType> ItT( this->m_TargetImage,
    this->m_TargetImage->GetLargestPossibleRegion() );

  unsigned int numberOfNeighbors = ItN.GetNeighborhood().Size();
  for( ItN.GoToBegin(), ItT.GoToBegin(); !ItN.IsAtEnd(); ++ItN, ++ItT )
    {
    PixelType centerPixel = ItN.GetCenterPixel();
    if( centerPixel != this->m_BackgroundValue ||
      ItT.Get() < this->m_ThresholdValue )
      {
      continue;
      }
    for( unsigned int m = 0; m < numberOfNeighbors; m++ )
      {
      bool mInBounds = false;
      PixelType mPixel = ItN.GetPixel( m, mInBounds );
      if( !mInBounds ||
        m == static_cast<unsigned int>( 0.5 * numberOfNeighbors ) ||
        mPixel == this->m_BackgroundValue )
        {
        continue;
        }
      typename LabelArrayType::iterator mloc
        = find( labels.begin(), labels.end(), mPixel );
      for( unsigned int n = 0; n < numberOfNeighbors; n++ )
        {
        bool nInBounds = false;
        PixelType nPixel = ItN.GetPixel( n, nInBounds );
        if( mPixel == nPixel || !nInBounds ||
          n == static_cast<unsigned int>( 0.5 * numberOfNeighbors ) ||
          nPixel == this->m_BackgroundValue )
          {
          continue;
          }
        typename LabelArrayType::iterator nloc
          = find( labels.begin(), labels.end(), nPixel );
        labelNeighborCount( mloc - labels.begin(), nloc - labels.begin() )++;
        labelNeighborCount( nloc - labels.begin(), mloc - labels.begin() )++;
        }
      }
    }

  unsigned int max_i = 0;
  unsigned int max_j = 0;
  unsigned long maxCount = 0;

  for( unsigned int i = 0; i < labelNeighborCount.Rows(); i++ )
    {
    for( unsigned int j = 0; j < labelNeighborCount.Cols(); j++ )
      {
      if( labelNeighborCount( i, j ) > maxCount )
        {
        maxCount = labelNeighborCount( i, j );
        max_i = i;
        max_j = j;
        }
      }
    }

  this->m_OrderedObjectLabels.clear();
  this->m_OrderedObjectLabels.push_back( labels[max_i] );
  this->m_OrderedObjectLabels.push_back( labels[max_j] );

  while( this->m_OrderedObjectLabels.size() < labels.size() )
    {
    maxCount = 0;
    typename LabelArrayType::iterator it;
    for( it = labels.begin(); it != labels.end(); ++it )
      {
      unsigned int j = it - labels.begin();
      unsigned long count = 0;
      if( find( this->m_OrderedObjectLabels.begin(),
        this->m_OrderedObjectLabels.end(), *it ) ==
        this->m_OrderedObjectLabels.end() )
        {
        for( unsigned int i = 0; i < labelNeighborCount.Rows(); i++ )
          {
          if( labelNeighborCount( i, j ) > 0 &&
            find( this->m_OrderedObjectLabels.begin(),
            this->m_OrderedObjectLabels.end(),
            labels[i] ) != this->m_OrderedObjectLabels.end() )
            {
            count += labelNeighborCount( i, j );
            }
          }
        }
      if( count > maxCount )
        {
        maxCount = count;
        max_j = j;
        }
      }
    this->m_OrderedObjectLabels.push_back( labels[max_j] );
    }
}

template <class TImage>
void
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::GlueObjects()
{
  this->m_GluingLabels.clear();

  typename ImageType::Pointer joinedImage = ImageType::New();

  typedef ImageDuplicator<ImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( this->GetOutput() );
  duplicator->Update();

  joinedImage = duplicator->GetOutput();

  typedef NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );

  typedef BinaryThresholdImageFilter<RealImageType, ImageType>
    ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( this->m_TargetImage );
  thresholder->SetInsideValue( NumericTraits<PixelType>::One );
  thresholder->SetOutsideValue( NumericTraits<PixelType>::Zero );
  thresholder->SetUpperThreshold( NumericTraits<RealType>::max() );
  thresholder->SetLowerThreshold( this->m_ThresholdValue );
  thresholder->Update();

  typedef SignedMaurerDistanceMapImageFilter<ImageType, RealImageType>
    DistancerType;
  typename DistancerType::Pointer distancer = DistancerType::New();
  distancer->SetInput( thresholder->GetOutput() );
  distancer->SetSquaredDistance( true );
  distancer->SetUseImageSpacing( true );
  distancer->SetInsideIsPositive( true );
  distancer->SetNumberOfThreads( 1 );
  distancer->Update();

  /**
   * Step through the ordered this->m_OrderedObjectLabels and glue
   * label i to label i+1
   */

  for( unsigned int i = 0; i < this->m_OrderedObjectLabels.size()-1; i++ )
    {
    std::cout << "Joining objects with labels "
      << this->m_OrderedObjectLabels[i] << " and "
      << this->m_OrderedObjectLabels[i+1]
      << "." << std::endl;

    PixelType currentGluingLabel = NumericTraits<PixelType>::Zero;
    while( currentGluingLabel++ < NumericTraits<PixelType>::max() )
      {
      if( find( this->m_OrderedObjectLabels.begin(),
        this->m_OrderedObjectLabels.end(), currentGluingLabel ) ==
        this->m_OrderedObjectLabels.end() &&
        find( this->m_GluingLabels.begin(),
        this->m_GluingLabels.end(), currentGluingLabel ) ==
        this->m_GluingLabels.end() )
        {
        this->m_GluingLabels.push_back( currentGluingLabel );
        break;
        }
      }

    typename RealImageType::Pointer targetImage = RealImageType::New();
    targetImage->SetRegions( joinedImage->GetRequestedRegion() );
    targetImage->SetOrigin( joinedImage->GetOrigin() );
    targetImage->SetSpacing( joinedImage->GetSpacing() );
    targetImage->SetDirection( joinedImage->GetDirection() );
    targetImage->Allocate();
    targetImage->FillBuffer( NumericTraits<RealType>::Zero );

    typename ImageType::Pointer sourceImage = ImageType::New();
    sourceImage->SetRegions( this->GetInput()->GetRequestedRegion() );
    sourceImage->SetOrigin( this->GetInput()->GetOrigin() );
    sourceImage->SetSpacing( this->GetInput()->GetSpacing() );
    sourceImage->SetDirection( this->GetInput()->GetDirection() );
    sourceImage->Allocate();
    sourceImage->FillBuffer( this->m_BackgroundValue );

    /**
     * Iterate through the image and find the 'candidate pixels'.  These are
     * defined as those pixels which have pixels with labels
     * his->m_OrderedObjectLabels[i] and his->m_OrderedObjectLabels[i+1]
     * and no other this->m_OrderedObjectLabels in their neighborhood.
     */
    unsigned long numberOfCandidatePixels = 0;

    NeighborhoodIterator<ImageType> ItN( radius, joinedImage,
      joinedImage->GetRequestedRegion() );
    unsigned int numberOfNeighbors = ItN.GetNeighborhood().Size();
    for( ItN.GoToBegin(); !ItN.IsAtEnd(); ++ItN )
      {
      if( ItN.GetCenterPixel() == this->m_OrderedObjectLabels[i] ||
        ItN.GetCenterPixel() == this->m_OrderedObjectLabels[i+1] )
        {
        sourceImage->SetPixel( ItN.GetIndex(),
          this->m_OrderedObjectLabels[i+1] );
        }
      else if( ItN.GetCenterPixel() == this->m_BackgroundValue &&
        this->m_TargetImage->GetPixel( ItN.GetIndex() )
        > this->m_ThresholdValue )
        {
        bool foundLabel_i = false;
        bool foundLabel_ip1 = false;
        bool foundOtherLabel = false;
        for( unsigned int n = 0; n < numberOfNeighbors; n++ )
          {
          if( n != static_cast<unsigned int>( 0.5 * numberOfNeighbors ) )
            {
            bool inBounds = false;
            PixelType pixel = ItN.GetPixel( n, inBounds );

            if( inBounds && pixel == this->m_OrderedObjectLabels[i] )
              {
              foundLabel_i = true;
              }
            else if( inBounds && pixel == this->m_OrderedObjectLabels[i+1] )
              {
              foundLabel_ip1 = true;
              }
            else if( inBounds && pixel != this->m_BackgroundValue )
              {
              foundOtherLabel = true;
              }
            }
          }
        if( foundLabel_i && foundLabel_ip1 && !foundOtherLabel )
          {
          targetImage->SetPixel( ItN.GetIndex(),
            NumericTraits<RealType>::One );
          numberOfCandidatePixels++;
          }
        }
      }

    /**
     * If there are 'candidate' pixels , we need to iterate through the image
     * again and find those 'candidate' pixels which violate the Jordan curve
     * theorem, i.e. those pixels which have more than two components in
     * their immediate neighborhood.
     */

    if( numberOfCandidatePixels > 0 )
      {
      IndexContainerType nonCandidateIndices;
      nonCandidateIndices.clear();

      typename ImageType::Pointer neighborhoodImage = ImageType::New();
      typename ImageType::SizeType neighborhoodSize;
      typename ImageType::IndexType neighborhoodIndex;
      neighborhoodSize.Fill( 3 );
      neighborhoodIndex.Fill( 1 );
      neighborhoodImage->SetRegions( neighborhoodSize );
      neighborhoodImage->Allocate();

      NeighborhoodIterator<RealImageType> ItN( radius, targetImage,
        targetImage->GetRequestedRegion() );
      unsigned int numberOfNeighbors = ItN.GetNeighborhood().Size();

      for( ItN.GoToBegin(); !ItN.IsAtEnd(); ++ItN )
        {
        if( ItN.GetCenterPixel() == NumericTraits<RealType>::One )
          {
          neighborhoodImage->FillBuffer( NumericTraits<PixelType>::Zero );
          for( unsigned int n = 0; n < numberOfNeighbors; n++ )
            {
            bool inBounds = false;
            RealType neighborhoodPixel = ItN.GetPixel( n, inBounds );
            if( inBounds )
              {
              if( neighborhoodPixel == NumericTraits<RealType>::One ||
                this->GetOutput()->GetPixel( ItN.GetIndex( n ) )
                == NumericTraits<PixelType>::Zero )
                {
                neighborhoodImage->SetPixel( neighborhoodIndex
                  + ItN.GetOffset( n ), NumericTraits<PixelType>::Zero );
                }
              else
                {
                neighborhoodImage->SetPixel( neighborhoodIndex
                  + ItN.GetOffset( n ), NumericTraits<PixelType>::One );
                }
              }
            else
              {
              neighborhoodImage->SetPixel( neighborhoodIndex
                + ItN.GetOffset( n ), NumericTraits<PixelType>::Zero );
              }
            }
          typedef ConnectedComponentImageFilter<ImageType, ImageType>
            ConnecterType;
          typename ConnecterType::Pointer connecter = ConnecterType::New();
          connecter->SetInput( neighborhoodImage );
          connecter->Update();

          typedef RelabelComponentImageFilter<ImageType, ImageType>
            RelabelerType;
          typename RelabelerType::Pointer relabeler = RelabelerType::New();
          relabeler->SetInput( connecter->GetOutput() );
          relabeler->Update();

          if( relabeler->GetNumberOfObjects() != 2 )
            {
            nonCandidateIndices.push_back( ItN.GetIndex() );
            }
          }
        }
      typename IndexContainerType::iterator it;
      for( it = nonCandidateIndices.begin(); it != nonCandidateIndices.end();
        ++it )
        {
        targetImage->SetPixel( *it, NumericTraits<RealType>::Zero );
        }

      /**
       * After removing those pixels from candidacy which violate the Jordan
       * curve theorem, we find the largest target component which will
       * comprise potential gluing region.
       */

      typedef BinaryThresholdImageFilter<RealImageType, ImageType>
        ThresholderType;
      typename ThresholderType::Pointer thresholder =
        ThresholderType::New();
      thresholder->SetInput( targetImage );
      thresholder->SetLowerThreshold( 0.5 );
      thresholder->SetUpperThreshold( 1.5 );
      thresholder->SetInsideValue( 1 );
      thresholder->SetOutsideValue( 0 );
      thresholder->Update();

      typedef ConnectedComponentImageFilter<ImageType, ImageType>
        ConnecterType;
      typename ConnecterType::Pointer connecter = ConnecterType::New();
      connecter->SetInput( thresholder->GetOutput() );
      connecter->Update();

      typedef RelabelComponentImageFilter<ImageType, ImageType>
        RelabelerType;
      typename RelabelerType::Pointer relabeler = RelabelerType::New();
      relabeler->SetInput( connecter->GetOutput() );
      relabeler->Update();

      ImageRegionIterator<ImageType> ItR( relabeler->GetOutput(),
        relabeler->GetOutput()->GetRequestedRegion() );
      ImageRegionIterator<RealImageType> ItT( targetImage,
        targetImage->GetRequestedRegion() );
      for( ItR.GoToBegin(), ItT.GoToBegin(); !ItR.IsAtEnd(); ++ItR, ++ItT )
        {
        if( ItR.Get() > NumericTraits<PixelType>::One )
          {
          ItT.Set( NumericTraits<RealType>::Zero );
          ItR.Set( NumericTraits<PixelType>::Zero );
          }
        }

      /**
       * Check for holes in a plane.
       */
      if( this->m_GluingStrategy == 0 )
        {
        RealType maxDistance = NumericTraits<RealType>::NonpositiveMin();
        IndexType bestIndex = sourceImage->GetRequestedRegion().GetIndex();
  
        ImageRegionIterator<ImageType> ItS( sourceImage,
          sourceImage->GetRequestedRegion() );
        for( ItN.GoToBegin(), ItS.GoToBegin(); !ItN.IsAtEnd(); ++ItN, ++ItS )
          {
          if( ItN.GetCenterPixel() == NumericTraits<RealType>::One )
            {
            bool wellComposed = false;
            if( ImageDimension == 2 )
              {
              wellComposed = this->IsChangeWellComposed2D( ItN.GetIndex() );
              }
            else
              {
              wellComposed = this->IsChangeWellComposed3D( ItN.GetIndex() );
              }
            if( wellComposed )
              {
              neighborhoodImage->FillBuffer( NumericTraits<PixelType>::Zero );
              for( unsigned int n = 0; n < numberOfNeighbors; n++ )
                {
                bool inBounds = false;
                RealType neighborhoodPixel = ItN.GetPixel( n, inBounds );
                if( inBounds )
                  {
                  if( neighborhoodPixel == NumericTraits<RealType>::One ||
                    this->GetOutput()->GetPixel( ItN.GetIndex( n ) )
                    == NumericTraits<PixelType>::Zero )
                    {
                    neighborhoodImage->SetPixel( neighborhoodIndex
                      + ItN.GetOffset( n ), NumericTraits<PixelType>::Zero );
                    }
                  else
                    {
                    neighborhoodImage->SetPixel( neighborhoodIndex
                      + ItN.GetOffset( n ), NumericTraits<PixelType>::One );
                    }
                  }
                else
                  {
                  neighborhoodImage->SetPixel( neighborhoodIndex
                    + ItN.GetOffset( n ), NumericTraits<PixelType>::Zero );
                  }
                }
              typedef ConnectedComponentImageFilter<ImageType, ImageType>
                ConnecterType;
              typename ConnecterType::Pointer connecter = ConnecterType::New();
              connecter->SetInput( neighborhoodImage );
              connecter->Update();
  
              typedef RelabelComponentImageFilter<ImageType, ImageType>
                RelabelerType;
              typename RelabelerType::Pointer relabeler = RelabelerType::New();
              relabeler->SetInput( connecter->GetOutput() );
              relabeler->Update();
  
              if( relabeler->GetNumberOfObjects() == 2 )
                {
                RealType distance
                  = distancer->GetOutput()->GetPixel( ItN.GetIndex() );
                if( maxDistance < distance )
                  {
                  maxDistance = distance;
                  bestIndex = ItN.GetIndex();
                  }
                }
              }
            }
          }
  
        sourceImage->SetPixel( bestIndex, this->m_OrderedObjectLabels[i+1] );
  
        typename Self::Pointer topologyFilter = Self::New();
        topologyFilter->SetSourceImage( sourceImage );
        topologyFilter->SetTargetImage( targetImage );
        topologyFilter->SetNumberOfIterations(
          NumericTraits<unsigned int>::max() );
        topologyFilter->SetGlueObjects( false );
        topologyFilter->SetGrowOnly( true );
        topologyFilter->Update();
  
        sourceImage = topologyFilter->GetOutput();
        }
      else if( this->m_GluingStrategy == 1 )
        {  
        typedef BinaryThresholdImageFilter<ImageType, ImageType>
          InverterType;
        typename InverterType::Pointer inverter = InverterType::New();
        inverter->SetInput( relabeler->GetOutput() );
        inverter->SetInsideValue( NumericTraits<PixelType>::Zero );
        inverter->SetOutsideValue( NumericTraits<PixelType>::One );
        inverter->SetLowerThreshold( 1 );
        inverter->SetUpperThreshold( 1 );
        inverter->Update();
  
        typedef WellComposedImageFilter<ImageType> WellComposerType;
        typename WellComposerType::Pointer wellComposer = WellComposerType::New();
        wellComposer->SetInput( inverter->GetOutput() );
        wellComposer->SetTotalNumberOfLabels( 2 );
  //      wellComposer->Update();
  
        typename InverterType::Pointer reinverter = InverterType::New();
        reinverter->SetInput( wellComposer->GetOutput() );
        reinverter->SetInsideValue( NumericTraits<PixelType>::Zero );
        reinverter->SetOutsideValue( this->m_OrderedObjectLabels[i+1] );
        reinverter->SetLowerThreshold( 1 );
        reinverter->SetUpperThreshold( 1 );
        reinverter->Update();
        
        sourceImage = reinverter->GetOutput();
        }
      else  
        {
        typedef BinaryThresholdImageFilter<ImageType, ImageType>
          ThresholderType;
        typename ThresholderType::Pointer thresholder = ThresholderType::New();
        thresholder->SetInput( relabeler->GetOutput() );
        thresholder->SetInsideValue( this->m_OrderedObjectLabels[i+1] );
        thresholder->SetOutsideValue( NumericTraits<PixelType>::Zero );
        thresholder->SetLowerThreshold( NumericTraits<PixelType>::One );
        thresholder->SetUpperThreshold( NumericTraits<PixelType>::One );
        thresholder->Update();

        sourceImage = thresholder->GetOutput(); 
        }
      }
 
    ImageRegionIterator<ImageType> ItJ( joinedImage,
      joinedImage->GetRequestedRegion() );
    ImageRegionIteratorWithIndex<ImageType> ItO( this->GetOutput(),
      this->GetOutput()->GetRequestedRegion() );
    ImageRegionIterator<ImageType> ItS( sourceImage,
      sourceImage->GetRequestedRegion() );

    ItJ.GoToBegin();
    ItO.GoToBegin();
    ItS.GoToBegin();
    while( !ItJ.IsAtEnd() )
      {
      if( ItS.Get() == this->m_OrderedObjectLabels[i+1] &&
        ItO.Get() == this->m_BackgroundValue )
        {
        ItO.Set( currentGluingLabel );
        ItJ.Set( this->m_OrderedObjectLabels[i+1] );
        }
      if( ItJ.Get() == this->m_OrderedObjectLabels[i] )
        {
        ItJ.Set( this->m_OrderedObjectLabels[i+1] );
        }
      ++ItJ;
      ++ItO;
      ++ItS;
      }
    }
}

template <class TImage>
void
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::FillCavities()
{

  std::cout << "Filling cavities." << std::endl;

  this->m_CavityLabels.clear();

  /**
   * Create cavity image
   */

  typename ImageType::Pointer cavityImage = ImageType::New();
  cavityImage->SetRegions( this->GetInput()->GetRequestedRegion() );
  cavityImage->SetOrigin( this->GetInput()->GetOrigin() );
  cavityImage->SetSpacing( this->GetInput()->GetSpacing() );
  cavityImage->SetDirection( this->GetInput()->GetDirection() );
  cavityImage->Allocate();

  ImageRegionIterator<ImageType> ItO( this->GetOutput(),
    this->GetOutput()->GetLargestPossibleRegion() );
  ImageRegionIterator<ImageType> ItC( cavityImage,
    cavityImage->GetRequestedRegion() );
  for( ItO.GoToBegin(), ItC.GoToBegin(); !ItO.IsAtEnd(); ++ItC, ++ItO )
    {
    if( ItO.Get() != this->m_BackgroundValue )
      {
      ItC.Set( this->m_ForegroundValue );
      }
    }

  typename Self::Pointer topologyFilter = Self::New();
  topologyFilter->SetSourceImage( cavityImage );
  topologyFilter->SetTargetImage( this->m_TargetImage );
  topologyFilter->SetNumberOfIterations(
    NumericTraits<unsigned int>::max() );
  topologyFilter->SetGrowOnly( true );
  topologyFilter->SetGlueObjects( false );
  topologyFilter->SetFillCavities( false );
  topologyFilter->Update();

  ImageRegionIterator<ImageType> ItT( topologyFilter->GetOutput(),
    topologyFilter->GetOutput()->GetRequestedRegion() );
  ItO.GoToBegin();
  ItT.GoToBegin();
  ItC.GoToBegin();
  while( !ItO.IsAtEnd() )
    {
    if( ItO.Get() == this->m_BackgroundValue
      && ItT.Get() == this->m_ForegroundValue )
      {
      ItC.Set( this->m_ForegroundValue );
      }
    else
      {
      ItC.Set( this->m_BackgroundValue );
      }
    ++ItO;
    ++ItT;
    ++ItC;
    }

  /**
   * Determine the number of cavities.
   */
  typedef ConnectedComponentImageFilter<ImageType, ImageType> ConnecterType;
  typename ConnecterType::Pointer connecter = ConnecterType::New();
  connecter->SetInput( cavityImage );
  connecter->Update();

  typedef RelabelComponentImageFilter<ImageType, ImageType> RelabelerType;
  typename RelabelerType::Pointer relabeler = RelabelerType::New();
  relabeler->SetInput( connecter->GetOutput() );
  relabeler->Update();

  for( unsigned int m = 0; m < relabeler->GetNumberOfObjects(); m++ )
    {

    PixelType currentCavityLabel = NumericTraits<PixelType>::Zero;
    while( currentCavityLabel++ < NumericTraits<PixelType>::max() )
      {
      if( find( this->m_OrderedObjectLabels.begin(),
        this->m_OrderedObjectLabels.end(), currentCavityLabel ) ==
        this->m_OrderedObjectLabels.end() &&
        find( this->m_GluingLabels.begin(),
        this->m_GluingLabels.end(), currentCavityLabel ) ==
        this->m_GluingLabels.end() &&
        find( this->m_CavityLabels.begin(),
        this->m_CavityLabels.end(), currentCavityLabel ) ==
        this->m_CavityLabels.end() )
        {
        this->m_CavityLabels.push_back( currentCavityLabel );
        break;
        }
      }

    ImageRegionIterator<ImageType> ItR( relabeler->GetOutput(),
      relabeler->GetOutput()->GetRequestedRegion() );
    for( ItO.GoToBegin(), ItR.GoToBegin(); !ItO.IsAtEnd(); ++ItO, ++ItR )
      {
      if( ItR.Get() == static_cast<PixelType>( m+1 ) )
        {
        ItO.Set( currentCavityLabel );
        }
      }
    }
}

template <class TImage>
void
TopologyPreservingDigitalSurfaceEvolutionImageFilter<TImage>
::PrintSelf(
  std::ostream& os,
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Target image: "
    << this->m_TargetImage << std::endl;
  os << indent << "Number of iterations: "
    << this->m_NumberOfIterations << std::endl;
  os << indent << "ThresholdValue: "
    << this->m_ThresholdValue << std::endl;
  if( this->m_GrowOnly )
    {
    std::cout << "Grow only." << std::endl;
    }
  if( this->m_GlueObjects )
    {
    os << indent << "Glue objects. " << std::endl;
    os << indent << "   Object label gluing order: ";
    for( unsigned int n = 0; n < this->m_OrderedObjectLabels.size(); n++ )
      {
      os << this->m_OrderedObjectLabels[n] << " ";
      }
    os << std::endl;
    os << indent << "   Glue labels: ";
    for( unsigned int n = 0; n < this->m_GluingLabels.size(); n++ )
      {
      os << this->m_GluingLabels[n] << " ";
      }
    os << std::endl;
    if( this->m_FillCavities )
      {
      os << indent << "Fill cavities. " << std::endl;
      os << indent << "   Cavity labels: ";
      for( unsigned int n = 0; n < this->m_CavityLabels.size(); n++ )
        {
        os << this->m_CavityLabels[n] << " ";
        }
      os << std::endl;
      }
    }
}

} // end namespace itk

#endif
