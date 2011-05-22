#ifndef __itkWellComposedImageFilter_txx
#define __itkWellComposedImageFilter_txx

#include "itkWellComposedImageFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkConstantPadImageFilter.h"
#include "itkExtractImageFilter.h"

namespace itk
{

template<class TImage>
WellComposedImageFilter<TImage>
::WellComposedImageFilter()
{
  this->m_TotalNumberOfLabels = 2;

  if ( ImageDimension == 2 )
    {
    this->m_FullInvariance = true;
    this->InitializeIndices2D();
    }
  else if ( ImageDimension == 3 )
    {
    this->InitializeIndices3D();
    }
  else
    {
    itkExceptionMacro( << "Image dimension must be equal to 2 or 3." );
    }
}

template<class TImage>
WellComposedImageFilter<TImage>
::~WellComposedImageFilter()
{
}

template<class TImage>
void
WellComposedImageFilter<TImage>
::GenerateData()
{
  this->AllocateOutputs();

  // Pad the boundary to eliminate boundary effects
  typedef ConstantPadImageFilter<ImageType, ImageType> PadderType;
  typename PadderType::Pointer padder = PadderType::New();
  padder->SetInput( this->GetOutput() );
  padder->SetConstant( 0 );
  unsigned long upperfactors[ImageDimension];
  unsigned long lowerfactors[ImageDimension];
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    upperfactors[i] = 1;
    lowerfactors[i] = 1;
    }
  padder->SetPadUpperBound( upperfactors );
  padder->SetPadLowerBound( lowerfactors );
  padder->Update();

  this->GraftOutput( padder->GetOutput() );

  if ( ImageDimension == 2 )
    {
    this->MakeImageWellComposed2D();
    }
  else if ( ImageDimension == 3 )
    {
    this->MakeImageWellComposed3D();
    }

  // Extract the output image region to match input image region
  typedef itk::ExtractImageFilter<ImageType, ImageType> CropperType;
  typename CropperType::Pointer cropper = CropperType::New();
  typename ImageType::RegionType region;
  region.SetSize( this->GetInput()->GetRequestedRegion().GetSize() );
  region.SetIndex( this->GetInput()->GetRequestedRegion().GetIndex() );
  cropper->SetInput( this->GetOutput() );
  cropper->SetExtractionRegion( region );
  cropper->SetDirectionCollapseToSubmatrix();
  cropper->Update();

  this->GraftOutput( cropper->GetOutput() );
}

/*
 * 2-D
 */
template<class TImage>
void
WellComposedImageFilter<TImage>
::MakeImageWellComposed2D()
{
  unsigned long NumberOfC1Configurations2D = 0;
  unsigned long NumberOfC2Configurations2D = 0;
  unsigned long NumberOfC3Configurations2D = 0;
  unsigned long NumberOfC4Configurations2D = 0;

  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
        this->GetOutput()->GetRequestedRegion() );

  Array<char> neighborhoodPixels( 9 );

  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( !It.InBounds() )
      {
      continue;
      }

    PixelType currentLabel = NumericTraits<PixelType>::Zero;

    while ( currentLabel <
      static_cast<PixelType>( this->m_TotalNumberOfLabels ) )
      {
      // Check for critical configurations: 4 90-degree rotations
      for ( unsigned int i = 0; i < 4; i++ )
        {
        for ( unsigned int j = 0; j < 9; j++ )
          {
          neighborhoodPixels[j] =
            ( It.GetPixel( this->m_RotationIndices[i][j] ) == currentLabel );
          }

        if ( this->IsCriticalC1Configuration2D( neighborhoodPixels ) )
          {
          NumberOfC1Configurations2D++;
          for ( int k = currentLabel; k >= 0; k-- )
            {
            if ( this->IsChangeSafe2D( k,
                   It.GetIndex( this->m_RotationIndices[i][4] ) ) )
              {
              It.SetPixel( this->m_RotationIndices[i][4],
                static_cast<PixelType>( k ) );
              currentLabel = k;
              break;
              }
            }
          break;
          }
        else if ( this->IsCriticalC2Configuration2D( neighborhoodPixels ) )
          {
          NumberOfC2Configurations2D++;
          for ( int k = currentLabel; k >= 0; k-- )
            {
            if ( this->IsChangeSafe2D( k,
                   It.GetIndex( this->m_RotationIndices[i][4] ) ) )
              {
              It.SetPixel( this->m_RotationIndices[i][4],
                static_cast<PixelType>( k ) );
              currentLabel = k;
              break;
              }
            }
          break;
          }
        else if ( this->IsCriticalC3Configuration2D( neighborhoodPixels ) )
          {
          NumberOfC3Configurations2D++;
          int k4;
          for ( k4 = currentLabel; k4 >= 0; k4-- )
            {
            if ( this->IsChangeSafe2D( k4,
                   It.GetIndex( this->m_RotationIndices[i][4] ) ) )
              {
              It.SetPixel( this->m_RotationIndices[i][4],
                static_cast<PixelType>( k4 ) );
              break;
              }
            }
          int k7;
          for ( k7 = currentLabel; k7 >= 0; k7-- )
            {
            if ( this->IsChangeSafe2D( k7,
                   It.GetIndex( this->m_RotationIndices[i][7] ) ) )
              {
              It.SetPixel( this->m_RotationIndices[i][7],
                static_cast<PixelType>( k7 ) );
              break;
              }
            }
          currentLabel = vnl_math_min( static_cast<PixelType>( k4 ),
                                       static_cast<PixelType>( k7 ) );
          break;
          }

        else if ( this->IsCriticalC4Configuration2D( neighborhoodPixels ) )
          {
          NumberOfC4Configurations2D++;

          int k4;
          for ( k4 = currentLabel; k4 >= 0; k4-- )
            {
            if ( this->IsChangeSafe2D( k4,
                   It.GetIndex( this->m_RotationIndices[i][4] ) ) )
              {
              It.SetPixel( this->m_RotationIndices[i][4],
                static_cast<PixelType>( k4 ) );
              break;
              }
            }
          int k7;
          for ( k7 = currentLabel; k7 >= 0; k7-- )
            {
            if ( this->IsChangeSafe2D( k7,
                   It.GetIndex( this->m_RotationIndices[i][7] ) ) )
              {
              It.SetPixel( this->m_RotationIndices[i][7],
                static_cast<PixelType>( k7 ) );
              break;
              }
            }
          int k6 = currentLabel;
          if ( this->IsSpecialCaseOfC4Configuration2D( currentLabel,
                 It.GetIndex(), It.GetIndex( this->m_RotationIndices[i][6] ),
                 It.GetIndex( this->m_RotationIndices[i][7] ) ) )
            {
            for ( k6 = currentLabel; k6 >= 0; k6-- )
              {
              if ( this->IsChangeSafe2D( k6,
                     It.GetIndex( this->m_RotationIndices[i][6] ) ) )
                {
                It.SetPixel( this->m_RotationIndices[i][6],
                  static_cast<PixelType>( k6 ) );
                break;
                }
              }
            }
          currentLabel = vnl_math_min( static_cast<PixelType>( k4 ),
                                       static_cast<PixelType>( k7 ) );
          currentLabel = vnl_math_min( currentLabel,
                                       static_cast<PixelType>( k6 ) );

          break;
          }
        if ( !this->m_FullInvariance )
          {
          break;
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
            ( It.GetPixel( this->m_ReflectionIndices[i][j] ) == currentLabel );
          }

        if ( !this->m_FullInvariance
             && this->IsCriticalC1Configuration2D( neighborhoodPixels ) )
          {
          NumberOfC1Configurations2D++;
          for ( int k = currentLabel; k >= 0; k-- )
            {
            if ( this->IsChangeSafe2D( k,
                   It.GetIndex( this->m_ReflectionIndices[i][4] ) ) )
              {
              It.SetPixel( this->m_ReflectionIndices[i][4],
                static_cast<PixelType>( k ) );
              currentLabel = k;
              break;
              }
            }
          break;
          }
        else if ( !this->m_FullInvariance
                  && this->IsCriticalC2Configuration2D( neighborhoodPixels ) )
          {
          NumberOfC2Configurations2D++;
          for ( int k = currentLabel; k >= 0; k-- )
            {
            if ( this->IsChangeSafe2D( k,
                   It.GetIndex( this->m_ReflectionIndices[i][4] ) ) )
              {
              It.SetPixel( this->m_ReflectionIndices[i][4],
                static_cast<PixelType>( k ) );
              currentLabel = k;
              break;
              }
            }
          break;
          }
        else if ( this->IsCriticalC3Configuration2D( neighborhoodPixels ) )
          {
          NumberOfC3Configurations2D++;
          int k4;
          for ( k4 = currentLabel; k4 >= 0; k4-- )
            {
            if ( this->IsChangeSafe2D( k4,
                   It.GetIndex( this->m_ReflectionIndices[i][4] ) ) )
              {
              It.SetPixel( this->m_ReflectionIndices[i][4],
                static_cast<PixelType>( k4 ) );
              break;
              }
            }
          int k7;
          for ( k7 = currentLabel; k7 >= 0; k7-- )
            {
            if ( this->IsChangeSafe2D( k7,
                   It.GetIndex( this->m_ReflectionIndices[i][7] ) ) )
              {
              It.SetPixel( this->m_ReflectionIndices[i][7],
                static_cast<PixelType>( k7 ) );
              break;
              }
            }
          currentLabel = vnl_math_min( static_cast<PixelType>( k4 ),
                                       static_cast<PixelType>( k7 ) );
          break;
          }
        else if ( this->IsCriticalC4Configuration2D( neighborhoodPixels ) )
          {
          NumberOfC4Configurations2D++;

          int k4;
          for ( k4 = currentLabel; k4 >= 0; k4-- )
            {
            if ( this->IsChangeSafe2D( k4,
                   It.GetIndex( this->m_ReflectionIndices[i][4] ) ) )
              {
              It.SetPixel( this->m_ReflectionIndices[i][4],
                static_cast<PixelType>( k4 ) );
              break;
              }
            }
          int k7;
          for ( k7 = currentLabel; k7 >= 0; k7-- )
            {
            if ( this->IsChangeSafe2D( k7,
                   It.GetIndex( this->m_ReflectionIndices[i][7] ) ) )
              {
              It.SetPixel( this->m_ReflectionIndices[i][7],
                static_cast<PixelType>( k7 ) );
              break;
              }
            }
          int k6 = currentLabel;
          if ( this->IsSpecialCaseOfC4Configuration2D( currentLabel,
                 It.GetIndex(), It.GetIndex( this->m_ReflectionIndices[i][6] ),
                 It.GetIndex( this->m_ReflectionIndices[i][7] ) ) )
            {
            for ( k6 = currentLabel; k6 >= 0; k6-- )
              {
              if ( this->IsChangeSafe2D( k6,
                     It.GetIndex( this->m_ReflectionIndices[i][6] ) ) )
                {
                It.SetPixel( this->m_ReflectionIndices[i][6],
                  static_cast<PixelType>( k6 ) );
                break;
                }
              }
            }
          currentLabel = vnl_math_min( static_cast<PixelType>( k4 ),
                                       static_cast<PixelType>( k7 ) );
          currentLabel = vnl_math_min( currentLabel,
                                       static_cast<PixelType>( k6 ) );

          break;
          }
        if ( !this->m_FullInvariance )
          {
          break;
          }
        }
      currentLabel++;
      }
    }

  itkDebugMacro( << NumberOfC1Configurations2D
                 << " C1 configurations were located and repaired." );
  itkDebugMacro( << NumberOfC2Configurations2D
                 << " C2 configurations were located and repaired." );
  itkDebugMacro( << NumberOfC3Configurations2D
                 << " C3 configurations were located and repaired." );
  itkDebugMacro( << NumberOfC4Configurations2D
                 << " C4 configurations were located and repaired." );
}

/*
 * 2-D
 */
template<class TImage>
bool
WellComposedImageFilter<TImage>
::IsChangeSafe2D( PixelType label, IndexType idx )
{
  Array<char> neighborhoodPixels( 9 );

  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
        this->GetOutput()->GetLargestPossibleRegion() );
  It.SetLocation( idx );

  for ( unsigned int i = 0; i < 4; i++ )
    {
    for ( unsigned int j = 0; j < 9; j++ )
      {
      neighborhoodPixels[j] =
        ( It.GetPixel( this->m_RotationIndices[i][j] ) == label );
      if ( this->m_RotationIndices[i][j] == 4 )
        {
        neighborhoodPixels[j] = !neighborhoodPixels[j];
        }
      }
    if ( this->IsCriticalC1Configuration2D( neighborhoodPixels ) ||
         this->IsCriticalC2Configuration2D( neighborhoodPixels ) ||
         this->IsCriticalC3Configuration2D( neighborhoodPixels ) ||
         this->IsCriticalC4Configuration2D( neighborhoodPixels ) )
      {
      return false;
      }
    }

  for ( unsigned int i = 0; i < 2; i++ )
    {
    for ( unsigned int j = 0; j < 9; j++ )
      {
      neighborhoodPixels[j] =
        ( It.GetPixel( this->m_ReflectionIndices[i][j] ) == label );
      if ( this->m_ReflectionIndices[i][j] == 4 )
        {
        neighborhoodPixels[j] = !neighborhoodPixels[j];
        }
      }
    if ( this->IsCriticalC3Configuration2D( neighborhoodPixels ) ||
         this->IsCriticalC4Configuration2D( neighborhoodPixels ) )
      {
      return false;
      }
    }
  return true;
}

/*
 * 2-D
 */
template<class TImage>
bool
WellComposedImageFilter<TImage>
::IsCriticalC1Configuration2D( Array<char> neighborhood )
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
WellComposedImageFilter<TImage>
::IsCriticalC2Configuration2D( Array<char> neighborhood )
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
WellComposedImageFilter<TImage>
::IsCriticalC3Configuration2D( Array<char> neighborhood )
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
WellComposedImageFilter<TImage>
::IsCriticalC4Configuration2D( Array<char> neighborhood )
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
WellComposedImageFilter<TImage>
::IsSpecialCaseOfC4Configuration2D( PixelType label,
                                    IndexType idx,
                                    IndexType idx6,
                                    IndexType idx7 )
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
           this->GetOutput()->GetPixel( idxa ) != label &&
           this->GetOutput()->GetPixel( idxb ) == label );
}

/*
 * 2-D
 */
template<class TImage>
void
WellComposedImageFilter<TImage>
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
void
WellComposedImageFilter<TImage>
::MakeImageWellComposed3D()
{
  if ( this->GetDebug() )
    {
    this->CountCriticalConfigurations3D();
    }

  /**
   * Find critical configurations for all the labels
   */
  this->m_CriticalConfigurationIndices.resize( this->m_TotalNumberOfLabels );
  for ( int i = this->m_TotalNumberOfLabels-1; i >= 0; i-- )
    {
    this->LocateCriticalConfigurations3D( static_cast<PixelType>( i ) );
    }

  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
    this->GetOutput()->GetLargestPossibleRegion() );

  /**
   * Remove critical configurations
   */
  bool foundCriticalConfiguration = true;
  while ( foundCriticalConfiguration )
    {
    foundCriticalConfiguration = false;
    int label;
    for ( label = this->m_TotalNumberOfLabels-1; label >= 0; label-- )
      {
      if ( !this->m_CriticalConfigurationIndices[label].empty() )
        {
        break;
        }
      }
    if ( label < 0 )
      {
      break;
      }
    itkDebugMacro( << "Removing the " << this->m_CriticalConfigurationIndices[label].size()
                   << " critical configurations of label " << label );
    while ( !this->m_CriticalConfigurationIndices[label].empty() )
      {
      foundCriticalConfiguration = true;

      It.SetLocation( this->m_CriticalConfigurationIndices[label].front() );
      this->m_CriticalConfigurationIndices[label].pop_front();

      Array<char> neighborhoodPixels( 8 );
      bool removedCriticalConfiguration = false;
      /**
       * Deal with the C1 configurations
       */
      for ( unsigned int i = 0; i < 3; i++ )
        {
        for ( unsigned int j = 0; j < 4; j++ )
          {
          neighborhoodPixels[j] = ( It.GetPixel( this->m_C1Indices[i][j] ) == label );
          }
        if ( this->IsCriticalC1Configuration3D( neighborhoodPixels ) )
          {
          this->RemoveCriticalC1Configuration3D( i, label, It.GetIndex() );
          removedCriticalConfiguration = true;
          }
        }

      /**
       * Deal with the C2 configurations
       */
      if ( !removedCriticalConfiguration )
        {
        for ( unsigned int j = 0; j < 8; j++ )
          {
          neighborhoodPixels[j]
            = ( It.GetPixel( this->m_C2Indices[0][j] ) == label );
          }
        if ( this->IsCriticalC2Configuration3D( neighborhoodPixels ) )
          {
          this->RemoveCriticalC2Configuration3D( 0, label, It.GetIndex() );
          }
        }
      }
    }
  if ( this->GetDebug() )
    {
    this->CountCriticalConfigurations3D();
    }
}

/*
 * 3-D
 */
template<class TImage>
void
WellComposedImageFilter<TImage>
::LocateCriticalConfigurations3D( PixelType label )
{
  this->m_CriticalConfigurationIndices[label].clear();

  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
    this->GetOutput()->GetLargestPossibleRegion() );

  Array<char> neighborhoodPixels( 8 );

  int countC1 = 0;
  int countC2 = 0;

  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( !It.InBounds() )
      {
      continue;
      }

    /**
     * Check for C1 critical configurations
     */
    bool foundC1 = false;
    for ( unsigned int i = 0; i < 3; i++ )
      {
      for ( unsigned int j = 0; j < 4; j++ )
        {
        neighborhoodPixels[j] = ( It.GetPixel( this->m_C1Indices[i][j] ) == label );
        }
      if ( this->IsCriticalC1Configuration3D( neighborhoodPixels ) )
        {
        foundC1 = true;
        countC1++;
        }
      }
    if ( foundC1 )
      {
      this->m_CriticalConfigurationIndices[label].push_back( It.GetIndex() );
      }

    /**
     * Check for C2 critical configurations
     */
    else
      {
      for ( unsigned int j = 0; j < 8; j++ )
        {
        neighborhoodPixels[j] = ( It.GetPixel( this->m_C2Indices[0][j] ) == label );
        }
      if ( this->IsCriticalC2Configuration3D( neighborhoodPixels ) )
        {
        this->m_CriticalConfigurationIndices[label].push_back( It.GetIndex() );
        countC2++;
        }
      }
    }
  itkDebugMacro( << "Label " << label << ": " << countC1 << " C1 and " << countC2
                 << " C2 critical configurations." );
}

template<class TImage>
void
WellComposedImageFilter<TImage>
::RemoveCriticalC1Configuration3D( int which, PixelType label, IndexType idx )
{
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
    this->GetOutput()->GetLargestPossibleRegion() );
  It.SetLocation( idx );

  IndexContainerType indices;
  PixelType newLabel = label;

  if ( It.GetPixel( this->m_C1Indices[which][0] ) == label )
    {
    if ( It.GetPixel( this->m_C1Indices[which][2] ) < label  )
      {
      indices.push_back( It.GetIndex( this->m_C1Indices[which][2] ) );
      }
    if ( It.GetPixel( this->m_C1Indices[which][3] ) < label  )
      {
      indices.push_back( It.GetIndex( this->m_C1Indices[which][3] ) );
      }

    if ( indices.empty() )
      {
      newLabel = vnl_math_max(
        It.GetPixel( this->m_C1Indices[which][2] ),
        It.GetPixel( this->m_C1Indices[which][3] ) );
      indices.push_back( It.GetIndex( this->m_C1Indices[which][0] ) );
      indices.push_back( It.GetIndex( this->m_C1Indices[which][1] ) );
      }
    }
  else
    {
    if ( It.GetPixel( this->m_C1Indices[which][0] ) < label  )
      {
      indices.push_back( It.GetIndex( this->m_C1Indices[which][0] ) );
      }
    if ( It.GetPixel( this->m_C1Indices[which][1] ) < label  )
      {
      indices.push_back( It.GetIndex( this->m_C1Indices[which][1] ) );
      }

    if ( indices.empty() )
      {
      newLabel = vnl_math_max(
        It.GetPixel( this->m_C1Indices[which][0] ),
        It.GetPixel( this->m_C1Indices[which][1] ) );
      indices.push_back( It.GetIndex( this->m_C1Indices[which][2] ) );
      indices.push_back( It.GetIndex( this->m_C1Indices[which][3] ) );
      }
    }

  IndexContainerType unsafeIndices;
  typename IndexContainerType::const_iterator it;
  for ( it = indices.begin(); it != indices.end(); ++it )
    {
    if ( this->IsChangeSafe3D( label, *it ) )
      {
      unsafeIndices.push_back( *it );
      }
    }
  typedef typename Statistics
     ::MersenneTwisterRandomVariateGenerator GeneratorType;
  typename GeneratorType::Pointer generator = GeneratorType::New();
  generator->SetSeed();

  if ( !unsafeIndices.empty() )
    {
    int n = generator->GetIntegerVariate( unsafeIndices.size()-1 );
    if ( newLabel == label )
      {
      this->InsertCriticalConfiguration3D(
        this->GetOutput()->GetPixel( unsafeIndices[n] ), unsafeIndices[n] );
      }
    else
      {
      this->InsertCriticalConfiguration3D( newLabel, unsafeIndices[n] );
      }
    this->GetOutput()->SetPixel( unsafeIndices[n], newLabel );
    }
  else
    {
    int n = generator->GetIntegerVariate( indices.size()-1 );
    this->InsertCriticalConfiguration3D( label, indices[n] );
    if ( newLabel == label )
      {
      this->InsertCriticalConfiguration3D(
        this->GetOutput()->GetPixel( indices[n] ), indices[n] );
      }
    else
      {
      this->InsertCriticalConfiguration3D( newLabel, indices[n] );
      }
    this->GetOutput()->SetPixel( indices[n], newLabel );
    }
}

template<class TImage>
void
WellComposedImageFilter<TImage>
::RemoveCriticalC2Configuration3D( int which, PixelType label, IndexType idx )
{
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
    this->GetOutput()->GetLargestPossibleRegion() );
  It.SetLocation( idx );

  PixelType newLabel = label;

  IndexContainerType indices;
  for ( unsigned int i = 0; i < 8; i++ )
    {
    if ( It.GetPixel( this->m_C2Indices[which][i] ) < label )
      {
      indices.push_back( It.GetIndex( this->m_C2Indices[which][i] ) );
      }
    }

  if ( indices.empty() )
    {
    newLabel = NumericTraits<PixelType>::Zero;
    for ( unsigned int i = 0; i < 8; i++ )
      {
      PixelType pix = It.GetPixel( this->m_C2Indices[which][i] );
      if ( pix == label )
        {
        indices.push_back( It.GetIndex( this->m_C2Indices[which][i] ) );
        }
      else if ( pix > newLabel )
        {
        newLabel = pix;
        }
      }
    }

  IndexContainerType unsafeIndices;
  typename IndexContainerType::const_iterator it;
  for ( it = indices.begin(); it != indices.end(); ++it )
    {
    if ( this->IsChangeSafe3D( label, *it ) )
      {
      unsafeIndices.push_back( *it );
      }
    }
  typedef typename Statistics
     ::MersenneTwisterRandomVariateGenerator GeneratorType;
  typename GeneratorType::Pointer generator = GeneratorType::New();
  generator->SetSeed();

  if ( !unsafeIndices.empty() )
    {
    int n = generator->GetIntegerVariate( unsafeIndices.size()-1 );
    if ( newLabel == label )
      {
      this->InsertCriticalConfiguration3D(
        this->GetOutput()->GetPixel( unsafeIndices[n] ), unsafeIndices[n] );
      }
    else
      {
      this->InsertCriticalConfiguration3D( newLabel, unsafeIndices[n] );
      }
    this->GetOutput()->SetPixel( unsafeIndices[n], newLabel );
    }
  else
    {
    int n = generator->GetIntegerVariate( indices.size()-1 );
    this->InsertCriticalConfiguration3D( label, indices[n] );
    if ( newLabel == label )
      {
      this->InsertCriticalConfiguration3D(
        this->GetOutput()->GetPixel( indices[n] ), indices[n] );
      }
    else
      {
      this->InsertCriticalConfiguration3D( newLabel, indices[n] );
      }
    this->GetOutput()->SetPixel( indices[n], newLabel );
    }
}

/*
 * 3-D
 */
template<class TImage>
bool
WellComposedImageFilter<TImage>
::IsChangeSafe3D( PixelType label, IndexType idx )
{
  Array<char> neighborhoodPixels( 8 );

  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
        this->GetOutput()->GetLargestPossibleRegion() );
  It.SetLocation( idx );

  // Check for C1 critical configurations
  for ( unsigned int i = 0; i < 12; i++ )
    {
    for ( unsigned int j = 0; j < 4; j++ )
      {
      neighborhoodPixels[j] = ( It.GetPixel( this->m_C1Indices[i][j] ) == label );
      if ( this->m_C1Indices[i][j] == 13 )
        {
        neighborhoodPixels[j] = !neighborhoodPixels[j];
        }
      }
    if ( this->IsCriticalC1Configuration3D( neighborhoodPixels ) )
      {
      return false;
      }
    }

  // Check for C2 critical configurations
  for ( unsigned int i = 0; i < 8; i++ )
    {
    for ( unsigned int j = 0; j < 8; j++ )
      {
      neighborhoodPixels[j] = ( It.GetPixel( this->m_C2Indices[i][j] ) == label );
      if ( this->m_C2Indices[i][j] == 13 )
        {
        neighborhoodPixels[j] = !neighborhoodPixels[j];
        }
      }
    if ( this->IsCriticalC2Configuration3D( neighborhoodPixels ) )
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
void
WellComposedImageFilter<TImage>
::InsertCriticalConfiguration3D( PixelType label,
                                 IndexType idx )
{
  Array<char> neighborhoodPixels( 8 );

  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
        this->GetOutput()->GetLargestPossibleRegion() );
  It.SetLocation( idx );

  // C1 configurations
  for ( unsigned int i = 0; i < 12; i++ )
    {
    for ( unsigned int j = 0; j < 4; j++ )
      {
      neighborhoodPixels[j] = ( It.GetPixel( this->m_C1Indices[i][j] ) == label );
      if ( this->m_C1Indices[i][j] == 13 )
        {
        neighborhoodPixels[j] = !neighborhoodPixels[j];
        }
      }
    if ( this->IsCriticalC1Configuration3D( neighborhoodPixels ) )
      {
      this->m_CriticalConfigurationIndices[label].push_back(
        It.GetIndex( this->m_C1Indices[i][0] ) );
      }
    }

  // C2 configurations
  for ( unsigned int i = 0; i < 8; i++ )
    {
    for ( unsigned int j = 0; j < 8; j++ )
      {
      neighborhoodPixels[j] = ( It.GetPixel( this->m_C2Indices[i][j] ) == label );
      if ( this->m_C2Indices[i][j] == 13 )
        {
        neighborhoodPixels[j] = !neighborhoodPixels[j];
        }
      }
    if ( this->IsCriticalC2Configuration3D( neighborhoodPixels ) )
      {
      this->m_CriticalConfigurationIndices[label].push_back(
        It.GetIndex( this->m_C2Indices[i][0] ) );
      }
    }
}

/*
 * 3-D
 */
template<class TImage>
bool
WellComposedImageFilter<TImage>
::IsCriticalC1Configuration3D( Array<char> neighborhood )
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
WellComposedImageFilter<TImage>
::IsCriticalC2Configuration3D( Array<char> neighborhood )
{
  // Check if Type 1 or Type 2
  for ( unsigned int i = 0; i < 4; i++ )
    {
    bool isC2 = false;
    if ( neighborhood[2*i] == neighborhood[2*i+1] )
      {
      isC2 = true;
      for ( unsigned int j = 0; j < 8; j++ )
        {
        if ( neighborhood[j] == neighborhood[2*i] &&
               j != 2*i && j != 2*i+1 )
          {
          isC2 = false;
          }
        }
      }
    if ( isC2 )
      {
      if ( neighborhood[2*i] )
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
WellComposedImageFilter<TImage>
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

  this->m_C1Indices[0][0] = 13;
  this->m_C1Indices[0][1] = 23;
  this->m_C1Indices[0][2] = 14;
  this->m_C1Indices[0][3] = 22;

  this->m_C1Indices[1][0] = 13;
  this->m_C1Indices[1][1] = 17;
  this->m_C1Indices[1][2] = 14;
  this->m_C1Indices[1][3] = 16;

  this->m_C1Indices[2][0] = 13;
  this->m_C1Indices[2][1] = 25;
  this->m_C1Indices[2][2] = 16;
  this->m_C1Indices[2][3] = 22;

  this->m_C1Indices[3][0] = 4;
  this->m_C1Indices[3][1] = 14;
  this->m_C1Indices[3][2] = 5;
  this->m_C1Indices[3][3] = 13;

  this->m_C1Indices[4][0] = 12;
  this->m_C1Indices[4][1] = 22;
  this->m_C1Indices[4][2] = 13;
  this->m_C1Indices[4][3] = 21;

  this->m_C1Indices[5][0] = 1;
  this->m_C1Indices[5][1] = 13;
  this->m_C1Indices[5][2] = 4;
  this->m_C1Indices[5][3] = 10;

  this->m_C1Indices[6][0] = 9;
  this->m_C1Indices[6][1] = 13;
  this->m_C1Indices[6][2] = 10;
  this->m_C1Indices[6][3] = 12;

  this->m_C1Indices[7][0] = 3;
  this->m_C1Indices[7][1] = 13;
  this->m_C1Indices[7][2] = 4;
  this->m_C1Indices[7][3] = 12;

  this->m_C1Indices[8][0] = 10;
  this->m_C1Indices[8][1] = 22;
  this->m_C1Indices[8][2] = 13;
  this->m_C1Indices[8][3] = 19;

  this->m_C1Indices[9][0] = 12;
  this->m_C1Indices[9][1] = 16;
  this->m_C1Indices[9][2] = 13;
  this->m_C1Indices[9][3] = 15;

  this->m_C1Indices[10][0] = 4;
  this->m_C1Indices[10][1] = 16;
  this->m_C1Indices[10][2] = 7;
  this->m_C1Indices[10][3] = 13;

  this->m_C1Indices[11][0] = 10;
  this->m_C1Indices[11][1] = 14;
  this->m_C1Indices[11][2] = 11;
  this->m_C1Indices[11][3] = 13;

  this->m_C2Indices[0][0] = 13;
  this->m_C2Indices[0][1] = 26;
  this->m_C2Indices[0][2] = 14;
  this->m_C2Indices[0][3] = 25;
  this->m_C2Indices[0][4] = 16;
  this->m_C2Indices[0][5] = 23;
  this->m_C2Indices[0][6] = 17;
  this->m_C2Indices[0][7] = 22;

  this->m_C2Indices[4][0] = 4;
  this->m_C2Indices[4][1] = 17;
  this->m_C2Indices[4][2] = 5;
  this->m_C2Indices[4][3] = 16;
  this->m_C2Indices[4][4] = 7;
  this->m_C2Indices[4][5] = 14;
  this->m_C2Indices[4][6] = 8;
  this->m_C2Indices[4][7] = 13;


  for ( unsigned int i = 1; i < 4; i++ )
    {
    int subtrahend;
    if ( i == 2 )
      {
      subtrahend = 2;
      }
    else
      {
      subtrahend = 1;
      }
    for ( unsigned int j = 0; j < 8; j++ )
      {
      this->m_C2Indices[i  ][j] = this->m_C2Indices[i-1][j] - subtrahend;
      this->m_C2Indices[i+4][j] = this->m_C2Indices[i+3][j] - subtrahend;
      }
    }
}

/*
 * 3-D: Debug
 */
template<class TImage>
bool
WellComposedImageFilter<TImage>
::IsCriticalC1Configuration3D( Array<PixelType> neighborhood )
{
  return ( ( ( neighborhood[1] == neighborhood[0] ) &&
             ( neighborhood[2] != neighborhood[0] ) &&
             ( neighborhood[3] != neighborhood[0] ) ) ||
           ( ( neighborhood[3] == neighborhood[2] ) &&
             ( neighborhood[0] != neighborhood[2] ) &&
             ( neighborhood[1] != neighborhood[2] ) ) );
}

/*
 * 3-D: Debug
 */
template<class TImage>
unsigned int
WellComposedImageFilter<TImage>
::IsCriticalC2Configuration3D( Array<PixelType> neighborhood )
{
  // Check for Subtype 1 of C.C. Type 2
  if ( ( ( neighborhood[0] == neighborhood[1] ) && ( neighborhood[0] != neighborhood[2] ) &&
         ( neighborhood[0] != neighborhood[5] ) && ( neighborhood[0] != neighborhood[7] ) &&
         ( neighborhood[0] != neighborhood[4] ) && ( neighborhood[0] != neighborhood[6] ) &&
         ( neighborhood[0] != neighborhood[3] ) ) ||
       ( ( neighborhood[5] == neighborhood[4] ) && ( neighborhood[5] != neighborhood[0] ) &&
         ( neighborhood[5] != neighborhood[2] ) && ( neighborhood[5] != neighborhood[7] ) &&
         ( neighborhood[5] != neighborhood[6] ) && ( neighborhood[5] != neighborhood[1] ) &&
         ( neighborhood[5] != neighborhood[3] ) ) ||
       ( ( neighborhood[2] == neighborhood[3] ) && ( neighborhood[2] != neighborhood[0] ) &&
         ( neighborhood[2] != neighborhood[5] ) && ( neighborhood[2] != neighborhood[7] ) &&
         ( neighborhood[2] != neighborhood[4] ) && ( neighborhood[2] != neighborhood[6] ) &&
         ( neighborhood[2] != neighborhood[1] ) ) ||
       ( ( neighborhood[7] == neighborhood[6] ) && ( neighborhood[7] != neighborhood[0] ) &&
         ( neighborhood[7] != neighborhood[2] ) && ( neighborhood[7] != neighborhood[5] ) &&
         ( neighborhood[7] != neighborhood[4] ) && ( neighborhood[7] != neighborhood[1] ) &&
         ( neighborhood[7] != neighborhood[3] ) ) )
    {
    return 1;
    }
  // Check for Subtype 2 of C.C. Type 2
  if ( ( ( neighborhood[2] == neighborhood[5] ) && ( neighborhood[2] == neighborhood[7] ) &&
         ( neighborhood[2] == neighborhood[4] ) && ( neighborhood[2] == neighborhood[6] ) &&
         ( neighborhood[2] == neighborhood[3] ) && ( neighborhood[0] != neighborhood[2] ) &&
         ( neighborhood[5] != neighborhood[2] ) ) ||
       ( ( neighborhood[0] == neighborhood[2] ) && ( neighborhood[0] == neighborhood[7] ) &&
         ( neighborhood[0] == neighborhood[6] ) && ( neighborhood[0] == neighborhood[1] ) &&
         ( neighborhood[0] == neighborhood[3] ) && ( neighborhood[0] != neighborhood[5] ) &&
         ( neighborhood[0] != neighborhood[4] ) ) ||
       ( ( neighborhood[0] == neighborhood[5] ) && ( neighborhood[0] == neighborhood[7] ) &&
         ( neighborhood[0] == neighborhood[4] ) && ( neighborhood[0] == neighborhood[6] ) &&
         ( neighborhood[0] == neighborhood[1] ) && ( neighborhood[0] != neighborhood[2] ) &&
         ( neighborhood[0] != neighborhood[3] ) ) ||
       ( ( neighborhood[0] == neighborhood[2] ) && ( neighborhood[0] == neighborhood[5] ) &&
         ( neighborhood[0] == neighborhood[4] ) && ( neighborhood[0] == neighborhood[1] ) &&
         ( neighborhood[0] == neighborhood[3] ) && ( neighborhood[0] != neighborhood[7] ) &&
         ( neighborhood[0] != neighborhood[6] ) ) )
    {
    return 2;
    }
  return 0;
}

/*
 * 3-D: Debug
 */
template<class TImage>
void
WellComposedImageFilter<TImage>
::CountCriticalConfigurations3D()
{
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
    this->GetOutput()->GetLargestPossibleRegion() );

  Array<PixelType> neighborhoodPixels( 8 );

  int countC1 = 0;
  int countC2 = 0;

  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( !It.InBounds() )
      {
      continue;
      }

    /**
     * Check for C1 critical configurations
     */
    bool foundC1 = false;
    for ( unsigned int i = 0; i < 3; i++ )
      {
      for ( unsigned int j = 0; j < 4; j++ )
        {
        neighborhoodPixels[j] = It.GetPixel( this->m_C1Indices[i][j] );
        }
      if ( this->IsCriticalC1Configuration3D( neighborhoodPixels ) )
        {
        foundC1 = true;
        countC1++;
        }
      }
    /**
     * Check for C2 critical configurations
     */
    if ( !foundC1 )
      {
      for ( unsigned int j = 0; j < 8; j++ )
        {
        neighborhoodPixels[j] = It.GetPixel( this->m_C2Indices[0][j] );
        }
      if ( this->IsCriticalC2Configuration3D( neighborhoodPixels ) )
        {
        countC2++;
        }
      }
    }
  itkDebugMacro( << "There are a total of " << countC1 << " C1 and " << countC2
                 << " C2 critical configurations." );
  std::cout << "There are a total of " << countC1 << " C1 and " << countC2
                 << " C2 critical configurations." << std::endl;
}

template <class TImage>
void
WellComposedImageFilter<TImage>
::PrintSelf(
  std::ostream& os,
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Total number of labels: "
     << this->m_TotalNumberOfLabels << std::endl;
}

} // end namespace itk

#endif
