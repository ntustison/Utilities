#ifndef __itkWellComposedImageFilter_hxx
#define __itkWellComposedImageFilter_hxx

#include "itkWellComposedImageFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkImageLinearIteratorWithIndex.h"

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
  this->ConvertBoundaryPixels( NumericTraits<PixelType>::Zero );
  
  if ( ImageDimension == 2 )
    {
    this->MakeImageWellComposed2D();
    }
  else // if ( ImageDimension == 3 )
    {
    this->MakeImageWellComposed3D();
    }
}

template<class TImage>
void
WellComposedImageFilter<TImage>
::ConvertBoundaryPixels( PixelType value )
{
  ImageLinearIteratorWithIndex<ImageType> It( this->GetOutput(),
    this->GetOutput()->GetRequestedRegion() );
  
  unsigned long N = 0;

  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    It.SetDirection( d );
    It.GoToBegin();
    while ( !It.IsAtEnd() )
    {
      It.GoToBeginOfLine();
      if ( It.Get() != value )
        { 
        It.Set( value );
        N++;
        }
      It.GoToEndOfLine();
      --It;
      if ( It.Get() != value )
        { 
        It.Set( value );
        N++;
        }
      It.NextLine();
      }
    }
  itkDebugMacro( << N << " boundary pixels were converted." );
}

/*
 * 2-D
 */
template<class TImage>
void
WellComposedImageFilter<TImage>
::MakeImageWellComposed2D()
{
 typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
        this->GetOutput()->GetRequestedRegion() );

  Array<char> neighborhoodPixels( 9 );

  for ( int label = this->m_TotalNumberOfLabels-1; label > 0 ; label-- )
    {
    unsigned long NumberOfC1Configurations2D = 0;
    unsigned long NumberOfC2Configurations2D = 0;
    unsigned long NumberOfC3Configurations2D = 0;
    unsigned long NumberOfC4Configurations2D = 0;
    
    for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if ( !It.InBounds())
        {
        continue;
        } 
      
      // Check for critical configurations: 4 90-degree rotations
      
      for ( unsigned int i = 0; i < 4; i++ )
        {
        for ( unsigned int j = 0; j < 9; j++ )
          {
          neighborhoodPixels[j] = 
            ( It.GetPixel( this->m_RotationIndices[i][j] ) >= label );
          }
        if ( this->IsCriticalC1Configuration2D( neighborhoodPixels ) )
          {
          It.SetPixel( this->m_RotationIndices[i][4], static_cast<PixelType>( label ) );
          NumberOfC1Configurations2D++;
          break;
          }
        else if ( this->IsCriticalC2Configuration2D( neighborhoodPixels ) )
          {
          It.SetPixel( this->m_RotationIndices[i][4], static_cast<PixelType>( label ) );
          NumberOfC2Configurations2D++;
          break;
          }
        else if ( this->IsCriticalC3Configuration2D( neighborhoodPixels ) )
          {
          It.SetPixel( this->m_RotationIndices[i][4], static_cast<PixelType>( label ) );
          It.SetPixel( this->m_RotationIndices[i][7], static_cast<PixelType>( label ) );
          NumberOfC3Configurations2D++;
          break;
          }
        else if ( this->IsCriticalC4Configuration2D( neighborhoodPixels ) )
          {
          It.SetPixel( this->m_RotationIndices[i][4], static_cast<PixelType>( label ) );
          It.SetPixel( this->m_RotationIndices[i][7], static_cast<PixelType>( label ) );
          if ( this->IsSpecialCaseOfC4Configuration2D( label,
                 It.GetIndex(), It.GetIndex( this->m_RotationIndices[i][6] ), 
                 It.GetIndex( this->m_RotationIndices[i][7] ) ) )
            {
            It.SetPixel( this->m_RotationIndices[i][6], static_cast<PixelType>( label ) );
            }
          NumberOfC4Configurations2D++;
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
            ( It.GetPixel( this->m_ReflectionIndices[i][j] ) >= label );
          }
        if ( !this->m_FullInvariance && this->IsCriticalC1Configuration2D( neighborhoodPixels ) )
          {
          It.SetPixel( this->m_ReflectionIndices[i][4], static_cast<PixelType>( label ) );
          NumberOfC1Configurations2D++;
          break;
          }
        else if ( !this->m_FullInvariance && this->IsCriticalC2Configuration2D( neighborhoodPixels ) )
          {
          It.SetPixel( this->m_ReflectionIndices[i][4], static_cast<PixelType>( label )  );
          NumberOfC2Configurations2D++;
          break;
          }
        else if ( this->IsCriticalC3Configuration2D( neighborhoodPixels ) )
          {
          It.SetPixel( this->m_ReflectionIndices[i][4], static_cast<PixelType>( label )  );
          It.SetPixel( this->m_ReflectionIndices[i][7], static_cast<PixelType>( label )  );
          NumberOfC3Configurations2D++;
          break;
          }
        else if ( this->IsCriticalC4Configuration2D( neighborhoodPixels ) )
          {
          It.SetPixel( this->m_ReflectionIndices[i][4], static_cast<PixelType>( label )  );
          It.SetPixel( this->m_ReflectionIndices[i][7], static_cast<PixelType>( label )  );
          if ( this->IsSpecialCaseOfC4Configuration2D( label,
                 It.GetIndex(), It.GetIndex( this->m_ReflectionIndices[i][6] ), 
                 It.GetIndex( this->m_ReflectionIndices[i][7] ) ) )
            {
            It.SetPixel( this->m_ReflectionIndices[i][6], 
                         static_cast<PixelType>( label )  );
            }
          NumberOfC4Configurations2D++;
          break;
          }
        if ( !this->m_FullInvariance )
          {
          break;
          }
        }  
      }

    itkDebugMacro( << NumberOfC1Configurations2D 
                   << " C1 configurations were located and repaired ("
                   << "label = " << label << ")." ); 
    itkDebugMacro( << NumberOfC2Configurations2D 
                   << " C2 configurations were located and repaired ("
                   << "label = " << label << ")." ); 
    itkDebugMacro( << NumberOfC3Configurations2D 
                   << " C3 configurations were located and repaired ("
                   << "label = " << label << ")." ); 
    itkDebugMacro( << NumberOfC4Configurations2D 
                   << " C4 configurations were located and repaired ("
                   << "label = " << label << ")." ); 
    }
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
           this->GetOutput()->GetPixel( idxa ) <  label &&
           this->GetOutput()->GetPixel( idxb ) >= label );
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
/**
 * 3-D
 */
template<class TImage>
void
WellComposedImageFilter<TImage>
::MakeImageWellComposed3D()
{
  /**
   * Find critical configurations for all the labels
   */
  this->m_CriticalConfigurationIndices.resize( this->m_TotalNumberOfLabels );

  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
    this->GetOutput()->GetLargestPossibleRegion() );

  /**
   * Remove critical configurations
   */
  for ( int label = this->m_TotalNumberOfLabels-1; label >= 0; label-- ) 
    {
    this->LocateCriticalConfigurations3D( static_cast<PixelType>( label ) );

    itkDebugMacro( << "Removing the " 
                   << this->m_CriticalConfigurationIndices.size()
                   << " critical configurations of label " << label );
    
    this->m_NewCriticalConfigurationIndices.clear();
    
    while ( !this->m_CriticalConfigurationIndices.empty() ) 
      {      
      IndexType idx = this->m_CriticalConfigurationIndices.front();
      this->m_CriticalConfigurationIndices.pop_front();
    
      It.SetLocation( idx );
      
      Array<char> neighborhoodPixels( 8 );
      /**
       * Deal with the C1 configurations
       */
      for ( unsigned int i = 0; i < 3; i++ )
        {
        for ( unsigned int j = 0; j < 4; j++ )
          {
          neighborhoodPixels[j] = ( It.GetPixel( this->m_C1Indices[i][j] ) >= label );
          }
        if ( this->IsCriticalC1Configuration3D( neighborhoodPixels ) )
          {
          this->RemoveCriticalC1Configuration3D( i, label, It.GetIndex() );
          }
        }
        
      /**
       * Deal with the C2 configurations
       */
      for ( unsigned int j = 0; j < 8; j++ )
        {
        neighborhoodPixels[j] = ( It.GetPixel( this->m_C2Indices[0][j] ) >= label );
        }
      if ( this->IsCriticalC2Configuration3D( neighborhoodPixels ) )
        {
        this->RemoveCriticalC2Configuration3D( 0, label, It.GetIndex() );
        }
      if ( this->m_CriticalConfigurationIndices.empty() )
        {
        itkDebugMacro( << "Label " << label << " created " 
                       << this->m_NewCriticalConfigurationIndices.size()
                       << " new configurations." );
        this->m_CriticalConfigurationIndices = this->m_NewCriticalConfigurationIndices;
        this->m_NewCriticalConfigurationIndices.clear();
        }        
      }
    }
  if ( this->GetDebug() )
    {
    for ( int label = this->m_TotalNumberOfLabels-1; label > 0; label-- )
      {
      this->LocateCriticalConfigurations3D( static_cast<PixelType>( label ) ); 
      }
    }  
}    

/*
 * 3-D
 */
template<class TImage>
bool
WellComposedImageFilter<TImage>
::RemoveCriticalC1Configuration3D( int which, PixelType label, IndexType idx )
{
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
    this->GetOutput()->GetLargestPossibleRegion() );
  It.SetLocation( idx );

  IndexContainerType indices;
  for ( unsigned int i = 0; i < 4; i++ )
    {
    if ( It.GetPixel( this->m_C1Indices[which][i] ) < label )
      {
      indices.push_back( It.GetIndex( this->m_C1Indices[which][i] ) );
      }
    }

  for ( unsigned int i = 0; i < indices.size(); i++ )
    {
    if ( this->IsChangeSafe3D( label, indices[i] ) )
      {
      this->GetOutput()->SetPixel( indices[i], label );
      return true;
      }
    }  

  typedef typename Statistics
     ::MersenneTwisterRandomVariateGenerator GeneratorType; 
  typename GeneratorType::Pointer generator = GeneratorType::New();
  int n = generator->GetIntegerVariate( 1 );
    
  this->InsertCriticalConfiguration3D( label, indices[n] );
  return false;
}

/*
 * 3-D
 */
template<class TImage>
bool
WellComposedImageFilter<TImage>
::RemoveCriticalC2Configuration3D( int which, PixelType label, IndexType idx )
{
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
    this->GetOutput()->GetLargestPossibleRegion() );
  It.SetLocation( idx );

  IndexContainerType indices;
  for ( unsigned int i = 0; i < 8; i++ )
    {
    if ( It.GetPixel( this->m_C2Indices[which][i] ) < label )
      {
      indices.push_back( It.GetIndex( this->m_C2Indices[which][i] ) );
      }
    }

  for ( unsigned int i = 0; i < indices.size(); i++ )
    {
    if ( this->IsChangeSafe3D( label, indices[i] ) )
      {
      this->GetOutput()->SetPixel( indices[i], label );
      return true;
      }
    }  

  typedef typename Statistics
     ::MersenneTwisterRandomVariateGenerator GeneratorType; 
  typename GeneratorType::Pointer generator = GeneratorType::New();
  int n = generator->GetIntegerVariate( indices.size()-1 );
    
  this->InsertCriticalConfiguration3D( label, indices[n] );
  return false;
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
      neighborhoodPixels[j] = ( It.GetPixel( this->m_C1Indices[i][j] ) >= label );
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
      neighborhoodPixels[j] = ( It.GetPixel( this->m_C2Indices[i][j] ) >= label );
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

  It.SetCenterPixel( label );
      
  // C1 configurations
  for ( unsigned int i = 0; i < 12; i++ )
    {
    for ( unsigned int j = 0; j < 4; j++ )
      {
      neighborhoodPixels[j] = ( It.GetPixel( this->m_C1Indices[i][j] ) >= label );
      }
    if ( this->IsCriticalC1Configuration3D( neighborhoodPixels ) )
      {
      this->m_NewCriticalConfigurationIndices.push_back( 
        It.GetIndex( this->m_C1Indices[i][1] ) ); 
      }
    }  

  // C2 configurations
  for ( unsigned int i = 0; i < 8; i++ )
    {
    for ( unsigned int j = 0; j < 8; j++ )
      {
      neighborhoodPixels[j] = ( It.GetPixel( this->m_C2Indices[i][j] ) >= label );
      }
    if ( this->IsCriticalC2Configuration3D( neighborhoodPixels ) )
      {
      this->m_NewCriticalConfigurationIndices.push_back( 
        It.GetIndex( this->m_C2Indices[i][1] ) );
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
::LocateCriticalConfigurations3D( PixelType label )
{
  this->m_NewCriticalConfigurationIndices.clear();

  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
    this->GetOutput()->GetLargestPossibleRegion() );

  unsigned long countC1 = 0;
  unsigned long countC2 = 0;

  Array<char> neighborhoodPixels( 8 );

  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( !It.InBounds() )
      {
      continue;
      }
      
    bool indexAdded = false;  
    /**
     * Check for C1 critical configurations
     */
    for ( unsigned int i = 0; i < 3; i++ )
      {
      for ( unsigned int j = 0; j < 4; j++ )
        {
        neighborhoodPixels[j] = ( It.GetPixel( this->m_C1Indices[i][j] ) >= label );
        }
      if ( this->IsCriticalC1Configuration3D( neighborhoodPixels ) )
        {
        countC1++;
        indexAdded = true;
        }
      }

    if ( indexAdded )
      {
      this->m_NewCriticalConfigurationIndices.push_back( It.GetIndex() );
      }
    else
      {
      /** 
       * Check for C2 critical configurations
       */
      for ( unsigned int j = 0; j < 8; j++ )
        {
        neighborhoodPixels[j] = ( It.GetPixel( this->m_C2Indices[0][j] ) >= label );
        }
      if ( this->IsCriticalC2Configuration3D( neighborhoodPixels ) )
        {
        this->m_NewCriticalConfigurationIndices.push_back( It.GetIndex() );
        countC2++; 
        }
      }
    }  
  itkDebugMacro( << "Label " << label << " contains " << countC1 << " C1 and "
                 << countC2 << " C2 critical configurations." ); 
    
  this->m_CriticalConfigurationIndices = this->m_NewCriticalConfigurationIndices;
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
    if ( i == 2 )
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