#ifndef __itkBinaryWellComposedImageFilter_txx
#define __itkBinaryWellComposedImageFilter_txx

#include "itkBinaryWellComposedImageFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkConstantPadImageFilter.h"
#include "itkExtractImageFilter.h"

namespace itk
{

template<class TImage>
BinaryWellComposedImageFilter<TImage>
::BinaryWellComposedImageFilter()
{
  this->m_BackgroundValue = NumericTraits<PixelType>::Zero;
  this->m_ForegroundValue = NumericTraits<PixelType>::One;

  this->m_FullInvariance = true;

  if ( ImageDimension == 2 )
    {
    this->InitializeIndices();
    }
  else if ( ImageDimension == 3 )
    {
    this->InitializeOffsetsAndIndices();
    }
  else
    {
    itkExceptionMacro( << "Image dimension must be equal to 2 or 3." );
    }
}

template<class TImage>
BinaryWellComposedImageFilter<TImage>
::~BinaryWellComposedImageFilter()
{
}

template<class TImage>
void
BinaryWellComposedImageFilter<TImage>
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
    this->Make2DImageWellComposed();
    }
  else if ( ImageDimension == 3 )
    {
    this->CountCriticalConfigurations();
    exit( 0 );
    if ( this->m_NumberOfC1Configurations3D > 0 
         || this->m_NumberOfC2Configurations3D > 0 )
      { 
      this->LocateCriticalConfigurations();
      this->Make3DImageWellComposed();
      }
    }
  // Extract the output image region to match input image region
  typedef itk::ExtractImageFilter<ImageType, ImageType> CropperType;
  typename CropperType::Pointer cropper = CropperType::New();
  typename ImageType::RegionType region;
  typename ImageType::RegionType::SizeType size;
  typename ImageType::RegionType::IndexType index;
  region.SetSize( this->GetInput()->GetRequestedRegion().GetSize() );
  region.SetIndex( this->GetInput()->GetRequestedRegion().GetIndex() );
  cropper->SetInput( this->GetOutput() );
  cropper->SetExtractionRegion( region );
  cropper->Update();  
  
  this->GraftOutput( cropper->GetOutput() );
}

/*
 * 2-D
 */
template<class TImage>
void
BinaryWellComposedImageFilter<TImage>
::Make2DImageWellComposed()
{
  this->m_NumberOfC1Configurations = 0;
  this->m_NumberOfC2Configurations = 0;
  this->m_NumberOfC3Configurations = 0;
  this->m_NumberOfC4Configurations = 0;

  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
        this->GetOutput()->GetRequestedRegion() );

  Array<PixelType> NeighborhoodPixels( 9 );

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
        NeighborhoodPixels[j] = It.GetPixel( this->m_RotationIndices[i][j] );
        }
      if ( this->IsCriticalC1Configuration( NeighborhoodPixels ) )
        {
        It.SetPixel( this->m_RotationIndices[i][4], this->m_ForegroundValue );
        this->m_NumberOfC1Configurations++;
        break;
        }
      else if ( this->IsCriticalC2Configuration( NeighborhoodPixels ) )
        {
        It.SetPixel( this->m_RotationIndices[i][4], this->m_ForegroundValue );
        this->m_NumberOfC2Configurations++;
        break;
        }
      else if ( this->IsCriticalC3Configuration( NeighborhoodPixels ) )
        {
        It.SetPixel( this->m_RotationIndices[i][4], this->m_ForegroundValue );
        It.SetPixel( this->m_RotationIndices[i][7], this->m_ForegroundValue );
        this->m_NumberOfC3Configurations++;
        break;
        }
      else if ( this->IsCriticalC4Configuration( NeighborhoodPixels ) )
        {
        It.SetPixel( this->m_RotationIndices[i][4], this->m_ForegroundValue );
        It.SetPixel( this->m_RotationIndices[i][7], this->m_ForegroundValue );
        if ( this->IsSpecialCaseOfC4Configuration( 
               It.GetIndex(), It.GetIndex( this->m_RotationIndices[i][6] ), 
               It.GetIndex( this->m_RotationIndices[i][7] ) ) )
          {
          It.SetPixel( this->m_RotationIndices[i][6], this->m_ForegroundValue );
          }
        this->m_NumberOfC4Configurations++;
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
        NeighborhoodPixels[j] = It.GetPixel( this->m_ReflectionIndices[i][j] );
        }
      if ( !this->m_FullInvariance && this->IsCriticalC1Configuration( NeighborhoodPixels ) )
        {
        It.SetPixel( this->m_ReflectionIndices[i][4], this->m_ForegroundValue );
        this->m_NumberOfC1Configurations++;
        break;
        }
      else if ( !this->m_FullInvariance && this->IsCriticalC2Configuration( NeighborhoodPixels ) )
        {
        It.SetPixel( this->m_ReflectionIndices[i][4], this->m_ForegroundValue );
        this->m_NumberOfC2Configurations++;
        break;
        }
      else if ( this->IsCriticalC3Configuration( NeighborhoodPixels ) )
        {
        It.SetPixel( this->m_ReflectionIndices[i][4], this->m_ForegroundValue );
        It.SetPixel( this->m_ReflectionIndices[i][7], this->m_ForegroundValue );
        this->m_NumberOfC3Configurations++;
        break;
        }
      else if ( this->IsCriticalC4Configuration( NeighborhoodPixels ) )
        {
        It.SetPixel( this->m_ReflectionIndices[i][4], this->m_ForegroundValue );
        It.SetPixel( this->m_ReflectionIndices[i][7], this->m_ForegroundValue );
        if ( this->IsSpecialCaseOfC4Configuration( 
               It.GetIndex(), It.GetIndex( this->m_ReflectionIndices[i][6] ), 
               It.GetIndex( this->m_ReflectionIndices[i][7] ) ) )
          {
          It.SetPixel( this->m_ReflectionIndices[i][6], 
                       this->m_ForegroundValue );
          }
        this->m_NumberOfC4Configurations++;
        break;
        }
      if ( !this->m_FullInvariance )
        {
        break;
        }
      }  
    }
    
  itkDebugMacro( << this->m_NumberOfC1Configurations 
                 << " C1 configurations were located and repaired." ); 
  itkDebugMacro( << this->m_NumberOfC2Configurations 
                 << " C2 configurations were located and repaired." ); 
  itkDebugMacro( << this->m_NumberOfC3Configurations 
                 << " C3 configurations were located and repaired." ); 
  itkDebugMacro( << this->m_NumberOfC4Configurations 
                 << " C4 configurations were located and repaired." ); 
}

/*
 * 2-D
 */
template<class TImage>
bool
BinaryWellComposedImageFilter<TImage>
::IsCriticalC1Configuration( Array<PixelType> neighborhood )
{
  return ( neighborhood[0] == this->m_BackgroundValue &&
           neighborhood[1] == this->m_ForegroundValue &&
           neighborhood[3] == this->m_ForegroundValue &&
           neighborhood[4] == this->m_BackgroundValue &&
           neighborhood[8] == this->m_BackgroundValue );
}

/*
 * 2-D
 */
template<class TImage>
bool
BinaryWellComposedImageFilter<TImage>
::IsCriticalC2Configuration( Array<PixelType> neighborhood )
{
  return ( neighborhood[0] == this->m_BackgroundValue &&
           neighborhood[1] == this->m_ForegroundValue &&
           neighborhood[3] == this->m_ForegroundValue &&
           neighborhood[4] == this->m_BackgroundValue &&
           neighborhood[8] == this->m_ForegroundValue &&
           ( neighborhood[5] == this->m_ForegroundValue ||
             neighborhood[7] == this->m_ForegroundValue ) );
}

/*
 * 2-D
 */
template<class TImage>
bool
BinaryWellComposedImageFilter<TImage>
::IsCriticalC3Configuration( Array<PixelType> neighborhood )
{
  return ( neighborhood[0] == this->m_BackgroundValue &&
           neighborhood[1] == this->m_ForegroundValue &&
           neighborhood[3] == this->m_ForegroundValue &&
           neighborhood[4] == this->m_BackgroundValue &&
           neighborhood[5] == this->m_BackgroundValue &&
           neighborhood[6] == this->m_ForegroundValue &&
           neighborhood[7] == this->m_BackgroundValue &&
           neighborhood[8] == this->m_ForegroundValue );
}

/*
 * 2-D
 */
template<class TImage>
bool
BinaryWellComposedImageFilter<TImage>
::IsCriticalC4Configuration( Array<PixelType> neighborhood )
{
  return ( neighborhood[0] == this->m_BackgroundValue &&
           neighborhood[1] == this->m_ForegroundValue &&
           neighborhood[3] == this->m_ForegroundValue &&
           neighborhood[4] == this->m_BackgroundValue &&
           neighborhood[5] == this->m_BackgroundValue &&
           neighborhood[6] == this->m_BackgroundValue &&
           neighborhood[7] == this->m_BackgroundValue &&
           neighborhood[8] == this->m_ForegroundValue );
}

/*
 * 2-D
 */
template<class TImage>
bool
BinaryWellComposedImageFilter<TImage>
::IsSpecialCaseOfC4Configuration( IndexType idx, 
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
           this->GetOutput()->GetPixel( idxa ) == this->m_BackgroundValue &&
           this->GetOutput()->GetPixel( idxb ) == this->m_ForegroundValue );
}

/*
 * 2-D
 */
template<class TImage>
void
BinaryWellComposedImageFilter<TImage>
::InitializeIndices()
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
BinaryWellComposedImageFilter<TImage>
::LocateCriticalConfigurations()
{
  this->m_CriticalConfigurationIndices.clear();
 
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
        this->GetOutput()->GetLargestPossibleRegion() );

  Array<unsigned char> p( 8 );

  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( !It.InBounds() )
      {
      continue;
      } 
   
    for ( unsigned int i = 0; i < 8; i++ )
      {
      if ( It.GetPixel( this->m_Offsets8[i] ) == this->m_ForegroundValue )
        {
        p[i] = NumericTraits<unsigned char>::One;
        }
      else
        {
        p[i] = NumericTraits<unsigned char>::Zero;
        }
      }

    IndexType idx = It.GetIndex();   
    
    // Find C1 configurations
    for ( unsigned int i = 0; i < 6; i++ )
      { 
      for ( unsigned int j = 0; j < 2; j++ )
        {  
        if ( p[this->m_C1IndicesI[i][0]] == static_cast<unsigned char>( j ) && 
             p[this->m_C1IndicesI[i][1]] == static_cast<unsigned char>( j ) &&
             p[this->m_C1IndicesI[i][2]] != static_cast<unsigned char>( j ) && 
             p[this->m_C1IndicesI[i][3]] != static_cast<unsigned char>( j ) )
          {
          if ( this->IsChangeSafe( idx 
                + this->m_Offsets8[this->m_C1IndicesI[i][2*j]] ) ) 
            {
            p[this->m_C1IndicesI[i][2*j]] = NumericTraits<unsigned char>::One;  
            It.SetPixel( this->m_Offsets8[this->m_C1IndicesI[i][2*j]], 
                         this->m_ForegroundValue );
            }
          else if ( this->IsChangeSafe( idx 
                     + this->m_Offsets8[this->m_C1IndicesI[i][2*j+1]] ) ) 
            {
            p[this->m_C1IndicesI[i][2*j+1]] = NumericTraits<unsigned char>::One;
            It.SetPixel( this->m_Offsets8[this->m_C1IndicesI[i][2*j+1]], 
                         this->m_ForegroundValue );
            }
          else 
            {
            this->m_CriticalConfigurationIndices.push_back( idx );
            }
          }
        }
      }  

    // Find C2 configurations  
    for ( unsigned int i = 0; i < 4; i++ )
      { 
      if ( p[this->m_C2IndicesI[i][0]] &&  p[this->m_C2IndicesI[i][1]] &&
          !p[this->m_C2IndicesI[i][2]] && !p[this->m_C2IndicesI[i][3]] &&
          !p[this->m_C2IndicesI[i][4]] && !p[this->m_C2IndicesI[i][5]] &&
          !p[this->m_C2IndicesI[i][6]] && !p[this->m_C2IndicesI[i][7]] ) 
        {
        bool SafeConfigurationFound = false;
        for ( unsigned int j = 2; j < 8; j++ )
          {
          if ( this->IsChangeSafe( idx 
                + this->m_Offsets8[this->m_C2IndicesI[i][j]] ) )
            {
            p[this->m_C2IndicesI[i][j]] = NumericTraits<unsigned char>::One;  
            It.SetPixel( this->m_Offsets8[this->m_C2IndicesI[i][j]], 
                         this->m_ForegroundValue );
            SafeConfigurationFound = true;
            break;
            }
          }  
        if ( !SafeConfigurationFound )
          {
          this->m_CriticalConfigurationIndices.push_back( idx );
          }
        }
      else if ( !p[this->m_C2IndicesI[i][0]] && !p[this->m_C2IndicesI[i][1]] &&
                 p[this->m_C2IndicesI[i][2]] &&  p[this->m_C2IndicesI[i][3]] &&
                 p[this->m_C2IndicesI[i][4]] &&  p[this->m_C2IndicesI[i][5]] &&
                 p[this->m_C2IndicesI[i][6]] &&  p[this->m_C2IndicesI[i][7]] ) 
        {
        bool SafeConfigurationFound = false;
        for ( unsigned int j = 0; j < 2; j++ )
          {
          if ( this->IsChangeSafe( idx 
                + this->m_Offsets8[this->m_C2IndicesI[i][j]] ) )
            {
            p[this->m_C2IndicesI[i][j]] = NumericTraits<unsigned char>::One;  
            It.SetPixel( this->m_Offsets8[this->m_C2IndicesI[i][j]], 
                         this->m_ForegroundValue );
            SafeConfigurationFound = true;
            break;
            }
          }  
        if ( !SafeConfigurationFound )
          {
          this->m_CriticalConfigurationIndices.push_back( idx );
          }
        }
      }
    }
}

/*
 * 3-D
 */
template<class TImage>
void
BinaryWellComposedImageFilter<TImage>
::Make3DImageWellComposed()
{
  Array<unsigned char> p( 8 );
  unsigned int iter = 0;
 
  this->m_NumberOfC1Configurations3D = 0;
  this->m_NumberOfC2Configurations3D = 0;

  while ( !this->m_CriticalConfigurationIndices.empty() ) 
    {
    
    itkDebugMacro( << "Iteration " << iter++ << " (" << 
                   this->m_CriticalConfigurationIndices.size() << " critical configurations). " );
    
    this->m_NewCriticalConfigurationIndices.clear();

    while ( !this->m_CriticalConfigurationIndices.empty() ) 
      {
      IndexType idx = this->m_CriticalConfigurationIndices.front();
      this->m_CriticalConfigurationIndices.pop_front();

      for (unsigned int i = 0; i < 8; i++)
        {
        if ( this->GetOutput()->GetPixel( idx  + this->m_Offsets8[i] ) 
               == this->m_ForegroundValue )
          {
          p[i] = NumericTraits<unsigned char>::One;
          }
        else
          {
          p[i] = NumericTraits<unsigned char>::Zero;
          }
        }  
          
      //  Deal with the C1 configurations
      for ( unsigned int j = 0; j < 6; j++ )
        {
        for ( unsigned int k = 0; k <= 1; k++ )
          {
          if ( p[this->m_C1IndicesI[j][0]]
                == static_cast<unsigned char>( k ) && 
               p[this->m_C1IndicesI[j][1]]
                == static_cast<unsigned char>( k ) &&
               p[this->m_C1IndicesI[j][2]]
                != static_cast<unsigned char>( k ) && 
               p[this->m_C1IndicesI[j][3]]
                != static_cast<unsigned char>( k ) )
            {
            this->MakeRandomChange(
              p[this->m_C1IndicesI[j][2*k  ]], 
              idx + this->m_Offsets8[this->m_C1IndicesI[j][2*k  ]],
              p[this->m_C1IndicesI[j][2*k+1]], 
              idx + this->m_Offsets8[this->m_C1IndicesI[j][2*k+1]] );
            }
          }
        }
      
      //  Deal with the C2 configurations
      for ( unsigned int j = 0; j < 4; j++ )
        { 
        if ( p[this->m_C2IndicesI[j][0]] &&  p[this->m_C2IndicesI[j][1]] &&
            !p[this->m_C2IndicesI[j][2]] && !p[this->m_C2IndicesI[j][3]] &&
            !p[this->m_C2IndicesI[j][4]] && !p[this->m_C2IndicesI[j][5]] &&
            !p[this->m_C2IndicesI[j][6]] && !p[this->m_C2IndicesI[j][7]] ) 
          {
          this->MakeRandomChange( 
            p[this->m_C2IndicesI[j][2]], 
            idx + this->m_Offsets8[this->m_C2IndicesI[j][2]], 
            p[this->m_C2IndicesI[j][3]], 
            idx + this->m_Offsets8[this->m_C2IndicesI[j][3]], 
            p[this->m_C2IndicesI[j][4]], 
            idx + this->m_Offsets8[this->m_C2IndicesI[j][4]],
            p[this->m_C2IndicesI[j][5]], 
            idx + this->m_Offsets8[this->m_C2IndicesI[j][5]], 
            p[this->m_C2IndicesI[j][6]], 
            idx + this->m_Offsets8[this->m_C2IndicesI[j][6]], 
            p[this->m_C2IndicesI[j][7]], 
            idx + this->m_Offsets8[this->m_C2IndicesI[j][7]] );
          break;
          }
        else if ( !p[this->m_C2IndicesI[j][0]] && 
                  !p[this->m_C2IndicesI[j][1]] &&
                   p[this->m_C2IndicesI[j][2]] &&  
                   p[this->m_C2IndicesI[j][3]] &&
                   p[this->m_C2IndicesI[j][4]] &&  
                   p[this->m_C2IndicesI[j][5]] &&
                   p[this->m_C2IndicesI[j][6]] &&  
                   p[this->m_C2IndicesI[j][7]] ) 
          {
          this->MakeRandomChange(
            p[this->m_C2IndicesI[j][0]], 
            idx + this->m_Offsets8[this->m_C2IndicesI[j][0]],
            p[this->m_C2IndicesI[j][1]], 
            idx + this->m_Offsets8[this->m_C2IndicesI[j][1]] );
          break;
          }
        }
      }
    this->m_CriticalConfigurationIndices 
      = this->m_NewCriticalConfigurationIndices;
    }
}

/*
 * 3-D
 */
template<class TImage>
bool
BinaryWellComposedImageFilter<TImage>
::IsChangeSafe( IndexType idx )
{
  Array<unsigned char> p( 4 );
  for ( unsigned int i = 0; i < 12; i++ )
    {
    for ( unsigned int j = 0; j < 4; j++ )
      {
      if ( this->m_C1IndicesII[i][j] == 13 )
        {   
        if ( this->GetOutput()->GetPixel( idx ) 
              == this->m_ForegroundValue )
          {
          p[j] = NumericTraits<unsigned char>::Zero;
          }
        else
          {
          p[j] = NumericTraits<unsigned char>::One;
          }
        }
      else
        {   
        if ( this->GetOutput()->GetPixel( idx + 
               this->m_Offsets27[this->m_C1IndicesII[i][j]] )
               == this->m_ForegroundValue )
          {
          p[j] = NumericTraits<unsigned char>::One;
          }
        else
          {
          p[j] = NumericTraits<unsigned char>::Zero;
          }
        }
      }
    if ( this->IsCriticalC1Configuration3D( p ) )
      {
      return false;
      }  
    } 

  p.SetSize( 8 );
  for ( unsigned int i = 0; i < 8; i++ )
    {
    for ( unsigned int j = 0; j < 8; j++ )
      {
      if ( this->m_C2IndicesII[i][j] == 13 )
        {   
        if ( this->GetOutput()->GetPixel( idx ) 
               == this->m_ForegroundValue ) 
          {
          p[j] = NumericTraits<unsigned char>::Zero;
          }
        else
          {
          p[j] = NumericTraits<unsigned char>::One;
          }
        }
      else
        {   
        if ( this->GetOutput()->GetPixel( idx 
               + this->m_Offsets27[this->m_C2IndicesII[i][j]] )
               == this->m_ForegroundValue )
          {
          p[j] = NumericTraits<unsigned char>::One;
          }
        else
          {  
          p[j] = NumericTraits<unsigned char>::Zero;
          }
        }
      }
    if ( this->IsCriticalC2Configuration3D( p ) )
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
BinaryWellComposedImageFilter<TImage>
::MakeRandomChange( unsigned char &p0, IndexType idx0, 
                    unsigned char &p1, IndexType idx1 )
{
  typedef typename Statistics
     ::MersenneTwisterRandomVariateGenerator GeneratorType; 
  typename GeneratorType::Pointer generator = GeneratorType::New();

  if ( this->IsChangeSafe( idx0 ) )
    {
    p0 = NumericTraits<unsigned char>::One;
    this->GetOutput()->SetPixel( idx0, this->m_ForegroundValue );
    } 
  else
    {
    if ( this->IsChangeSafe( idx1 ) )
      {
      p1 = NumericTraits<unsigned char>::One;
      this->GetOutput()->SetPixel( idx1, this->m_ForegroundValue );
      } 
    else
      { 
      if ( !generator->GetIntegerVariate( 1 ) )
        {  
        p0 = NumericTraits<unsigned char>::One;
        this->InsertCriticalConfiguration( idx0 );
        this->GetOutput()->SetPixel( idx0, this->m_ForegroundValue );
        }
      else
        {
        p1 = NumericTraits<unsigned char>::One;
        this->InsertCriticalConfiguration( idx1 );
        this->GetOutput()->SetPixel( idx1, this->m_ForegroundValue );
        }  
      }
    }
}

template<class TImage>
void
BinaryWellComposedImageFilter<TImage>
::MakeRandomChange( unsigned char &p0, IndexType idx0, 
                    unsigned char &p1, IndexType idx1,
                    unsigned char &p2, IndexType idx2, 
                    unsigned char &p3, IndexType idx3,  
                    unsigned char &p4, IndexType idx4, 
                    unsigned char &p5, IndexType idx5)
{
  if ( this->IsChangeSafe( idx0 ) )
    {
    p0 = NumericTraits<unsigned char>::One;
    this->GetOutput()->SetPixel( idx0, 
      this->m_ForegroundValue );
    return;
    } 
  if ( this->IsChangeSafe( idx1 ) )
    {
    p1 = NumericTraits<unsigned char>::One;
    this->GetOutput()->SetPixel( idx1, 
      this->m_ForegroundValue );
    return;
    } 
  if ( this->IsChangeSafe( idx2 ) )
    {
    p2 = NumericTraits<unsigned char>::One;
    this->GetOutput()->SetPixel( idx2, this->m_ForegroundValue );
    return;
    } 
  if ( this->IsChangeSafe( idx3 ) )
    {
    p3 = NumericTraits<unsigned char>::One;
    this->GetOutput()->SetPixel( idx3, this->m_ForegroundValue );
    return;
    } 
  if ( this->IsChangeSafe( idx4 ) )
    {
    p4 = NumericTraits<unsigned char>::One;
    this->GetOutput()->SetPixel( idx4, this->m_ForegroundValue );
    return;
    } 
  if ( this->IsChangeSafe( idx5 ) )
    {
    p5 = NumericTraits<unsigned char>::One;
    this->GetOutput()->SetPixel( idx5, this->m_ForegroundValue );
    return;
    } 

  typedef typename Statistics
    ::MersenneTwisterRandomVariateGenerator GeneratorType; 
  typename GeneratorType::Pointer generator = GeneratorType::New();

  switch ( generator->GetIntegerVariate( 5 ) )
    {
    case 0: 
      p0 = NumericTraits<unsigned char>::One;
      this->InsertCriticalConfiguration( idx0 );
      this->GetOutput()->SetPixel( idx0, this->m_ForegroundValue );
      break;
    case 1: 
      p1 = NumericTraits<unsigned char>::One;
      this->InsertCriticalConfiguration( idx1 );
      this->GetOutput()->SetPixel( idx1, this->m_ForegroundValue );
      break;
    case 2: 
      p2 = NumericTraits<unsigned char>::One;
      this->InsertCriticalConfiguration( idx2 );
      this->GetOutput()->SetPixel( idx2, this->m_ForegroundValue );
      break;
    case 3: 
      p3 = NumericTraits<unsigned char>::One;
      this->InsertCriticalConfiguration( idx3 );
      this->GetOutput()->SetPixel( idx3, this->m_ForegroundValue );
      break;
    case 4: 
      p4 = NumericTraits<unsigned char>::One;
      this->InsertCriticalConfiguration( idx4 );
      this->GetOutput()->SetPixel( idx4, this->m_ForegroundValue );
      break;
    case 5: 
      p5 = NumericTraits<unsigned char>::One;
      this->InsertCriticalConfiguration( idx5 );
      this->GetOutput()->SetPixel( idx5, this->m_ForegroundValue );
      break;
    }
}


/*
 * 3-D
 */
template<class TImage>
void
BinaryWellComposedImageFilter<TImage>
::InsertCriticalConfiguration( IndexType idx )
{
  Array<unsigned char> p( 4 );  

  for ( unsigned int i = 0; i < 12; i++ )
    {
    for ( unsigned int j = 0; j < 4; j++ )
      {
      if ( this->m_C1IndicesII[i][j] == 13 )
        {
        if ( this->GetOutput()->GetPixel(idx) 
             == this->m_ForegroundValue )
          {
          p[j] = NumericTraits<unsigned char>::Zero;
          }
        else
          {
          p[j] = NumericTraits<unsigned char>::One;
          }
        }
      else
        {
        if ( this->GetOutput()->GetPixel( idx 
              + this->m_Offsets27[this->m_C1IndicesII[i][j]] )
                 == this->m_ForegroundValue ) 
          {
          p[j] = NumericTraits<unsigned char>::One;
          }
        else
          {
          p[j] = NumericTraits<unsigned char>::Zero;
          }
        }
      }
    if ( this->IsCriticalC1Configuration3D( p ) )
      {
      this->m_NewCriticalConfigurationIndices.push_back( 
        idx + this->m_Offsets27[this->m_C1IndicesII[i][0]] );
      this->m_NumberOfC1Configurations3D++;
      }  
    } 

  p.SetSize( 8 );
  for ( unsigned int i = 0; i < 8; i++ )
    {
    for ( unsigned int j = 0; j < 8; j++ )
      {
      if ( this->m_C2IndicesII[i][j] == 13 )
        {
        if ( this->GetOutput()->GetPixel( idx )
             == this->m_ForegroundValue )
          {
          p[j] = NumericTraits<unsigned char>::Zero;
          }
        else
          {
          p[j] = NumericTraits<unsigned char>::One;
          }
        }
      else
        {
        if ( this->GetOutput()->GetPixel( idx 
              + this->m_Offsets27[this->m_C2IndicesII[i][j]] )
                 == this->m_ForegroundValue )
          {
          p[j] = NumericTraits<unsigned char>::One;
          }
        else
          {
          p[j] = NumericTraits<unsigned char>::Zero;
          }
        } 
      }
    if ( this->IsCriticalC2Configuration3D( p ) )
      {
      this->m_NewCriticalConfigurationIndices.push_back( 
        idx + this->m_Offsets27[this->m_C2IndicesII[i][0]] );
      this->m_NumberOfC2Configurations3D++;
      }  
    } 
}

/*
 * 3-D
 */
template<class TImage>
bool
BinaryWellComposedImageFilter<TImage>
::IsCriticalC1Configuration3D(Array<unsigned char> p)
{
  return (p[0] == p[2]) && (p[1] == p[3]) && (p[0] != p[1]);
}

/*
 * 3-D
 */
template<class TImage>
bool
BinaryWellComposedImageFilter<TImage>
::IsCriticalC2Configuration3D(Array<unsigned char> p)
{
  return ( (  p[0] && !p[1] && !p[2] && !p[3] && 
             !p[4] && !p[5] &&  p[6] && !p[7] ) ||
           ( !p[0] &&  p[1] &&  p[2] &&  p[3] &&  
              p[4] &&  p[5] && !p[6] &&  p[7] ) ||
           ( !p[0] && !p[1] &&  p[2] && !p[3] &&  
              p[4] && !p[5] && !p[6] && !p[7] ) ||
           (  p[0] &&  p[1] && !p[2] &&  p[3] && 
             !p[4] &&  p[5] &&  p[6] &&  p[7] ) ||
           ( !p[0] &&  p[1] && !p[2] && !p[3] && 
             !p[4] && !p[5] && !p[6] &&  p[7] ) ||
           (  p[0] && !p[1] &&  p[2] &&  p[3] &&  
              p[4] &&  p[5] &&  p[6] && !p[7] ) ||
           ( !p[0] && !p[1] && !p[2] &&  p[3] && 
             !p[4] &&  p[5] && !p[6] && !p[7] ) ||
           (  p[0] &&  p[1] &&  p[2] && !p[3] &&  
              p[4] && !p[5] &&  p[6] &&  p[7] ) );
}

/*
 * 3-D
 */
template<class TImage>
void
BinaryWellComposedImageFilter<TImage>
::CountCriticalConfigurations()
{
  this->m_NumberOfC1Configurations3D = 0;
  this->m_NumberOfC2Configurations3D = 0;

  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->GetOutput(),
    this->GetOutput()->GetRequestedRegion() );

  Array<unsigned char> p( 8 );

  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( !It.InBounds() )
      {
      continue;
      } 
   
    for ( unsigned int i = 0; i < 8; i++ )
      {
      if ( It.GetPixel( this->m_Offsets8[i] ) == this->m_ForegroundValue )
        {
        p[i] = NumericTraits<unsigned char>::One;
        }
      else
        {
        p[i] = NumericTraits<unsigned char>::Zero;
        }
      }
    IndexType idx = It.GetIndex();   
    
    // Find C1 configurations
    for ( unsigned int i = 0; i < 6; i++ )
      { 
      for ( unsigned int j = 0; j <= 1; j++ )
        {  
        if ( p[this->m_C1IndicesI[i][0]] != static_cast<unsigned char>( j ) && 
             p[this->m_C1IndicesI[i][1]] != static_cast<unsigned char>( j ) &&
             p[this->m_C1IndicesI[i][2]] == static_cast<unsigned char>( j ) && 
             p[this->m_C1IndicesI[i][3]] == static_cast<unsigned char>( j ) )
          {
          this->m_NumberOfC1Configurations3D++; 
          break;
          }
        }
      }  
    
    // Find C2 configurations
    for ( unsigned int i = 0; i < 4; i++ )
      { 
      bool foundC2 = false;
      for ( unsigned int j = 0; j <= 1; j++ )
        {  
        if ( p[this->m_C2IndicesI[i][0]] != static_cast<unsigned char>( j ) &&
             p[this->m_C2IndicesI[i][1]] != static_cast<unsigned char>( j ) &&
             p[this->m_C2IndicesI[i][2]] == static_cast<unsigned char>( j ) &&
             p[this->m_C2IndicesI[i][3]] == static_cast<unsigned char>( j ) &&
             p[this->m_C2IndicesI[i][4]] == static_cast<unsigned char>( j ) &&
             p[this->m_C2IndicesI[i][5]] == static_cast<unsigned char>( j ) &&
             p[this->m_C2IndicesI[i][6]] == static_cast<unsigned char>( j ) &&
             p[this->m_C2IndicesI[i][7]] == static_cast<unsigned char>( j ))
          {
          this->m_NumberOfC2Configurations3D++;
          foundC2 = true;
          break;
          }
        }
      if ( foundC2 )
        {
        break;
        }    
      }
    }
  itkDebugMacro( << "Number of C1 configurations = "
                 << this->m_NumberOfC1Configurations3D );
  itkDebugMacro( << "Number of C2 configurations = " 
                 << this->m_NumberOfC2Configurations3D );
  std::cout << "There are a total of " << this->m_NumberOfC1Configurations3D << " C1 and " << this->m_NumberOfC2Configurations3D
                 << " C2 critical configurations." << std::endl; 
                 
}

/*
 * 3-D
 */
template<class TImage>
void
BinaryWellComposedImageFilter<TImage>
::InitializeOffsetsAndIndices()
{
  for ( unsigned int i = 0; i <  6; i++ )
    {
    this->m_C1IndicesI[i].SetSize( 4 );
    }
  for ( unsigned int i = 0; i <  4; i++ )
    {
    this->m_C2IndicesI[i].SetSize( 8 );
    }
  for ( unsigned int i = 0; i < 12; i++ )  
    {
    this->m_C1IndicesII[i].SetSize( 4 );
    }
  for ( unsigned int i = 0; i <  8; i++ )  
    {
    this->m_C2IndicesII[i].SetSize( 8 );
    }

  this->m_Offsets8[0][0] = 0;  
  this->m_Offsets8[0][1] = 0;   
  this->m_Offsets8[0][2] = 0;   

  this->m_Offsets8[1][0] = 1;  
  this->m_Offsets8[1][1] = 0;   
  this->m_Offsets8[1][2] = 0;   

  this->m_Offsets8[2][0] = 1;  
  this->m_Offsets8[2][1] = 0;   
  this->m_Offsets8[2][2] = 1;   

  this->m_Offsets8[3][0] = 0;  
  this->m_Offsets8[3][1] = 0;   
  this->m_Offsets8[3][2] = 1;   

  this->m_Offsets8[4][0] = 0;  
  this->m_Offsets8[4][1] = 1;   
  this->m_Offsets8[4][2] = 0;   

  this->m_Offsets8[5][0] = 1;  
  this->m_Offsets8[5][1] = 1;   
  this->m_Offsets8[5][2] = 0;   

  this->m_Offsets8[6][0] = 1;  
  this->m_Offsets8[6][1] = 1;   
  this->m_Offsets8[6][2] = 1;   

  this->m_Offsets8[7][0] = 0;  
  this->m_Offsets8[7][1] = 1;   
  this->m_Offsets8[7][2] = 1;   

  unsigned int n = 0;
  for ( int i = -1; i <= 1; i++ )
    {  
    for ( int j = -1; j <= 1; j++ )
      {
      for ( int k = -1; k <= 1; k++ )
        {
        this->m_Offsets27[n][0] = i;
        this->m_Offsets27[n][1] = j;
        this->m_Offsets27[n][2] = k;
        n++;
        }
      }
    } 

  this->m_C1IndicesI[0][0] = 0;
  this->m_C1IndicesI[0][1] = 2;
  this->m_C1IndicesI[0][2] = 1;
  this->m_C1IndicesI[0][3] = 3;

  this->m_C1IndicesI[1][0] = 4;
  this->m_C1IndicesI[1][1] = 6;
  this->m_C1IndicesI[1][2] = 5;
  this->m_C1IndicesI[1][3] = 7;

  this->m_C1IndicesI[2][0] = 0;  
  this->m_C1IndicesI[2][1] = 5;  
  this->m_C1IndicesI[2][2] = 1;  
  this->m_C1IndicesI[2][3] = 4;  

  this->m_C1IndicesI[3][0] = 3;
  this->m_C1IndicesI[3][1] = 6;
  this->m_C1IndicesI[3][2] = 2;
  this->m_C1IndicesI[3][3] = 7;

  this->m_C1IndicesI[4][0] = 0;
  this->m_C1IndicesI[4][1] = 7;
  this->m_C1IndicesI[4][2] = 3;
  this->m_C1IndicesI[4][3] = 4;

  this->m_C1IndicesI[5][0] = 2;
  this->m_C1IndicesI[5][1] = 5;
  this->m_C1IndicesI[5][2] = 1;
  this->m_C1IndicesI[5][3] = 6;

  this->m_C2IndicesI[0][0] = 0;
  this->m_C2IndicesI[0][1] = 6;
  this->m_C2IndicesI[0][2] = 1;
  this->m_C2IndicesI[0][3] = 2;
  this->m_C2IndicesI[0][4] = 3;
  this->m_C2IndicesI[0][5] = 4;
  this->m_C2IndicesI[0][6] = 5;
  this->m_C2IndicesI[0][7] = 7;

  this->m_C2IndicesI[1][0] = 2;
  this->m_C2IndicesI[1][1] = 4;
  this->m_C2IndicesI[1][2] = 0;
  this->m_C2IndicesI[1][3] = 1;
  this->m_C2IndicesI[1][4] = 3;
  this->m_C2IndicesI[1][5] = 5;
  this->m_C2IndicesI[1][6] = 6;
  this->m_C2IndicesI[1][7] = 7;

  this->m_C2IndicesI[2][0] = 1;
  this->m_C2IndicesI[2][1] = 7;
  this->m_C2IndicesI[2][2] = 0;
  this->m_C2IndicesI[2][3] = 2;
  this->m_C2IndicesI[2][4] = 3;
  this->m_C2IndicesI[2][5] = 4;
  this->m_C2IndicesI[2][6] = 5;
  this->m_C2IndicesI[2][7] = 6;

  this->m_C2IndicesI[3][0] = 3;
  this->m_C2IndicesI[3][1] = 5;
  this->m_C2IndicesI[3][2] = 0;
  this->m_C2IndicesI[3][3] = 1;
  this->m_C2IndicesI[3][4] = 2;
  this->m_C2IndicesI[3][5] = 4;
  this->m_C2IndicesI[3][6] = 6;
  this->m_C2IndicesI[3][7] = 7;

  this->m_C1IndicesII[ 0][0] =  3;
  this->m_C1IndicesII[ 0][1] = 12;
  this->m_C1IndicesII[ 0][2] = 13;
  this->m_C1IndicesII[ 0][3] =  4;

  this->m_C1IndicesII[ 1][0] = 12;
  this->m_C1IndicesII[ 1][1] = 21;
  this->m_C1IndicesII[ 1][2] = 22;
  this->m_C1IndicesII[ 1][3] = 13;

  this->m_C1IndicesII[ 2][0] =  4;
  this->m_C1IndicesII[ 2][1] = 13;
  this->m_C1IndicesII[ 2][2] = 14;
  this->m_C1IndicesII[ 2][3] =  5;

  this->m_C1IndicesII[ 3][0] = 13;
  this->m_C1IndicesII[ 3][1] = 22;
  this->m_C1IndicesII[ 3][2] = 23;
  this->m_C1IndicesII[ 3][3] = 14;

  this->m_C1IndicesII[ 4][0] =  1;
  this->m_C1IndicesII[ 4][1] = 10;
  this->m_C1IndicesII[ 4][2] = 13;
  this->m_C1IndicesII[ 4][3] =  4;

  this->m_C1IndicesII[ 5][0] = 10;
  this->m_C1IndicesII[ 5][1] = 19;
  this->m_C1IndicesII[ 5][2] = 22;
  this->m_C1IndicesII[ 5][3] = 13;

  this->m_C1IndicesII[ 6][0] =  9;
  this->m_C1IndicesII[ 6][1] = 10;
  this->m_C1IndicesII[ 6][2] = 13;
  this->m_C1IndicesII[ 6][3] = 12;

  this->m_C1IndicesII[ 7][0] = 10;
  this->m_C1IndicesII[ 7][1] = 11;
  this->m_C1IndicesII[ 7][2] = 14;
  this->m_C1IndicesII[ 7][3] = 13;

  this->m_C1IndicesII[ 8][0] =  4;
  this->m_C1IndicesII[ 8][1] = 13;
  this->m_C1IndicesII[ 8][2] = 16;
  this->m_C1IndicesII[ 8][3] =  7;

  this->m_C1IndicesII[ 9][0] = 13;
  this->m_C1IndicesII[ 9][1] = 22;
  this->m_C1IndicesII[ 9][2] = 25;
  this->m_C1IndicesII[ 9][3] = 16;

  this->m_C1IndicesII[10][0] = 12;
  this->m_C1IndicesII[10][1] = 13;
  this->m_C1IndicesII[10][2] = 16;
  this->m_C1IndicesII[10][3] = 15;

  this->m_C1IndicesII[11][0] = 13;
  this->m_C1IndicesII[11][1] = 14;
  this->m_C1IndicesII[11][2] = 17;
  this->m_C1IndicesII[11][3] = 16;

  this->m_C2IndicesII[0][0] =  0;
  this->m_C2IndicesII[0][1] =  9;
  this->m_C2IndicesII[0][2] = 10;
  this->m_C2IndicesII[0][3] =  1;
  this->m_C2IndicesII[0][4] =  3;
  this->m_C2IndicesII[0][5] = 12;
  this->m_C2IndicesII[0][6] = 13;
  this->m_C2IndicesII[0][7] =  4;

  this->m_C2IndicesII[1][0] =  9;
  this->m_C2IndicesII[1][1] = 18;
  this->m_C2IndicesII[1][2] = 19;
  this->m_C2IndicesII[1][3] = 10;
  this->m_C2IndicesII[1][4] = 12;
  this->m_C2IndicesII[1][5] = 21;
  this->m_C2IndicesII[1][6] = 22;
  this->m_C2IndicesII[1][7] = 13;

  this->m_C2IndicesII[2][0] =  1;
  this->m_C2IndicesII[2][1] = 10;
  this->m_C2IndicesII[2][2] = 11;
  this->m_C2IndicesII[2][3] =  2;
  this->m_C2IndicesII[2][4] =  4;
  this->m_C2IndicesII[2][5] = 13;
  this->m_C2IndicesII[2][6] = 14;
  this->m_C2IndicesII[2][7] =  5;

  this->m_C2IndicesII[3][0] = 10;
  this->m_C2IndicesII[3][1] = 19;
  this->m_C2IndicesII[3][2] = 20;
  this->m_C2IndicesII[3][3] = 11;
  this->m_C2IndicesII[3][4] = 13;
  this->m_C2IndicesII[3][5] = 22;
  this->m_C2IndicesII[3][6] = 23;
  this->m_C2IndicesII[3][7] = 14;

  this->m_C2IndicesII[4][0] =  3;
  this->m_C2IndicesII[4][1] = 12;
  this->m_C2IndicesII[4][2] = 13;
  this->m_C2IndicesII[4][3] =  4;
  this->m_C2IndicesII[4][4] =  6;
  this->m_C2IndicesII[4][5] = 15;
  this->m_C2IndicesII[4][6] = 16;
  this->m_C2IndicesII[4][7] =  7;

  this->m_C2IndicesII[5][0] = 12;
  this->m_C2IndicesII[5][1] = 21;
  this->m_C2IndicesII[5][2] = 22;
  this->m_C2IndicesII[5][3] = 13;
  this->m_C2IndicesII[5][4] = 15;
  this->m_C2IndicesII[5][5] = 24;
  this->m_C2IndicesII[5][6] = 25;
  this->m_C2IndicesII[5][7] = 16;

  this->m_C2IndicesII[6][0] =  4;
  this->m_C2IndicesII[6][1] = 13;
  this->m_C2IndicesII[6][2] = 14;
  this->m_C2IndicesII[6][3] =  5;
  this->m_C2IndicesII[6][4] =  7;
  this->m_C2IndicesII[6][5] = 16;
  this->m_C2IndicesII[6][6] = 17;
  this->m_C2IndicesII[6][7] =  8;

  this->m_C2IndicesII[7][0] = 13;
  this->m_C2IndicesII[7][1] = 22;
  this->m_C2IndicesII[7][2] = 23;
  this->m_C2IndicesII[7][3] = 14;
  this->m_C2IndicesII[7][4] = 16;
  this->m_C2IndicesII[7][5] = 25;
  this->m_C2IndicesII[7][6] = 26;
  this->m_C2IndicesII[7][7] = 17;
}

template <class TImage>
void
BinaryWellComposedImageFilter<TImage>
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Background Value: " << this->m_BackgroundValue << std::endl;
  os << indent << "Foreground Value: " << this->m_ForegroundValue << std::endl;
}

} // end namespace itk

#endif
