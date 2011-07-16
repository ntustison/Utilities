#ifndef __itkBinaryWellComposed2DImageFilter_hxx
#define __itkBinaryWellComposed2DImageFilter_hxx

#include "itkBinaryWellComposed2DImageFilter.h"
#include "itkImageLinearIteratorWithIndex.h"

namespace itk
{

template<class TImage>
BinaryWellComposed2DImageFilter<TImage>
::BinaryWellComposed2DImageFilter()
{
  this->m_BackgroundValue = NumericTraits<PixelType>::Zero;
  this->m_ForegroundValue = NumericTraits<PixelType>::One;

  this->m_FullInvariance = true;
  this->InitializeIndices();
}

template<class TImage>
BinaryWellComposed2DImageFilter<TImage>
::~BinaryWellComposed2DImageFilter()
{
}

template<class TImage>
void
BinaryWellComposed2DImageFilter<TImage>
::GenerateData()
{
  if ( ImageDimension != 2 )
    {
    itkExceptionMacro( << "Image dimension must be equal to 2." );
    }
  this->AllocateOutputs();  
  this->ConvertBoundaryForegroundPixelsToBackgroundPixels();
  this->MakeImageWellComposed();
}

template<class TImage>
void
BinaryWellComposed2DImageFilter<TImage>
::ConvertBoundaryForegroundPixelsToBackgroundPixels()
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
      if ( It.Get() == this->m_ForegroundValue )
        { 
        It.Set(this->m_BackgroundValue );
        N++;
        }
      It.GoToEndOfLine();
      --It;
      if ( It.Get() == this->m_ForegroundValue )
        { 
        It.Set( this->m_BackgroundValue );
        N++;
        }
      It.NextLine();
      }
    }
  itkDebugMacro( << N << " boundary foreground pixels were"
                 << " converted to background pixels." );
}

template<class TImage>
void
BinaryWellComposed2DImageFilter<TImage>
::MakeImageWellComposed()
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


template<class TImage>
bool
BinaryWellComposed2DImageFilter<TImage>
::IsCriticalC1Configuration( Array<PixelType> neighborhood )
{
  return ( neighborhood[0] == this->m_BackgroundValue &&
           neighborhood[1] == this->m_ForegroundValue &&
           neighborhood[3] == this->m_ForegroundValue &&
           neighborhood[4] == this->m_BackgroundValue &&
           neighborhood[8] == this->m_BackgroundValue );
}

template<class TImage>
bool
BinaryWellComposed2DImageFilter<TImage>
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

template<class TImage>
bool
BinaryWellComposed2DImageFilter<TImage>
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

template<class TImage>
bool
BinaryWellComposed2DImageFilter<TImage>
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

template<class TImage>
bool
BinaryWellComposed2DImageFilter<TImage>
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

template<class TImage>
void
BinaryWellComposed2DImageFilter<TImage>
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

template <class TImage>
void
BinaryWellComposed2DImageFilter<TImage>
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Background Value: " << this->m_BackgroundValue << std::endl;
  os << indent << "Foreground Value: " << this->m_ForegroundValue << std::endl;
  os << indent << "Number of C1 configurations: " 
               << this->m_NumberOfC1Configurations << std::endl;
  os << indent << "Number of C2 configurations: " 
               << this->m_NumberOfC2Configurations << std::endl;
  os << indent << "Number of C3 configurations: " 
               << this->m_NumberOfC3Configurations << std::endl;
  os << indent << "Number of C4 configurations: " 
               << this->m_NumberOfC4Configurations << std::endl;
}

} // end namespace itk

#endif
