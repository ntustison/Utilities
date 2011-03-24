#ifndef __itkBinaryWellComposed3DImageFilter_txx
#define __itkBinaryWellComposed3DImageFilter_txx

#include "itkBinaryWellComposed3DImageFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkImageLinearIteratorWithIndex.h"

namespace itk
{

template<class TImage>
BinaryWellComposed3DImageFilter<TImage>
::BinaryWellComposed3DImageFilter()
{
  this->m_BackgroundValue = NumericTraits<PixelType>::Zero;
  this->m_ForegroundValue = NumericTraits<PixelType>::One;

  this->InitializeOffsetsAndIndices();
}

template<class TImage>
BinaryWellComposed3DImageFilter<TImage>
::~BinaryWellComposed3DImageFilter() {}

template<class TImage>
void
BinaryWellComposed3DImageFilter<TImage>
::GenerateData()
{
  if ( ImageDimension != 3 )
    {
    itkExceptionMacro( << "Image dimension must be equal to 3." );
    }

  this->AllocateOutputs();
  this->ConvertBoundaryForegroundPixelsToBackgroundPixels();
  this->CountCriticalConfigurations();
   
  if ( this->m_NumberOfC1Configurations > 0 
       || this->m_NumberOfC2Configurations > 0 )
    { 
    this->LocateCriticalConfigurations();
    this->MakeImageWellComposed();
    }
  this->CountCriticalConfigurations();
}

template<class TImage>
void
BinaryWellComposed3DImageFilter<TImage>
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

template<class TImage>
void
BinaryWellComposed3DImageFilter<TImage>
::MakeImageWellComposed()
{
  Array<unsigned char> p(8);
  unsigned int iter = 0;
 
  itkDebugMacro( << "Iteration " << iter++ 
                 << " (Current number of critical configurations = " 
                 << this->m_CriticalConfigurationIndices.size() << ")" );

  this->m_NewCriticalConfigurationIndices.clear();
  this->m_NumberOfC1Configurations = this->m_NumberOfC2Configurations = 0;

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

template<class TImage>
bool
BinaryWellComposed3DImageFilter<TImage>
::IsChangeSafe(IndexType idx)
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
    if ( this->IsCriticalC1Configuration( p ) )
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
    if ( this->IsCriticalC2Configuration( p ) )
      {
      return false;
      }  
    } 
  return true;
}

template<class TImage>
void
BinaryWellComposed3DImageFilter<TImage>
::MakeRandomChange(unsigned char &p0, IndexType idx0, 
                   unsigned char &p1, IndexType idx1)
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
BinaryWellComposed3DImageFilter<TImage>
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


template<class TImage>
void
BinaryWellComposed3DImageFilter<TImage>
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
    if ( this->IsCriticalC1Configuration( p ) )
      {
      this->m_NewCriticalConfigurationIndices.push_back( 
        idx + this->m_Offsets27[this->m_C1IndicesII[i][0]] );
      this->m_NumberOfC1Configurations++;
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
    if (this->IsCriticalC2Configuration(p))
      {
      this->m_NewCriticalConfigurationIndices.push_back( 
        idx + this->m_Offsets27[this->m_C2IndicesII[i][0]] );
      this->m_NumberOfC2Configurations++;
      }  
    } 
}

template<class TImage>
bool
BinaryWellComposed3DImageFilter<TImage>
::IsCriticalC1Configuration(Array<unsigned char> p)
{
  return (p[0] == p[2]) && (p[1] == p[3]) && (p[0] != p[1]);
}

template<class TImage>
bool
BinaryWellComposed3DImageFilter<TImage>
::IsCriticalC2Configuration(Array<unsigned char> p)
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

template<class TImage>
void
BinaryWellComposed3DImageFilter<TImage>
::CountCriticalConfigurations()
{
  this->m_NumberOfC1Configurations = 0;
  this->m_NumberOfC2Configurations = 0;

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
          this->m_NumberOfC1Configurations++; 
          }
        }
      }  
    
    // Find C2 configurations
    for ( unsigned int i = 0; i < 4; i++ )
      { 
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
          this->m_NumberOfC2Configurations++;
          break;
          }
        }
      }
    }
  itkDebugMacro( << "Initial number of C1 configurations = "
                 << this->m_NumberOfC1Configurations );
  itkDebugMacro( << "Initial number of C2 configurations = " 
                 << this->m_NumberOfC2Configurations );
}

template<class TImage>
void
BinaryWellComposed3DImageFilter<TImage>
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
BinaryWellComposed3DImageFilter<TImage>
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

/**
 * Standard "PrintSelf" method
 */
template <class TImage>
void
BinaryWellComposed3DImageFilter<TImage>
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
