#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryReinhardtMorphologicalImageFilter.h"
#include "itkBinaryThresholdImageFilter.h" 
#include "itkConnectedComponentImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRelabelComponentImageFilter.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " labelImage outputImage [Rf] [Rsa] [Rba] [v]" << std::endl;
    exit( 1 );
    }

  /**
   * This routine implements the optional smoothing step of the lung segmentation algorithm discussed in
   * Hu, et al., "Automatic Lung Segmentation for Accurate Quantitation of Volumetric
   * X-Ray CT Images", IEEE-TMI 20(6):490-498, 2001.  
   *
   * Input:  the labeled segmentation results from the previous three steps
   *
   * Output: one *smoothed* label image with the lung and main airways separated from the background
   *    of the image and the body.  The following labeling is given as
   *    1. The body has a label of '1'.
   *    2. The lungs have a label of '2' (left) and '3' (right).
   *    3. The airways have a label of '4'.
   *
   * This routine is meant to be used in the pipeline
   *
   * inputImage --> LungExtraction --> SegmentAirways --> SeparateLungs --> SmoothLungs ---> initialLabeling    
   *
   */

  typedef float RealType;
  typedef int PixelType;
  const unsigned int ImageDimension = 3;
  
  typedef itk::Image<unsigned int, 3> LabelImageType;
  typedef itk::Image<unsigned int, 2> LabelSliceType;

  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
  LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( argv[1] );
  labelReader->Update();

  typedef itk::ImageDuplicator<LabelImageType> DuplicatorType;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( labelReader->GetOutput() );
  duplicator->Update();

  LabelImageType::Pointer output = duplicator->GetOutput();
  LabelImageType::Pointer labelImage = LabelImageType::New();
  labelImage = labelReader->GetOutput();


  /**
   * Designate the parameters which influence the amount of smoothing.
   */

  unsigned int Rf = 2;
  if ( argc > 3 )
    {
    Rf = static_cast<unsigned int>( atoi( argv[3] ) );
    } 
  unsigned int Rba = 8;
  if ( argc > 4 )
    {
    Rba = static_cast<unsigned int>( atoi( argv[4] ) );
    } 
  unsigned int Rsa = 3;
  if ( argc > 5 )
    {
    Rsa = static_cast<unsigned int>( atoi( argv[5] ) );
    } 
  unsigned int v = 300;
  if ( argc > 6 )
    {
    v = static_cast<unsigned int>( atoi( argv[6] ) );
    } 


  /**
   * Do the smoothing slice by slice
   */

  LabelImageType::RegionType region; 
  LabelImageType::SizeType size = labelImage->GetLargestPossibleRegion().GetSize();
  size[2] = 0; 
  region.SetSize( size );

  LabelSliceType::Pointer labelSlice = LabelSliceType::New();

  for ( int s = labelImage->GetLargestPossibleRegion().GetSize()[2] - 1; 
          s >= 0; s-- )
    {
    LabelImageType::IndexType index;
    index.Fill( 0 );
    index[2] = s;
    region.SetIndex( index );
    
    typedef itk::ExtractImageFilter<LabelImageType, LabelSliceType> LabelExtracterType;
    LabelExtracterType::Pointer labelExtracter = LabelExtracterType::New();
    labelExtracter->SetInput( labelImage );
    labelExtracter->SetExtractionRegion( region );
    labelExtracter->Update();

    /**
     * Remove pulmonary vessels
     */
    if ( Rf > 0 )
      {
      for ( unsigned int lung = 2; lung <= 3; lung++ )
        {
        typedef itk::BinaryBallStructuringElement<LabelSliceType::PixelType,
          ImageDimension-1> BallStructuringElementType;
        BallStructuringElementType ballStructuringElement;
        ballStructuringElement.SetRadius( Rf ); 
        ballStructuringElement.CreateStructuringElement();
        
        typedef itk::BinaryMorphologicalClosingImageFilter<LabelSliceType, LabelSliceType, 
          BallStructuringElementType> CloserType;
        CloserType::Pointer closer = CloserType::New();
        closer->SetInput( labelExtracter->GetOutput() );
        closer->SetKernel( ballStructuringElement );
        closer->SetForegroundValue( static_cast<LabelImageType::PixelType>( lung ) );
        closer->Update();
  
        itk::ImageRegionIteratorWithIndex<LabelSliceType> It( closer->GetOutput(), 
          closer->GetOutput()->GetLargestPossibleRegion() );
        for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
          { 
          if ( It.Get() == static_cast<LabelImageType::PixelType>( lung ) )
            {
            LabelImageType::IndexType idx;
            for ( unsigned int d = 0; d < ImageDimension-1; d++ )
              {
              idx[d] = It.GetIndex()[d];
              }
            idx[ImageDimension-1] = s;
            output->SetPixel( idx, static_cast<LabelImageType::PixelType>( lung ) );
            }
          }    
        }
      }  

    /**
     * Remove small airways
     */

    if ( Rsa > 0 )
      {
      for ( unsigned int lung = 2; lung <= 3; lung++ )
        {
        typedef itk::BinaryThresholdImageFilter<LabelSliceType, LabelSliceType> ThresholderType;
        ThresholderType::Pointer thresholder = ThresholderType::New();
        thresholder->SetInput( labelExtracter->GetOutput() );
        thresholder->SetOutsideValue( itk::NumericTraits<LabelSliceType::PixelType>::Zero );
        thresholder->SetInsideValue( itk::NumericTraits<LabelSliceType::PixelType>::One );
        thresholder->SetLowerThreshold( static_cast<LabelImageType::PixelType>( lung ) );
        thresholder->SetUpperThreshold( static_cast<LabelImageType::PixelType>( lung ) );
        thresholder->Update();

        typedef itk::BinaryBallStructuringElement<LabelSliceType::PixelType,
          ImageDimension-1> BallStructuringElementType;
        BallStructuringElementType ballStructuringElement;
        ballStructuringElement.SetRadius( Rsa ); 
        ballStructuringElement.CreateStructuringElement();
        
        typedef itk::BinaryErodeImageFilter<LabelSliceType, LabelSliceType, 
          BallStructuringElementType> EroderType;
        EroderType::Pointer eroder = EroderType::New();
        eroder->SetInput( thresholder->GetOutput() );
        eroder->SetKernel( ballStructuringElement );
        eroder->SetForegroundValue( itk::NumericTraits<LabelImageType::PixelType>::One );
        eroder->SetBackgroundValue( itk::NumericTraits<LabelImageType::PixelType>::Zero );
        eroder->Update();
  
        typedef itk::ConnectedComponentImageFilter<LabelSliceType, LabelSliceType> ConnectedComponentType;
        ConnectedComponentType::Pointer connecter = ConnectedComponentType::New();
        connecter->SetInput( eroder->GetOutput() );
        connecter->FullyConnectedOff();
        connecter->Update();
  
        typedef itk::RelabelComponentImageFilter<LabelSliceType, LabelSliceType> LabelerType;
        LabelerType::Pointer labeler = LabelerType::New();
        labeler->SetInput( connecter->GetOutput() );
        labeler->Update();
  
        typedef itk::BinaryErodeImageFilter<LabelSliceType, LabelSliceType, 
          BallStructuringElementType> DilaterType;
        DilaterType::Pointer dilater = DilaterType::New();
        dilater->SetInput( labeler->GetOutput() );
        dilater->SetKernel( ballStructuringElement );
        dilater->SetForegroundValue( itk::NumericTraits<LabelImageType::PixelType>::One );
        dilater->SetBackgroundValue( itk::NumericTraits<LabelImageType::PixelType>::Zero );
        dilater->Update();
  
        itk::ImageRegionIteratorWithIndex<LabelSliceType> It( dilater->GetOutput(), 
          dilater->GetOutput()->GetLargestPossibleRegion() );
        for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
          { 
          if ( It.Get() == 1 )
            {
            LabelImageType::IndexType idx;
            for ( unsigned int d = 0; d < ImageDimension-1; d++ )
              {
              idx[d] = It.GetIndex()[d];
              }
            idx[ImageDimension-1] = s;
            output->SetPixel( idx, static_cast<LabelImageType::PixelType>( lung ) );
            }
          }    
        }
      }
   
    /**
     * Remove large airways
     */
    if ( Rba > 0 )
      {
      for ( unsigned int lung = 2; lung <= 3; lung++ )
        {
        typedef itk::BinaryBallStructuringElement<LabelSliceType::PixelType,
          ImageDimension-1> BallStructuringElementType;
        BallStructuringElementType ballStructuringElement;
        ballStructuringElement.SetRadius( Rba ); 
        ballStructuringElement.CreateStructuringElement();
        
        typedef itk::BinaryErodeImageFilter<LabelSliceType, LabelSliceType, 
          BallStructuringElementType> EroderType;
        EroderType::Pointer eroder = EroderType::New();
        eroder->SetInput( labelExtracter->GetOutput() );
        eroder->SetKernel( ballStructuringElement );
        eroder->SetForegroundValue( static_cast<LabelImageType::PixelType>( lung ) );
        eroder->SetBackgroundValue( itk::NumericTraits<LabelImageType::PixelType>::Zero );
        eroder->Update();
  
        typedef itk::ConnectedComponentImageFilter<LabelSliceType, LabelSliceType> ConnectedComponentType;
        ConnectedComponentType::Pointer connecter = ConnectedComponentType::New();
        connecter->SetInput( eroder->GetOutput() );
        connecter->FullyConnectedOff();
        connecter->Update();
  
        typedef itk::RelabelComponentImageFilter<LabelSliceType, LabelSliceType> LabelerType;
        LabelerType::Pointer labeler = LabelerType::New();
        labeler->SetInput( connecter->GetOutput() );
        labeler->Update();
  
        typedef itk::BinaryErodeImageFilter<LabelSliceType, LabelSliceType, 
          BallStructuringElementType> DilaterType;
        DilaterType::Pointer dilater = DilaterType::New();
        dilater->SetInput( labeler->GetOutput() );
        dilater->SetKernel( ballStructuringElement );
        dilater->SetForegroundValue( itk::NumericTraits<LabelImageType::PixelType>::One );
        dilater->SetBackgroundValue( itk::NumericTraits<LabelImageType::PixelType>::Zero );
        dilater->Update();
  
        itk::ImageRegionIterator<LabelSliceType> It( dilater->GetOutput(), 
          dilater->GetOutput()->GetLargestPossibleRegion() );
        itk::ImageRegionIterator<LabelSliceType> ItL( labelExtracter->GetOutput(), 
          labelExtracter->GetOutput()->GetLargestPossibleRegion() );
        for ( It.GoToBegin(), ItL.GoToBegin(); !It.IsAtEnd(); ++It, ++ItL )
          { 
          if ( It.Get() != 1 || 
               ItL.Get() == static_cast<LabelImageType::PixelType>( lung ) )
            {
            It.Set( 0 );
            }
          }    
  
        ConnectedComponentType::Pointer connecter2 = ConnectedComponentType::New();
        connecter2->SetInput( dilater->GetOutput() );
        connecter2->FullyConnectedOff();
        connecter2->Update();
  
        LabelerType::Pointer labeler2 = LabelerType::New();
        labeler2->SetInput( connecter2->GetOutput() );
        labeler2->SetMinimumObjectSize( v );
        labeler2->Update();
  
        itk::ImageRegionIteratorWithIndex<LabelSliceType> ItL2( labeler2->GetOutput(), 
          labeler2->GetOutput()->GetLargestPossibleRegion() );
        for ( ItL.GoToBegin(), ItL2.GoToBegin(); !ItL.IsAtEnd(); ++ItL, ++ItL2 )
          { 
          if ( ItL.Get() == static_cast<LabelImageType::PixelType>( lung ) && 
                 ItL2.Get() == 0 )
            {
            LabelImageType::IndexType idx;
            for ( unsigned int d = 0; d < ImageDimension-1; d++ )
              {
              idx[d] = ItL2.GetIndex()[d];
              }
            idx[ImageDimension-1] = s;
            output->SetPixel( idx, static_cast<LabelImageType::PixelType>( lung ) );
            }
          }    
        }
      }  

    /**
     * Additional corrections
     */
    typedef itk::ConstantPadImageFilter<LabelSliceType, LabelSliceType> PadderType;
    PadderType::Pointer padder = PadderType::New();
    padder->SetInput( labelExtracter->GetOutput() );
    padder->SetConstant( itk::NumericTraits<LabelImageType::PixelType>::Zero );
    unsigned long upperFactors[ImageDimension-1];
    unsigned long lowerFactors[ImageDimension-1];
    for ( unsigned int d = 0; d < ImageDimension-1; d++ )
      {
      upperFactors[d] = 1;
      lowerFactors[d] = 1;
      } 
    padder->SetPadLowerBound( lowerFactors );
    padder->SetPadUpperBound( upperFactors );
    padder->Update();
 
    typedef itk::BinaryBallStructuringElement< 
                        LabelSliceType::PixelType,
                        ImageDimension-1>             StructuringElementType;
  
    typedef itk::BinaryReinhardtMorphologicalImageFilter<
      LabelSliceType, LabelSliceType, StructuringElementType>  FilterType;
    FilterType::Pointer reinhardt = FilterType::New();
    reinhardt->SetInput( padder->GetOutput() );
    reinhardt->SetForegroundValue( static_cast<LabelImageType::PixelType>( 1 ) );
    reinhardt->SetEmploySaltAndPepperRepair( true );
    reinhardt->SetSaltAndPepperMinimumSizeInPixels( 
      static_cast<unsigned int>( 10.0 / padder->GetOutput()->GetSpacing()[0] ) );
    reinhardt->SetEmployMinimumDiameterFilter( false );
    reinhardt->SetEmployUnwantedCavityDeletion( false );
    reinhardt->SetEmployMinimumSizeFilter( false );
    reinhardt->SetEmployMaximumDiameterFilter( false );
    reinhardt->SetEmployConnectivityFilter( false );
    reinhardt->SetEmployBoundarySmoother( false );
    reinhardt->SetEmployUnclassifiedPixelProcessing( false );
    reinhardt->Update();

    itk::ImageRegionIteratorWithIndex<LabelSliceType> It( labelExtracter->GetOutput(), 
      labelExtracter->GetOutput()->GetLargestPossibleRegion() );
    for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
      { 
      LabelSliceType::PixelType rlabel = reinhardt->GetOutput()->GetPixel( It.GetIndex() );
      if ( It.Get() != static_cast<LabelImageType::PixelType>( 1 ) &&
             rlabel == static_cast<LabelImageType::PixelType>( 1 ) )
        {
        LabelImageType::IndexType idx;
        for ( unsigned int d = 0; d < ImageDimension-1; d++ )
          {
          idx[d] = It.GetIndex()[d];
          }
        idx[ImageDimension-1] = s;
        output->SetPixel( idx, static_cast<LabelImageType::PixelType>( 1 ) );
        }
      }  

    }

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( output );
  writer->Update();

  return 0;
}
