#include "itkBinaryReinhardtMorphologicalImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkPadImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include <fstream.h>

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputImageFile outputImageFile" << std::endl;     
    exit( 0 );
    }   

  /**
   * This routine implements the first step of the lung segmentation algorithm discussed in
   * Hu, et al., "Automatic Lung Segmentation for Accurate Quantitation of Volumetric
   * X-Ray CT Images", IEEE-TMI 20(6):490-498, 2001.  
   *
   * Input:  one CT image of the lung.  The assumptions on the input are as follows:
   *    1. Background has the largest volume
   *    2. The image is read as InputPixelType int and has dimension 3.
   *    3. The sagittal, coronal, and axial directions corresponds with image dimension 1, 2, and 3
   *        respectively.
   *    4. The start index is [0, 0, 0].
   *    5. Superior slices have higher index values than inferior slices.
   *
   * Output: one label image with the lung and main airways separated from the background
   *    of the image and the body.  The following labeling is given as
   *    1. The body has a label of '1'.
   *    2. The lungs and airways have a label value of '2'.
   *
   * The steps we employ are as follows:
   *    1. An optimal threshold value is obtained on this anisotropic diffusion image
   *       using an Otsu threshold filter (which takes half the time as the iterative
   *       procedure of Hu et al. 
   *  
   * This routine is meant to be used in the pipeline
   *
   * inputImage --> LungExtraction --> SegmentAirways --> SeparateLungs --> initialLabeling    
   *
   */

  typedef int PixelType;
  const unsigned int ImageDimension = 3;
  typedef float RealType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<int, ImageDimension> LabelImageType;
  typedef itk::Image<int, ImageDimension-1> LabelSliceType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  /**
   * Threshold selection
   */

  typedef itk::OtsuThresholdImageFilter<ImageType, LabelImageType> FilterType;
  FilterType::Pointer binaryFilter = FilterType::New();
  binaryFilter->SetInput( reader->GetOutput() );
  binaryFilter->SetInsideValue( 0 );
  binaryFilter->SetOutsideValue( 1 );
  binaryFilter->SetNumberOfHistogramBins( 100 );
  binaryFilter->Update();

  typedef itk::ConnectedComponentImageFilter<LabelImageType, LabelImageType> ConnectedComponentType;
  ConnectedComponentType::Pointer connecter = ConnectedComponentType::New();
  connecter->SetInput( binaryFilter->GetOutput() );
  connecter->FullyConnectedOff();
  connecter->Update();

  typedef itk::RelabelComponentImageFilter<LabelImageType, LabelImageType> RelabelerType;
  RelabelerType::Pointer relabeler = RelabelerType::New();
  relabeler->SetInput( connecter->GetOutput() );
  relabeler->InPlaceOff();
  relabeler->Update();

  typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType> ThresholderType;
  ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( relabeler->GetOutput() );
  thresholder->SetInsideValue( 0 );
  thresholder->SetOutsideValue( 1 );
  thresholder->SetLowerThreshold( 1 );
  thresholder->SetUpperThreshold( 1 );
  thresholder->Update();

  ConnectedComponentType::Pointer connecter3 = ConnectedComponentType::New();
  connecter3->SetInput( thresholder->GetOutput() );
  connecter3->FullyConnectedOff();
  connecter3->Update();

  RelabelerType::Pointer relabeler3 = RelabelerType::New();
  relabeler3->SetInput( connecter3->GetOutput() );
  relabeler3->InPlaceOff();
  relabeler3->Update();

  /** 
   * At this point, given the assumption that the background has the largest volume, 
   * the background has a label value of 1
   * and the body has a value of 0.  We want to switch this labeleing so that
   * the body has a value of 1 and the background has a value of 0.  We also 
   * get rid of the labels that have small corresponding volumes. After this step
   * we assume the following labeling.
   *   0 -> Background
   *   1 -> Body
   *   2 -> both lungs
   */

  bool needToSeparateLungs = false;
  if ( relabeler3->GetSizeOfObjectInPhysicalUnits( 3 ) < 
        0.75 * relabeler3->GetSizeOfObjectInPhysicalUnits( 2 ) )
    {
    needToSeparateLungs = true;
    }   

  itk::ImageRegionIterator<LabelImageType> It( relabeler3->GetOutput(), 
    relabeler3->GetOutput()->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( It.Get() == 0 )
      {
      It.Set( 1 );
      } 
    else if ( It.Get() == 1 || It.Get() > 3 )
      {
      It.Set( 0 );
      } 
    else if ( It.Get() > 2 && needToSeparateLungs )
      {
      It.Set( 1 );
      } 
    }  

  /** 
   * Because of the inversion, there might be some spurious pixels labeled as
   * Body pixels.  We apply the connected components again to get rid of these
   * spurious labels 
   */

  ConnectedComponentType::Pointer connecter2 = ConnectedComponentType::New();
  connecter2->SetInput( relabeler3->GetOutput() );
  connecter2->FullyConnectedOff();
  connecter2->Update();

  RelabelerType::Pointer relabeler2 = RelabelerType::New();
  relabeler2->SetInput( connecter2->GetOutput() );
  relabeler2->InPlaceOff();
  relabeler2->Update();

  itk::ImageRegionIterator<LabelImageType> It2( relabeler2->GetOutput(), 
    relabeler2->GetOutput()->GetLargestPossibleRegion() );
  for ( It.GoToBegin(), It2.GoToBegin(); !It2.IsAtEnd(); ++It, ++It2 )
    {
    if ( It2.Get() > 1 )
      {
      It2.Set( 0 );
      } 
    else if ( It.Get() == 2 || ( It.Get() > 2 && needToSeparateLungs ) )
      {
      It2.Set( It.Get() );
      }
    }  

  /**
   * Fill the unwanted cavities in the lungs and body
   */
  LabelImageType::RegionType region;
  LabelImageType::SizeType size = relabeler2->GetOutput()->GetLargestPossibleRegion().GetSize();
  LabelImageType::IndexType index;
  index.Fill( 0 );
  size[2] = 0;
  region.SetSize( size );

  for ( int s = relabeler2->GetOutput()->GetLargestPossibleRegion().GetSize()[2] - 1; 
          s >= 0; s-- )
    {
    index[2] = s;
    region.SetIndex( index );

    typedef itk::ExtractImageFilter<LabelImageType, LabelSliceType> LabelExtracterType;
    LabelExtracterType::Pointer labelExtracter = LabelExtracterType::New();
    labelExtracter->SetInput( relabeler2->GetOutput() );
    labelExtracter->SetExtractionRegion( region );
    labelExtracter->Update();

    /**
     * Fill in the unwanted cavities and repair the salt and pepper noise in the lungs
     */
      { 
      typedef itk::BinaryBallStructuringElement< 
                          LabelImageType::PixelType,
                          ImageDimension>             StructuringElementType;
     
      typedef itk::BinaryReinhardtMorphologicalImageFilter<
        LabelSliceType, LabelSliceType, StructuringElementType>  FilterType;
  
      FilterType::Pointer filter = FilterType::New();
      filter->SetInput( labelExtracter->GetOutput() );
  
      filter->SetForegroundValue( 2 );

      filter->SetEmploySaltAndPepperRepair( true );
      filter->SetEmployMinimumDiameterFilter( false );
      filter->SetEmployUnwantedCavityDeletion( true );
      filter->SetEmployMinimumSizeFilter( false );
      filter->SetEmployMaximumDiameterFilter( false );
      filter->SetEmployConnectivityFilter( false );
      filter->SetEmployBoundarySmoother( false );
      filter->SetEmployUnclassifiedPixelProcessing( false );
      filter->Update();
      
      itk::ImageRegionIteratorWithIndex<LabelSliceType> It( filter->GetOutput(),
        filter->GetOutput()->GetLargestPossibleRegion() );
      for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
        {
        LabelImageType::IndexType index;
        for ( unsigned int d = 0; d < ImageDimension-1; d++ )
          {
          index[d] = It.GetIndex()[d];
          }
        index[ImageDimension-1] = s;
        if ( It.Get() == filter->GetForegroundValue() ) 
          {  
          relabeler2->GetOutput()->SetPixel( index, filter->GetForegroundValue() );
          }
        }   
      }   
    }


  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( relabeler2->GetOutput() );
  writer->Update();

}
