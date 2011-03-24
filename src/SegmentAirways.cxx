#include "itkBinaryDiamondStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBoundedSpaceDilateImageFilter.h"
#include "itkBinaryBoundedSpaceErodeImageFilter.h"
#include "itkBinaryReinhardtMorphologicalImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkHoughTransform2DCirclesImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkNumericTraits.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkTimeProbe.h"


int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " inputImageFile labelImageFile outputImageFile" << std::endl;
    exit( 0 );
    }

  /**
   * This routine implements the second step of the lung segmentation algorithm discussed in
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
   *    6. The trachea can be seen in the first axial slice (may need to be changed).
   * Input: Also, an initial segmentation is provided in the form of a label image where the
   *   body has a label of '1', background has a label of '0', and the lungs and main airways
   *   have a label of '2'.

   * Output: one label image with the main airways separated from the background, lung and
   *    body.  The following labeling is given as
   *    1. The body has a label of '1'.
   *    2. The lungs have a label value of '2'.
   *    3. The airways have a label value of '4'.
   *
   * The steps we employ are as follows:
   *    1. Find the trachea in the first slice using a Hough transform.
   *    2. Until the carina is reached, we proceed slice by slice where the result from
   *       the previous slice is intersected with the segmentation of the current slice.
   *       A bounded space dilation is used to grow the region in the current slice.
   *    3. Once two regions result, we know that we have proceeded past the carina.
   *    4. We continue to grow on each slice for each of the two branches until the
   *       area from the region growing indicates that we've entered into the lung.
   *
   * This routine is meant to be used in the pipeline
   *
   * inputImage --> LungExtraction --> SegmentAirways --> SeparateLungs --> initialLabeling
   */

  /**
   * Set up the initial typedefs and read in the images.
   */
  const unsigned int ImageDimension = 3;

  typedef itk::Image<int, ImageDimension> LabelImageType;
  typedef itk::Image<LabelImageType::PixelType, ImageDimension-1> LabelSliceType;

  typedef itk::Image<int, ImageDimension> ImageType;
  typedef itk::Image<ImageType::PixelType, ImageDimension-1> SliceType;

  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Image<RealImageType::PixelType, ImageDimension-1> RealSliceType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
  LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( argv[2] );
  labelReader->Update();

  LabelImageType::Pointer labelImage = labelReader->GetOutput();

  ImageType::Pointer image = ImageType::New();
  image = reader->GetOutput();

  /**
   * Now apply the Hough transform to the first slice to find the trachea.
   */

  SliceType::Pointer slice = SliceType::New();
  ImageType::RegionType region;
  ImageType::RegionType::SizeType size = image->GetLargestPossibleRegion().GetSize();
  size[2] = 0;
  ImageType::IndexType index;
  index.Fill( 0 );
  region.SetSize( size );

  unsigned int seedSliceIndex = image->GetLargestPossibleRegion().GetSize()[2] - 1;

  unsigned int voxelMinRadius 
    = static_cast<unsigned int>( 2.0 / reader->GetOutput()->GetSpacing()[0] ); 
  unsigned int voxelMaxRadius 
    = static_cast<unsigned int>( 15.0 / reader->GetOutput()->GetSpacing()[0] ); 

  for ( int n = seedSliceIndex; n >= 0; n-- )
    {
    index[2] = n;
    region.SetIndex( index );

    typedef itk::ExtractImageFilter<LabelImageType, LabelSliceType> LabelExtracterType;
    LabelExtracterType::Pointer labelextracter = LabelExtracterType::New();
    labelextracter->SetInput( labelImage );
    labelextracter->SetExtractionRegion( region );
    labelextracter->Update();

    unsigned int numberOfAirwayVoxels = 0;
    itk::ImageRegionIteratorWithIndex<LabelSliceType> It( labelextracter->GetOutput(),
      labelextracter->GetOutput()->GetLargestPossibleRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if( It.Get() == 2 )
        {
        numberOfAirwayVoxels++;
        } 
      }

    if( vcl_sqrt( static_cast<float>( numberOfAirwayVoxels ) / vnl_math::pi ) >=
        voxelMinRadius+2 )
      {
      typedef itk::ExtractImageFilter<ImageType, SliceType> ExtracterType;
      ExtracterType::Pointer extracter = ExtracterType::New();
      extracter->SetInput( image );
      extracter->SetExtractionRegion( region );
      extracter->Update();

      slice = extracter->GetOutput();
      seedSliceIndex = index[2];
      break;
      }
      
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if( It.Get() == 2 )
        {
        LabelImageType::IndexType idx;
        idx[0] = It.GetIndex()[0];
        idx[1] = It.GetIndex()[1];
        idx[2] = index[2];
        labelImage->SetPixel( idx, 1 );
        } 
      }
    }

  typedef itk::ExtractImageFilter<SliceType, RealSliceType> SliceExtracterType;
  SliceExtracterType::Pointer sliceExtracter = SliceExtracterType::New();
  SliceType::RegionType portion;
  SliceType::RegionType::SizeType portionSize;
  portionSize[0] = static_cast<unsigned int>( 0.5 * size[0] );
  portionSize[1] = static_cast<unsigned int>( 0.5 * size[1] );
  SliceType::IndexType portionIndex;
  portionIndex[0] = static_cast<unsigned int>( 0.25 * size[0] );
  portionIndex[1] = static_cast<unsigned int>( 0.25 * size[1] );
  portion.SetSize( portionSize );
  portion.SetIndex( portionIndex );
  sliceExtracter->SetInput( slice );
  sliceExtracter->SetExtractionRegion( portion );
  sliceExtracter->Update();

  typedef itk::ExtractImageFilter<LabelImageType, LabelSliceType> LabelExtracterType;
  LabelExtracterType::Pointer labelExtracter = LabelExtracterType::New();
  labelExtracter->SetInput( labelImage );
  labelExtracter->SetExtractionRegion( region );
  labelExtracter->Update();

  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<RealSliceType, RealSliceType> GradientFilterType;
  GradientFilterType::Pointer gradFilter =  GradientFilterType::New();
  gradFilter->SetInput( sliceExtracter->GetOutput() );
  gradFilter->SetSigma( 0.1 );
  gradFilter->Update();


  typedef itk::HoughTransform2DCirclesImageFilter<
    RealType, RealType> HoughFilterType;
  HoughFilterType::Pointer hough = HoughFilterType::New();
  hough->SetInput( gradFilter->GetOutput() );
  hough->SetMinimumRadius( voxelMinRadius );
  hough->SetMaximumRadius( voxelMaxRadius );
  hough->SetThreshold( itk::NumericTraits<SliceType::PixelType>::NonpositiveMin() );
  hough->SetSigmaGradient( 1.0 );
  hough->SetNumberOfCircles( 15 );
  hough->SetSweepAngle( 0.3 );
  hough->SetVariance( 1.0 );
  hough->Update();

  HoughFilterType::CirclesListType circles = hough->GetCircles( hough->GetNumberOfCircles() );
  HoughFilterType::CirclesListType::iterator iter;
  RealType minDistance = itk::NumericTraits<RealType>::max();
  RealType minCenterX = 0.0;
  RealType minCenterY = 0.0;
  RealType minRadius = itk::NumericTraits<RealType>::max();

  RealType centerAnteriorY = 0.0;
  RealType centerAnteriorX = 0.5*
    static_cast<RealType>( image->GetLargestPossibleRegion().GetSize()[0] - 1 );

  for ( iter = circles.begin(); iter != circles.end(); ++iter )
    {
    RealType centerX = (*iter)->GetObjectToParentTransform()->GetOffset()[0];
    RealType centerY = (*iter)->GetObjectToParentTransform()->GetOffset()[1];

    RealType distance = vnl_math_sqr( centerX - centerAnteriorX ) +
      0*vnl_math_sqr( centerY - centerAnteriorY );

    LabelSliceType::IndexType index;
    index[0] = static_cast<unsigned int>( centerX + 0.5 );
    index[1] = static_cast<unsigned int>( centerY + 0.5 );

    bool tracheaFound = false;

    LabelSliceType::PixelType centerValue 
      = labelExtracter->GetOutput()->GetPixel( index );
    if( centerValue == static_cast<LabelSliceType::PixelType>( 2 ) )
      {
      tracheaFound = true; 
      }  
    if( !tracheaFound )
      {
      LabelSliceType::IndexType localIndex; 
      for ( int i = -static_cast<int>( (*iter)->GetRadius()[0] ); 
        i <= static_cast<int>( (*iter)->GetRadius()[0] ); i++ )
        {
        localIndex[0] = static_cast<int>( centerX ) + i;
        for ( int j = -static_cast<int>( (*iter)->GetRadius()[0] ); 
          j <= static_cast<int>( (*iter)->GetRadius()[0] ); j++ )
          {
          localIndex[1] = static_cast<int>( centerY ) + j;
          if( labelExtracter->GetOutput()->GetPixel( localIndex ) == 
            static_cast<LabelSliceType::PixelType>( 2 ) )
            {
            tracheaFound = true;
            break; 
            }
          }
        if( tracheaFound )
          {
          break; 
          }  
        }  
      }

    if ( distance  < minDistance && tracheaFound )
      {
      minDistance = distance;
      minRadius = (*iter)->GetRadius()[0];
      minCenterX = centerX;
      minCenterY = centerY;
      }
    }
  
  if( minRadius > 1e10 )
    {
    std::cerr << "Radius not found." << std::endl; 
    return EXIT_FAILURE;
    }

  LabelSliceType::Pointer seedSlice = LabelSliceType::New();
  seedSlice->SetOrigin( labelExtracter->GetOutput()->GetOrigin() );
  seedSlice->SetSpacing( labelExtracter->GetOutput()->GetSpacing() );
  seedSlice->SetRegions( labelExtracter->GetOutput()->GetLargestPossibleRegion() );
  seedSlice->Allocate();
  seedSlice->FillBuffer( 0 );

  for ( int i = -static_cast<int>( minRadius ); i <= static_cast<int>( minRadius ); i++ )
    {
    index[0] = static_cast<int>( minCenterX ) + i;
    for ( int j = -static_cast<int>( minRadius ); j <= static_cast<int>( minRadius ); j++ )
      {
      index[1] = static_cast<int>( minCenterY ) + j;
      RealType distance = vcl_sqrt( vnl_math_sqr( minCenterX - static_cast<RealType>( index[0] ) ) +
        vnl_math_sqr( minCenterY - static_cast<RealType>( index[1] ) ) );

      if ( distance <= minRadius )
        {
        LabelImageType::PixelType label = labelImage->GetPixel( index );
        LabelSliceType::IndexType idx;
        idx[0] = index[0];
        idx[1] = index[1];
        if ( label == 2 || label == 3 )
          {
          seedSlice->SetPixel( idx, 4 );
          }
        }
      }
    }

  /**
   * After finding the trachea in the first slice we use bounded space dilation (region
   * growing) to fill in the connected region in that first slice.
   */

  typedef itk::BinaryDiamondStructuringElement<
                      LabelImageType::PixelType,
                      ImageDimension-1>             StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius( 1 );
  structuringElement.CreateStructuringElement();

  typedef itk::BinaryBoundedSpaceDilateImageFilter<
    LabelSliceType, LabelSliceType, StructuringElementType>  DilaterType;
  DilaterType::Pointer dilater = DilaterType::New();
  dilater->SetInput( seedSlice );
  dilater->SetScaling( 15 );
  dilater->SetKernel( structuringElement );
  dilater->SetBoundedSpaceImage( labelExtracter->GetOutput() );
  dilater->SetBoundedSpaceValue( 2 );
  dilater->SetForegroundValue( 4 );
  dilater->Update();

  itk::ImageRegionIteratorWithIndex<LabelSliceType> It( seedSlice,
    seedSlice->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    for ( unsigned int d = 0; d < ImageDimension-1; d++ )
      {
      index[d] = It.GetIndex()[d];
      }
//    LabelSliceType::PixelType seedPixel = It.Get();

    if ( dilater->GetOutput()->GetPixel( It.GetIndex() ) == 4 )
      {
      labelImage->SetPixel( index, 4 );
      }
    }

  /**
   * Now we proceed processing each axial slice individually.  The previous
   * axial slice provides a seed for the current slice.  Bounded space dilation
   * is used to fill out the regions.  Also, if the result of the bounded space
   * dilation is an area greater than 3 times the standard deviation of the
   * current average area, then we don't use that result for that slice.
   */

  RealType N = 0;
  RealType var = 0.0;
  RealType meanArea = 0.0;
 
  seedSlice = dilater->GetOutput();

  typedef itk::LabelStatisticsImageFilter<LabelSliceType, LabelSliceType> StatsFilterType;
  StatsFilterType::Pointer stats = StatsFilterType::New();
  stats->SetInput( seedSlice );
  stats->SetLabelInput( seedSlice );
  stats->Update();

  meanArea = stats->GetCount( 4 );
  for ( unsigned int d = 0; d < ImageDimension-1; d++ )
    {
    meanArea *= seedSlice->GetSpacing()[d];
    }
  N++;

  size[2] = 0;
  index.Fill( 0 );
  region.SetSize( size );

  for ( int s = seedSliceIndex-1; s >= 0; s-- )
    {

    index.Fill( 0 );
    index[2] = s;
    region.SetIndex( index );

    typedef itk::ExtractImageFilter<LabelImageType, LabelSliceType> LabelExtracterType;
    LabelExtracterType::Pointer labelExtracter = LabelExtracterType::New();
    labelExtracter->SetInput( labelImage );
    labelExtracter->SetExtractionRegion( region );
    labelExtracter->Update();

    itk::ImageRegionIterator<LabelSliceType> ItL( labelExtracter->GetOutput(),
      labelExtracter->GetOutput()->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<LabelSliceType> ItS( seedSlice,
      seedSlice->GetLargestPossibleRegion() );
    for ( ItL.GoToBegin(), ItS.GoToBegin(); !ItL.IsAtEnd(); ++ItL, ++ItS )
      {
      if ( ItL.Get() != 2 || ItS.Get() != 4 )
        {
        ItS.Set( 0 );
        }
      }

    typedef itk::ConnectedComponentImageFilter<LabelSliceType, LabelSliceType> ConnectedComponentType;
    ConnectedComponentType::Pointer connecter = ConnectedComponentType::New();
    connecter->SetInput( seedSlice );
    connecter->FullyConnectedOn();
    connecter->Update();

    typedef itk::RelabelComponentImageFilter<LabelSliceType, LabelSliceType> LabelerType;
    LabelerType::Pointer labeler = LabelerType::New();
    labeler->SetInput( connecter->GetOutput() );
    labeler->Update();

    for ( unsigned int n = 1; n <= labeler->GetNumberOfObjects(); n++ )
      {
      typedef itk::BinaryDiamondStructuringElement<
                          LabelImageType::PixelType,
                          ImageDimension-1>             StructuringElementType;
      StructuringElementType structuringElement;
      structuringElement.SetRadius( 1 );
      structuringElement.CreateStructuringElement();

      typedef itk::BinaryBoundedSpaceDilateImageFilter<
        LabelSliceType, LabelSliceType, StructuringElementType>  DilaterType;
      DilaterType::Pointer dilater = DilaterType::New();
      dilater->SetInput( labeler->GetOutput() );
      dilater->SetScaling( 100 );
      dilater->SetKernel( structuringElement );
      dilater->SetBoundedSpaceImage( labelExtracter->GetOutput() );
      dilater->SetBoundedSpaceValue( 2 );
      dilater->SetForegroundValue( n );
      dilater->SetBackgroundValue( 0 );
      dilater->Update();

      typedef itk::LabelStatisticsImageFilter<LabelSliceType, LabelSliceType> StatsFilterType;
      StatsFilterType::Pointer stats = StatsFilterType::New();
      stats->SetInput( dilater->GetOutput() );
      stats->SetLabelInput( dilater->GetOutput() );
      stats->Update();

      RealType newArea = stats->GetCount( n );
      for ( unsigned int d = 0; d < ImageDimension-1; d++ )
        {
        newArea *= dilater->GetOutput()->GetSpacing()[d];
        }
      if ( newArea < meanArea + 3 * vcl_sqrt( var ) || N < 10 )
        {

        N++;
        meanArea = meanArea *( N - 1.0 )/N + newArea/N;
        if ( N > 1.0 )
          {
          var = var*( N - 1.0 )/N + meanArea*meanArea/( N - 1.0 );
          }

        itk::ImageRegionIteratorWithIndex<LabelSliceType> It( dilater->GetOutput(),
          dilater->GetOutput()->GetLargestPossibleRegion() );
        for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
          {
          LabelImageType::IndexType index;
          for ( unsigned int d = 0; d < ImageDimension-1; d++ )
            {
            index[d] = It.GetIndex()[d];
            }
          index[ImageDimension-1] = s;

          if ( It.Get() == static_cast<int>( n ) )
            {
            labelImage->SetPixel( index, 4 );
            It.Set( 4 );
            }
          }
        }
      }

    if ( s >= 0 )
      {
      index.Fill( 0 );
      index[2] = s;
      region.SetIndex( index );

      typedef itk::ExtractImageFilter<LabelImageType, LabelSliceType> LabelExtracterType;
      LabelExtracterType::Pointer labelExtracter = LabelExtracterType::New();
      labelExtracter->SetInput( labelImage );
      labelExtracter->SetExtractionRegion( region );
      labelExtracter->Update();

      seedSlice = labelExtracter->GetOutput();
      }
    }

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( labelImage );
  writer->Update();

};
