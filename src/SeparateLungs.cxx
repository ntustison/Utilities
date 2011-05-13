#include "itkAddImageFilter.h"
#include "itkBinaryBoundedSpaceDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkGrayscaleFillholeImageFilter.h"
#include "itkLiveWireImageFunction.h"
#include "itkPathIterator.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkStatisticsImageFilter.h"

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage labelImage outputImage" << std::endl;
    exit( 1 );
    }

  /**
   * This routine implements the third and final step of the lung segmentation algorithm discussed in
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
   *    6. The user must also provide the initial estimates for the body and lung values.
   *
   * Output: one label image with the lung and main airways separated from the background
   *    of the image and the body.  The following labeling is given as
   *    1. The body has a label of '1'.
   *    2. The lungs have a label of '2' (left) and '3' (right).
   *    3. The airways have a label of '4'.
   *
   * The steps we employ are as follows:
   *    1. Proceed slice by slice (axial).  Using connected components erode the
   *       objects in the current slice until there are 0 or 2 objects.
   *    2. If there are two objects *after erosion* we then assume that those two
   *       objects are the left and right lungs.  We then use bounded space dilation
   *       to grow the two regions until they overlap.  Any regions that overlap are
   *       then severed using livewire segmentation.
   *    3. Clean-up then uses the location of the airways as a center divider to take
   *       care of problematic regions.
   *
   * This routine is meant to be used in the pipeline
   *
   * inputImage --> LungExtraction --> SegmentAirways --> SeparateLungs --> initialLabeling
   *
   */

  typedef float RealType;
  typedef int PixelType;
  const unsigned int ImageDimension = 3;

  const unsigned lungLabel = 2;

  typedef itk::Image<PixelType, 3> ImageType;
  typedef itk::Image<unsigned int, 3> LabelImageType;

  typedef itk::Image<PixelType, 2> SliceType;
  typedef itk::Image<unsigned int, 2> LabelSliceType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
  LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( argv[2] );
  labelReader->Update();

  typedef itk::ImageDuplicator<LabelImageType> DuplicatorType;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( labelReader->GetOutput() );
  duplicator->Update();

  LabelImageType::Pointer output = duplicator->GetOutput();

  LabelImageType::Pointer labelImage = LabelImageType::New();
  labelImage = labelReader->GetOutput();

  RealType centerAirwayX = 0.0;
  RealType N = 0.0;

  itk::ImageRegionIteratorWithIndex<LabelImageType> It( labelImage,
    labelImage->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( It.Get() == 4 )
      {
      centerAirwayX += static_cast<RealType>( It.GetIndex()[0] );
      N++;
      }
    if ( It.Get() != lungLabel )
      {
      It.Set( 0 );
      }
    }
  centerAirwayX /= N;


  LabelImageType::RegionType region;
  LabelImageType::SizeType size = labelImage->GetLargestPossibleRegion().GetSize();
  size[2] = 0;
  region.SetSize( size );

  LabelSliceType::Pointer labelSlice = LabelSliceType::New();

  for ( int s = labelImage->GetLargestPossibleRegion().GetSize()[2] - 1;
          s >= 0; s-- )
    {
    ImageType::IndexType index;
    index.Fill( 0 );
    index[2] = s;
    region.SetIndex( index );

    typedef itk::ExtractImageFilter<LabelImageType, LabelSliceType> LabelExtracterType;
    LabelExtracterType::Pointer labelExtracter = LabelExtracterType::New();
    labelExtracter->SetInput( labelImage );
    labelExtracter->SetExtractionRegion( region );
    labelExtracter->SetDirectionCollapseToIdentity();
    labelExtracter->Update();

    typedef itk::GrayscaleFillholeImageFilter<LabelSliceType, LabelSliceType> HoleFillerType;
    HoleFillerType::Pointer holeFiller = HoleFillerType::New();
    holeFiller->SetInput( labelExtracter->GetOutput() );
    holeFiller->FullyConnectedOff();
    holeFiller->Update();

    typedef itk::ConnectedComponentImageFilter<LabelSliceType, LabelSliceType> ConnectedComponentType;
    ConnectedComponentType::Pointer connecter = ConnectedComponentType::New();
    connecter->SetInput( holeFiller->GetOutput() );
    connecter->FullyConnectedOff();
    connecter->Update();

    typedef itk::RelabelComponentImageFilter<LabelSliceType, LabelSliceType> LabelerType;
    LabelerType::Pointer labeler = LabelerType::New();
    labeler->SetInput( connecter->GetOutput() );
    labeler->SetMinimumObjectSize( 1000 );
    labeler->Update();

    unsigned int numberOfObjects = labeler->GetNumberOfObjects();

    labelSlice = labeler->GetOutput();

    unsigned int radius = 1;
    bool eroded = false;
    while ( numberOfObjects != 2 && numberOfObjects != 0 )
      {
      typedef itk::BinaryBallStructuringElement<
                          LabelSliceType::PixelType,
                          ImageDimension-1>             StructuringElementType;
      StructuringElementType structuringElement;
      structuringElement.SetRadius( radius++ );
      structuringElement.CreateStructuringElement();

      typedef itk::BinaryErodeImageFilter<
        LabelSliceType, LabelSliceType, StructuringElementType>  EroderType;
      EroderType::Pointer eroder = EroderType::New();
      eroder->SetInput( holeFiller->GetOutput() );
      eroder->SetKernel( structuringElement );
      eroder->SetForegroundValue( 2 );
      eroder->Update();

      typedef itk::ConnectedComponentImageFilter<LabelSliceType, LabelSliceType> ConnectedComponentType;
      ConnectedComponentType::Pointer connecter = ConnectedComponentType::New();
      connecter->SetInput( eroder->GetOutput() );
      connecter->FullyConnectedOff();
      connecter->Update();

      typedef itk::RelabelComponentImageFilter<LabelSliceType, LabelSliceType> LabelerType;
      LabelerType::Pointer labeler2 = LabelerType::New();
      labeler2->SetInput( connecter->GetOutput() );
      labeler2->SetMinimumObjectSize( 1000 );
      labeler2->Update();

      numberOfObjects = labeler2->GetNumberOfObjects();

      if ( numberOfObjects == 2 )
        {

        RealType x1 = 0.0;
        RealType N1 = 0.0;
        RealType x2 = 0.0;
        RealType N2 = 0.0;

        itk::ImageRegionIteratorWithIndex<LabelSliceType> It( labeler2->GetOutput(),
          labeler2->GetOutput()->GetLargestPossibleRegion() );
        for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
          {
          if ( It.Get() == 1 )
            {
            x1 += It.GetIndex()[0];
            N1++;
            }
          else if ( It.Get() == 2 )
            {
            x2 += It.GetIndex()[0];
            N2++;
            }
          }

        if ( ( x1/N1 < centerAirwayX && x2/N2 > centerAirwayX )
             || ( x1/N1 > centerAirwayX && x2/N2 < centerAirwayX ) )
          {
          labelSlice = eroder->GetOutput();
          eroded = true;
          }
        else
          {
          eroded = false;
          }
        }
      }

    if ( numberOfObjects == 2 && eroded )
      {

      typedef itk::ConnectedComponentImageFilter<LabelSliceType, LabelSliceType> ConnectedComponentType;
      ConnectedComponentType::Pointer connecter = ConnectedComponentType::New();
      connecter->SetInput( labelSlice );
      connecter->FullyConnectedOff();
      connecter->Update();

      typedef itk::RelabelComponentImageFilter<LabelSliceType, LabelSliceType> LabelerType;
      LabelerType::Pointer labeler3 = LabelerType::New();
      labeler3->SetInput( connecter->GetOutput() );
      labeler3->Update();

      typedef itk::BinaryBallStructuringElement<
                          LabelImageType::PixelType,
                          ImageDimension-1>             StructuringElementType;
      StructuringElementType structuringElement;
      structuringElement.SetRadius( 1 );
      structuringElement.CreateStructuringElement();

      typedef itk::BinaryBoundedSpaceDilateImageFilter<
        LabelSliceType, LabelSliceType, StructuringElementType>  DilaterType;

      DilaterType::Pointer dilater1 = DilaterType::New();
      dilater1->SetInput( labeler3->GetOutput() );
      dilater1->SetScaling( radius+1 );
      dilater1->SetKernel( structuringElement );
      dilater1->SetBoundedSpaceImage( holeFiller->GetOutput() );
      dilater1->SetBoundedSpaceValue( 2 );
      dilater1->SetForegroundValue( 1 );
      dilater1->Update();

      RealType x1 = 0.0;
      RealType N1 = 0.0;
      itk::ImageRegionIteratorWithIndex<LabelSliceType> It1( dilater1->GetOutput(),
        dilater1->GetOutput()->GetLargestPossibleRegion() );
      for ( It1.GoToBegin(); !It1.IsAtEnd(); ++It1 )
        {
        if ( It1.Get() != dilater1->GetForegroundValue() )
          {
          It1.Set( 0 );
          }
        else
          {
          x1 += It1.GetIndex()[0];
          N1++;
          }
        }

      DilaterType::Pointer dilater2 = DilaterType::New();
      dilater2->SetInput( labeler3->GetOutput() );
      dilater2->SetScaling( radius+1 );
      dilater2->SetKernel( structuringElement );
      dilater2->SetBoundedSpaceImage( holeFiller->GetOutput() );
      dilater2->SetBoundedSpaceValue( 2 );
      dilater2->SetForegroundValue( 2 );
      dilater2->Update();

      RealType x2 = 0.0;
      RealType N2 = 0.0;
      itk::ImageRegionIteratorWithIndex<LabelSliceType> It2( dilater2->GetOutput(),
        dilater2->GetOutput()->GetLargestPossibleRegion() );
      for ( It2.GoToBegin(); !It2.IsAtEnd(); ++It2 )
        {
        if ( It2.Get() != dilater2->GetForegroundValue() )
          {
          It2.Set( 0 );
          }
        else
          {
          x2 += It2.GetIndex()[0];
          N2++;
          }
        }

     for ( It1.GoToBegin(); !It1.IsAtEnd(); ++It1 )
       {
       if ( It1.Get() == dilater1->GetForegroundValue()
            && x1/N1 > centerAirwayX )
         {
         It1.Set( dilater2->GetForegroundValue() );
         }
       }

     for ( It2.GoToBegin(); !It2.IsAtEnd(); ++It2 )
       {
       if ( It2.Get() == dilater2->GetForegroundValue()
            && x2/N2 < centerAirwayX )
         {
         It2.Set( dilater1->GetForegroundValue() );
         }
       }

      typedef itk::AddImageFilter<LabelSliceType, LabelSliceType, LabelSliceType> AdderType;
      AdderType::Pointer adder = AdderType::New();
      adder->SetInput1( dilater1->GetOutput() );
      adder->SetInput2( dilater2->GetOutput() );
      adder->Update();

      labelSlice = adder->GetOutput();

      typedef itk::StatisticsImageFilter<LabelSliceType> StatsType;
      StatsType::Pointer stats = StatsType::New();
      stats->SetInput( adder->GetOutput() );
      stats->Update();

      /**
       * Now we want to use the live wire function to sever the lungs.
       */
      if ( stats->GetMaximum() == 3 )
        {

        typedef itk::BinaryThresholdImageFilter
          <LabelSliceType, LabelSliceType> BinaryFilterType;
        BinaryFilterType::Pointer binary = BinaryFilterType::New();
        binary->SetInput( adder->GetOutput() );
        binary->SetUpperThreshold( 3 );
        binary->SetLowerThreshold( 3 );
        binary->SetInsideValue( 1 );
        binary->SetOutsideValue( 0 );
        binary->Update();

        typedef itk::ConnectedComponentImageFilter<LabelSliceType, LabelSliceType> ConnectedComponentType;
        ConnectedComponentType::Pointer connecter = ConnectedComponentType::New();
        connecter->SetInput( binary->GetOutput() );
        connecter->FullyConnectedOff();
        connecter->Update();

        typedef itk::RelabelComponentImageFilter<LabelSliceType, LabelSliceType> LabelerType;
        LabelerType::Pointer labeler4 = LabelerType::New();
        labeler4->SetInput( connecter->GetOutput() );
        labeler4->Update();

        typedef itk::ExtractImageFilter<ImageType, SliceType> ExtracterType;
        ExtracterType::Pointer extracter = ExtracterType::New();
        extracter->SetInput( reader->GetOutput() );
        extracter->SetExtractionRegion( region );
        extracter->SetDirectionCollapseToIdentity();
        extracter->Update();

        typedef itk::LiveWireImageFunction<SliceType> LiveWireType;
        LiveWireType::Pointer livewire = LiveWireType::New();
        livewire->SetInputImage( extracter->GetOutput() );

        livewire->SetGradientMagnitudeWeight( 0.5 );
        livewire->SetZeroCrossingWeight( 0.3 );
        livewire->SetGradientDirectionWeight( 0.2 );

        for ( unsigned int i = 0; i < labeler4->GetNumberOfObjects(); i++ )
          {
          LabelSliceType::IndexType idxSeed1;
          LabelSliceType::IndexType idxSeed2;

          idxSeed1[1] = itk::NumericTraits<int>::max();
          idxSeed2[1] = itk::NumericTraits<int>::Zero;

          itk::ImageRegionIteratorWithIndex<LabelSliceType> It( labeler4->GetOutput(),
            labeler4->GetOutput()->GetLargestPossibleRegion() );
          for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
            {
            if ( It.Get() == static_cast<LabelSliceType::PixelType>( i+1 ) )
              {
              LabelSliceType::IndexType idx = It.GetIndex();
              if ( idx[1] < idxSeed1[1] )
                {
                idxSeed1 = idx;
                }
              if ( idx[1] > idxSeed2[1] )
                {
                idxSeed2 = idx;
                }
              }
            }

          livewire->SetAnchorSeed( idxSeed1 );

          LiveWireType::OutputType::Pointer path = livewire->EvaluateAtIndex( idxSeed2 );
          typedef itk::PathIterator<LabelSliceType, LiveWireType::OutputType> IteratorType;

          IteratorType ItP( labelSlice, path );
          ItP.GoToBegin();
          while ( !ItP.IsAtEnd() )
            {
            labelSlice->SetPixel( ItP.GetIndex(), 5 );
            ++ItP;
            }

          }
        itk::ImageLinearIteratorWithIndex<LabelSliceType> It( labelSlice,
          labelSlice->GetLargestPossibleRegion() );
        itk::ImageLinearIteratorWithIndex<LabelSliceType> ItL( labelExtracter->GetOutput(),
          labelExtracter->GetOutput()->GetLargestPossibleRegion() );
        for ( It.GoToBegin(), ItL.GoToBegin(); !It.IsAtEnd(); It.NextLine(), ItL.NextLine() )
          {
          It.GoToBeginOfLine();
          ItL.GoToBeginOfLine();
          bool penDown = false;
          while( !It.IsAtEndOfLine() )
            {
            LabelImageType::IndexType idx;
            for ( unsigned int d = 0; d < ImageDimension-1; d++ )
              {
              idx[d] = It.GetIndex()[d];
              }
            idx[ImageDimension-1] = s;

            if ( It.Get() == 1 && ItL.Get() == 2 )
              {
              output->SetPixel( idx, 2 );
              }
            else if ( It.Get() == 2 && ItL.Get() == 2 )
              {
              output->SetPixel( idx, 3 );
              }
            else if ( It.Get() == 3 && ItL.Get() == 2 && !penDown )
              {
              output->SetPixel( idx, 2 );
              }
            else if ( It.Get() == 5 && ItL.Get() == 2 )
              {
              penDown = true;
              output->SetPixel( idx, 3 );
              }
            else if ( It.Get() == 3 && ItL.Get() == 2 && penDown )
              {
              output->SetPixel( idx, 3 );
              }
            ++It;
            ++ItL;
            }
          }
        }
      else
        {
        labelSlice = adder->GetOutput();
        eroded = false;
        }
      }

    if ( !eroded )
      {
      RealType x1 = 0.0;
      RealType N1 = 0.0;
      RealType x2 = 0.0;
      RealType N2 = 0.0;

      itk::ImageRegionIteratorWithIndex<LabelSliceType> It( labelSlice,
        labelSlice->GetLargestPossibleRegion() );
      itk::ImageRegionIteratorWithIndex<LabelSliceType> ItL( labelExtracter->GetOutput(),
        labelExtracter->GetOutput()->GetLargestPossibleRegion() );
      for ( It.GoToBegin(), ItL.GoToBegin(); !It.IsAtEnd(); ++It, ++ItL )
        {
        if ( It.Get() == 1 && ItL.Get() == 2 )
          {
          x1 += It.GetIndex()[0];
          N1++;
          }
        else if ( It.Get() == 2 && ItL.Get() == 2 )
          {
          x2 += It.GetIndex()[0];
          N2++;
          }
        }

      for ( It.GoToBegin(), ItL.GoToBegin(); !It.IsAtEnd(); ++It, ++ItL )
        {
        LabelImageType::IndexType idx;
        for ( unsigned int d = 0; d < ImageDimension-1; d++ )
          {
          idx[d] = It.GetIndex()[d];
          }
        idx[ImageDimension-1] = s;

        if ( It.Get() == 1 && ItL.Get() == 2 )
          {
          if ( x1/N1 < centerAirwayX )
            {
            output->SetPixel( idx, 2 );
            }
          else
            {
            output->SetPixel( idx, 3 );
            }
          }
        else if ( It.Get() == 2 && ItL.Get() == 2 )
          {
          if ( x2/N2 < centerAirwayX )
            {
            output->SetPixel( idx, 2 );
            }
          else
            {
            output->SetPixel( idx, 3 );
            }
          }
        }
      }
    }

  // Now perform any clean-up
  for ( int s = labelImage->GetLargestPossibleRegion().GetSize()[2] - 1;
          s >= 0; s-- )
    {
    ImageType::IndexType index;
    index.Fill( 0 );
    index[2] = s;
    region.SetIndex( index );

    typedef itk::ExtractImageFilter<LabelImageType, LabelSliceType> LabelExtracterType;
    LabelExtracterType::Pointer labelExtracter = LabelExtracterType::New();
    labelExtracter->SetInput( output );
    labelExtracter->SetExtractionRegion( region );
    labelExtracter->SetDirectionCollapseToIdentity()
    labelExtracter->Update();

    typedef itk::BinaryThresholdImageFilter<LabelSliceType, LabelSliceType> BinaryFilterType;
    BinaryFilterType::Pointer binary = BinaryFilterType::New();
    binary->SetInput( labelExtracter->GetOutput() );
    binary->SetUpperThreshold( 2 );
    binary->SetLowerThreshold( 2 );
    binary->SetInsideValue( 2 );
    binary->SetOutsideValue( 0 );
    binary->Update();

    typedef itk::ConnectedComponentImageFilter<LabelSliceType, LabelSliceType> ConnectedComponentType;
    ConnectedComponentType::Pointer connecter = ConnectedComponentType::New();
    connecter->SetInput( binary->GetOutput() );
    connecter->FullyConnectedOff();
    connecter->Update();

    typedef itk::RelabelComponentImageFilter<LabelSliceType, LabelSliceType> LabelerType;
    LabelerType::Pointer labeler = LabelerType::New();
    labeler->SetInput( connecter->GetOutput() );
    labeler->Update();

    for ( unsigned int n = 0; n < labeler->GetNumberOfObjects(); n++ )
      {
      itk::ImageRegionIteratorWithIndex<LabelSliceType> ItL( labeler->GetOutput(),
        labeler->GetOutput()->GetLargestPossibleRegion() );

      RealType x = 0;
      RealType N = 0;
      for ( ItL.GoToBegin(); !ItL.IsAtEnd(); ++ItL )
        {
        if ( ItL.Get() == n+1 )
          {
          x += static_cast<RealType>( ItL.GetIndex()[0] );
          N++;
          }
        }

      if ( x/N > centerAirwayX )
        {
        for ( ItL.GoToBegin(); !ItL.IsAtEnd(); ++ItL )
          {
          LabelImageType::IndexType idx;
          for ( unsigned int d = 0; d < ImageDimension-1; d++ )
            {
            idx[d] = ItL.GetIndex()[d];
            }
          idx[ImageDimension-1] = s;
          if ( ItL.Get() == n+1 )
            {
            output->SetPixel( idx, 3 );
            }
          }
        }
      }
    }

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( output );
  writer->Update();

  return 0;
}
