#include "itkBinaryThresholdImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkTimeProbe.h"

#include <string>
#include <fstream.h>

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cout << "Usage: " << argv[0] << " maximalResponseImage maskImage outputLabelImage direction thresholdPercentage" << std::endl;     
    exit( 0 );
    }   

  itk::TimeProbe timer;
  timer.Start();

  typedef float RealType;

  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Image<PixelType, ImageDimension-1> SliceType;
  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;
  typedef itk::Image<unsigned int, ImageDimension-1> LabelSliceType;

  typedef itk::ImageFileReader<LabelImageType> MaskReaderType;
  MaskReaderType::Pointer maskreader = MaskReaderType::New();
  maskreader->SetFileName( argv[2] );
  maskreader->Update();

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  unsigned int dim = atoi( argv[4] );
  RealType t = atof( argv[5] );

  LabelImageType::Pointer labels = LabelImageType::New();
  labels->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  labels->SetOrigin( reader->GetOutput()->GetOrigin() );
  labels->SetSpacing( reader->GetOutput()->GetSpacing() );
  labels->Allocate();
  labels->FillBuffer( 0 );

  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    if ( d == dim )
      {
      continue;
      } 
    for ( unsigned int i = 0; i < reader->GetOutput()->GetLargestPossibleRegion().GetSize()[d]; i++ )
      {
      RealImageType::RegionType region;
      RealImageType::RegionType::SizeType size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
      size[d] = 0;
      RealImageType::IndexType index;
      index.Fill( 0 );
      index[d] = i;
      region.SetIndex( index );
      region.SetSize( size );

      typedef itk::ExtractImageFilter<RealImageType, SliceType> ExtracterType;
      ExtracterType::Pointer extracter = ExtracterType::New();
      extracter->SetInput( reader->GetOutput() );
      extracter->SetExtractionRegion( region );
      extracter->Update();

      typedef itk::ExtractImageFilter<LabelImageType, LabelSliceType> MaskExtracterType;
      MaskExtracterType::Pointer maskExtracter = MaskExtracterType::New();
      maskExtracter->SetInput( maskreader->GetOutput() );
      maskExtracter->SetExtractionRegion( region );
      maskExtracter->Update();

      typedef itk::LabelStatisticsImageFilter<SliceType, LabelSliceType> HistogramGeneratorType;
      HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
      stats->SetInput( extracter->GetOutput() );
      stats->SetLabelInput( maskExtracter->GetOutput() );
      stats->UseHistogramsOff();
      stats->Update();
  
      unsigned int numberOfBins = 100;

      for ( unsigned int n = 1; n <= 2; n++ )
        {
        if ( stats->GetCount( n ) == 0 )
          {
          continue;
          }   

        RealType delta = ( static_cast<RealType>( stats->GetMaximum( n ) ) 
          - static_cast<RealType>( stats->GetMinimum( n ) ) )
          / static_cast<RealType>( numberOfBins );
        RealType lowerBound = static_cast<RealType>( stats->GetMinimum( n ) ) - delta;
        RealType upperBound = static_cast<RealType>( stats->GetMaximum( n ) ) + delta;
      
        stats->UseHistogramsOn();
        stats->SetHistogramParameters( numberOfBins, lowerBound, upperBound );
        stats->Update();
    
        typedef HistogramGeneratorType::HistogramType  HistogramType;
        const HistogramType *histogram = stats->GetHistogram( n );


/*
        unsigned long cumulativeSum[numberOfBins];
        cumulativeSum[0] = histogram->GetFrequency( 0, 0 );
        RealType percentage = static_cast<RealType>( cumulativeSum[0] ) 
          / static_cast<RealType>( histogram->GetTotalFrequency() );        
      
        unsigned int bin = 1;

        while ( percentage < t && bin < numberOfBins )
          {
          std::cout << "   " << percentage << ", " << t << ", " << bin << ", " << numberOfBins << std::endl;
          cumulativeSum[bin] = cumulativeSum[bin-1] + histogram->GetFrequency( bin, 0 );
          percentage = static_cast<RealType>( cumulativeSum[bin] ) 
            / static_cast<RealType>( histogram->GetTotalFrequency() );
          bin++;        
          }

  
        RealType m2 = static_cast<RealType>( histogram->GetMeasurement( bin, 0 ) );
        RealType m1 = static_cast<RealType>( histogram->GetMeasurement( bin-1, 0 ) );
        RealType p2 = static_cast<RealType>( cumulativeSum[bin] ) 
          / static_cast<RealType>( histogram->GetTotalFrequency() );        
        RealType p1 = static_cast<RealType>( cumulativeSum[bin-1] ) 
          / static_cast<RealType>( histogram->GetTotalFrequency() );        
        if ( p1 == p2 )
          {
          continue;
          } 
        RealType lowerThreshold = m2 
          - ( p2 - t ) * ( m2 - m1 ) / ( p2 - p1 );     
        */          

        RealType lower = histogram->Quantile( 0, t );
        RealType upper = stats->GetMaximum( n );
        if ( lower >= upper )
          {
          continue;
          }

        typedef itk::BinaryThresholdImageFilter<SliceType, SliceType> ThresholderType;
        ThresholderType::Pointer thresholder = ThresholderType::New();
        thresholder->SetInput( extracter->GetOutput() );
        thresholder->SetInsideValue( 255 );
        thresholder->SetOutsideValue( 0 );
        thresholder->SetLowerThreshold( lower ); 
        thresholder->SetUpperThreshold( upper ); 
        thresholder->Update();

        itk::ImageRegionIteratorWithIndex<LabelSliceType> ItM( maskExtracter->GetOutput(), 
          maskExtracter->GetOutput()->GetLargestPossibleRegion() );
        itk::ImageRegionIteratorWithIndex<SliceType> ItT( thresholder->GetOutput(), 
          thresholder->GetOutput()->GetLargestPossibleRegion() );
         
        ItM.GoToBegin();
        ItT.GoToBegin();
        while ( !ItM.IsAtEnd() )
          {
          if ( ItT.Get() > 0 && ItM.Get() == n )
            {
            RealImageType::IndexType index;
            unsigned int count = 0;
            for ( unsigned int dd = 0; dd < ImageDimension; dd++ )
              {
              if ( d == dd )
                {
                index[dd] = i;
                } 
              else
                { 
                index[dd] = ItM.GetIndex()[count++];
                }  
              }  
            labels->SetPixel( index, labels->GetPixel( index ) + 1 );
            }  
          ++ItM; 
          ++ItT;
          }
        }  
      }  
    }

  itk::ImageRegionIteratorWithIndex<LabelImageType> ItM( maskreader->GetOutput(), 
    maskreader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<LabelImageType> ItT( labels, 
    labels->GetLargestPossibleRegion() );
    
  for ( ItM.GoToBegin(), ItT.GoToBegin(); !ItM.IsAtEnd(); ++ItM, ++ItT )
    {
    if ( ItT.Get() == ImageDimension-1 )
      {
      ItT.Set( ItM.Get() );
      }  
    else
      {
      ItT.Set( 0 );
      } 
    }

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( labels );
  writer->Update(); 


}
