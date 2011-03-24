#include "itkArray.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkTimeProbe.h"


#include "global.h"

int main( unsigned int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage probabilityImagePrefix labelImage intensityPercentage neighborhoodSigma numberOfSamples" << std::endl;
    exit( 1 );
    }
  itk::TimeProbe timer;
  timer.Start();  


  typedef float RealType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<int, ImageDimension> LabelImageType;
  typedef itk::Array<RealType> ArrayType;
  typedef itk::Image<ArrayType, ImageDimension> ArrayImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
  LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( argv[3] );
  labelReader->Update();

  ArrayImageType::Pointer probabilityImage = ArrayImageType::New();
  probabilityImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  probabilityImage->SetOrigin( reader->GetOutput()->GetOrigin() );
  probabilityImage->SetSpacing( reader->GetOutput()->GetSpacing() );
  probabilityImage->Allocate();

  typedef itk::LabelStatisticsImageFilter<ImageType, LabelImageType> StatsFilterType;
  StatsFilterType::Pointer stats = StatsFilterType::New();
  stats->SetInput( reader->GetOutput() );
  stats->SetLabelInput( labelReader->GetOutput() );
  stats->Update(); 

  PixelType maxIntensity = itk::NumericTraits<PixelType>::NonpositiveMin(); 
  PixelType minIntensity = itk::NumericTraits<PixelType>::max();
  for ( unsigned int i = 0; i < stats->GetNumberOfLabels(); i++ )
    {
    if ( stats->GetMaximum( i ) > maxIntensity )
      {
      maxIntensity = stats->GetMaximum( i );
      } 
    if ( stats->GetMinimum( i ) < minIntensity )
      {
      minIntensity = stats->GetMinimum( i );
      } 
    }   

  RealType percentage = atof( argv[4] );
  if ( percentage > 1 )
    {
    percentage /= 100;
    } 
  RealType sigma = ( maxIntensity - minIntensity ) * percentage;
  RealType neighborhoodSigma = atof( argv[5] );
  unsigned int numberOfSamples = atoi( argv[6] );

  RealType prefactor = 1.0 / vcl_sqrt( 2.0*vnl_math::pi*sigma*sigma );

  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
  GeneratorType::Pointer generator = GeneratorType::New();
  generator->SetSeed();

  unsigned long totalNumberOfPixels = 1;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    totalNumberOfPixels *= reader->GetOutput()->GetLargestPossibleRegion().GetSize()[i];
    } 
  unsigned int iter = 0;
  unsigned long numberOfPixelsChanged = itk::NumericTraits<unsigned long>::max();
  while ( numberOfPixelsChanged > 0.1*totalNumberOfPixels && iter++ < 10 )
    {     
    numberOfPixelsChanged = 0;

    itk::ImageRegionIteratorWithIndex<ImageType> It( reader->GetOutput(),
      reader->GetOutput()->GetLargestPossibleRegion() ); 

    for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
      { 
      ImageType::IndexType index = It.GetIndex();  
      ImageType::PixelType intensity = It.Get();

      unsigned long count = 0;
      ArrayType probability( stats->GetNumberOfLabels() );
      probability.Fill( 0.0 );
      ArrayType N( stats->GetNumberOfLabels() );
      N.Fill( 0.0 );
      while ( count++ < numberOfSamples )
        {
        ImageType::IndexType sampleIndex;
        for ( unsigned int d = 0; d < ImageDimension; d++ )
          {
          sampleIndex[d] = index[d] 
            + static_cast<int>( generator->GetNormalVariate( 0, neighborhoodSigma ) + 0.5 );
          }  
        if ( reader->GetOutput()->GetLargestPossibleRegion().IsInside( sampleIndex ) )            
          {
          ImageType::PixelType sampleIntensity = reader->GetOutput()->GetPixel( sampleIndex );
          LabelImageType::PixelType sampleLabel = labelReader->GetOutput()->GetPixel( sampleIndex );
          probability[sampleLabel] += ( prefactor * 
            vcl_exp( 0.5 * ( sampleIntensity - intensity ) / ( sigma * sigma ) ) );
          N[sampleLabel]++;
          }  
        }
      RealType totalProbability = 0;
      RealType maxProbability = 0;
      int maxLabel;
      for ( unsigned int n = 0; n < stats->GetNumberOfLabels(); n++ )
        {
        if ( N[n] == 0 )
          {
          continue;
          }  
        probability[n] /= N[n];
        if ( probability[n] > maxProbability )
          {
          maxLabel = n;
          maxProbability = probability[n];
          } 
        totalProbability += probability[n];
        }
      for ( unsigned int n = 0; n < stats->GetNumberOfLabels(); n++ )
        {
        probability[n] /= totalProbability;
        }  
      LabelImageType::PixelType label = labelReader->GetOutput()->GetPixel( index );
      if ( maxLabel != label ) 
        {
        numberOfPixelsChanged++;
        labelReader->GetOutput()->SetPixel( index, maxLabel );
        } 
      probabilityImage->SetPixel( index, probability ); 
      }
    std::cout << "Iteration " << iter << ": " 
      << static_cast<RealType>( numberOfPixelsChanged ) 
         / static_cast<RealType>( totalNumberOfPixels )*100 << "% pixels changed." << std::endl;
    } 

  {
  std::string filename = std::string( argv[2] ) + std::string( "Label.nii" );

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput( labelReader->GetOutput() );
  writer->Update();

  }

  for ( unsigned int i = 0; i < stats->GetNumberOfLabels(); i++ )
    {
    ImageType::Pointer output = ImageType::New();
    output->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
    output->SetOrigin( reader->GetOutput()->GetOrigin() );
    output->SetSpacing( reader->GetOutput()->GetSpacing() );
    output->Allocate();

    itk::ImageRegionIteratorWithIndex<ArrayImageType> It( probabilityImage,
      probabilityImage->GetLargestPossibleRegion() ); 
    itk::ImageRegionIteratorWithIndex<ImageType> ItO( output,
      output->GetLargestPossibleRegion() ); 
   
    for ( It.GoToBegin(), ItO.GoToBegin(); !It.IsAtEnd(); ++It, ++ItO )
      {
      ItO.Set( It.Get()[i] );  
      } 
   
    itk::OStringStream buf;
    buf << i;

    std::string filename = std::string( argv[2] ) + std::string( "." ) 
      + buf.str() + ".nii";

    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( filename.c_str() );
    writer->SetInput( output );
    writer->Update();
    }
  timer.Stop();
  std::cout << "Time elapsed:  " << timer.GetMeanTime() << std::endl;

  return 0;
}
