#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include "itkListSampleToHistogramGenerator.h"
#include "itkListSample.h"
#include "itkVector.h"
#include "itkMeanCalculator.h"
#include "itkCovarianceCalculator.h"

#include <fstream.h>

int main(int argc, char *argv[])        
{
  if ( argc != 4 )
    {
    std::cout << "Usage: " << argv[0] << " control_file_list subject_file_list number_of_bins" << std::endl;
    exit( 1 );
    }
  
  const unsigned int ImageDimension = 3;
  typedef int PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageRegionIterator<ImageType> IteratorType;

  typedef itk::ImageFileReader<ImageType> ReaderType;

  ImageType::Pointer mask = ImageType::New();  

  std::string fn1 = std::string( argv[1] );
  std::string fn2 = std::string( argv[2] );

  unsigned long ct1 = 0;
  unsigned long ct2 = 0;
  const unsigned int maxChar = 512;
  char lineBuffer[maxChar]; 
  char filenm[maxChar];

  unsigned int filecount1 = 0;
  unsigned int filecount2 = 0;
  std::ifstream inputStreamA( fn1.c_str(), std::ios::in );
  if ( !inputStreamA.is_open() )
    {
    std::cout << "Can't open file: " << argv[1] << std::endl;  
    return -1;
    }
  while ( !inputStreamA.eof() )
    {
    inputStreamA.getline( lineBuffer, maxChar, '\n' ); 
    if ( sscanf( lineBuffer, "%s ",filenm) != 1 )
      {
      continue;
      }
    else
      {
      filecount1++;
      }
    }
  inputStreamA.close();  
  
  std::ifstream inputStreamB( fn2.c_str(), std::ios::in );
  if ( !inputStreamB.is_open() )
    {
    std::cout << "Can't open file: " << argv[2] << std::endl;  
    return -1;
    }
  while ( !inputStreamB.eof() )
    {
    inputStreamB.getline( lineBuffer, maxChar, '\n' ); 
    if ( sscanf( lineBuffer, "%s ",filenm) != 1 )
      {
      continue;
      }
    else
      {
      filecount2++;
      }
    }
  inputStreamB.close();
 
  std::cout << " NFiles1 " << filecount1 << " NFiles2 " << filecount2 << std::endl;

  std::vector<bool> controlbool( filecount1 + filecount2 );
  std::vector<std::string> filenames( filecount1 + filecount2 );

  unsigned int ct = 0;
  inputStreamA.open( fn1.c_str(), std::ios::in );
  while ( !inputStreamA.eof() )
    {
    inputStreamA.getline( lineBuffer, maxChar, '\n' ); 
    if ( sscanf( lineBuffer, "%s ",filenm) != 1 )
      {
      continue;
      }
      else
      {
      filenames[ct] = filenm;
      controlbool[ct] = true;
      ct++;
      }
    }
  inputStreamA.close(); 
 
  inputStreamB.open( fn2.c_str(), std::ios::in );
  if ( !inputStreamB.is_open() )
    {
    std::cout << "Can't open parameter file: " << argv[1] << std::endl;  
    return -1;
    }
  while ( !inputStreamB.eof() )
    {
    inputStreamB.getline( lineBuffer, maxChar, '\n' ); 
    if ( sscanf( lineBuffer, "%s ",filenm) != 1 )
      {
      continue;
      }
    else
      {
      filenames[ct] = filenm;
      controlbool[ct] = false;
      ct++;  
      }
    }
  inputStreamB.close();
  
  for ( unsigned int i = 0; i< filecount1 + filecount2; i++) 
    {
    std::cout << " n1 " << filenames[i] << " is " << controlbool[i] << std::endl;
    }

  // Calculation of statistics corresponding to first set of files 
  typedef itk::Vector<PixelType, 1> VectorType;

  typedef itk::Statistics::ListSample<VectorType> SampleType;
  SampleType::Pointer sample1 = SampleType::New();
  SampleType::Pointer sample2 = SampleType::New();

  std::ifstream inputStream( fn1.c_str(), std::ios::in );
  if ( !inputStream.is_open() )
    {
    std::cout << "Can't open file: " << argv[1] << std::endl;  
    return -1;
    }
  
  while ( !inputStream.eof() )
    {
    inputStream.getline( lineBuffer, maxChar, '\n' ); 
    if ( sscanf( lineBuffer, "%s ",filenm) != 1 )
      {
      continue;
      }
    else
      {
      ImageType::Pointer image2 = ImageType::New();
      
      ReaderType::Pointer reader2 = ReaderType::New();
      reader2->SetFileName( ( filenm + std::string( "cropped.hdr" ) ).c_str() ); 
      reader2->Update(); 
      ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( ( filenm + std::string( "segmentation_cropped.hdr" ) ).c_str() ); 
      reader->Update(); 
      
      try
        {   
        image2 = reader2->GetOutput(); 
        mask = reader->GetOutput(); 
        }
      catch(...)
        {
        std::cout << " Error reading " << std::string(filenm) << std::endl;
        return 0;
        } 

      IteratorType It1( image2, image2->GetLargestPossibleRegion() );
      IteratorType It2( mask, mask->GetLargestPossibleRegion() );
      
      It1.GoToBegin();
      It2.GoToBegin();
      while ( !It1.IsAtEnd() )
        {
        if ( It2.Get() != itk::NumericTraits<PixelType>::Zero )
          {
          VectorType V;
          V[0] = It1.Get();
          sample1->PushBack( V );
          }  
        ++It1;
        ++It2;
        }
      }
    }
  inputStream.close();

  typedef itk::Statistics::ListSampleToHistogramGenerator<SampleType, float> GeneratorType;
  GeneratorType::Pointer generator1 = GeneratorType::New();
  GeneratorType::HistogramType::SizeType size;
  size.Fill( atoi( argv[3] ) );
  
  VectorType minH;
  minH.Fill( 1 );
  VectorType maxH;
  maxH.Fill( 2100 );
  
  generator1->SetListSample( sample1 );
  generator1->SetMarginalScale( 10.0 );
  generator1->SetNumberOfBins( size );
  generator1->SetHistogramMin( minH );
  generator1->SetHistogramMax( maxH );
  generator1->Update();

  // Calculation of statistics corresponding to second set of files 

  std::ifstream inputStream2( fn2.c_str(), std::ios::in );
  if ( !inputStream2.is_open() )
    {
    std::cout << "Can't open parameter file: " << argv[1] << std::endl;  
    return -1;
    }
  while ( !inputStream2.eof() )
    {
    inputStream2.getline( lineBuffer, maxChar, '\n' ); 
  
    if ( sscanf( lineBuffer, "%s ",filenm) != 1 )
      {
      continue;
      }
    else
      {
      ImageType::Pointer image2 = ImageType::New();
      
      ReaderType::Pointer reader2 = ReaderType::New();
      reader2->SetFileName( ( filenm + std::string( "cropped.hdr" ) ).c_str() ); 
      reader2->Update(); 
      ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( ( filenm + std::string( "segmentation_cropped.hdr" ) ).c_str() ); 
      reader->Update(); 
      
      try
        {   
        image2 = reader2->GetOutput(); 
        mask = reader->GetOutput(); 
        }
      catch(...)
        {
        std::cout << " Error reading " << std::string(filenm) << std::endl;
        return 0;
        } 

      IteratorType It1( image2, image2->GetLargestPossibleRegion() );
      IteratorType It2( mask, mask->GetLargestPossibleRegion() );
      
      It1.GoToBegin();
      It2.GoToBegin();
      while ( !It1.IsAtEnd() )
        {
        if ( It2.Get() != itk::NumericTraits<PixelType>::Zero )
          {
          VectorType V;
          V[0] = It1.Get();
          sample2->PushBack( V );
          }  
        ++It1;
        ++It2;
        }
      }
    }
  inputStream2.close();
  
  GeneratorType::Pointer generator2 = GeneratorType::New();
  generator2->SetListSample( sample2 );
  generator2->SetMarginalScale( 10.0 );
  generator2->SetNumberOfBins( size );
  generator2->SetHistogramMin( minH );
  generator2->SetHistogramMax( maxH );
  generator2->Update();
  
  GeneratorType::HistogramType::ConstPointer histogram1 = generator1->GetOutput();
  GeneratorType::HistogramType::ConstPointer histogram2 = generator2->GetOutput();
  
  ofstream str1( "histogram1.dat" );
  ofstream str2( "histogram2.dat" );

  typedef itk::Statistics::MeanCalculator<SampleType> MeanCalculatorType;
  MeanCalculatorType::Pointer mean1 = MeanCalculatorType::New();
  mean1->SetInputSample( sample1 );
  mean1->Update();
  MeanCalculatorType::Pointer mean2 = MeanCalculatorType::New();
  mean2->SetInputSample( sample2 );
  mean2->Update();

  typedef itk::Statistics::CovarianceCalculator<SampleType> CovarianceCalculatorType;
  CovarianceCalculatorType::Pointer covar1 = CovarianceCalculatorType::New();
  covar1->SetInputSample( sample1 );
  covar1->SetMean( mean1->GetOutput() );
  covar1->Update();
  CovarianceCalculatorType::Pointer covar2 = CovarianceCalculatorType::New();
  covar2->SetInputSample( sample2 );
  covar2->SetMean( mean2->GetOutput() );
  covar2->Update();
  

  unsigned int i;
  for( i = 0; i < histogram1->Size(); i++ )
    {
    str1 << histogram1->GetMeasurement( i, 0 ) << " "
        << histogram1->GetFrequency( i, 0 ) << " "
        << histogram1->GetTotalFrequency() << " " 
        << mean1->GetOutput()[0][0] << " " 
        << (*covar1->GetOutput())(0, 0) << std::endl;
    }
    
  for( i = 0; i < histogram2->Size(); i++ )
    {
    str2 << histogram2->GetMeasurement( i, 0 ) << " "
        << histogram2->GetFrequency( i, 0 ) << " "
        << histogram2->GetTotalFrequency() << " " 
        << mean2->GetOutput()[0][0] << " " 
        << (*covar2->GetOutput())(0, 0) << std::endl;
    }
  
  float var1 = (*covar1->GetOutput())(0, 0);
  float n1 = static_cast<float>( histogram1->GetTotalFrequency() );
  float var2 = (*covar2->GetOutput())(0, 0);
  float n2 = static_cast<float>( histogram2->GetTotalFrequency() );
  
  float tvalue = ( mean1->GetOutput()[0][0] - mean2->GetOutput()[0][0] ) 
                 / sqrt( var1/n1 + var2/n2 ); 
  float lambda = ( var1/n1 + var2/n2 ) * ( var1/n1 + var2/n2 ) 
                 / ( var1*var1/(n1*n1*(n1-1)) 
                   + var2*var2/(n2*n2*(n2-1)) );

  std::cout << tvalue << ", " << lambda << std::endl;

  return 0;
 
}     


      

