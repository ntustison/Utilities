/*=========================================================================
  
  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: PermutationTests.cxx,v $
  Language:  C++      
  Date:      $Date: 2008/05/03 01:52:26 $
  Version:   $Revision: 1.1 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
  
=========================================================================*/


#include <vector>
#include <cstdlib> 
#include <ctime> 
#include <iostream>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include "itkMersenneTwisterRandomVariateGenerator.h"            

#include "itkMinimumMaximumImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include "itkDiscreteGaussianImageFilter.h"

template <class TImage>
typename TImage::Pointer 
MakeNewImage( typename TImage::Pointer image, typename TImage::PixelType initval )
{
  typename TImage::Pointer newimage = TImage::New();
  newimage->SetLargestPossibleRegion( image->GetLargestPossibleRegion() );
  newimage->SetBufferedRegion( image->GetLargestPossibleRegion() );
  newimage->SetLargestPossibleRegion( image->GetLargestPossibleRegion() );
  newimage->Allocate(); 
  newimage->SetSpacing( image->GetSpacing() );
  newimage->SetOrigin( image->GetOrigin() );
  newimage->FillBuffer( initval );
  return newimage;
}

template <class TImage>
typename TImage::Pointer 
SmoothImage( typename TImage::Pointer image, float sig )
{
  typedef itk::DiscreteGaussianImageFilter<TImage, TImage> dgf;
  typename dgf::Pointer filter = dgf::New();
  filter->SetVariance( sig );
  filter->SetUseImageSpacingOn();
  filter->SetMaximumError( 0.01f );
  filter->SetInput( image );
  filter->Update();
  return filter->GetOutput();
}

template <class TImage>
unsigned int
GetClusterStat( typename TImage::Pointer image, 
                float Tthreshold, 
                unsigned int minSize, 
                unsigned int whichstat,
                std::string outfn, 
                bool TRUTH)
{
  typedef float RealPixelType;
  typedef TImage ImageType;
  
  typedef itk::Image<int, TImage::ImageDimension> InternalImageType;

  typedef itk::BinaryThresholdImageFilter<ImageType, InternalImageType> ThresholdFilterType;  
  typename ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
  threshold->SetInput( image );
  threshold->SetInsideValue( itk::NumericTraits<int>::One );
  threshold->SetOutsideValue( itk::NumericTraits<int>::Zero );
  threshold->SetLowerThreshold( Tthreshold );
  threshold->SetUpperThreshold( itk::NumericTraits<RealPixelType>::max() );
  threshold->Update();

  typedef itk::ConnectedComponentImageFilter<InternalImageType, InternalImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( threshold->GetOutput() );

  typedef itk::RelabelComponentImageFilter<InternalImageType, InternalImageType> RelabelType;
  typename RelabelType::Pointer relabel = RelabelType::New();
  filter->SetFullyConnected( true );
  relabel->SetInput( filter->GetOutput() );
  relabel->SetMinimumObjectSize( minSize );
  try
    {
    relabel->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Relabel: exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

  std::vector<unsigned int> histogram( relabel->GetNumberOfObjects() + 1, 0 );
  
  itk::ImageRegionIteratorWithIndex<InternalImageType> It( 
      relabel->GetOutput(), relabel->GetOutput()->GetLargestPossibleRegion() );  
  
  for (  It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
      if ( It.Get() > 0 ) 
        {
        histogram[It.Get()] = histogram[It.Get()]+1;
        }
    }

  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( It.Get() > 0 ) 
      {
      It.Set( histogram[It.Get()] );
      }
    }

  if (TRUTH)
    {
    typedef itk::ImageFileWriter<InternalImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( ( outfn + std::string( "Clusters.hdr" ) ).c_str() );
    writer->SetInput( relabel->GetOutput() ); 
    writer->Write();   
    }

  return histogram[whichstat];
}  

int main(int argc, char *argv[])        
{
  typedef float RealPixelType;
  const unsigned int ImageDimension = 3;
  typedef itk::Vector<float, ImageDimension>         VectorType;
  typedef itk::Image<VectorType,ImageDimension>      FieldType;
  typedef itk::Image<RealPixelType,ImageDimension>   ImageType;
  typedef itk::ImageFileReader<ImageType>            ReaderType;
  typedef itk::ImageFileWriter<ImageType>            WriterType;
  typedef ImageType::IndexType IndexType;
  typedef ImageType::SizeType SizeType;
  typedef ImageType::SpacingType SpacingType;
  typedef itk::ImageRegionIterator<ImageType> IteratorType;

   
  if ( argc < 5 )     
  { 
    std::cout << "Useage ex:  "<< std::endl; 
    std::cout << argv[0] << " controlslist.txt subjectslist.txt uselog outfn smoothvarCont smoothvarSubj whichstat NPermutations Tthreshold {ClustThresh} {roiimage.hdr}  " << std::endl; 
    std::cout << " if uselog then we take the log of the input image ( for jacobians) " << std::endl;
    std::cout << " DER defines output filename prefix. " << std::endl;
    std::cout << " smoothvar  - this entry gives the amount of smoothing applied to variance estimates.  Helpful when sample size is small. " << std::endl;
    std::cout << " uselog  - bool, says if you want to use logs. " << std::endl;
    std::cout << " cont.txt  -  a list of control filenames, 1 per line " << std::endl;
    std::cout << " subj.txt  -  a list of subject filenames, 1 per line " << std::endl;
    std::cout << " whichstat -- 0 = size,  1 = sum,  2 = mean " << std::endl;
    return 1;
  }           

  std::string fn1 = std::string( argv[1] );
  std::string fn2 = std::string( argv[2] );
  bool uselog = atoi( argv[3] );
  std::string outfn = std::string( argv[4] );
  float smoothvar1 = atof( argv[5] );
  float smoothvar2 = atof( argv[6] );
  unsigned int whichstat = atoi( argv[7] );
  unsigned int NPermutations = atoi( argv[8] );
  float Tthreshold = atof( argv[9] );
  unsigned int ClustThresh = 10;
  if ( argc > 10 ) ClustThresh = atoi( argv[10] );
  std::cout << " params : uselog " << uselog << " smooth? " << smoothvar1 << std::endl;
  std::string roifn = "";
  if (argc > 11) 
    {
    roifn = std::string( argv[11] );
    }

  ImageType::Pointer image2 = NULL; 
  ImageType::Pointer avgimage1 = NULL; 
  ImageType::Pointer varimage1 = NULL; 
  ImageType::Pointer avgimage2 = NULL; 
  ImageType::Pointer varimage2 = NULL; 
  ImageType::Pointer ttestimg = NULL; 
  ImageType::Pointer ROIimg = NULL; 
  ImageType::Pointer weightImage = NULL;
  ImageType::Pointer weight1Image = NULL; 
  ImageType::Pointer weight2Image = NULL; 
  ImageType::Pointer ones = NULL;
  

  std::cout << "NPermutations = " << NPermutations << std::endl;


  if  ( argc > 11 )
    {
    std::cout <<" reading roi image " << roifn << std::endl;
    ReaderType::Pointer reader2 = ReaderType::New();
    reader2->SetFileName(roifn.c_str()); 
    reader2->UpdateLargestPossibleRegion();
    try
      {   
      ROIimg = reader2->GetOutput(); 
      }
    catch(...)
      {
      std::cout << " Error reading ROI image " << std::endl;
      return 0;
      }
    }

  // now do the recursive average
  unsigned long ct1 = 0;
  unsigned long ct2 = 0;
  const unsigned int maxChar = 512;
  char lineBuffer[maxChar]; 
  char filenm[maxChar];

  unsigned int clustersizes;
  unsigned int fakeclustersizes;

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
      ReaderType::Pointer reader2 = ReaderType::New();
      reader2->SetFileName( filenm ); 
      reader2->UpdateLargestPossibleRegion();
      try
        {   
        image2 = reader2->GetOutput(); 
        }
      catch(...)
        {
        std::cout << " Error reading " << std::string(filenm) << std::endl;
        return 0;
        } 
      if ( !avgimage1 )
        {
        avgimage1 = MakeNewImage<ImageType>( reader2->GetOutput(), 0.0 );
        avgimage2 = MakeNewImage<ImageType>( reader2->GetOutput(), 0.0 );
        varimage1 = MakeNewImage<ImageType>( reader2->GetOutput(), 0.0 );
        varimage2 = MakeNewImage<ImageType>( reader2->GetOutput(), 0.0 );
        ttestimg = MakeNewImage<ImageType>( reader2->GetOutput(), 0.0 );
        if ( !ROIimg )
          {
          ROIimg = MakeNewImage<ImageType>( reader2->GetOutput(), 1.0 );
          }
        }
      float weight = static_cast<float>( ct1 + 1 );  
      float wt2 = 1.0/weight;
      float wt1 = 1.0 - wt2;
      float wt3;

      if ( ct1 > 0 )
        {
        wt3 = 1.0/( weight - 1.0 );
        }
      
      IteratorType It1( avgimage1, avgimage1->GetLargestPossibleRegion() );
      IteratorType It2( image2, image2->GetLargestPossibleRegion() );
      IteratorType It3( varimage1, varimage1->GetLargestPossibleRegion() );
      IteratorType It4( ROIimg, ROIimg->GetLargestPossibleRegion() );
      
      It1.GoToBegin();
      It2.GoToBegin();
      It3.GoToBegin();
      It4.GoToBegin();
      while ( !It1.IsAtEnd() )
        {
        if ( It4.Get() != itk::NumericTraits<RealPixelType>::Zero )
          {
          float pix1 = static_cast<float>( It1.Get() );
          float pix2 = static_cast<float>( It2.Get() );
          if ( uselog )
            {
            if ( pix2 > 1.0e-11 )
              {
              pix2 = log( pix2 );
              }
            else
              {
              pix2 = 0.0;
              }
            }
          It1.Set( pix1*wt1 + pix2*wt2 );

          if ( ct1 > 0 )
            {
            float pix3 = static_cast<float>( It3.Get() );
            It3.Set( wt1*pix3 + wt3*( pix2 - pix1 )*( pix2 - pix1 ) );
            }
          }  
        ++It1;
        ++It2;
        ++It3;  
        ++It4;  
        }
      ct1++; 
      }
    }
  inputStream.close();

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
      ReaderType::Pointer reader2 = ReaderType::New();
      reader2->SetFileName( filenm ); 
      reader2->UpdateLargestPossibleRegion();
      try
        {   
        image2 = reader2->GetOutput(); 
        }
      catch(...)
        {
        std::cout << " Error reading " << std::string( filenm ) << std::endl;
        return 0;
        } 
        
      float weight = static_cast<float>( ct2 + 1 );  
      float wt2 = 1.0/weight;
      float wt1 = 1.0 - wt2;
      float wt3;

      if ( ct2 > 0 )
        {
        wt3 = 1.0/( weight - 1.0 );
        }
      
      IteratorType It1( avgimage2, avgimage2->GetLargestPossibleRegion() );
      IteratorType It2( image2, image2->GetLargestPossibleRegion() );
      IteratorType It3( varimage2, varimage2->GetLargestPossibleRegion() );
      IteratorType It4( ROIimg, ROIimg->GetLargestPossibleRegion() );
      
      It1.GoToBegin();
      It2.GoToBegin();
      It3.GoToBegin();
      It4.GoToBegin();
      while ( !It1.IsAtEnd() )
        {
        if ( It4.Get() != itk::NumericTraits<RealPixelType>::Zero )
          {
          float pix1 = static_cast<float>( It1.Get() );
          float pix2 = static_cast<float>( It2.Get() );
          if ( uselog )
            {
            if ( pix2 > 1.0e-11 )
              {
              pix2 = log( pix2 );
              }
            else
              {
              pix2 = 0.0;
              }
            }
          It1.Set( pix1*wt1 + pix2*wt2 );

          if ( ct2 > 0 )
            {
            float pix3 = static_cast<float>( It3.Get() );
            It3.Set( wt1*pix3 + wt3*( pix2 - pix1 )*( pix2 - pix1 ) );
            }
          }  
        ++It1;
        ++It2;
        ++It3;  
        ++It4;  
        }
      ct2++; 
      }
    }
  inputStream2.close();

  // now compute the timages 
  //   first, scale variances & smooth, if so desired

  if ( smoothvar1 > 0.0 ) 
    {
    varimage1 = SmoothImage<ImageType>( varimage1, smoothvar1 );
    }
  if ( smoothvar2 > 0.0 ) 
    {
    varimage2 = SmoothImage<ImageType>( varimage2, smoothvar2 );
    }
  
  // Create the t-test image     
  
  std::cout << " t-test begin " << filecount1 << " & " << filecount2 << std::endl;

  RealPixelType n1 = static_cast<RealPixelType>( filecount1 );
  RealPixelType n2 = static_cast<RealPixelType>( filecount2 );

  IteratorType It1( avgimage1, avgimage1->GetLargestPossibleRegion() );
  IteratorType It2( avgimage2, avgimage2->GetLargestPossibleRegion() );
  IteratorType It3( varimage1, varimage1->GetLargestPossibleRegion() );
  IteratorType It4( varimage2, varimage2->GetLargestPossibleRegion() );
  IteratorType It5( ROIimg, ROIimg->GetLargestPossibleRegion() );
  IteratorType It6( ttestimg, ttestimg->GetLargestPossibleRegion() );
  
  It1.GoToBegin();
  It2.GoToBegin();
  It3.GoToBegin();
  It4.GoToBegin();
  It5.GoToBegin();
  It6.GoToBegin();
  while ( !It1.IsAtEnd() )
    {
    if ( It5.Get() != itk::NumericTraits<RealPixelType>::Zero )
      {
      RealPixelType den = sqrt( It3.Get()/n1 + It4.Get()/n2 );
      if ( den > 1e-6 )
        {
        It6.Set( ( It1.Get() - It2.Get() ) / den);
        }
      else
        {
        It6.Set( static_cast<RealPixelType>( 0 ) );
        }
      }  
    else
      {
      It6.Set( static_cast<RealPixelType>( 0 ) );
      }
    ++It1;
    ++It2;
    ++It3;  
    ++It4;  
    ++It5;  
    ++It6;  
    }
      
  std::cout << " t-test end " << std::endl;

  WriterType::Pointer writer = WriterType::New();


  
  writer->SetFileName(  (outfn+std::string("ttest.hdr")).c_str());
  writer->SetInput( ttestimg ); 
  writer->Write();   
  writer->SetFileName(  (outfn+std::string("avg1.hdr")).c_str());
  writer->SetInput( avgimage1 ); 
  writer->Write();   
  writer->SetFileName(  (outfn+std::string("avg2.hdr")).c_str());
  writer->SetInput( avgimage2 ); 
  writer->Write();   
  writer->SetFileName(  (outfn+std::string("var1.hdr")).c_str());
  writer->SetInput( varimage1 ); 
  writer->Write();   
  writer->SetFileName(  (outfn+std::string("var2.hdr")).c_str());
  writer->SetInput( varimage2 ); 
  writer->Write();   

  return 0;

  std::cout << " Thresh " << Tthreshold << std::endl;
  clustersizes = GetClusterStat<ImageType>( ttestimg, Tthreshold, ClustThresh, whichstat, outfn, true);
  
  std::cout << " writing output " << outfn << " TRUE maxclust " << clustersizes << std::endl;

  srand( time( NULL ) );
  int randind;  
  
  //itk::MersenneTwisterRandomVariateGenerator::Pointer rand = MersenneTwisterRandomVariateGenerator::New();

  // set up the histogram of clustersizes
  // the histogram length is of maximum cluster size 

  std::vector<unsigned int> histogramofsizes( clustersizes + 1, 0 );

  for ( unsigned int permct = 0; permct < NPermutations; permct++ )
    {
    for ( unsigned int i = 0; i< (filecount1+filecount2); i++ )  
      {
      controlbool[i] = false;
      }
    int ncont = 0;
    while ( ncont != filecount1 )
      {
      randind = ( rand() % ( filecount1+filecount2 ) ); 
      if ( controlbool[randind] == false ) 
        {
        controlbool[randind]=true;
        ncont++;
        }
      }
      
    avgimage1 = NULL;
    
    ct1 = 0;
    for ( unsigned int qq = 0; qq < ( filecount1+filecount2 ); qq++ )
      {
      if ( controlbool[qq] == true )
        {
        std::string filenma=filenames[qq];
        ReaderType::Pointer reader2 = ReaderType::New();
        reader2->SetFileName(filenma.c_str()); 
        reader2->UpdateLargestPossibleRegion();
        try
          {   
          image2 = reader2->GetOutput(); 
          }
        catch(...)
          {
          std::cout << " Error reading " << std::string(filenm) << std::endl;
          return 0;
          } 
        if ( !avgimage1 )
          {
          avgimage1 = MakeNewImage<ImageType>( reader2->GetOutput(), 0.0 );
          avgimage2 = MakeNewImage<ImageType>( reader2->GetOutput(), 0.0 );
          varimage1 = MakeNewImage<ImageType>( reader2->GetOutput(), 0.0 );
          varimage2 = MakeNewImage<ImageType>( reader2->GetOutput(), 0.0 );
          ttestimg = MakeNewImage<ImageType>( reader2->GetOutput(), 0.0 );
          if ( !ROIimg )
            {
            ROIimg = MakeNewImage<ImageType>( reader2->GetOutput(), 1.0 );
            }
          }

        float weight = static_cast<float>( ct1 + 1 );  
        float wt2 = 1.0/weight;
        float wt1 = 1.0 - wt2;
        float wt3;

        if ( ct1 > 0 )
          {
          wt3 = 1.0/( weight - 1.0 );
          }
        
        IteratorType It1( avgimage1, avgimage1->GetLargestPossibleRegion() );
        IteratorType It2( image2, image2->GetLargestPossibleRegion() );
        IteratorType It3( varimage1, varimage1->GetLargestPossibleRegion() );
        IteratorType It4( ROIimg, ROIimg->GetLargestPossibleRegion() );
        
        It1.GoToBegin();
        It2.GoToBegin();
        It3.GoToBegin();
        It4.GoToBegin();
        while ( !It1.IsAtEnd() )
          {
          if ( It4.Get() != itk::NumericTraits<RealPixelType>::Zero )
            {
            float pix1 = static_cast<float>( It1.Get() );
            float pix2 = static_cast<float>( It2.Get() );
            if ( uselog )
              {
              if ( pix2 > 1.0e-11 )
                {
                pix2 = log( pix2 );
                }
              else
                {
                pix2 = 0.0;
                }
              }
            It1.Set( pix1*wt1 + pix2*wt2 );

            if ( ct1 > 0 )
              {
              float pix3 = static_cast<float>( It3.Get() );
              It3.Set( wt1*pix3 + wt3*( pix2 - pix1 )*( pix2 - pix1 ) );
              }
            }  
          ++It1;
          ++It2;
          ++It3;  
          ++It4;  
          }
        ct1++; 
        }
      }
  
    ct2 = 0;
    for ( unsigned int qq = 0; qq < ( filecount1+filecount2 ); qq++ )
      {
      if ( controlbool[qq] == true )
        {
        std::string filenma=filenames[qq];
        ReaderType::Pointer reader2 = ReaderType::New();
        reader2->SetFileName(filenma.c_str()); 
        reader2->UpdateLargestPossibleRegion();
        try
          {   
          image2 = reader2->GetOutput(); 
          }
        catch(...)
          {
          std::cout << " Error reading " << std::string(filenm) << std::endl;
          return 0;
          } 

        float weight = static_cast<float>( ct2 + 1 );  
        float wt2 = 1.0/weight;
        float wt1 = 1.0 - wt2;
        float wt3;

        if ( ct2 > 0 )
          {
          wt3 = 1.0/( weight - 1.0 );
          }
        
        IteratorType It1( avgimage2, avgimage2->GetLargestPossibleRegion() );
        IteratorType It2( image2, image2->GetLargestPossibleRegion() );
        IteratorType It3( varimage2, varimage2->GetLargestPossibleRegion() );
        IteratorType It4( ROIimg, ROIimg->GetLargestPossibleRegion() );
        
        It1.GoToBegin();
        It2.GoToBegin();
        It3.GoToBegin();
        It4.GoToBegin();
        while ( !It1.IsAtEnd() )
          {
          if ( It4.Get() != itk::NumericTraits<RealPixelType>::Zero )
            {
            float pix1 = static_cast<float>( It1.Get() );
            float pix2 = static_cast<float>( It2.Get() );
            if ( uselog )
              {
              if ( pix2 > 1.0e-11 )
                {
                pix2 = log( pix2 );
                }
              else
                {
                pix2 = 0.0;
                }
              }
            It1.Set( pix1*wt1 + pix2*wt2 );

            if ( ct1 > 0 )
              {
              float pix3 = static_cast<float>( It3.Get() );
              It3.Set( wt1*pix3 + wt3*( pix2 - pix1 )*( pix2 - pix1 ) );
              }
            }  
          ++It1;
          ++It2;
          ++It3;  
          ++It4;  
          }
        ct2++; 
        }
      }  

     // now compute the timages 
    //   first, scale variances & smooth, if so desired
      //  std::cout << " smooth begin " << std::endl;
    if ( smoothvar1 > 0.0 ) 
      {
      varimage1 = SmoothImage<ImageType>( varimage1, smoothvar1 );
      }
    if ( smoothvar2 > 0.0 ) 
      {
      varimage2 = SmoothImage<ImageType>( varimage2, smoothvar2 );
      }

    RealPixelType n1 = static_cast<RealPixelType>( filecount1 );
    RealPixelType n2 = static_cast<RealPixelType>( filecount2 );

    IteratorType It1( avgimage1, avgimage1->GetLargestPossibleRegion() );
    IteratorType It2( avgimage2, avgimage2->GetLargestPossibleRegion() );
    IteratorType It3( varimage1, varimage1->GetLargestPossibleRegion() );
    IteratorType It4( varimage2, varimage2->GetLargestPossibleRegion() );
    IteratorType It5( ROIimg, ROIimg->GetLargestPossibleRegion() );
    IteratorType It6( ttestimg, ttestimg->GetLargestPossibleRegion() );
    
    It1.GoToBegin();
    It2.GoToBegin();
    It3.GoToBegin();
    It4.GoToBegin();
    It5.GoToBegin();
    It6.GoToBegin();
    while ( !It1.IsAtEnd() )
      {
      if ( It5.Get() != itk::NumericTraits<RealPixelType>::Zero )
        {
        RealPixelType den = sqrt( It3.Get()/n1 + It4.Get()/n2 );
        if ( den != itk::NumericTraits<RealPixelType>::Zero )
          {
          It6.Set( ( It1.Get() - It2.Get() ) / den);
          }
        else
          {
          It6.Set( static_cast<RealPixelType>( 0 ) );
          }
        }  
      else
        {
        It6.Set( static_cast<RealPixelType>( 0 ) );
        }
      ++It1;
      ++It2;
      ++It3;  
      ++It4;  
      ++It5;  
      ++It6;  
      }
      
    std::cout << " t-test end " << std::endl;

    fakeclustersizes = GetClusterStat<ImageType>( ttestimg, Tthreshold, ClustThresh, whichstat, outfn, false );
     
    unsigned int csz = fakeclustersizes;
    if ( csz > histogramofsizes.size() - 1 ) 
      {
      csz = histogramofsizes.size() - 1;
      }
    for ( int qq = 0; qq <= csz; qq++ ) 
      {
      histogramofsizes[qq] += 1;
      }
    std::cout  << "P#: " << permct << " prob " << static_cast<float>( histogramofsizes[clustersizes] )/static_cast<float>( permct ) << std::endl;
    }  
    
  if ( NPermutations > 0 )
    {
    std::cout << std::endl;
    std::cout << " PERMUTATIONS DONE " << std::endl;
    std::cout << " Permutation Results: " << std::endl;
    std::cout << std::endl;

    for (unsigned int qq = static_cast<unsigned int>( ClustThresh ); qq < histogramofsizes.size(); qq++ )
      {
      float prob = static_cast<float>( histogramofsizes[qq] )/static_cast<float>( NPermutations );
      std::cout << " size " << qq << " ct " << histogramofsizes[qq] << " prob " << prob <<  std::endl;
      }

    ImageType::Pointer clusts;
    // now read the Cluster image and relabel it as a probability image.
    std::string tfn=outfn+"Clusters.hdr";
    ReaderType::Pointer reader2 = ReaderType::New();
    reader2->SetFileName(tfn.c_str()); 
    reader2->UpdateLargestPossibleRegion();
    try
      {   
      clusts = reader2->GetOutput(); 
      }
    catch(...)
      {
      std::cout << " Error reading ROI image " << std::endl;
      return 0;
      } 

    IteratorType It( clusts,  clusts->GetLargestPossibleRegion() ); 
    IteratorType Itt( ROIimg, ROIimg->GetLargestPossibleRegion() );
    for(  It.GoToBegin(), Itt.GoToBegin(); !It.IsAtEnd(); ++It, ++Itt)
      {
      if ( It.Get() > 0 && Itt.Get() != itk::NumericTraits<RealPixelType>::Zero )
        {
        float prob = static_cast<float>( histogramofsizes[static_cast<int>( It.Get() )] )
                   / static_cast<float> ( NPermutations );
        It.Set( 1.0 - prob );
        }
      }
    writer->SetFileName(  (outfn+std::string("OneMinusPval.hdr")).c_str());
    writer->SetInput( clusts ); 
    writer->Write();   
    }
  return 0;
 
}     


      

