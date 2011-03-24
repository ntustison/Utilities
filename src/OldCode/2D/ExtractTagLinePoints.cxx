
#include "itkGaborFilterBankTagging2DImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkTimeProbe.h"

#include <string>
#include <fstream.h>

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cout << "Usage: " << argv[0] << " imageFile maskFile outputFilePrefix tagSpacing "
              << "[thresholdPercentage] [numberOfAngleSteps] [numberOfTagSpacingSteps] [angleOffset] [tagSpacingFactor]" << std::endl;     
    exit( 0 );
    }   

  itk::TimeProbe timer;
  timer.Start();

  typedef float RealType;

  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<LabelImageType> MaskReaderType;
  MaskReaderType::Pointer maskreader = MaskReaderType::New();
  maskreader->SetFileName( argv[2] );
  maskreader->Update();

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  /** 
   * do each of the three sets of tag planes 
   */
  for ( unsigned int d = 0; d < ImageDimension; d++ )  
    {  
    std::cout << "Analyzing tag plane set " << d << std::endl; 

    /**
     * Calculate initial segmentation with Gabor filter bank
     */
    
    itk::TimeProbe timer1;
    timer1.Start();
    std::cout << "      Calculating initial segmentation with Gabor filter bank.   " << std::flush;    
    typedef itk::GaborFilterBankTagging2DImageFilter<RealImageType, LabelImageType> GaborFilterType;
    GaborFilterType::Pointer gaborFilter = GaborFilterType::New();
//    gaborFilter->DebugOn();

    gaborFilter->SetInput( reader->GetOutput() );
    gaborFilter->SetMaskImage( maskreader->GetOutput() );
    if ( argc > 5 )
      {
      if ( atof( argv[5] ) < 0.0 )
        {
        gaborFilter->SetUseAutomaticThresholding( true );
        } 
      else
        {
        gaborFilter->SetThresholdPercentage( atof( argv[5] ) );
        gaborFilter->SetUseAutomaticThresholding( false );
        } 
      }  
    else
      { 
      gaborFilter->SetUseAutomaticThresholding( true );
//      gaborFilter->SetThresholdPercentage( 0.9 );
      }  

    RealType tagSpacingFactor = 0.25;
    if ( argc > 9 )
      {
      tagSpacingFactor = atof( argv[9] );
      } 

    gaborFilter->SetTagSpacingMinimum( ( 1.0 - tagSpacingFactor ) * atof( argv[4] ) );
    gaborFilter->SetTagSpacingMaximum( ( 1.0 + tagSpacingFactor ) * atof( argv[4] ) );

    if ( argc > 7 )
      {
      gaborFilter->SetNumberOfTagSpacingSteps( atoi( argv[7] ) );
      }  
    else
      {
      gaborFilter->SetNumberOfTagSpacingSteps( 3 );
      }  
  
    RealType angleOffset = 10.0 * vnl_math::pi/180.0;
    if ( argc > 8 )
      {
      angleOffset = atof( argv[8] ) * vnl_math::pi/180.0;;
      }  

    /**
     * Rotation[0, 1, 2] = [theta, psi, phi] = [rot_x, rot_y, rot_z]
     * where imaging plane is parallel to x-y plane
     */  

    GaborFilterType::ArrayType beginRotation, endRotation;
    GaborFilterType::UnsignedIntArrayType samples;
    unsigned int numberOfSamples = 3;  
    if ( argc > 6 )
      {
      numberOfSamples = atoi( argv[6] );
      }  

    switch ( d )
      {
      case 0:
        beginRotation[0] = 0 - angleOffset;
        endRotation[0] = 0 + angleOffset;
        beginRotation[1] = 0;
        endRotation[1] = 0;
        samples[0] = numberOfSamples;
        samples[1] = 1;
        break;
      case 1:
        beginRotation[0] = 0.5 * vnl_math::pi - angleOffset;
        endRotation[0] = 0.5 * vnl_math::pi + angleOffset;
        beginRotation[1] = 0;
        endRotation[1] = 0;
        samples[0] = numberOfSamples;
        samples[1] = 1;
        break;
      } 

    gaborFilter->SetRotationAngleMinimum( beginRotation );
    gaborFilter->SetRotationAngleMaximum( endRotation );
    gaborFilter->SetNumberOfRotationAngleSteps( samples );
    
    gaborFilter->Update();
    timer1.Stop();
    std::cout << "Done. (" << timer1.GetMeanTime() << ")" << std::endl;

    itk::OStringStream buf_d;
    buf_d << d;
    
    std::string filename;

    filename = std::string( argv[3] ) + std::string( "GaborPoints." )
      + buf_d.str() + std::string( ".txt" );
 
    ofstream str( filename.c_str() );
    str << "0 0 0 0" << std::endl;

    itk::ImageRegionIteratorWithIndex<LabelImageType> ItL( gaborFilter->GetOutput(),
      gaborFilter->GetOutput()->GetLargestPossibleRegion() );
    for ( ItL.GoToBegin(); !ItL.IsAtEnd(); ++ItL )
      {
      if ( ItL.Get() > 0 )
        {
        LabelImageType::PointType point;
        gaborFilter->GetOutput()->TransformIndexToPhysicalPoint( ItL.GetIndex(), point );
        if ( ImageDimension == 2 )
          {
          str << point[0] << " " << point[1] << " 0 " << ItL.Get() << std::endl;   
          } 
        else if ( ImageDimension == 3 )
          {
          str << point[0] << " " << point[1] << " " << point[2] << " " << ItL.Get() << std::endl;   
          } 
        } 
      }  

    str << "0 0 0 0" << std::endl;
    str.close();


    filename = std::string( argv[3] ) + std::string( "GaborMaximalResponse." )
      + buf_d.str() + std::string( ".nii.gz" );
  
    typedef itk::ImageFileWriter<GaborFilterType::RealImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( filename.c_str() );
    writer->SetInput( gaborFilter->GetMaximumResponseImage() );
    writer->Update();

    filename = std::string( argv[3] ) + std::string( "GaborLabel." )
      + buf_d.str() + std::string( ".nii.gz" );

    typedef itk::ImageFileWriter<LabelImageType> LabelWriterType;
    LabelWriterType::Pointer labelwriter = LabelWriterType::New();
    labelwriter->SetFileName( filename.c_str() );
    labelwriter->SetInput( gaborFilter->GetOutput() );
    labelwriter->Update();
  

    }       

}
