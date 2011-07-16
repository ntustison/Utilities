#include "itkGaborFilterBankTaggingImageFilter.h"
#include "itkGaborFilterBankTagging2DImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkTimeProbe.h"

#include <string>
#include <fstream.h>

template <unsigned int ImageDimension>
int ExtractTagLinePoints3D( int argc, char *argv[] )
{
  itk::TimeProbe timer;
  timer.Start();

  typedef float RealType;

  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<LabelImageType> MaskReaderType;
  typename MaskReaderType::Pointer maskreader = MaskReaderType::New();
  maskreader->SetFileName( argv[3] );
  maskreader->Update();

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
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
    typedef itk::GaborFilterBankTaggingImageFilter<RealImageType, LabelImageType> GaborFilterType;
    typename GaborFilterType::Pointer gaborFilter = GaborFilterType::New();
//    gaborFilter->DebugOn();

    gaborFilter->SetInput( reader->GetOutput() );
    gaborFilter->SetMaskImage( maskreader->GetOutput() );
    if ( argc > 6 )
      {
      if ( atof( argv[6] ) < 0.0 )
        {
        gaborFilter->SetUseAutomaticThresholding( true );
        } 
      else
        {
        gaborFilter->SetThresholdPercentage( atof( argv[6] ) );
        gaborFilter->SetUseAutomaticThresholding( false );
        } 
      }  
    else
      { 
      gaborFilter->SetUseAutomaticThresholding( true );
//      gaborFilter->SetThresholdPercentage( 0.9 );
      }  

    RealType angleOffset = 10.0 * vnl_math::pi/180.0;
    if ( argc > 9 )
      {
      angleOffset = atof( argv[9] ) * vnl_math::pi/180.0;;
      }  

    RealType tagSpacingFactor = 0.25;
    if ( argc > 10 )
      {
      tagSpacingFactor = atof( argv[10] );
      } 

    gaborFilter->SetTagSpacingMinimum( ( 1.0 - tagSpacingFactor ) * atof( argv[5] ) );
    gaborFilter->SetTagSpacingMaximum( ( 1.0 + tagSpacingFactor ) * atof( argv[5] ) );

    if ( argc > 8 )
      {
      gaborFilter->SetNumberOfTagSpacingSteps( atoi( argv[8] ) );
      }  
    else
      {
      gaborFilter->SetNumberOfTagSpacingSteps( 3 );
      }  
  
    /**
     * Rotation[0, 1, 2] = [theta, psi, phi] = [rot_x, rot_y, rot_z]
     * where imaging plane is parallel to x-y plane
     */  

    typename GaborFilterType::ArrayType beginRotation;
    beginRotation.Fill( 0.0 );
    typename GaborFilterType::ArrayType endRotation;
    endRotation.Fill( 0.0 );
    typename GaborFilterType::UnsignedIntArrayType samples;
    samples.Fill( 0 );
    unsigned int numberOfSamples = 3;  
    if ( argc > 7 )
      {
      numberOfSamples = atoi( argv[7] );
      }  

    switch ( d )
      {
      case 0:
        beginRotation[0] = 0;
        endRotation[0] = 0;
        beginRotation[1] = 0 - angleOffset;
        endRotation[1] = 0 + angleOffset;
        beginRotation[2] = 0 - angleOffset;
        endRotation[2] = 0 + angleOffset;
        samples[0] = 1;
        samples[1] = numberOfSamples;
        samples[2] = numberOfSamples;
        break;
      case 1:
        beginRotation[0] = 0 - angleOffset;
        endRotation[0] = 0 + angleOffset;
        beginRotation[1] = 0;
        endRotation[1] = 0;
        beginRotation[2] = 0.5 * vnl_math::pi - angleOffset;
        endRotation[2] = 0.5 * vnl_math::pi + angleOffset;
        samples[0] = numberOfSamples;
        samples[1] = 1;
        samples[2] = numberOfSamples;
        break;
      case 2:
        beginRotation[0] = 0 - angleOffset;
        endRotation[0] = 0 + angleOffset;
        beginRotation[1] = 0.5 * vnl_math::pi - angleOffset;
        endRotation[1] = 0.5 * vnl_math::pi + angleOffset;
        beginRotation[2] = 0;
        endRotation[2] = 0;
        samples[0] = numberOfSamples;
        samples[1] = numberOfSamples;
        samples[2] = 1;
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

//    filename = std::string( argv[4] ) + std::string( "GaborPoints." )
//      + buf_d.str() + std::string( ".txt" );
// 
//    ofstream str( filename.c_str() );
//    str << "0 0 0 0" << std::endl;
//
//    itk::ImageRegionIteratorWithIndex<LabelImageType> ItL( gaborFilter->GetOutput(),
//      gaborFilter->GetOutput()->GetLargestPossibleRegion() );
//    for ( ItL.GoToBegin(); !ItL.IsAtEnd(); ++ItL )
//      {
//      if ( ItL.Get() > 0 )
//        {
//        typename LabelImageType::PointType point;
//        gaborFilter->GetOutput()->TransformIndexToPhysicalPoint( ItL.GetIndex(), point );
//        if ( ImageDimension == 2 )
//          {
//          str << point[0] << " " << point[1] << " 0 " << ItL.Get() << std::endl;   
//          } 
//        else if ( ImageDimension == 3 )
//          {
//          str << point[0] << " " << point[1] << " " << point[2] << " " << ItL.Get() << std::endl;   
//          } 
//        } 
//      }  
//
//    str << "0 0 0 0" << std::endl;
//    str.close();

    filename = std::string( argv[4] ) + std::string( "GaborMaximalResponse." )
      + buf_d.str() + std::string( ".nii.gz" );
  
    typedef itk::ImageFileWriter<typename GaborFilterType::RealImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( filename.c_str() );
    writer->SetInput( gaborFilter->GetMaximumResponseImage() );
    writer->Update();

//    filename = std::string( argv[4] ) + std::string( "GaborLabel." )
//      + buf_d.str() + std::string( ".nii.gz" );
//
//    typedef itk::ImageFileWriter<LabelImageType> LabelWriterType;
//    typename LabelWriterType::Pointer labelwriter = LabelWriterType::New();
//    labelwriter->SetFileName( filename.c_str() );
//    labelwriter->SetInput( gaborFilter->GetOutput() );
//    labelwriter->Update();
//
    }       
  return 0;
}

template <unsigned int ImageDimension>
int ExtractTagLinePoints2D( int argc, char *argv[] )
{
  itk::TimeProbe timer;
  timer.Start();

  typedef float RealType;

  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<LabelImageType> MaskReaderType;
  typename MaskReaderType::Pointer maskreader = MaskReaderType::New();
  maskreader->SetFileName( argv[3] );
  maskreader->Update();

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
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
    typename GaborFilterType::Pointer gaborFilter = GaborFilterType::New();
//    gaborFilter->DebugOn();

    gaborFilter->SetInput( reader->GetOutput() );
    gaborFilter->SetMaskImage( maskreader->GetOutput() );
    if ( argc > 6 )
      {
      if ( atof( argv[6] ) < 0.0 )
        {
        gaborFilter->SetUseAutomaticThresholding( true );
        } 
      else
        {
        gaborFilter->SetThresholdPercentage( atof( argv[6] ) );
        gaborFilter->SetUseAutomaticThresholding( false );
        } 
      }  
    else
      { 
      gaborFilter->SetUseAutomaticThresholding( true );
//      gaborFilter->SetThresholdPercentage( 0.9 );
      }  

    RealType tagSpacingFactor = 0.25;
    if ( argc > 10 )
      {
      tagSpacingFactor = atof( argv[10] );
      } 

    gaborFilter->SetTagSpacingMinimum( ( 1.0 - tagSpacingFactor ) * atof( argv[5] ) );
    gaborFilter->SetTagSpacingMaximum( ( 1.0 + tagSpacingFactor ) * atof( argv[5] ) );

    if ( argc > 8 )
      {
      gaborFilter->SetNumberOfTagSpacingSteps( atoi( argv[8] ) );
      }  
    else
      {
      gaborFilter->SetNumberOfTagSpacingSteps( 3 );
      }  
  
    RealType angleOffset = 10.0 * vnl_math::pi/180.0;
    if ( argc > 9 )
      {
      angleOffset = atof( argv[9] ) * vnl_math::pi/180.0;;
      }  

    /**
     * Rotation[0, 1, 2] = [theta, psi, phi] = [rot_x, rot_y, rot_z]
     * where imaging plane is parallel to x-y plane
     */  

    typename GaborFilterType::ArrayType beginRotation, endRotation;
    typename GaborFilterType::UnsignedIntArrayType samples;
    unsigned int numberOfSamples = 3;  
    if ( argc > 7 )
      {
      numberOfSamples = atoi( argv[7] );
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

//    filename = std::string( argv[4] ) + std::string( "GaborPoints." )
//      + buf_d.str() + std::string( ".txt" );
// 
//    ofstream str( filename.c_str() );
//    str << "0 0 0 0" << std::endl;
//
//    itk::ImageRegionIteratorWithIndex<LabelImageType> ItL( gaborFilter->GetOutput(),
//      gaborFilter->GetOutput()->GetLargestPossibleRegion() );
//    for ( ItL.GoToBegin(); !ItL.IsAtEnd(); ++ItL )
//      {
//      if ( ItL.Get() > 0 )
//        {
//        typename LabelImageType::PointType point;
//        gaborFilter->GetOutput()->TransformIndexToPhysicalPoint( ItL.GetIndex(), point );
//        if ( ImageDimension == 2 )
//          {
//          str << point[0] << " " << point[1] << " 0 " << ItL.Get() << std::endl;   
//          } 
//        else if ( ImageDimension == 3 )
//          {
//          str << point[0] << " " << point[1] << " " << point[2] << " " << ItL.Get() << std::endl;   
//          } 
//        } 
//      }  
//
//    str << "0 0 0 0" << std::endl;
//    str.close();


    filename = std::string( argv[4] ) + std::string( "GaborMaximalResponse." )
      + buf_d.str() + std::string( ".nii.gz" );
  
    typedef itk::ImageFileWriter<typename GaborFilterType::RealImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( filename.c_str() );
    writer->SetInput( gaborFilter->GetMaximumResponseImage() );
    writer->Update();

//    filename = std::string( argv[4] ) + std::string( "GaborLabel." )
//      + buf_d.str() + std::string( ".nii.gz" );
//
//    typedef itk::ImageFileWriter<LabelImageType> LabelWriterType;
//    typename LabelWriterType::Pointer labelwriter = LabelWriterType::New();
//    labelwriter->SetFileName( filename.c_str() );
//    labelwriter->SetInput( gaborFilter->GetOutput() );
//    labelwriter->Update();
  

    }       
  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 7 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage maskImage outputPrefix tagSpacing "
              << "[thresholdPercentage] [numberOfAngleSteps] [numberOfTagSpacingSteps] [angleOffset] [tagSpacingFactor]" << std::endl;     
    exit( 0 );
    }   

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     ExtractTagLinePoints2D<2>( argc, argv );
     break;
   case 3:
     ExtractTagLinePoints3D<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

