// #include "itkConstNeighborhoodIterator.h"
// #include "itkCrossCorrelationRegistrationFunction.h"
// #include "itkDisplacementFieldTransform.h"
// #include "itkFiniteDifferenceFunction.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImportImageFilter.h"
#include "itkVector.h"
//
// #include "itkANTSNeighborhoodCorrelationImageToImageMetricv4.h"

#include <fstream>


template <unsigned int ImageDimension>
int Test( unsigned int argc, char *argv[] )
{
  typedef double RealType;

  typedef itk::Image<RealType, ImageDimension> ImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;

  std::ofstream str( "helix.txt" );

  for( float t = 0.0; t < 5 * 2 * vnl_math::pi; t+=0.1 )
    {
//     float x = vcl_cos( t );
//     float y = vcl_sin( t );
//     float z = t;
    float x = 5*t;
    float y = 3*t;
    float z = t;
    str << x << " " << y << " " << z << " 1" << std::endl;
    }
  str.close();


//
//   typedef itk::DisplacementFieldTransform<RealType, ImageDimension> DisplacementFieldTransformType;
//
//   typedef itk::ImageFileReader<ImageType> ReaderType;
//
//   typename ReaderType::Pointer fixedreader = ReaderType::New();
//   fixedreader->SetFileName( argv[2] );
//   fixedreader->Update();
//
//   typename ReaderType::Pointer movingreader = ReaderType::New();
//   movingreader->SetFileName( argv[3] );
//   movingreader->Update();
//
//   VectorType zeroVector( 0.0 );
//
//   //////////////////////////////////////////////////////////////////
//   // ITKv4 metric
//   //////////////////////////////////////////////////////////////////
//
//   typename DisplacementFieldType::Pointer identityField = DisplacementFieldType::New();
//   identityField->CopyInformation( fixedreader->GetOutput() );
//   identityField->SetRegions( fixedreader->GetOutput()->GetRequestedRegion() );
//   identityField->Allocate();
//   identityField->FillBuffer( zeroVector );
//
//   typename DisplacementFieldTransformType::Pointer identityDisplacementFieldTransform = DisplacementFieldTransformType::New();
//   identityDisplacementFieldTransform->SetDisplacementField( identityField );
//
//   typedef itk::ANTSNeighborhoodCorrelationImageToImageMetricv4<ImageType, ImageType> CorrelationMetricType;
//   typename CorrelationMetricType::Pointer correlationMetric = CorrelationMetricType::New();
//   typename CorrelationMetricType::RadiusType radius;
//   radius.Fill( 4 );
//   correlationMetric->SetRadius( radius );
//   correlationMetric->SetFixedImage( fixedreader->GetOutput() );
//   correlationMetric->SetFixedTransform( identityDisplacementFieldTransform );
//   correlationMetric->SetMovingImage( movingreader->GetOutput() );
//   correlationMetric->SetMovingTransform( identityDisplacementFieldTransform );
//   correlationMetric->SetVirtualDomainImage( fixedreader->GetOutput() );
//   correlationMetric->SetUseMovingImageGradientFilter( false );
//   correlationMetric->SetUseFixedImageGradientFilter( false );
//   correlationMetric->Initialize();
//
//   RealType value;
//   typename CorrelationMetricType::DerivativeType metricDerivative;
//   correlationMetric->GetValueAndDerivative( value, metricDerivative );
//
//   const unsigned long numberOfPixels = static_cast<unsigned long>( metricDerivative.Size() / ImageDimension );
//   const bool importFilterWillReleaseMemory = false;
//
//   VectorType *metricDerivativeFieldPointer = reinterpret_cast<VectorType *>( metricDerivative.data_block() );
//
//   typedef itk::ImportImageFilter<VectorType, ImageDimension> ImporterType;
//   typename ImporterType::Pointer importer = ImporterType::New();
//   importer->SetImportPointer( metricDerivativeFieldPointer, numberOfPixels, importFilterWillReleaseMemory );
//   importer->SetRegion( fixedreader->GetOutput()->GetBufferedRegion() );
//   importer->SetOrigin( fixedreader->GetOutput()->GetOrigin() );
//   importer->SetSpacing( fixedreader->GetOutput()->GetSpacing() );
//   importer->SetDirection( fixedreader->GetOutput()->GetDirection() );
//   importer->Update();
//
//   std::cout << "ITKv4 value: " << value << std::endl;
//
//   {
//   typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
//   typename WriterType::Pointer writer = WriterType::New();
//   writer->SetFileName( "itkv4.nii.gz" );
//   writer->SetInput( importer->GetOutput() );
//   writer->Update();
//   }
//
//   //////////////////////////////////////////////////////////////////
//   // ants metric
//   //////////////////////////////////////////////////////////////////
//
//   typedef itk::CrossCorrelationRegistrationFunction<ImageType, ImageType,
//                                                     DisplacementFieldType>         CCMetricType;
//   typename CCMetricType::Pointer ccmet = CCMetricType::New();
//
//   typename DisplacementFieldType::Pointer updateField = DisplacementFieldType::New();
//   updateField->CopyInformation( fixedreader->GetOutput() );
//   updateField->SetRegions( fixedreader->GetOutput()->GetLargestPossibleRegion() );
//   updateField->Allocate();
//   updateField->FillBuffer( zeroVector );
//
//   void * globalData = ccmet->GetGlobalDataPointer();
//
//   ccmet->SetFixedImage( fixedreader->GetOutput() );
//   ccmet->SetMovingImage( movingreader->GetOutput() );
//   ccmet->SetRadius( radius );
//   ccmet->SetGradientStep( 1.e2 );
//   ccmet->SetNormalizeGradient( false );
//   ccmet->InitializeIteration();
//
//
//   double ccval = 0;
//   typedef itk::ConstNeighborhoodIterator<ImageType> ScanIteratorType;
//   typedef itk::NeighborhoodIterator<DisplacementFieldType> NIteratorType;
//   typename ImageType::RegionType region = fixedreader->GetOutput()->GetLargestPossibleRegion();
//   ScanIteratorType asamIt( radius, fixedreader->GetOutput(), region);
//   unsigned long    ct = 0;
//
//   typedef itk::FiniteDifferenceFunction<DisplacementFieldType>     FiniteDifferenceFunctionType;
//   typedef typename FiniteDifferenceFunctionType::NeighborhoodType
//   NeighborhoodIteratorType;
//   NeighborhoodIteratorType nD(radius, updateField, region);
//   NIteratorType nD2(radius, updateField, region);
//
//
//   double metricvalue = 0.0;
//
//   itk::ImageRegionIteratorWithIndex<ImageType> iter( fixedreader->GetOutput(), region );
//   for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
//     {
//     typename ImageType::IndexType index = iter.GetIndex();
//     double    val = 0;
//     double    fmip = 0, mmip = 0, ffip = 0;
//     asamIt.SetLocation(index);
//     nD.SetLocation( index );
//     VectorType temp = ccmet->ComputeUpdate(nD, globalData);
//     nD2.SetLocation( index );
//     nD2.SetCenterPixel( temp );
//
//     typename ImageType::PointType point;
//     fixedreader->GetOutput()->TransformIndexToPhysicalPoint( index, point );
//
//     for( unsigned int i = 0; i < asamIt.Size(); i++ )
//       {
//       typename ImageType::IndexType locind = asamIt.GetIndex(i);
//       if( region.IsInside( locind ) )
//         {
//         double f = fixedreader->GetOutput()->GetPixel(locind);
//         double m = movingreader->GetOutput()->GetPixel(locind);
//         fmip += (f * m);  mmip += (m * m);  ffip += (f * f);
//         }
//       }
//     double denom = mmip * ffip;
//     if( denom == 0 )
//       {
//       val = 1;
//       }
//     else
//       {
//       val = fmip / vcl_sqrt(denom);
//       }
//     ccval += val;
//     ct++;
//     }
//   if( ct >  0 )
//     {
//     metricvalue = ccval / ct;
//     }
//
//   std::cout << "ants value: " << ccmet->GetEnergy() / ct << std::endl;
//
//   {
//   typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
//   typename WriterType::Pointer writer = WriterType::New();
//   writer->SetFileName( "ants.nii.gz" );
//   writer->SetInput( updateField );
//   writer->Update();
//   }

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
//   if ( argc < 4 )
//     {
//     std::cout << argv[0] << " imageDimension fixedImage movingImage" << std::endl;
//     exit( 1 );
//     }

  switch( atoi( argv[1] ) )
   {
   case 2:
     Test<2>( argc, argv );
     break;
   case 3:
     Test<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


