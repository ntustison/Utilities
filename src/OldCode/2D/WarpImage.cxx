#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorImageFileReader.h"
#include "itkVectorImageFileWriter.h"
#include "itkVectorLinearInterpolateImageFunction.h"

#include "itkWarpImageFilter.h"
#include "global.h"


int main( unsigned int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: WarpImage input_filename output_filename deformationfield_1_filename deformationfield_2_filename ... deformationfield_n_filename " << std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  DeformationFieldType::Pointer totalField = DeformationFieldType::New();
  totalField->SetOrigin( reader->GetOutput()->GetOrigin() );
  totalField->SetSpacing( reader->GetOutput()->GetSpacing() );
  totalField->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  totalField->Allocate();
  VectorType V;
  V.Fill( 0 );
  totalField->FillBuffer( V );

  for ( unsigned int i = 3; i < argc; i++ )
    {
    typedef itk::VectorImageFileReader<RealImageType,
        DeformationFieldType> FieldReaderType;
    FieldReaderType::Pointer fieldreader = FieldReaderType::New();
    fieldreader->SetFileName( argv[i] );
    fieldreader->SetUseAvantsNamingConvention( true );
    fieldreader->Update();

    typedef itk::VectorLinearInterpolateImageFunction<DeformationFieldType, RealType>
      DeformationFieldInterpolatorType;
    DeformationFieldInterpolatorType::Pointer interpolator = DeformationFieldInterpolatorType::New();
    interpolator->SetInputImage( fieldreader->GetOutput() );

    itk::ImageRegionIteratorWithIndex<DeformationFieldType> It
      ( totalField, totalField->GetLargestPossibleRegion() );
    for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      VectorType V = It.Get();
      DeformationFieldType::PointType point;
      totalField->TransformIndexToPhysicalPoint( It.GetIndex(), point );

      DeformationFieldInterpolatorType::PointType pt;
      for ( unsigned int d = 0; d < ImageDimension; d++ )
        {
        pt[d] = point[d] + V[d];
        }

      bool isOutside = false;
      for ( unsigned int d = 0; d < ImageDimension; d++ )
        {
        if ( pt[d] < totalField->GetOrigin()[d] ||
             pt[d] > totalField->GetOrigin()[d] +
               ( totalField->GetLargestPossibleRegion().GetSize()[d]-1 )
                 *totalField->GetSpacing()[d] )
          {
          isOutside = true;
          break;
          }
        }
      if ( !isOutside )
        {
        DeformationFieldInterpolatorType::OutputType disp;
        disp = interpolator->Evaluate( pt );
        for ( unsigned int d = 0; d < ImageDimension; d++ )
          {
          V[d] += disp[d];
          }
        }
      It.Set( V );
      }
    }

  typedef itk::WarpImageFilter<ImageType, ImageType, DeformationFieldType> WarperType;
  WarperType::Pointer warper = WarperType::New();
  warper->SetInput( reader->GetOutput() );
  warper->SetDeformationField( totalField );
  warper->SetOutputSpacing( totalField->GetSpacing() );
  warper->SetOutputOrigin( totalField->GetOrigin() );

  warper->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( warper->GetOutput() );
  writer->Update();

  typedef itk::VectorImageFileWriter<
      DeformationFieldType, RealImageType> FieldWriterType;
  FieldWriterType::Pointer fieldwriter = FieldWriterType::New();
  fieldwriter->SetInput( totalField );
  fieldwriter->SetFileName( argv[2] );
  fieldwriter->SetUseAvantsNamingConvention( true );
  fieldwriter->Update();

  return 0;
}
