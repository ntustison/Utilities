#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkContinuousIndex.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkPointSet.h"
#include "itkTimeProbe.h"
#include "itkVector.h"
#include "itkVectorImageFileReader.h"
#include "itkVectorImageFileWriter.h"
#include "itkVectorLinearInterpolateImageFunction.h"

#include <stdio.h>
#include <vector>

#include "vnl/vnl_vector.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] 
              << " inputDeformationField outputDeformationField [order] [nlevels] [maskImage]" << std::endl;
    exit( 1 );
    } 

  typedef float RealType;

  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Image<unsigned int, ImageDimension> MaskImageType;

  MaskImageType::Pointer mask = MaskImageType::New();
  mask = NULL;

  if ( argc > 5 )
    {
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName( argv[5] );
    maskReader->Update();
    mask = maskReader->GetOutput();
    }

  /**
   * Read in vector field
   */
  typedef itk::VectorImageFileReader<RealImageType, DeformationFieldType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::PointSet<VectorType, ImageDimension> DeformationFieldPointSetType;
  DeformationFieldPointSetType::Pointer fieldPoints = 
    DeformationFieldPointSetType::New();    
  unsigned long count = 0;

  itk::ImageRegionIteratorWithIndex<DeformationFieldType>
    It( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    VectorType vector = It.Get();
    DeformationFieldType::PointType point;
    reader->GetOutput()->TransformIndexToPhysicalPoint( It.GetIndex(), point );

    point += vector;
    itk::ContinuousIndex<double, ImageDimension> cidx;
    reader->GetOutput()->TransformPhysicalPointToContinuousIndex( point, cidx );

    if ( mask )
      {
      typedef itk::LinearInterpolateImageFunction<MaskImageType, double> InterpolatorType;
      InterpolatorType::Pointer interpolator = InterpolatorType::New();
      interpolator->SetInputImage( mask );
      if ( interpolator->EvaluateAtContinuousIndex( cidx ) > 0 )
        {
        continue;  
        }   
      }  

    if ( !reader->GetOutput()->GetLargestPossibleRegion().IsInside( cidx ) )
      {
      continue;
      } 

    DeformationFieldPointSetType::PointType fieldPoint;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      fieldPoint[i] = point[i];  
      }  
    fieldPoints->SetPoint( count, fieldPoint );
    fieldPoints->SetPointData( count, -vector );   
    count++; 
    }   
 
  unsigned int order = 3;
  if ( argc > 3 )
    {
    order = atoi( argv[3] );
    }
  unsigned int nlevels = 5;
  if ( argc > 4 )
    {
    nlevels = atoi( argv[4] );
    }
  
  typedef itk::BSplineScatteredDataPointSetToImageFilter
    <DeformationFieldPointSetType, DeformationFieldType> BSplineFilterType;
  BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

  BSplineFilterType::ArrayType close;  
  close.Fill( false );
  BSplineFilterType::ArrayType ncps;
  ncps.Fill( order+1 );

  //bspliner->DebugOn();
  bspliner->SetOrigin( reader->GetOutput()->GetOrigin() );
  bspliner->SetSpacing( reader->GetOutput()->GetSpacing() );
  bspliner->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  bspliner->SetCloseDimension( close );
  bspliner->SetGenerateOutputImage( true );
  bspliner->SetNumberOfLevels( nlevels );
  bspliner->SetSplineOrder( order );
  bspliner->SetNumberOfControlPoints( ncps );
  bspliner->SetInput( fieldPoints );

  itk::TimeProbe timer;
  timer.Start();  
  bspliner->Update();
  timer.Stop();
  std::cout << "  Timing results: " << timer.GetMeanTime() << std::endl;

  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::VectorImageFileWriter<DeformationFieldType, RealImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetUseAvantsNamingConvention( true );
  writer->SetInput( bspliner->GetOutput() );
  writer->Update();

  return 0;
}
