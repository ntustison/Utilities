#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPointSet.h"
#include "itkTimeProbe.h"
#include "itkVector.h"
#include "itkVectorImageFileWriter.h"
#include "itkVectorLinearInterpolateImageFunction.h"

#include <stdio.h>
#include <vector>
#include <fstream.h>

#include "vnl/vnl_vector.h"

#include "global.h"

int main( unsigned int argc, char *argv[] )
{
  if ( argc < 3*ImageDimension + 7 )
    {
    std::cout << "Usage: " << argv[0] << "  inputMovingLandmarkFile inputFixedLandmarkFile outputField "
              << "origin[0] ... origin[dimension-1] " << std::endl
              << "spacing[0] ... spacing[dimension-1] " << std::endl
              << "size[0] ... size[dimension-1] " << std::endl
              << "nlevels bsplineOrder ncps" << std::endl
              << "(optional: outputControlPointImage)" << std::endl
              << std::endl;
    exit( 1 );
    }

  typedef float RealType;

  RealType x;
  RealType y;
  RealType z;
  int l, lf_max, lm_max;

  std::vector<RealType> Xm;
  std::vector<RealType> Ym;
  std::vector<RealType> Zm;
  std::vector<int> Lm;
  ifstream strM( argv[1] );

  lm_max = 0;
  while ( strM >> x >> y >> z >> l )
    {
    Xm.push_back( x );
    Ym.push_back( y );
    Zm.push_back( z );
    Lm.push_back( l );
    if ( lm_max < l )
      {
      lm_max = l;
      } 
    }
  strM.close();

  std::vector<RealType> Xf;
  std::vector<RealType> Yf;
  std::vector<RealType> Zf;
  std::vector<int> Lf;
  ifstream strF( argv[2] );

  while ( strF >> x >> y >> z >> l )
    {
    Xf.push_back( x );
    Yf.push_back( y );
    Zf.push_back( z );
    Lf.push_back( l );
    if ( lf_max < l )
      {
      lf_max = l;
      } 
    }

  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;
  typedef itk::PointSet<VectorType, ImageDimension> DeformationFieldPointSetType;
  DeformationFieldPointSetType::Pointer fieldPoints = 
    DeformationFieldPointSetType::New();    
  typedef itk::BSplineScatteredDataPointSetToImageFilter
    <DeformationFieldPointSetType, DeformationFieldType> BSplineFilterType;
  BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

  unsigned int Npoints = 0;
  for ( unsigned int i = 1; i < Lm.size()-1; i++ )
    {
    RealType min_distance = itk::NumericTraits<RealType>::max();  
    int min_idx = -1;

    for ( unsigned int j = 1; j < Lf.size()-1; j++ )
      {
      if ( Lm[j] == Lm[i] )
        {
        RealType distance = ( Xm[i] - Xf[j] ) * ( Xm[i] - Xf[j] )
                          + ( Ym[i] - Yf[j] ) * ( Ym[i] - Yf[j] )
                          + ( Zm[i] - Zf[j] ) * ( Zm[i] - Zf[j] );
        if ( distance < min_distance )
          {
          min_distance = distance;
          min_idx = j;
          }
        }
      }
    if ( min_idx != -1 )
      {
      DeformationFieldPointSetType::PointType point;
      point[0] = Xm[i];
      point[1] = Ym[i];
      if ( ImageDimension == 3 )
        {
        point[2] = Zm[i];
        } 
      VectorType vector;
      vector[0] = Xf[min_idx] - Xm[i];
      vector[1] = Yf[min_idx] - Ym[i];
      if ( ImageDimension == 3 )
        {
        vector[2] = Zf[min_idx] - Zm[i];
        } 
      fieldPoints->SetPoint( Npoints, point );
      fieldPoints->SetPointData( Npoints, vector );
      Npoints++;
      }
    }             
  DeformationFieldType::PointType origin;
  DeformationFieldType::SpacingType spacing;
  DeformationFieldType::SizeType size;

  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin[i] = atof( argv[4+i] ); 
    spacing[i] = atof( argv[4+ImageDimension+i] );
    size[i] = atoi( argv[4+2*ImageDimension+i] );
    } 

  BSplineFilterType::ArrayType close;  
  close.Fill( false );
  BSplineFilterType::ArrayType ncps;
  ncps.Fill( atoi( argv[4+ImageDimension*3+2] ) );

  bspliner->SetOrigin( origin );
  bspliner->SetSpacing( spacing );
  bspliner->SetSize( size );
  bspliner->SetCloseDimension( close );
  bspliner->SetGenerateOutputImage( true );
  bspliner->SetNumberOfLevels( atoi( argv[4+ImageDimension*3] ) );
  bspliner->SetSplineOrder( atoi( argv[4+ImageDimension*3+1] ) );
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
  writer->SetFileName( argv[3] );
  writer->SetUseAvantsNamingConvention( true );
  writer->SetInput( bspliner->GetOutput() );
  writer->Update();

  return 0;
}
