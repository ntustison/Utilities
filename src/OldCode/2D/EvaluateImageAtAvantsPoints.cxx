#include "itkImageFileReader.h"
#include "itkLinearInterpolateImageFunction.h"

#include <stdio.h>
#include <fstream.h>

#include "vnl/vnl_vector.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << "  inputImage inputLandmarkFile outputFile radius " << std::endl;
    exit( 1 );
    }

  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();  

  typedef itk::LinearInterpolateImageFunction<ImageType> InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( reader->GetOutput() );
  
  InterpolatorType::PointType delta;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    delta[i] = reader->GetOutput()->GetSpacing()[i];
    } 

  unsigned int radius = atoi( argv[4] );

  ifstream str( argv[2] );
  ofstream strt( argv[3] );

  strt << "0 0 0 0 0" << std::endl;

  RealType totalMetric = 0.0;
  RealType NN = 0.0;

  RealType x[4];
  while ( str >> x[0] >> x[1] >> x[2] >> x[3] )
    {
    InterpolatorType::PointType start;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      start[i] = x[i] - radius*delta[i];
      }   
  
    RealType N = 0.0;
    RealType metric = 0.0;

    InterpolatorType::PointType point;
    
    if ( ImageDimension == 2 )
      {
      for ( unsigned int i = 0; i < 2*radius+1; i++ )
        { 
        for ( unsigned int j = 0; j < 2*radius+1; j++ )
          {
          point[0] = start[0] + i*delta[0];
          point[1] = start[1] + j*delta[1];

          if ( interpolator->IsInsideBuffer( point ) )
            {
            N++;
            metric += interpolator->Evaluate( point );   
            } 
          }
        }
      }        
    else if ( ImageDimension == 3 )
      {
      for ( unsigned int i = 0; i < 2*radius+1; i++ )
        { 
        for ( unsigned int j = 0; j < 2*radius+1; j++ )
          {
          for ( unsigned int k = 0; k < 2*radius+1; k++ )
            {
            point[0] = start[0] + i*delta[0];
            point[1] = start[1] + j*delta[1];
            point[2] = start[2] + k*delta[2];

            if ( interpolator->IsInsideBuffer( point ) )
              {
              N++;
              metric += interpolator->Evaluate( point );   
              }
            }   
          }
        }
      }
    strt << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << ' ' << ( metric/N ) << std::endl;
    totalMetric += ( metric/N );            
    NN++;     
    }
  strt << "0 0 0 0 0" << std::endl;

  std::cout << "Average metric value = " << totalMetric/NN << std::endl;


  return 0;
}
