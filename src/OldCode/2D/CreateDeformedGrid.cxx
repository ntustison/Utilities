#include "itkVector.h"
#include "itkVectorImageFileReader.h"
#include "itkVectorLinearInterpolateImageFunction.h"

#include "global.h"
#include "fstream.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " deformationField outputGridFile.txt " 
              << "[numberOfGridLinesX] [numberOfGridLinesY] [numberOfGridLinesZ]" << std::endl;
    exit( 1 );
    }

  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> VectorImageType;

  /**
   * Read in vector field
   */
  typedef itk::VectorImageFileReader<RealImageType, VectorImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::VectorLinearInterpolateImageFunction<VectorImageType, RealType> InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( reader->GetOutput() );
 
  ofstream str( argv[2] );

  str << "0 0 0 0" << std::endl;

  VectorImageType::PointType maxBoundary;
  VectorImageType::SpacingType gridSpacing;
  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    maxBoundary[d] = ( reader->GetOutput()->GetLargestPossibleRegion().GetSize()[d] - 1 ) 
      * reader->GetOutput()->GetSpacing()[d] + reader->GetOutput()->GetOrigin()[d];

    unsigned int numberOfGridLines = 10;
    if ( static_cast<unsigned>( argc ) > 3 + d )
      {
      numberOfGridLines = static_cast<unsigned int>( atoi( argv[d+3] ) );
      }  

    gridSpacing[d] = ( maxBoundary[d] - reader->GetOutput()->GetOrigin()[d] ) 
      / static_cast<RealType>( numberOfGridLines - 1 ) - 0.0001;
    }

  unsigned int lineCount = 0;
  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    RealType delta = reader->GetOutput()->GetSpacing()[d];

    InterpolatorType::PointType point; 

    for ( unsigned int c = 0; c < ImageDimension; c++ )
      {
      point[c] = reader->GetOutput()->GetOrigin()[c];
      }

    while ( true )
      {      
      InterpolatorType::OutputType vector = interpolator->Evaluate( point );  

      str << point[0] + vector[0] << " " << point[1] + vector[1];    
      if ( ImageDimension == 2 )
        {
        str << " 0 " << lineCount+1 << std::endl;
        }
      else
        {  
        str << " " << point[2] + vector[2] << " " << lineCount+1 << std::endl;
        }
        
      point[d] += delta;

      unsigned int doneFlag = 0;  
      if ( point[d] > maxBoundary[d] )
        {

        lineCount++;
        point[d] = reader->GetOutput()->GetOrigin()[d];
        doneFlag++;

        for ( unsigned int c = 0; c < ImageDimension-1; c++ )
          {
          point[(d+c+1)%ImageDimension] += gridSpacing[(d+c+1)%ImageDimension];
          if ( point[(d+c+1)%ImageDimension] > maxBoundary[(d+c+1)%ImageDimension] )
            {
            point[(d+c+1)%ImageDimension] = reader->GetOutput()->GetOrigin()[(d+c+1)%ImageDimension];  
            doneFlag++;
            }    
          else
            {
            break;
            }
          }  
        }
      if ( doneFlag == ImageDimension )
        {
        break;
        }     
      }   
    }    
  str << "0 0 0 0" << std::endl;
  str.close();  

  return 0;
}
