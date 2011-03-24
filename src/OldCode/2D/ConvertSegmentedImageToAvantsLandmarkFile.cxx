#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "global.h"
#include "string.h"
#include <fstream.h>

#define isnan(x) ((x) != (x))

int main( unsigned int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " segmentedImage outputFile [label]" << std::endl;
    exit( 1 );
    }
  
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  unsigned long count = 1;

  itk::ImageRegionIteratorWithIndex<ImageType> It
    ( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );

  ofstream str( argv[2] );
  
  str << "0 0 0 0" << std::endl;
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    int label = It.Get();
    if ( ( argc >= 4 && label == atoi( argv[3] ) ) || ( argc == 3 && label > 0 ) )
      {
      ImageType::PointType point;
      reader->GetOutput()->TransformIndexToPhysicalPoint( It.GetIndex(), point );
      str << point[0] << " " << point[1] << " ";
      if ( ImageDimension == 3 )
        {
        str << point[2] << " ";
        }
      if ( argc == 3 )
        {
        str << label << std::endl;
        } 
      else
        {
        str << count << std::endl;
        } 
      count++;
      }   
    }    
  str << "0 0 0 0" << std::endl;


/*  
  ImageType::PointType point;
  point.Fill( 0.0 );
  std::vector<ImageType::PointType> points( stats->GetMaximum()+1, point );
  std::vector<unsigned int> N( stats->GetMaximum()+1, 0 );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( It.Get() > 0 )
      {
      ImageType::PointType pt1 = points[It.Get()];
      ImageType::PointType pt2;
      reader->GetOutput()->TransformIndexToPhysicalPoint( It.GetIndex(), pt2 );
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        pt1[i] += pt2[i];
        } 
      points[It.Get()] = pt1;
      N[It.Get()]++;          
      }
    }      
    
  std::string filename = std::string( argv[2] ) + std::string( ".txt" );
  ofstream str( filename.c_str() );
  
  str << "0 0 0 0" << std::endl;
  for ( unsigned int i = 1; i <= stats->GetMaximum(); i++ )
    {
    ImageType::PointType pt = points[i];
    for ( unsigned int j = 0; j < ImageDimension; j++ )
      {
      pt[j] /= static_cast<RealType>( N[i] );
      } 
    str << pt[0] << " " << pt[1] << " " << std::flush;
    if ( ImageDimension == 2 )
      {
      str << 0 << " " << i << std::endl;
      }
    if ( ImageDimension == 3 )
      {
      str << pt[2] << " " << i << std::endl;
      }
    }
*/  
  return 0;
}
