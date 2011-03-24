#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include "global.h"
#include "string.h"
#include "fstream.h"

int main( int argc, char *argv[] )
{
  if ( argc != 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage outputImage" << std::endl;
    exit( 1 );
    }
  
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  
  char extension[4] = "";
  int n = strlen( argv[1] );
  int idx = 0;
  for ( int i = n-3; i < n; i++ )
    {
    extension[idx++] = argv[1][i];
    }

  // Convert Marcelo's image file format to standard format
  if ( strcmp( extension, "img" ) == 0 ) 
    {
    
    fstream str( argv[1], std::ios::in );

    char line[101];
    for ( unsigned int i = 0; i < 4; i++ )
      {
      str.getline( line, 100 );
      }
    ImageType::SizeType size;  
    ImageType::SpacingType spacing;  
    ImageType::PointType origin;
    str >> line >> size[0] >> size[1] >> size[2];
    str >> line >> spacing[0] >> spacing[1] >> spacing[2];
    str >> line >> origin[0] >> origin[1] >> origin[2];
      
    str.getline( line, 100 );
    str.getline( line, 100 );
    char garbage[10] = "";
    char basura[10] = "";
    char type[20] = "";
    str >> garbage >> basura >> type;
    str.getline( line, 100 );
    
    ImageType::Pointer image = ImageType::New();
    ImageType::RegionType region;
    region.SetSize( size );
    image->SetOrigin( origin );
    image->SetSpacing( spacing );
    image->SetRegions( region );
    image->Allocate();
 
    itk::ImageRegionIterator<ImageType> It
      ( image, image->GetLargestPossibleRegion() );   

    for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if ( strcmp( type, "character" ) == 0 )
        {
        char value;
        str.read( (char*)&value, sizeof(char) );
        It.Set( static_cast<PixelType>( value ) ); 
        }
      if ( strcmp( type, "double" ) == 0 )
        {
        double value;
        str.read( (char*)&value, sizeof(double) );
        It.Set( static_cast<PixelType>( value ) ); 
        }
      }  
    str.close();  


    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[2] );
    writer->SetInput( image );
    writer->Update();
    }
  else
    {   
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[1] );
    reader->Update();
  
    fstream str( argv[2], std::ios::out );
    str << "# vtk DataFile Version 2.0" << std::endl;
    str << "LITTLE ENDIAN" << std::endl;
    str << "BINARY" << std::endl;
    str << "DATASET STRUCTURED_POINTS" << std::endl;
    str << "DIMENSIONS " 
        << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0] << " "    
        << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1] << " "    
        << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[2] << std::endl;
    str << "SPACING "    
        << reader->GetOutput()->GetSpacing()[0] << " "    
        << reader->GetOutput()->GetSpacing()[1] << " "    
        << reader->GetOutput()->GetSpacing()[2] << std::endl;
    str << "ORIGIN "    
        << reader->GetOutput()->GetOrigin()[0] << " "    
        << reader->GetOutput()->GetOrigin()[1] << " "    
        << reader->GetOutput()->GetOrigin()[2] << std::endl;
    str << "POINT_DATA "    
        << reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0] *   
           reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1] *    
           reader->GetOutput()->GetLargestPossibleRegion().GetSize()[2] << std::endl;
    str << "SCALARS data ";

    PixelType p;

    if ( strcmp( typeid( p ).name(), "c" ) == 0 )
      {
      str << "character" << std::endl;
      }
    else if ( strcmp( typeid( p ).name(), "d" ) == 0 )
      {
      str << "double" << std::endl;
      }

    itk::ImageRegionIterator<ImageType> It
      ( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );   

    for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      PixelType p = It.Get();
      str.write( (char*)&p, sizeof(PixelType) );   
      }
    str.close();   
    }
  return 0;
}
