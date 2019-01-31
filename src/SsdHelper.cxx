#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkExtractImageFilter.h"
#include "itkLabelGeometryImageFilter.h"

#include <iostream>

int SsdHelper( int argc, char * argv[] )
{

  typedef int LabelType;
  typedef itk::Image<LabelType, 3> LabelImageType;
  typedef itk::Image<LabelType, 2> LabelSliceType;
  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;

  LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( argv[1] );
  
  LabelImageType::Pointer labelImage = labelReader->GetOutput();
  labelImage->Update();
  labelImage->DisconnectPipeline();
  
  unsigned int whichAxis = 0;
  if( argc > 2 )
    {
    whichAxis = static_cast<unsigned int>( atoi( argv[2] ) );
    }

  LabelImageType::RegionType::SizeType size = 
    labelImage->GetRequestedRegion().GetSize();
  unsigned int numberOfSlices = size[whichAxis];
  LabelImageType::IndexType startIndex = 
    labelImage->GetRequestedRegion().GetIndex();

  LabelImageType::RegionType region;
  size[whichAxis] = 0;

  LabelImageType::IndexType index;
  index.Fill( 0 );

  std::cout << "slice,label,xmin,xmax,ymin,ymax" << std::endl;

  for( unsigned int i = 0; i < numberOfSlices; i++ )
    {
    index[whichAxis] = startIndex[whichAxis] + i;
    region.SetIndex( index );
    region.SetSize( size );

    typedef itk::ExtractImageFilter<LabelImageType, LabelSliceType> ExtracterType;
    ExtracterType::Pointer extracter = ExtracterType::New();
    extracter->SetInput( labelImage );
    extracter->SetExtractionRegion( region );
    extracter->SetDirectionCollapseToIdentity();

    typedef itk::LabelGeometryImageFilter<LabelSliceType> FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput( extracter->GetOutput() );
    filter->Update();

    FilterType::LabelsType allLabels = filter->GetLabels();
    FilterType::LabelsType::iterator allLabelsIt;

    for( allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end(); allLabelsIt++ )
      {
      if( *allLabelsIt == 0 )
        {
        continue;
        }
        
      FilterType::BoundingBoxType boundingBox = filter->GetBoundingBox( *allLabelsIt );

      if( boundingBox[0] == boundingBox[1] || boundingBox[2] == boundingBox[3] )
        {
        continue;
        }  

      std::cout << i << ',' << *allLabelsIt << ',' <<
        boundingBox[0] << ',' << boundingBox[1] << ',' <<
        boundingBox[2] << ',' << boundingBox[3] << std::endl;
      }
    }

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << argv[0] << " labelImage [whichAxis=0]"
      << std::endl;
    return EXIT_FAILURE;
    }
  SsdHelper( argc, argv );
}

