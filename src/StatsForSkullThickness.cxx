/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: StatsForSkullThickness.cxx,v $
  Language:  C++
  Date:      $Date: 2010/02/22 15:34:48 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkBresenhamLine.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVector.h"


#include <vector>
#include <algorithm>
#include <fstream>

template <unsigned int ImageDimension>
int StatsForSkullThickness( int argc, char* argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef typename ImageType::IndexType IndexType;

  typedef unsigned int LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::ImageFileReader<LabelImageType> MaskReaderType;
  typename MaskReaderType::Pointer maskReader = MaskReaderType::New();
  maskReader->SetFileName( argv[3] );
  maskReader->Update();

  float priorDistance = atof( argv[5] );

  std::vector<LabelType> labels;

  itk::ImageRegionIteratorWithIndex<LabelImageType> It( maskReader->GetOutput(),
    maskReader->GetOutput()->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( It.Get() > 2 &&
      std::find( labels.begin(), labels.end(), It.Get() ) == labels.end() )
      {
      labels.push_back( It.Get() );
      }
    }

  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( It.Get() > 2 )
      {
      typename std::vector<LabelType>::iterator it =
        std::find( labels.begin(), labels.end(), It.Get() );
      It.Set( it - labels.begin() );
      }
    }

  std::vector<IndexType> indices1;
  std::vector<IndexType> indices2;  indices2.resize( labels.size() );

  labels.clear();

  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( It.Get() > 2 )
      {
      std::vector<LabelType>::iterator it = std::find( labels.begin(), labels.end(), It.Get() );
      if( it == labels.end() )
        {
        labels.push_back( It.Get() );
        indices1.push_back( It.GetIndex() );
        }
      else
        {
        indices2[it - labels.begin()] = It.GetIndex();
        }
      }
    }

  std::string distanceFileName = std::string( argv[4] ) + std::string( "distance.csv" );
  std::string intensityFileName = std::string( argv[4] ) + std::string( "intensity.csv" );

  std::ofstream strD( distanceFileName.c_str() );
  std::ofstream strI( intensityFileName.c_str() );

  typedef itk::Vector<float, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;

  typename DisplacementFieldType::Pointer displacementField = DisplacementFieldType::New();
  displacementField = NULL;
  if( argc > 6 )
    {
    typedef itk::ImageFileReader<DisplacementFieldType> FieldReaderType;
    typename FieldReaderType::Pointer fieldreader = FieldReaderType::New();
    fieldreader->SetFileName( argv[6] );

    displacementField = fieldreader->GetOutput();
    displacementField->Update();
    displacementField->DisconnectPipeline();
    }

  typename std::vector<IndexType>::const_iterator it1;
  for( it1 = indices1.begin(); it1 != indices1.end(); ++it1 )
    {
    IndexType index1 = *it1;
    IndexType index2 = indices2[it1 - indices1.begin()];

    typename ImageType::PointType point1;
    maskReader->GetOutput()->TransformIndexToPhysicalPoint( index1, point1 );

    typename ImageType::PointType point2;
    maskReader->GetOutput()->TransformIndexToPhysicalPoint( index2, point2 );

    if( displacementField )
      {
      point1 += displacementField->GetPixel( index1 );
      point2 += displacementField->GetPixel( index2 );

      reader->GetOutput()->TransformPhysicalPointToIndex( point1, index1 );
      reader->GetOutput()->TransformPhysicalPointToIndex( point2, index2 );
      }

    float distance = point1.EuclideanDistanceTo( point2 );
    if( distance > priorDistance )
      {
      continue;
      }

    strD << labels[it1 - indices1.begin()] - 2 << "," << distance << std::endl;

    typedef itk::BresenhamLine<ImageDimension> LinerType;
    LinerType liner;
    typename LinerType::IndexArray indices = liner.BuildLine( index1, index2 );

    float averageProfile = 0.0;
    float N = 0.0;

    typename LinerType::IndexArray::const_iterator it;
    for( it = indices.begin(); it != indices.end(); ++it )
      {
      if( *it == index2 )
        {
        break;
        }
      PixelType intensity = reader->GetOutput()->GetPixel( *it );
      averageProfile += intensity;
      N++;
      }
    averageProfile /= N;

    strI << labels[it1 - indices1.begin()] - 2 << "," << averageProfile << std::endl;
    }

  strD.close();
  strI.close();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension inputImage maskImage outputPrefix priorDistance [displacementField]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     StatsForSkullThickness<2>( argc, argv );
     break;
   case 3:
     StatsForSkullThickness<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}
