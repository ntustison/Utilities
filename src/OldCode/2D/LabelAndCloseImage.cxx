#include <stdio.h>

#include "itkAndImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkGridImageSource.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkVector.h"
#include "itkVectorContainer.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0] << " inputImage outputImage" << std::endl;
    exit( 0 ); 
    } 

  typedef itk::Image<unsigned int, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  typedef itk::ImageFileWriter<ImageType> WriterType;

  // Try to separate the lobes.
  
  typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentType;
  ConnectedComponentType::Pointer connecter = ConnectedComponentType::New();
  connecter->SetInput( reader->GetOutput() );
  typedef itk::RelabelComponentImageFilter<ImageType, ImageType> RelabelerType;
  RelabelerType::Pointer relabeler = RelabelerType::New();
  relabeler->SetInput( connecter->GetOutput() );
  relabeler->Update();

  typedef itk::Vector<double, ImageDimension> VectorType;
  typedef itk::VectorContainer<unsigned int, VectorType> VectorContainerType;
  VectorContainerType::Pointer vectors = VectorContainerType::New();
  vectors->Initialize();
  for ( unsigned int i = 0; i < relabeler->GetNumberOfObjects(); i++ )
    {
    vectors->CreateElementAt( i );
    VectorType v;
    v.Fill( 0 );
    vectors->SetElement( i, v );
    }  

  itk::ImageRegionIteratorWithIndex<ImageType> It( relabeler->GetOutput(), 
    relabeler->GetOutput()->GetLargestPossibleRegion() );
  
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( It.Get() > 0 )
      {
      VectorType v = vectors->GetElement( static_cast<unsigned int>( It.Get()-1 ) );
      ImageType::PointType point;
      relabeler->GetOutput()->TransformIndexToPhysicalPoint( It.GetIndex(), point );
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        v[i] += point[i];
        }
      vectors->SetElement( static_cast<unsigned int>( It.Get() )-1, v );
      }     
    }

  for ( unsigned int i = 0; i < relabeler->GetNumberOfObjects(); i++ )
    {
    VectorType v = vectors->GetElement( i );
    v /= relabeler->GetSizeOfObjectsInPixels()[i];
    vectors->SetElement( i, v );
    }

  itk::ContinuousIndex<double, ImageDimension> idx;
  idx[0] = 0.5*relabeler->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
  idx[1] = 0;
  idx[2] = 0;
  ImageType::PointType mid;
  relabeler->GetOutput()->TransformContinuousIndexToPhysicalPoint( idx, mid );
 
  vnl_vector<unsigned int> labels( relabeler->GetNumberOfObjects() );
  for ( unsigned int i = 0; i < relabeler->GetNumberOfObjects(); i++ )
    {
    VectorType v = vectors->GetElement( i );
    if ( vectors->GetElement( i )[0] < mid[0] )
      {
      labels[i] = 1;
      }
    else
      {
      labels[i] = 2;
      }
    }

  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( It.Get() > 0 )
      {
      It.Set( labels[static_cast<unsigned int>( It.Get() )-1] );
      }     
    }

  typedef itk::BinaryBallStructuringElement<unsigned int, ImageDimension> KernelType;
  KernelType ball;
  KernelType::SizeType size;
  size.Fill( 20 );
  ball.SetRadius( size );
  ball.CreateStructuringElement();

  typedef itk::BinaryMorphologicalClosingImageFilter<ImageType, ImageType, KernelType> CloserType;
  CloserType::Pointer closer1 = CloserType::New();
  closer1->SetInput( relabeler->GetOutput() );
  closer1->SetKernel( ball );
  closer1->SetSafeBorder( true );
  closer1->SetForegroundValue( 1 );
  closer1->Update();

  CloserType::Pointer closer2 = CloserType::New();
  closer2->SetInput( closer1->GetOutput() );
  closer2->SetKernel( ball );
  closer2->SetSafeBorder( true );
  closer2->SetForegroundValue( 2 );
  closer2->Update();

  WriterType::Pointer writercloser = WriterType::New();
  writercloser->SetFileName( argv[2] );
  writercloser->SetInput( closer2->GetOutput() );
  writercloser->Update();

}
