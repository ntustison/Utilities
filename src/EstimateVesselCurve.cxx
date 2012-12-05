
#include "itkBinaryThresholdImageFilter.h"
#include "itkBSplineControlPointImageFunction.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkCollidingFrontsImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkOtsuMultipleThresholdsImageFilter.h"
#include "itkPointSet.h"
#include "itkVector.h"

#include "vnl/vnl_cross.h"

#include <string>
#include <vector>
#include <fstream>

template<class TValue>
TValue Convert( std::string optionString )
			{
			TValue value;
			std::istringstream iss( optionString );
			iss >> value;
			return value;
			}

template<class TValue>
std::vector<TValue> ConvertVector( std::string optionString )
			{
			std::vector<TValue> values;
			std::string::size_type crosspos = optionString.find( 'x', 0 );

			if ( crosspos == std::string::npos )
					{
					values.push_back( Convert<TValue>( optionString ) );
					}
			else
					{
					std::string element = optionString.substr( 0, crosspos ) ;
					TValue value;
					std::istringstream iss( element );
					iss >> value;
					values.push_back( value );
					while ( crosspos != std::string::npos )
							{
							std::string::size_type crossposfrom = crosspos;
							crosspos = optionString.find( 'x', crossposfrom + 1 );
							if ( crosspos == std::string::npos )
									{
									element = optionString.substr( crossposfrom + 1, optionString.length() );
									}
							else
									{
									element = optionString.substr( crossposfrom + 1, crosspos ) ;
									}
							std::istringstream iss2( element );
							iss2 >> value;
							values.push_back( value );
							}
					}
			return values;
			}


int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " inputImage inputSeedImage outputSampledSplineAvantsFile" << std::endl;
    std::cout << "Note:  The seed image needs to have two labels (1 and 2) designating the two sets of seed points." << std::endl;
    exit( 1 );
    }

  typedef float RealType;
  typedef float PixelType;

  const unsigned int ImageDimension = 3;
  const unsigned int numberOfSamples = 100;
  const unsigned int numberOfIterations = 2;
  RealType gradientStep = 0.1;

  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::ImageFileReader<LabelImageType> SeedsReaderType;
  SeedsReaderType::Pointer seedsReader = SeedsReaderType::New();
  seedsReader->SetFileName( argv[2] );
  seedsReader->Update();

  // Otsu threshold the image and select the regions
  // corresponding to the MRA regions

  typedef itk::OtsuMultipleThresholdsImageFilter<ImageType, LabelImageType> OtsuType;
  OtsuType::Pointer otsu = OtsuType::New();
  otsu->SetInput( reader->GetOutput() );
  otsu->SetNumberOfHistogramBins( 200 );
  otsu->SetNumberOfThresholds( 4 );

  typedef itk::BinaryThresholdImageFilter<LabelImageType, ImageType> ThresholderType;
  ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( otsu->GetOutput() );
  thresholder->SetInsideValue( 1 );
  thresholder->SetOutsideValue( 0 );
  thresholder->SetLowerThreshold( 3 );
  thresholder->SetUpperThreshold( 3 );
  thresholder->Update();

  // Use colliding fronts to connect the seeds.  Also, create an
  // initial curve based on the seed points and the intermediate
  // point.  We'll use this to parameterize a better guess.

  typedef itk::Image<VectorType, 1> CurveImageType;

  typedef itk::PointSet<VectorType, 1> PointSetType;
  PointSetType::Pointer pointSet = PointSetType::New();
  pointSet->Initialize();

  typedef itk::BSplineScatteredDataPointSetToImageFilter<PointSetType, CurveImageType> BSplinerType;

  typedef itk::CollidingFrontsImageFilter<ImageType, ImageType> FilterType;

  typedef FilterType::NodeContainer NodeContainerType;
  typedef FilterType::NodeType NodeType;
  NodeContainerType::Pointer seeds1 = NodeContainerType::New();
  NodeContainerType::Pointer seeds2 = NodeContainerType::New();
  NodeContainerType::Pointer seeds3 = NodeContainerType::New();

  unsigned long seeds1Counter = 0;
  unsigned long seeds2Counter = 0;
  unsigned long seeds3Counter = 0;

  unsigned long count = 0;
  VectorType offset;
  offset.Fill( 0.0 );

  itk::ImageRegionConstIteratorWithIndex<LabelImageType> It( seedsReader->GetOutput(),
    seedsReader->GetOutput()->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( It.Get() == 1 )
      {
      LabelImageType::IndexType position = It.GetIndex();

      LabelImageType::PointType point;
      seedsReader->GetOutput()->TransformIndexToPhysicalPoint( position, point );
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        offset[d] += point[d];
        }
      count++;

      NodeType node;
      const double value = 0.0;

      node.SetValue( value );
      node.SetIndex( position );

      seeds1->InsertElement( seeds1Counter++, node );
      }
    else if( It.Get() == 2 )
      {
      LabelImageType::IndexType position = It.GetIndex();

      LabelImageType::PointType point;
      seedsReader->GetOutput()->TransformIndexToPhysicalPoint( position, point );
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        offset[d] += point[d];
        }
      count++;

      NodeType node;
      const double value = 0.0;

      node.SetValue( value );
      node.SetIndex( position );

      seeds2->InsertElement( seeds2Counter++, node );
      }
    else if( It.Get() == 3 )
      {
      LabelImageType::IndexType position = It.GetIndex();

      LabelImageType::PointType point;
      seedsReader->GetOutput()->TransformIndexToPhysicalPoint( position, point );
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        offset[d] += point[d];
        }
      count++;

      NodeType node;
      const double value = 0.0;

      node.SetValue( value );
      node.SetIndex( position );

      seeds3->InsertElement( seeds3Counter++, node );
      }
    }

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( thresholder->GetOutput() );
  filter->SetSeedPoints1( seeds1 );
  filter->SetSeedPoints2( seeds2 );
  filter->ApplyConnectivityOn();
  filter->Update();

  PixelType minimumValue = itk::NumericTraits<PixelType>::max();
  PixelType maximumValue = itk::NumericTraits<PixelType>::NonpositiveMin();

  itk::ImageRegionIterator<ImageType> It2( filter->GetOutput(),
    filter->GetOutput()->GetRequestedRegion() );
  for( It2.GoToBegin(); !It2.IsAtEnd(); ++It2 )
    {
    if( It2.Get() > maximumValue )
      {
      maximumValue = It2.Get();
      }
    if( It2.Get() < minimumValue )
      {
      minimumValue = It2.Get();
      }
    }

  PixelType slope = 1.0 / ( minimumValue - maximumValue );
  for( It2.GoToBegin(), It.GoToBegin(); !It2.IsAtEnd(); ++It2, ++It )
    {
    It2.Set( slope * ( It2.Get() - maximumValue ) );
    }


  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    offset[d] /= count;
    }

  BSplinerType::WeightsContainerType::Pointer weights = BSplinerType::WeightsContainerType::New();
  weights->Initialize();

  count = 0;
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( It.Get() == 1 )
      {
      LabelImageType::IndexType position = It.GetIndex();

      LabelImageType::PointType point;
      seedsReader->GetOutput()->TransformIndexToPhysicalPoint( position, point );
      VectorType bpoint;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        bpoint[d] = point[d] - offset[d];
        }
      PointSetType::PointType param;
      param[0] = 0.0;
      pointSet->SetPoint( count, param );
      pointSet->SetPointData( count, bpoint );
      weights->InsertElement( count, 1 );
      count++;
      }
    else if( It.Get() == 2 )
      {
      LabelImageType::IndexType position = It.GetIndex();

      LabelImageType::PointType point;
      seedsReader->GetOutput()->TransformIndexToPhysicalPoint( position, point );
      VectorType bpoint;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        bpoint[d] = point[d] - offset[d];
        }
      PointSetType::PointType param;
      param[0] = 1.0;
      pointSet->SetPoint( count, param );
      pointSet->SetPointData( count, bpoint );
      weights->InsertElement( count, 1 );
      count++;
      }
    else if( It.Get() == 3 )
      {
      LabelImageType::IndexType position = It.GetIndex();

      LabelImageType::PointType point;
      seedsReader->GetOutput()->TransformIndexToPhysicalPoint( position, point );
      VectorType bpoint;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        bpoint[d] = point[d] - offset[d];
        }
      PointSetType::PointType param;
      param[0] = 0.5;
      pointSet->SetPoint( count, param );
      pointSet->SetPointData( count, bpoint );
      weights->InsertElement( count, 1 );
      count++;
      }
    }

  BSplinerType::Pointer bspliner = BSplinerType::New();

  unsigned int iteration = 0;
  while( iteration++ < numberOfIterations )
    {
    bspliner->SetInput( pointSet );
    bspliner->SetGenerateOutputImage( false );

    CurveImageType::PointType origin;
    origin.Fill( 0.0 );
    bspliner->SetOrigin( origin );

    CurveImageType::SpacingType spacing;
    spacing[0] = 0.001;
    bspliner->SetSpacing( spacing );

    CurveImageType::SizeType size;
    size[0] = static_cast<unsigned int>( 1.0 / spacing[0] + 1 );
    bspliner->SetSize( size );

    BSplinerType::ArrayType order;
    order[0] = 3;
    bspliner->SetSplineOrder( order );

    BSplinerType::ArrayType ncps;
    ncps[0] = order[0] + 1;
    bspliner->SetNumberOfControlPoints( ncps );

    BSplinerType::ArrayType nlevels;
    nlevels[0] = 3;
    bspliner->SetNumberOfLevels( nlevels );

    BSplinerType::ArrayType close;
    close[0] = false;
    bspliner->SetCloseDimension( close );

    bspliner->SetPointWeights( weights );
    bspliner->Update();

    typedef itk::BSplineControlPointImageFunction<CurveImageType> BSplineFunctionType;
    BSplineFunctionType::Pointer bsplineFunction = BSplineFunctionType::New();
    bsplineFunction->SetSplineOrder( bspliner->GetSplineOrder() );
    bsplineFunction->SetSpacing( bspliner->GetSpacing() );
    bsplineFunction->SetSize( bspliner->GetSize() );
    bsplineFunction->SetOrigin( bspliner->GetOrigin() );
    bsplineFunction->SetInputImage( bspliner->GetPhiLattice() );

    typedef BSplineFunctionType::GradientType GradientType;

    std::vector<VectorType> gradients;
    std::vector<VectorType> normals;
    std::vector<VectorType> binormals;
    std::vector<VectorType> points;

    for( unsigned int n = 0; n < numberOfSamples+1; n++ )
      {
      PointSetType::PointType param;
      param[0] = static_cast<float>( n ) / ( numberOfSamples + 0.0001 );

      VectorType spatialPoint = bsplineFunction->Evaluate( param );
      GradientType gradientMatrix = bsplineFunction->EvaluateGradient( param );

      points.push_back( spatialPoint );
      points[n] += offset;

      VectorType grad;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        grad[d] = gradientMatrix[0][d];
        }
      grad.Normalize();
      gradients.push_back( grad );
      }

    for( unsigned int n = 0; n < gradients.size(); n++ )
      {
      VectorType normal;

      PointSetType::PointType param;
      param[0] = static_cast<float>( n ) / ( numberOfSamples + 0.0001 );

      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        normal[d] = bsplineFunction->EvaluateHessian( param, d )[0][0];
        }
      normal.Normalize();
      normals.push_back( normal );
      }

    for( unsigned int n = 0; n < gradients.size(); n++ )
      {
      vnl_vector<RealType> tmp = vnl_cross_3d( normals[n].GetVnlVector(), gradients[n].GetVnlVector() );
      tmp.normalize();
      VectorType binormal;
      binormal.SetVnlVector( tmp );

      binormals.push_back( binormal );
      }

    normals.clear();
    for( unsigned int n = 0; n < binormals.size(); n++ )
      {
      vnl_vector<RealType> tmp = vnl_cross_3d( binormals[n].GetVnlVector(), gradients[n].GetVnlVector() );
      tmp.normalize();
      VectorType normal;
      normal.SetVnlVector( tmp );

      normals.push_back( normal );
      }

    // Now traverse the curve and sample the normal plane

    typedef itk::LinearInterpolateImageFunction<ImageType> InterpolatorType;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    interpolator->SetInputImage( filter->GetOutput() );

    weights->Initialize();

    RealType voxelSpacingFactor = 0.0;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      voxelSpacingFactor += vnl_math_sqr( filter->GetOutput()->GetSpacing()[d] );
      }
    voxelSpacingFactor = vcl_sqrt( voxelSpacingFactor );

    pointSet->Initialize();

    count = 0;
    for( unsigned int n = 0; n < points.size(); n++ )
      {
      PointSetType::PointType param;
      param[0] = static_cast<float>( n ) / ( points.size() + 0.0001 );

      for( RealType u = -gradientStep; u <= gradientStep; u += ( gradientStep * 0.2 ) )
        {
        for( RealType v = -gradientStep; v <= gradientStep; v += ( gradientStep * 0.2 ) )
          {
          if( vcl_sqrt( vnl_math_sqr( u ) + vnl_math_sqr( v ) ) > 3.0 )
            {
            continue;
            }

          VectorType tmp = points[n] + normals[n] * u * voxelSpacingFactor + binormals[n] * v * voxelSpacingFactor;

          InterpolatorType::PointType samplePoint;
          for( unsigned int d = 0; d < ImageDimension; d++ )
            {
            samplePoint[d] = tmp[d];
            }

          if( ! interpolator->IsInsideBuffer( samplePoint ) )
            {
            continue;
            }

          RealType weight = interpolator->Evaluate( samplePoint );

          if( weight <= 0 )
            {
            continue;
            }

          weights->InsertElement( count, weight );

          pointSet->SetPoint( count, param );
          pointSet->SetPointData( count, tmp );
          count++;
          }
        }
      }
    }

  std::ofstream str( argv[3] );

  str << "0 0 0 0" << std::endl;
  itk::ImageRegionIterator<CurveImageType> ItB( bspliner->GetOutput(),
    bspliner->GetOutput()->GetLargestPossibleRegion() );
  for( ItB.GoToBegin(); !ItB.IsAtEnd(); ++ItB )
    {
    VectorType curvePoint = ItB.Get();
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      str << curvePoint[d] << " ";
      }
    str << "1" << std::endl;
    }
  str << "0 0 0 0" << std::endl;

//   typedef itk::ImageFileWriter<ImageType> WriterType;
//   WriterType::Pointer writer = WriterType::New();
//   writer->SetFileName( argv[3] );
//   writer->SetInput( filter->GetOutput() );
//   writer->Update();

  return 0;
}


