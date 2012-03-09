#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkInvertDisplacementFieldImageFilter.h"
//#include "itkBSplineInvertDisplacementFieldImageFilter.h"
#include "itkVector.h"

#include <string>
#include <vector>

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


template <unsigned int ImageDimension>
int Invert( int argc, char *argv[] )
{
  typedef itk::Vector<double, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;

  typedef itk::ImageFileReader<DisplacementFieldType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  unsigned int numberOfIterations = 20;
  double meanTolerance = 0.001;
  double maxTolerance = 0.1;
  if( argc > 4 )
    {
    numberOfIterations = atoi( argv[4] );
    }
  if( argc > 5 )
    {
    meanTolerance = atof( argv[5] );
    }
  if( argc > 6 )
    {
    maxTolerance = atof( argv[6] );
    }

  typedef itk::InvertDisplacementFieldImageFilter<DisplacementFieldType> InverterType;
  typename InverterType::Pointer inverter = InverterType::New();
  inverter->SetInput( reader->GetOutput() );
  inverter->SetMaximumNumberOfIterations( numberOfIterations );
  inverter->SetMeanErrorToleranceThreshold( meanTolerance );
  inverter->SetMaxErrorToleranceThreshold( maxTolerance );

  if( argc > 7 )
    {
    typename ReaderType::Pointer reader2 = ReaderType::New();
    reader2->SetFileName( argv[7] );
    reader2->Update();

    inverter->SetInverseFieldInitialEstimate( reader2->GetOutput() );
    }

  inverter->Update();

  typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( inverter->GetOutput() );
  writer->Update();

  return 0;
}

//template <unsigned int ImageDimension>
//int BSplineInvertField( int argc, char *argv[] )
//{
//  typedef itk::Vector<double, ImageDimension> VectorType;
//  typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;
//
//  typedef itk::ImageFileReader<DisplacementFieldType> ReaderType;
//  typename ReaderType::Pointer reader = ReaderType::New();
//  reader->SetFileName( argv[2] );
//  reader->Update();
//
//  typedef itk::BSplineInvertDisplacementFieldImageFilter<DisplacementFieldType> InverterType;
//
//  unsigned int splineOrder = 3;
//  unsigned int numberOfFittingLevels = 1;
//  typename InverterType::ArrayType numberOfControlPoints;
//  numberOfControlPoints.Fill( 4 );
//
//  if( argc > 4 )
//    {
//    numberOfControlPoints.Fill( atoi( argv[4] ) );
//    }
//  if( argc > 5 )
//    {
//    numberOfFittingLevels = atoi( argv[5] );
//    }
//  if( argc > 6 )
//    {
//    splineOrder = atoi( argv[6] );
//    }
//
//  typename InverterType::Pointer inverter = InverterType::New();
//  inverter->SetInput( reader->GetOutput() );
//  inverter->SetNumberOfFittingLevels( numberOfFittingLevels );
//  inverter->SetSplineOrder( splineOrder );
//  inverter->SetNumberOfControlPoints( numberOfControlPoints );
//  inverter->Update();
//
//  typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
//  typename WriterType::Pointer writer = WriterType::New();
//  writer->SetFileName( argv[3] );
//  writer->SetInput( inverter->GetOutput() );
//  writer->Update();
//
//  return 0;
//}

template <class DisplacementFieldType>
void
ComposeDiffs( typename DisplacementFieldType::Pointer fieldtowarpby,
  typename DisplacementFieldType::Pointer field,
  typename DisplacementFieldType::Pointer fieldout, float timesign)
{
  typedef double TReal;
  const unsigned int ImageDimension = DisplacementFieldType::ImageDimension;

  typedef itk::Image<TReal, ImageDimension> ImageType;


  typedef typename DisplacementFieldType::PixelType VectorType;

  typedef itk::Point<TReal,ImageDimension> VPointType;

  if (!fieldout)
    {
    fieldout=DisplacementFieldType::New();
    fieldout->SetSpacing( fieldtowarpby->GetSpacing() );
    fieldout->SetOrigin( fieldtowarpby->GetOrigin() );
    fieldout->SetDirection( fieldtowarpby->GetDirection() );
    fieldout->SetLargestPossibleRegion(fieldtowarpby->GetLargestPossibleRegion()  );
    fieldout->SetRequestedRegion( fieldtowarpby->GetLargestPossibleRegion()   );
    fieldout->SetBufferedRegion( fieldtowarpby->GetLargestPossibleRegion()  );
    fieldout->Allocate();
    VectorType zero;  zero.Fill(0);
    fieldout->FillBuffer(zero);
    }
    typedef typename DisplacementFieldType::PixelType VectorType;

    typedef DisplacementFieldType FieldType;
    typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType>         FieldIterator;
    typedef ImageType TRealImageType;

    typedef itk::ImageFileWriter<ImageType> writertype;

    typename ImageType::SpacingType oldspace = field->GetSpacing();
    typename ImageType::SpacingType newspace = fieldtowarpby->GetSpacing();


    typedef typename DisplacementFieldType::IndexType IndexType;
    typedef typename DisplacementFieldType::PointType PointType;


    typedef itk::VectorLinearInterpolateImageFunction<DisplacementFieldType,TReal> DefaultInterpolatorType;
    typename DefaultInterpolatorType::Pointer vinterp =  DefaultInterpolatorType::New();
    vinterp->SetInputImage(field);
    //    vinterp->SetParameters(NULL,1);


    VPointType pointIn1;
    VPointType pointIn2;
    typename DefaultInterpolatorType::ContinuousIndexType  contind; // married to pointIn2
    VPointType pointIn3;
    unsigned int ct=0;
    // iterate through fieldtowarpby finding the points that it maps to via field.
    // then take the difference from the original point and put it in the output field.
    //      std::cout << " begin iteration " << std::endl;
    FieldIterator m_FieldIter( fieldtowarpby, fieldtowarpby->GetLargestPossibleRegion());
    for(  m_FieldIter.GoToBegin(); !m_FieldIter.IsAtEnd(); ++m_FieldIter )
      {
      IndexType index = m_FieldIter.GetIndex();
      bool dosample = true;
      //      if (sub && m_TRealImage->GetPixel(index) < 0.5) dosample=false;
      if (dosample)
        {

    fieldtowarpby->TransformIndexToPhysicalPoint( index, pointIn1 );
    VectorType disp=m_FieldIter.Get();
    for (int jj=0; jj<ImageDimension; jj++)
      {
      pointIn2[jj]=disp[jj]+pointIn1[jj];
      }
    typename DefaultInterpolatorType::OutputType disp2;
    if (vinterp->IsInsideBuffer(pointIn2)) disp2 = vinterp->Evaluate( pointIn2 );
    else disp2.Fill(0);
    for (int jj=0; jj<ImageDimension; jj++) pointIn3[jj]=disp2[jj]*timesign+pointIn2[jj];

    VectorType out;
    for (int jj=0; jj<ImageDimension; jj++) out[jj]=pointIn3[jj]-pointIn1[jj];

    fieldout->SetPixel(m_FieldIter.GetIndex(),out);
    ct++;

        }//endif
      }//end iteration
}


template <unsigned int ImageDimension>
int Invert2( int argc, char *argv[] )
{

  typedef double TReal;

  unsigned int numberOfIterations = 20;
  TReal meanTolerance = 0.001;
  TReal maxTolerance = 0.1;
  if( argc > 4 )
    {
    numberOfIterations = atoi( argv[4] );
    }
  if( argc > 5 )
    {
    meanTolerance = atof( argv[5] );
    }
  if( argc > 6 )
    {
    maxTolerance = atoi( argv[6] );
    }

  TReal weight = 1.0;

  TReal mytoler=maxTolerance;
  unsigned int mymaxiter=numberOfIterations;

  typedef itk::Vector<TReal, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;
  typedef itk::Image<TReal, ImageDimension> ImageType;

  typedef itk::ImageFileReader<DisplacementFieldType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename DisplacementFieldType::Pointer field = reader->GetOutput();
  typedef typename DisplacementFieldType::Pointer DisplacementFieldPointer;

  VectorType zero; zero.Fill(0);

  typename ImageType::Pointer TRealImage = ImageType::New();
  TRealImage->SetLargestPossibleRegion( field->GetLargestPossibleRegion() );
  TRealImage->SetBufferedRegion( field->GetLargestPossibleRegion().GetSize() );
  TRealImage->SetSpacing(field->GetSpacing());
  TRealImage->SetOrigin(field->GetOrigin());
  TRealImage->SetDirection(field->GetDirection());
  TRealImage->Allocate();

  typedef typename DisplacementFieldType::IndexType IndexType;
  typedef typename VectorType::ValueType           ScalarType;
  typedef itk::ImageRegionIteratorWithIndex<DisplacementFieldType> Iterator;

  DisplacementFieldPointer lagrangianInitCond=DisplacementFieldType::New();
  lagrangianInitCond->SetSpacing( field->GetSpacing() );
  lagrangianInitCond->SetOrigin( field->GetOrigin() );
  lagrangianInitCond->SetDirection( field->GetDirection() );
  lagrangianInitCond->SetLargestPossibleRegion( field->GetLargestPossibleRegion() );
  lagrangianInitCond->SetRequestedRegion(field->GetRequestedRegion() );
  lagrangianInitCond->SetBufferedRegion( field->GetLargestPossibleRegion() );
  lagrangianInitCond->Allocate();
  DisplacementFieldPointer eulerianInitCond=DisplacementFieldType::New();
  eulerianInitCond->SetSpacing( field->GetSpacing() );
  eulerianInitCond->SetOrigin( field->GetOrigin() );
  eulerianInitCond->SetDirection( field->GetDirection() );
  eulerianInitCond->SetLargestPossibleRegion( field->GetLargestPossibleRegion() );
  eulerianInitCond->SetRequestedRegion(field->GetRequestedRegion() );
  eulerianInitCond->SetBufferedRegion( field->GetLargestPossibleRegion() );
  eulerianInitCond->Allocate();

  DisplacementFieldPointer inverseField=DisplacementFieldType::New();
  inverseField->SetSpacing( field->GetSpacing() );
  inverseField->SetOrigin( field->GetOrigin() );
  inverseField->SetDirection( field->GetDirection() );
  inverseField->SetLargestPossibleRegion( field->GetLargestPossibleRegion() );
  inverseField->SetRequestedRegion(field->GetRequestedRegion() );
  inverseField->SetBufferedRegion( field->GetLargestPossibleRegion() );
  inverseField->Allocate();
  inverseField->FillBuffer( zero );


  typedef typename DisplacementFieldType::SizeType SizeType;
  SizeType size=field->GetLargestPossibleRegion().GetSize();


  typename ImageType::SpacingType spacing = field->GetSpacing();
  TReal subpix=0.0;
  unsigned long npix=1;
  for (int j=0; j<ImageDimension; j++)  // only use in-plane spacing
  {
    npix*=field->GetLargestPossibleRegion().GetSize()[j];
  }
  subpix=pow((TReal)ImageDimension,(TReal)ImageDimension)*0.5;

  TReal max=0;
    Iterator iter( field, field->GetLargestPossibleRegion() );
    for(  iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
      IndexType  index=iter.GetIndex();
      VectorType vec1=iter.Get();
      VectorType newvec=vec1*weight;
      lagrangianInitCond->SetPixel(index,newvec);
      TReal mag=0;
      for (unsigned int jj=0; jj<ImageDimension; jj++) mag+=newvec[jj]*newvec[jj];
      mag=sqrt(mag);
      if (mag > max) max=mag;
    }

    eulerianInitCond->FillBuffer(zero);

    TReal scale=(1.)/max;
    if (scale > 1.) scale=1.0;
//    TReal initscale=scale;
    Iterator vfIter( inverseField, inverseField->GetLargestPossibleRegion() );

//  int num=10;
//  for (int its=0; its<num; its++)
    TReal difmag=10.0;
  unsigned int ct=0;
    TReal meandif=1.e8;
//    int badct=0;
//  while (difmag > subpix && meandif > subpix*0.1 && badct < 2 )//&& ct < 20 && denergy > 0)
//    TReal length=0.0;
    TReal stepl=2.;
    TReal lastdifmag=0;

    TReal epsilon = (TReal)size[0]/256;
    if (epsilon > 1) epsilon = 1;

    while ( ct < mymaxiter && difmag > mytoler && meandif > 0.001)
  {
    meandif=0.0;

    //this field says what position the eulerian field should contain in the E domain
    ComposeDiffs<DisplacementFieldType>(inverseField,lagrangianInitCond,    eulerianInitCond, 1);

    difmag=0.0;
    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
      {
    IndexType  index=vfIter.GetIndex();
    VectorType  update=eulerianInitCond->GetPixel(index);
    TReal mag=0;
    for (int j=0; j<ImageDimension;j++)
      {
        update[j]*=(-1.0);
        mag+=(update[j]/spacing[j])*(update[j]/spacing[j]);
                    }
    mag=sqrt(mag);
    meandif+=mag;
    if (mag > difmag) {difmag=mag; }
    //      if (mag < 1.e-2) update.Fill(0);

    eulerianInitCond->SetPixel(index,update);
    TRealImage->SetPixel(index,mag);
    }
    meandif/=(TReal)npix;
    if (ct == 0) epsilon = 0.75;
    else epsilon=0.5;
    stepl=difmag*epsilon;

    for(  vfIter.GoToBegin(); !vfIter.IsAtEnd(); ++vfIter )
    {
      TReal val = TRealImage->GetPixel(vfIter.GetIndex());
      VectorType update=eulerianInitCond->GetPixel(vfIter.GetIndex());
      if (val > stepl) update = update * (stepl/val);
      VectorType upd=vfIter.Get()+update * (epsilon);
      vfIter.Set(upd);
    }
    ct++;
    std::cout << "Iteration " << ct << ": mean error norm = " << meandif
      << ", max error norm = " << difmag <<
      "( " << mymaxiter << ", " << mytoler << " )" <<
      std::endl;
  }

  typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( inverseField );
  writer->Update();



}


int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension inputField outputField [numberOfIterations=20] [meanTolerance=0.001] [maxTolerance=0.1]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     Invert<2>( argc, argv );
     break;
   case 3:
     Invert<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

