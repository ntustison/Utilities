#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryFillholeImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryBoxStructuringElement.h"
#include "itkBinaryDiamondStructuringElement.h"
#include "itkBinaryThinning3DImageFilter.h"
#include "itkBinaryThinningImageFilter.h"

#include "itkCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageDuplicator.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkSliceBySliceImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

#include "vnl/vnl_math.h"

template <unsigned int ImageDimension>
int BinaryMorphology( int argc, char * argv[] )
{
  typedef short PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  unsigned int radius = 1;
  if ( argc > 5 )
    {
    radius = atoi( argv[5] );
    }

  PixelType foreground = itk::NumericTraits<PixelType>::One;
  PixelType background = itk::NumericTraits<PixelType>::Zero;
  if ( argc > 7 )
    {
    foreground = static_cast<PixelType>( atof( argv[7] ) );
    }
  if ( argc > 8 )
    {
    background = static_cast<PixelType>( atof( argv[8] ) );
    }

  unsigned int operation = static_cast<unsigned int>( atoi( argv[4] ) );

  if( operation == 5 )
    {
    typedef itk::BinaryFillholeImageFilter<ImageType>  FilterType;
    typename FilterType::Pointer  filter = FilterType::New();
    filter->SetInput( reader->GetOutput() );
    filter->SetForegroundValue( foreground );
    filter->SetFullyConnected( true );
    filter->Update();

    typedef itk::ImageFileWriter<ImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( filter->GetOutput() );
    writer->SetFileName( argv[3] );
    writer->Update();

    return EXIT_SUCCESS;
    }
  else if( operation == 6 )
    {
    typedef itk::VotingBinaryIterativeHoleFillingImageFilter<ImageType>  FilterType;
    typename FilterType::InputSizeType radii;
    radii.Fill( radius );

    typename FilterType::Pointer  filter = FilterType::New();
    filter->SetInput( reader->GetOutput() );
    filter->SetForegroundValue( foreground );
    filter->SetMajorityThreshold( 1 );  // 1 == default
    filter->SetRadius( radii );
    filter->SetMaximumNumberOfIterations( 100000000 );
    filter->Update();

    typedef itk::ImageFileWriter<ImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( filter->GetOutput() );
    writer->SetFileName( argv[3] );
    writer->Update();

    return EXIT_SUCCESS;
    }
  else if( operation == 7 )
    {
    typedef itk::Image<unsigned short, ImageDimension> ShortImageType;
    typedef itk::CastImageFilter<ImageType, ShortImageType> CasterType;
    typename CasterType::Pointer caster = CasterType::New();
    caster->SetInput( reader->GetOutput() );
    caster->Update();

    typedef itk::LabelStatisticsImageFilter<ShortImageType, ShortImageType>
      StatsFilterType;
    typename StatsFilterType::Pointer stats = StatsFilterType::New();
    stats->SetLabelInput( caster->GetOutput() );
    stats->SetInput( caster->GetOutput() );
    stats->Update();

    typename ShortImageType::RegionType region = stats->GetRegion( foreground );
    typename ShortImageType::IndexType lowerIndex = region.GetIndex();
    typename ShortImageType::IndexType upperIndex = region.GetUpperIndex();

    unsigned int direction = 0;
    if ( argc > 8 )
      {
      direction = static_cast<unsigned int>( atof( argv[8] ) );
      }

    itk::ImageRegionIteratorWithIndex<ImageType> It( reader->GetOutput(),
      reader->GetOutput()->GetLargestPossibleRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      typename ImageType::IndexType index = It.GetIndex();
      if( index[direction] == lowerIndex[direction] || index[direction] == upperIndex[direction] )
        {
        It.Set( 0 );
        }
      }

    typedef itk::ImageFileWriter<ImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( reader->GetOutput() );
    writer->SetFileName( argv[3] );
    writer->Update();

    return EXIT_SUCCESS;
    }


  if ( argc < 6 || atoi( argv[6] ) == 1 )
    {
    typedef itk::BinaryBallStructuringElement<
                        PixelType,
                        ImageDimension> StructuringElementType;
    StructuringElementType  element;
    element.SetRadius( radius );
    element.CreateStructuringElement();

    switch ( operation )
      {
      case 0:
        {
        typedef itk::BinaryDilateImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
        filter->Update();
        typedef itk::ImageFileWriter<ImageType>  WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( filter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 1:
        {
        typedef itk::BinaryErodeImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
        filter->Update();
        typedef itk::ImageFileWriter<ImageType>  WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( filter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 2:
        {
        typedef itk::BinaryMorphologicalClosingImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetForegroundValue( foreground );
        filter->Update();
        typedef itk::ImageFileWriter<ImageType>  WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( filter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 3:
        {
        typedef itk::BinaryMorphologicalOpeningImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
        filter->Update();
        typedef itk::ImageFileWriter<ImageType>  WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( filter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      default:
        {
        std::cerr << "Invalid operation choice." << std::endl;
        return EXIT_FAILURE;
        }
      }
    }
  else if ( atoi( argv[6] ) == 0 )
    {
    typedef itk::BinaryBoxStructuringElement<
                        PixelType,
                        ImageDimension>  StructuringElementType;
    StructuringElementType element;
    element.SetRadius( radius );
    element.CreateStructuringElement();

    switch ( operation )
      {
      case 0:
        {
        typedef itk::BinaryDilateImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
        filter->Update();
        typedef itk::ImageFileWriter<ImageType>  WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( filter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 1:
        {
        typedef itk::BinaryErodeImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
        filter->Update();
        typedef itk::ImageFileWriter<ImageType>  WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( filter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 2:
        {
        typedef itk::BinaryMorphologicalClosingImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetForegroundValue( foreground );
        filter->Update();
        typedef itk::ImageFileWriter<ImageType>  WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( filter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 3:
        {
        typedef itk::BinaryMorphologicalOpeningImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
        filter->Update();
        typedef itk::ImageFileWriter<ImageType>  WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( filter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      default:
        {
        std::cerr << "Invalid operation choice." << std::endl;
        return EXIT_FAILURE;
        }
      }
    }
  else
    {
    typedef itk::BinaryDiamondStructuringElement<
                        PixelType,
                        ImageDimension> StructuringElementType;
    StructuringElementType  element;
    element.SetRadius( radius );
    element.CreateStructuringElement();

    switch ( operation )
      {
      case 0:
        {
        typedef itk::BinaryDilateImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
        filter->Update();
        typedef itk::ImageFileWriter<ImageType>  WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( filter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 1:
        {
        typedef itk::BinaryErodeImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
        filter->Update();
        typedef itk::ImageFileWriter<ImageType>  WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( filter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 2:
        {
        typedef itk::BinaryMorphologicalClosingImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetForegroundValue( foreground );
        filter->Update();
        typedef itk::ImageFileWriter<ImageType>  WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( filter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 3:
        {
        typedef itk::BinaryMorphologicalOpeningImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
        filter->Update();
        typedef itk::ImageFileWriter<ImageType>  WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( filter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      default:
        {
        std::cerr << "Invalid operation choice." << std::endl;
        return EXIT_FAILURE;
        }
      }
    }
  return EXIT_SUCCESS;
}


int myDilate( int argc, char * argv[] )
{

//  ImageType::Pointer im, ImageType::Pointer dmask, int Tlabel, int R

    const unsigned int ImageDimension = 3;



    typedef short PixelType;
    typedef itk::Image<PixelType, ImageDimension> ImageType;

    typedef itk::ImageRegionIteratorWithIndex< ImageType > IndexIteratorType;

    typedef itk::ImageFileReader<ImageType>  ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[2] );
    reader->Update();

    int R = 1;
    if ( argc > 5 )
      {
      R = atoi( argv[5] );
      }

    PixelType Tlabel = itk::NumericTraits<PixelType>::One;
    if ( argc > 7 )
      {
      Tlabel = static_cast<PixelType>( atof( argv[7] ) );
      }

    ImageType::Pointer im = reader->GetOutput();
    im->DisconnectPipeline();

    typedef itk::ImageDuplicator<ImageType> DuplicatorType;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage( im );
    duplicator->Update();

    ImageType::Pointer dmask = duplicator->GetOutput();
    dmask->DisconnectPipeline();

    int x, y, z, j, k, l;
    ImageType::IndexType idx;
    ImageType::Pointer tmask = ImageType::New();
    tmask->CopyInformation( im );
    tmask->SetRegions( im->GetRequestedRegion() );
    tmask->Allocate();
    tmask->FillBuffer( 0 );

    int Dx=im->GetRequestedRegion().GetSize()[0];
    int Dy=im->GetRequestedRegion().GetSize()[1];
    int Dz=im->GetRequestedRegion().GetSize()[2];

    IndexIteratorType dmaskit(dmask, dmask->GetRequestedRegion());
    IndexIteratorType tmaskit(tmask, tmask->GetRequestedRegion());
    IndexIteratorType imit(im, im->GetRequestedRegion());
    for (dmaskit.GoToBegin(),tmaskit.GoToBegin(); !dmaskit.IsAtEnd(); ++dmaskit,++tmaskit){
      dmaskit.Set(0);
      tmaskit.Set(0);
    }
    int tc=0;
    for (imit.GoToBegin(); !imit.IsAtEnd(); ++imit){
      idx = imit.GetIndex();
      if ((Tlabel>=0 && imit.Value()==Tlabel) || (Tlabel<0 && imit.Value()>0)){
        tc++;
        x=idx[0];
        for (j=-R;j<R+1;j++){
          idx[0] = x+j;
          if (idx[0]<0 || idx[0]>=Dx)
            continue;

          dmaskit.SetIndex(idx);
          dmaskit.Set(1);
        }
      }
    }
//    cout<<"tc: "<<tc<<endl;
    for (dmaskit.GoToBegin(); !dmaskit.IsAtEnd(); ++dmaskit){
      idx = dmaskit.GetIndex();
      if (dmaskit.Value()==1){
        y=idx[1];
        for (k=-R;k<R+1;k++){
          idx[1] = y+k;
          if (idx[1]<0 || idx[1]>=Dy)
            continue;

          tmaskit.SetIndex(idx);
          tmaskit.Set(1);
        }
      }
    }

    for (tmaskit.GoToBegin(); !tmaskit.IsAtEnd(); ++tmaskit){
      idx = tmaskit.GetIndex();
      if (tmaskit.Value()==1){
        z=idx[2];
        for (l=-R;l<R+1;l++){
          idx[2] = z+l;
          if (idx[2]<0 || idx[2]>=Dz)
            continue;

          dmaskit.SetIndex(idx);
          dmaskit.Set(1);
        }
      }
    }

    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( dmask );
    writer->SetFileName( argv[3] );
    writer->Update();

    return EXIT_SUCCESS;
}



int BinaryMorphologySliceBySlice( int argc, char * argv[] )
{
  const unsigned int ImageDimension = 3;

  typedef short PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::SliceBySliceImageFilter<ImageType, ImageType> SliceFilterType;
  SliceFilterType::Pointer sliceFilter = SliceFilterType::New();
  sliceFilter->SetInput( reader->GetOutput() );


  unsigned int radius = 1;
  if ( argc > 5 )
    {
    radius = atoi( argv[5] );
    }

  PixelType foreground = itk::NumericTraits<PixelType>::One;
  PixelType background = itk::NumericTraits<PixelType>::Zero;
  if ( argc > 7 )
    {
    foreground = static_cast<PixelType>( atof( argv[7] ) );
    }
  if ( argc > 8 )
    {
    background = static_cast<PixelType>( atof( argv[8] ) );
    }

  unsigned int operation = static_cast<unsigned int>( atoi( argv[4] ) );

  if( operation == 5 )
    {
    typedef itk::BinaryFillholeImageFilter<SliceFilterType::InternalInputImageType>  FilterType;
    FilterType::Pointer  filter = FilterType::New();
    filter->SetForegroundValue( foreground );
    filter->SetFullyConnected( true );

    sliceFilter->SetFilter( filter );

    typedef itk::ImageFileWriter<ImageType>  WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( sliceFilter->GetOutput() );
    writer->SetFileName( argv[3] );
    writer->Update();

    return EXIT_SUCCESS;
    }
  else if( operation == 6 )
    {
    typedef itk::VotingBinaryIterativeHoleFillingImageFilter<SliceFilterType::InternalInputImageType>  FilterType;
    FilterType::InputSizeType radii;
    radii.Fill( radius );

    FilterType::Pointer  filter = FilterType::New();
    filter->SetForegroundValue( foreground );
    filter->SetMajorityThreshold( 1 );  // 1 == default
    filter->SetRadius( radii );
    filter->SetMaximumNumberOfIterations( 100000000 );
    filter->Update();

    sliceFilter->SetFilter( filter );

    typedef itk::ImageFileWriter<ImageType>  WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( sliceFilter->GetOutput() );
    writer->SetFileName( argv[3] );
    writer->Update();

    return EXIT_SUCCESS;
    }


  if ( argc < 6 || atoi( argv[6] ) == 1 )
    {
    typedef itk::BinaryBallStructuringElement<
                        PixelType,
                        ImageDimension-1> StructuringElementType;
    StructuringElementType  element;
    element.SetRadius( radius );
    element.CreateStructuringElement();

    switch ( operation )
      {
      case 0:
        {
        typedef itk::BinaryDilateImageFilter<SliceFilterType::InternalInputImageType, SliceFilterType::InternalOutputImageType,
          StructuringElementType >  FilterType;
        FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );

								sliceFilter->SetFilter( filter );

								typedef itk::ImageFileWriter<ImageType>  WriterType;
								WriterType::Pointer writer = WriterType::New();
								writer->SetInput( sliceFilter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 1:
        {
        typedef itk::BinaryErodeImageFilter<SliceFilterType::InternalInputImageType, SliceFilterType::InternalOutputImageType,
          StructuringElementType >  FilterType;
        FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );

								sliceFilter->SetFilter( filter );

								typedef itk::ImageFileWriter<ImageType>  WriterType;
								WriterType::Pointer writer = WriterType::New();
								writer->SetInput( sliceFilter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 2:
        {
        typedef itk::BinaryMorphologicalClosingImageFilter<SliceFilterType::InternalInputImageType, SliceFilterType::InternalOutputImageType,
          StructuringElementType >  FilterType;
        FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetForegroundValue( foreground );
								sliceFilter->SetFilter( filter );

								typedef itk::ImageFileWriter<ImageType>  WriterType;
								WriterType::Pointer writer = WriterType::New();
								writer->SetInput( sliceFilter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 3:
        {
        typedef itk::BinaryMorphologicalOpeningImageFilter<SliceFilterType::InternalInputImageType, SliceFilterType::InternalOutputImageType,
          StructuringElementType >  FilterType;
        FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );

								sliceFilter->SetFilter( filter );

								typedef itk::ImageFileWriter<ImageType>  WriterType;
								WriterType::Pointer writer = WriterType::New();
								writer->SetInput( sliceFilter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      default:
        {
        std::cerr << "Invalid operation choice." << std::endl;
        return EXIT_FAILURE;
        }
      }
    }
  else if ( atoi( argv[6] ) == 0 )
    {
    typedef itk::BinaryBoxStructuringElement<
                        PixelType,
                        ImageDimension-1>  StructuringElementType;
    StructuringElementType element;
    element.SetRadius( radius );
    element.CreateStructuringElement();

    switch ( operation )
      {
      case 0:
        {
        typedef itk::BinaryDilateImageFilter<SliceFilterType::InternalInputImageType, SliceFilterType::InternalOutputImageType,
          StructuringElementType >  FilterType;
        FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
        filter->Update();

								sliceFilter->SetFilter( filter );

								typedef itk::ImageFileWriter<ImageType>  WriterType;
								WriterType::Pointer writer = WriterType::New();
								writer->SetInput( sliceFilter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 1:
        {
        typedef itk::BinaryErodeImageFilter<SliceFilterType::InternalInputImageType, SliceFilterType::InternalOutputImageType,
          StructuringElementType >  FilterType;
        FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );

								sliceFilter->SetFilter( filter );

								typedef itk::ImageFileWriter<ImageType>  WriterType;
								WriterType::Pointer writer = WriterType::New();
								writer->SetInput( sliceFilter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 2:
        {
        typedef itk::BinaryMorphologicalClosingImageFilter<SliceFilterType::InternalInputImageType, SliceFilterType::InternalOutputImageType,
          StructuringElementType >  FilterType;
        FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetForegroundValue( foreground );
        filter->Update();

								sliceFilter->SetFilter( filter );

								typedef itk::ImageFileWriter<ImageType>  WriterType;
								WriterType::Pointer writer = WriterType::New();
								writer->SetInput( sliceFilter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 3:
        {
        typedef itk::BinaryMorphologicalOpeningImageFilter<SliceFilterType::InternalInputImageType, SliceFilterType::InternalOutputImageType,
          StructuringElementType >  FilterType;
        FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );

								sliceFilter->SetFilter( filter );

								typedef itk::ImageFileWriter<ImageType>  WriterType;
								WriterType::Pointer writer = WriterType::New();
								writer->SetInput( sliceFilter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      default:
        {
        std::cerr << "Invalid operation choice." << std::endl;
        return EXIT_FAILURE;
        }
      }
    }
  else
    {
    typedef itk::BinaryDiamondStructuringElement<
                        PixelType,
                        ImageDimension-1> StructuringElementType;
    StructuringElementType  element;
    element.SetRadius( radius );
    element.CreateStructuringElement();

    switch ( operation )
      {
      case 0:
        {
        typedef itk::BinaryDilateImageFilter<SliceFilterType::InternalInputImageType, SliceFilterType::InternalOutputImageType,
          StructuringElementType >  FilterType;
        FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
								sliceFilter->SetFilter( filter );

								typedef itk::ImageFileWriter<ImageType>  WriterType;
								WriterType::Pointer writer = WriterType::New();
								writer->SetInput( sliceFilter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 1:
        {
        typedef itk::BinaryErodeImageFilter<SliceFilterType::InternalInputImageType, SliceFilterType::InternalOutputImageType,
          StructuringElementType >  FilterType;
        FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );

								sliceFilter->SetFilter( filter );

								typedef itk::ImageFileWriter<ImageType>  WriterType;
								WriterType::Pointer writer = WriterType::New();
								writer->SetInput( sliceFilter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 2:
        {
        typedef itk::BinaryMorphologicalClosingImageFilter<SliceFilterType::InternalInputImageType, SliceFilterType::InternalOutputImageType,
          StructuringElementType >  FilterType;
        FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetForegroundValue( foreground );

								sliceFilter->SetFilter( filter );

								typedef itk::ImageFileWriter<ImageType>  WriterType;
								WriterType::Pointer writer = WriterType::New();
								writer->SetInput( sliceFilter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      case 3:
        {
        typedef itk::BinaryMorphologicalOpeningImageFilter<SliceFilterType::InternalInputImageType, SliceFilterType::InternalOutputImageType,
          StructuringElementType >  FilterType;
        FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );

								sliceFilter->SetFilter( filter );

								typedef itk::ImageFileWriter<ImageType>  WriterType;
								WriterType::Pointer writer = WriterType::New();
								writer->SetInput( sliceFilter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        }
      default:
        {
        std::cerr << "Invalid operation choice." << std::endl;
        return EXIT_FAILURE;
        }
      }
    }


  return EXIT_SUCCESS;
}

int BinaryThin2D( int argc, char * argv[] )
{
  typedef short PixelType;
  typedef itk::Image<PixelType, 2> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::BinaryThinningImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer  filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->Update();

  typedef itk::ImageFileWriter<ImageType>  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();

  return EXIT_SUCCESS;
}

int BinaryThin3D( int argc, char * argv[] )
{
  typedef short PixelType;
  typedef itk::Image<PixelType, 3> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::BinaryThinning3DImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer  filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->Update();

  typedef itk::ImageFileWriter<ImageType>  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " imageDimension inputImage outputImage operation "
      << "[radius] [type: box == 0, ball = 1, diamond = 2] [label]" << std::endl;
    std::cerr << "  operation: " << std::endl;
    std::cerr << "    0. dilate" << std::endl;
    std::cerr << "    1. erode " << std::endl;
    std::cerr << "    2. close " << std::endl;
    std::cerr << "    3. open " << std::endl;
    std::cerr << "    4. thin " << std::endl;
    std::cerr << "    5. fill holes" << std::endl;
    std::cerr << "    6. iterative hole filling (approx. convex hull)" << std::endl;
    std::cerr << "    7. remove first and last slice in 3-D (must also specify direction)" << std::endl;
    return EXIT_FAILURE;
    }

  if( *argv[1] == 'X' )
    {
    BinaryMorphologySliceBySlice( argc, argv );
    }
  else if( *argv[1] == 'H' )
    {
    myDilate( argc, argv );
    }
  else
    {
				if( atoi( argv[4] ) == 4 )
						{
						switch( atoi( argv[1] ) )
							{
							case 2:
									BinaryThin2D( argc, argv );
									break;
							case 3:
									BinaryThin3D( argc, argv );
									break;
							default:
										std::cerr << "Unsupported dimension" << std::endl;
										exit( EXIT_FAILURE );
							}
						}
				else
						{
						switch( atoi( argv[1] ) )
							{
							case 2:
									BinaryMorphology<2>( argc, argv );
									break;
							case 3:
									BinaryMorphology<3>( argc, argv );
									break;
							default:
										std::cerr << "Unsupported dimension" << std::endl;
										exit( EXIT_FAILURE );
							}
					}
			}
}

