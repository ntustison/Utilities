
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkComposeDisplacementFieldsImageFilter.h"
#include "itkDisplacementFieldTransform.h"
#include "itkGaussianOperator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImportImageFilter.h"
#include "itkInvertDisplacementFieldImageFilter.h"
#include "itkIterationReporter.h"
#include "itkMultiplyImageFilter.h"
#include "itkVectorNeighborhoodOperatorImageFilter.h"
#include "itkWindowConvergenceMonitoringFunction.h"


template<class DisplacementFieldType>
typename DisplacementFieldType::Pointer
InvertDisplacementField(const DisplacementFieldType * field, const DisplacementFieldType * inverseFieldEstimate)
{
  using InverterType = itk::InvertDisplacementFieldImageFilter<DisplacementFieldType>;

  typename InverterType::Pointer inverter = InverterType::New();
  inverter->SetInput(field);
  inverter->SetInverseFieldInitialEstimate(inverseFieldEstimate);
  inverter->SetMaximumNumberOfIterations(20);
  inverter->SetMeanErrorToleranceThreshold(0.001);
  inverter->SetMaxErrorToleranceThreshold(0.1);
  inverter->Update();

  typename DisplacementFieldType::Pointer inverseField = inverter->GetOutput();

  return inverseField;
}

template<class DisplacementFieldType, class RealType = float>
float GetSimilarityErrorValue(const DisplacementFieldType * field1, const DisplacementFieldType * field2)
{
  itk::ImageRegionConstIterator<DisplacementFieldType> It1(field1, field1->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<DisplacementFieldType> It2(field2, field2->GetLargestPossibleRegion());

  RealType value = 0.0;
  RealType N = 0.0;
  for(It1.GoToBegin(), It2.GoToBegin(); !It1.IsAtEnd(); ++It1, ++It2)
    {
    value += (It1.Get() - It2.Get()).GetNorm();
    N += 1.0;
    }
  return(value / N);
}

template<class DisplacementFieldType, class RealType = float>
float GetInverseErrorValue(const DisplacementFieldType * field, const DisplacementFieldType * inverseField)
{
  itk::ImageRegionConstIterator<DisplacementFieldType> ItF(field, field->GetLargestPossibleRegion());

  using InterpolatorType = itk::VectorLinearInterpolateImageFunction<DisplacementFieldType, RealType>;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage(inverseField);

  RealType value = 0.0;
  RealType N = 0.0;
  for(ItF.GoToBegin(); !ItF.IsAtEnd(); ++ItF)
    {
    typename DisplacementFieldType::PointType imagePoint;
    field->TransformIndexToPhysicalPoint(ItF.GetIndex(), imagePoint);
    typename DisplacementFieldType::PixelType displacement = ItF.Get();

    typename InterpolatorType::PointType inputPoint;
    for(unsigned int d = 0; d < DisplacementFieldType::ImageDimension; d++)
      {
      inputPoint[d] = imagePoint[d] + displacement[d];
      }

    typename InterpolatorType::OutputType inverseDisplacement;
    if(interpolator->IsInsideBuffer(inputPoint))
      {
      inverseDisplacement = interpolator->Evaluate(inputPoint);
      }
    else
      {
      inverseDisplacement.Fill(0.0);
      }

    value += ((inputPoint + inverseDisplacement) - imagePoint).GetNorm();
    N += 1.0;
    }

  return(value / N);
}

template<class DisplacementFieldType, class RealType = float>
typename DisplacementFieldType::Pointer
GetUpdateField(const DisplacementFieldType * currentField, const DisplacementFieldType * inputField, RealType learningRate)
{
  const typename DisplacementFieldType::PixelType zeroVector(0.0);

  using RealImageType = itk::Image<RealType, DisplacementFieldType::ImageDimension>;

  typename DisplacementFieldType::Pointer updateField = DisplacementFieldType::New();
  updateField->CopyInformation(inputField);
  updateField->SetRegions(inputField->GetLargestPossibleRegion());
  updateField->Allocate();
  updateField->FillBuffer(zeroVector);

  using InterpolatorType = itk::VectorLinearInterpolateImageFunction<DisplacementFieldType, RealType>;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage(inputField);

  itk::ImageRegionConstIteratorWithIndex<DisplacementFieldType> It(currentField, currentField->GetLargestPossibleRegion());
  itk::ImageRegionIterator<DisplacementFieldType> ItU(updateField, updateField->GetLargestPossibleRegion());

  for(It.GoToBegin(), ItU.GoToBegin(); !It.IsAtEnd(); ++It, ++ItU)
    {
    typename DisplacementFieldType::PointType imagePoint;
    currentField->TransformIndexToPhysicalPoint(It.GetIndex(), imagePoint);
    typename DisplacementFieldType::PixelType displacement = It.Get();

    typename InterpolatorType::PointType inputPoint;
    for(unsigned int d = 0; d < DisplacementFieldType::ImageDimension; d++)
      {
      inputPoint[d] = imagePoint[d] + displacement[d];
      }

    typename InterpolatorType::OutputType inputDisplacement;
    if(interpolator->IsInsideBuffer(inputPoint))
      {
      inputDisplacement = interpolator->Evaluate(inputPoint);
      }
    else
      {
      inputDisplacement.Fill(0.0);
      }

    typename DisplacementFieldType::PixelType updateDisplacement;
    for(unsigned int d = 0; d < DisplacementFieldType::ImageDimension; d++)
      {
      updateDisplacement[d] = inputDisplacement[d];
      }
    ItU.Set(updateDisplacement);
    }

  // Scale the update field

  typename DisplacementFieldType::SpacingType spacing = updateField->GetSpacing();
  itk::ImageRegionConstIterator<DisplacementFieldType> ItF(updateField, updateField->GetLargestPossibleRegion());

  RealType maxNorm = itk::NumericTraits<RealType>::NonpositiveMin();
  for(ItF.GoToBegin(); !ItF.IsAtEnd(); ++ItF)
    {
    typename DisplacementFieldType::PixelType vector = ItF.Get();

    RealType localNorm = 0;
    for(itk::SizeValueType d = 0; d < DisplacementFieldType::ImageDimension; d++)
      {
      localNorm += itk::Math::sqr(vector[d] / spacing[d]);
      }
    localNorm = std::sqrt(localNorm);

    if(localNorm > maxNorm)
      {
      maxNorm = localNorm;
      }
    }

  RealType scale = learningRate;
  if (maxNorm > itk::NumericTraits<RealType>::ZeroValue())
    {
    scale /= maxNorm;
    }

  using RealImageType = itk::Image<RealType, DisplacementFieldType::ImageDimension>;

  using MultiplierType = itk::MultiplyImageFilter<DisplacementFieldType, RealImageType, DisplacementFieldType>;
  typename MultiplierType::Pointer multiplier = MultiplierType::New();
  multiplier->SetInput(updateField);
  multiplier->SetConstant(scale);

  typename DisplacementFieldType::Pointer scaledUpdateField = multiplier->GetOutput();
  scaledUpdateField->Update();
  scaledUpdateField->DisconnectPipeline();

  return scaledUpdateField;
}

template<class DisplacementFieldType, class RealType = float>
typename DisplacementFieldType::Pointer
GaussianSmoothDisplacementField(const DisplacementFieldType * field, const RealType variance)
{
  using DuplicatorType = itk::ImageDuplicator<DisplacementFieldType>;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(field);
  duplicator->Update();

 typename DisplacementFieldType::Pointer smoothField = duplicator->GetOutput();

  if(variance <= 0.0)
    {
    return smoothField;
    }

  using GaussianSmoothingOperatorType = itk::GaussianOperator<RealType, DisplacementFieldType::ImageDimension>;
  GaussianSmoothingOperatorType gaussianSmoothingOperator;

  using GaussianSmoothingSmootherType =
    itk::VectorNeighborhoodOperatorImageFilter<DisplacementFieldType, DisplacementFieldType>;
  typename GaussianSmoothingSmootherType::Pointer smoother = GaussianSmoothingSmootherType::New();

  for(itk::SizeValueType d = 0; d < DisplacementFieldType::ImageDimension; d++)
    {
    // smooth along this dimension
    gaussianSmoothingOperator.SetDirection(d);
    gaussianSmoothingOperator.SetVariance(variance);
    gaussianSmoothingOperator.SetMaximumError(0.001);
    gaussianSmoothingOperator.SetMaximumKernelWidth(smoothField->GetLargestPossibleRegion().GetSize()[d]);
    gaussianSmoothingOperator.CreateDirectional();

    // todo: make sure we only smooth within the buffered region
    smoother->SetOperator(gaussianSmoothingOperator);
    smoother->SetInput(smoothField);
    try
      {
      smoother->Update();
      }
    catch(itk::ExceptionObject & exc)
      {
      std::string msg("Caught exception: ");
      msg += exc.what();
      exit(1);
      }

    smoothField = smoother->GetOutput();
    smoothField->Update();
    smoothField->DisconnectPipeline();
    }

  const typename DisplacementFieldType::PixelType zeroVector(0.0);

  // make sure boundary does not move
  RealType weight1 = 1.0;
  if(variance < 0.5)
    {
    weight1 = 1.0 - 1.0 * (variance / 0.5);
    }
  RealType weight2 = 1.0 - weight1;

  const typename DisplacementFieldType::RegionType region = field->GetLargestPossibleRegion();
  const typename DisplacementFieldType::SizeType   size = region.GetSize();
  const typename DisplacementFieldType::IndexType  startIndex = region.GetIndex();

  itk::ImageRegionConstIteratorWithIndex<DisplacementFieldType> ItF(field, field->GetLargestPossibleRegion());
  itk::ImageRegionIteratorWithIndex<DisplacementFieldType>      ItS(smoothField, smoothField->GetLargestPossibleRegion());
  for(ItF.GoToBegin(), ItS.GoToBegin(); !ItF.IsAtEnd(); ++ItF, ++ItS)
    {
    typename DisplacementFieldType::IndexType index = ItF.GetIndex();
    bool                                      isOnBoundary = false;
    for(unsigned int d = 0; d < DisplacementFieldType::ImageDimension; d++)
      {
      if(index[d] == startIndex[d] || index[d] ==
        static_cast<typename DisplacementFieldType::IndexType::IndexValueType>(size[d]) - startIndex[d] - 1)
        {
        isOnBoundary = true;
        break;
        }
      }
    if(isOnBoundary)
      {
      ItS.Set(zeroVector);
      }
    else
      {
      ItS.Set(ItS.Get() * weight1 + ItF.Get() * weight2);
      }
    }

  return smoothField;
}

template<unsigned int ImageDimension>
int EstimateDiffeomorphism(int argc, char *argv[])
{
  typedef double RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  using DisplacementFieldTransformType = itk::DisplacementFieldTransform<RealType, ImageDimension>;
  using DisplacementFieldType = typename DisplacementFieldTransformType::DisplacementFieldType;
  using DisplacementFieldPointer = typename DisplacementFieldType::Pointer;
  using DisplacementVectorType = typename DisplacementFieldType::PixelType;

  using ReaderType = itk::ImageFileReader<DisplacementFieldType>;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[2]);

  DisplacementFieldPointer inputField = reader->GetOutput();
  inputField->Update();
  inputField->DisconnectPipeline();

  const DisplacementVectorType zeroVector(0.0);

  typename DisplacementFieldType::Pointer forwardDisplacementField = DisplacementFieldType::New();
  forwardDisplacementField->CopyInformation(inputField);
  forwardDisplacementField->SetRegions(inputField->GetLargestPossibleRegion());
  forwardDisplacementField->Allocate();
  forwardDisplacementField->FillBuffer(zeroVector);

  typename DisplacementFieldType::Pointer inverseDisplacementField = DisplacementFieldType::New();
  inverseDisplacementField->CopyInformation(inputField);
  inverseDisplacementField->SetRegions(inputField->GetLargestPossibleRegion());
  inverseDisplacementField->Allocate();
  inverseDisplacementField->FillBuffer(zeroVector);

  typename DisplacementFieldTransformType::Pointer transform = DisplacementFieldTransformType::New();
  transform->SetDisplacementField(forwardDisplacementField);
  transform->SetInverseDisplacementField(inverseDisplacementField) ;

  // Monitor the convergence
  using ConvergenceMonitoringType = itk::Function::WindowConvergenceMonitoringFunction<RealType>;
  typename ConvergenceMonitoringType::Pointer convergenceMonitoring = ConvergenceMonitoringType::New();
  convergenceMonitoring->SetWindowSize(10);

  bool isConverged = false;
  unsigned int currentIteration = 0;

  /** **/
  unsigned int totalNumberOfIterations = 100;
  RealType convergenceThreshold = 0.0001;
  RealType totalSmoothingGaussianVariance = 0.0;
  RealType learningRate = 0.1;
  /** **/

  std::cout << std::setw(10) << "Iter";
  std::cout << std::setw(20) << "Similarity error";
  std::cout << std::setw(20) << "(Inverse error)";
  std::cout << std::setw(20) << "Convergence value" << std::endl;

  while(currentIteration++ < totalNumberOfIterations && isConverged == false)
    {
    DisplacementFieldPointer updateField = GetUpdateField<DisplacementFieldType>(forwardDisplacementField, inputField, learningRate);

    using ComposerType = itk::ComposeDisplacementFieldsImageFilter<DisplacementFieldType>;

    typename ComposerType::Pointer forwardComposer = ComposerType::New();
    forwardComposer->SetDisplacementField(updateField);
    forwardComposer->SetWarpingField(transform->GetDisplacementField());
    forwardComposer->Update();

    DisplacementFieldPointer forwardSmoothTotalField = GaussianSmoothDisplacementField<DisplacementFieldType>(
      forwardComposer->GetOutput(), totalSmoothingGaussianVariance);

    DisplacementFieldPointer inverseField = InvertDisplacementField<DisplacementFieldType>(
      forwardSmoothTotalField, transform->GetInverseDisplacementField());
    DisplacementFieldPointer forwardField = InvertDisplacementField<DisplacementFieldType>(
      inverseField, transform->GetDisplacementField());

    transform->SetDisplacementField(forwardField);
    transform->SetInverseDisplacementField(inverseField);

    float currentSimilarityErrorValue = GetSimilarityErrorValue<DisplacementFieldType>(forwardField, inputField);
    float currentInverseErrorValue = GetInverseErrorValue<DisplacementFieldType>(forwardField, inverseField);

    convergenceMonitoring->AddEnergyValue(currentSimilarityErrorValue);
    RealType currentConvergenceValue = convergenceMonitoring->GetConvergenceValue();

    std::cout << std::setw(10) << currentIteration;
    std::cout << std::setw(20) << currentSimilarityErrorValue;
    std::cout << std::setw(20) << currentInverseErrorValue;
    std::cout << std::setw(20) << currentConvergenceValue << std::endl;

    if(currentConvergenceValue < convergenceThreshold)
      {
      isConverged = true;
      }
    }

  // Write the output
  std::string outputFilePrefix(argv[3]);

  using WriterType = itk::ImageFileWriter<DisplacementFieldType>;

  {
  std::string filename = outputFilePrefix + std::string("Warp.nii.gz");
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename.c_str());
  writer->SetInput(transform->GetDisplacementField());
  writer->Update();
  }

  {
  std::string filename = outputFilePrefix + std::string("InverseWarp.nii.gz");
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename.c_str());
  writer->SetInput(transform->GetInverseDisplacementField());
  writer->Update();
  }

  return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
  if (argc < 7)
    {
    std::cout << argv[0] << " imageDimension inputDisplacementField outputPrefix" << std::endl;
    exit(0);
    }

  switch(atoi(argv[1]))
   {
   case 2:
     EstimateDiffeomorphism<2>(argc, argv);
     break;
   case 3:
     EstimateDiffeomorphism<3>(argc, argv);
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit(EXIT_FAILURE);
   }
}

