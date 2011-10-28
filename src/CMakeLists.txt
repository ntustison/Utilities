project( Utilities )

cmake_minimum_required(VERSION 2.4)


if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#
# Find ITK
#
include (${CMAKE_ROOT}/Modules/FindITK.cmake)
if (USE_ITK_FILE)
  include(${USE_ITK_FILE})
endif(USE_ITK_FILE)

set( CMAKE_INCLUDE_DIRECTORIES_BEFORE ON )

include_directories(
  ${CMAKE_SOURCE_DIR}/../ExperimentalITK/Algorithms/
  ${CMAKE_SOURCE_DIR}/../ExperimentalITK/BasicFilters/
  ${CMAKE_SOURCE_DIR}/../ExperimentalITK/Common/
  ${CMAKE_SOURCE_DIR}/../ExperimentalITK/IO/
  ${CMAKE_SOURCE_DIR}/../ExperimentalITK/Numerics/
  ${CMAKE_SOURCE_DIR}/../ExperimentalITK/Numerics/Statistics/
  )

set(CMAKE_MODULE_PATH ${ITK_SOURCE_DIR}/CMake)

# Manage FFT
#
# USE_FFTWD -- use double precision fftw
option(USE_FFTWD "Use double precision fftw if found" OFF)
mark_as_advanced(USE_FFTWD)
#
# USE_FFTWF -- use single precision fftw
option(USE_FFTWF "Use single precision fftw if found" OFF)
mark_as_advanced(USE_FFTWF)

# USE_SYSTEM_FFTW -- locate a readybuilt fftw installation
option(USE_SYSTEM_FFTW "Use an installed version of fftw" OFF)
mark_as_advanced(USE_SYSTEM_FFTW)

include (${CMAKE_ROOT}/Modules/FindVTK.cmake)
if (USE_VTK_FILE)
  include(${USE_VTK_FILE})
endif (USE_VTK_FILE)

link_directories(
  ${ITK_LIBRARY_PATH}
)

set(SOURCES "${CMAKE_SOURCE_DIR}/../ExperimentalITK/Common/itkCommandLineParser" "${CMAKE_SOURCE_DIR}/../ExperimentalITK/Common/itkCommandLineOption" )

set(UI_SOURCES "${CMAKE_SOURCE_DIR}/../ExperimentalITK/Common/antsCommandLineParser" "${CMAKE_SOURCE_DIR}/../ExperimentalITK/Common/antsCommandLineOption" )

add_executable( itkTimeAndMemoryProbeTest itkTimeAndMemoryProbeTest.cxx )
target_link_libraries( itkTimeAndMemoryProbeTest ${ITK_LIBRARIES})

#add_executable( BSplineExample BSplineExample.cxx )
#target_link_libraries( BSplineExample ${ITK_LIBRARIES})

add_executable( BSplineTest BSplineTest.cxx )
target_link_libraries( BSplineTest ${ITK_LIBRARIES})

add_executable(AdaptiveHistogramEqualizeImage AdaptiveHistogramEqualizeImage.cxx )
target_link_libraries(AdaptiveHistogramEqualizeImage ${ITK_LIBRARIES})

add_executable(AddGaussianNoiseToImage AddGaussianNoiseToImage.cxx )
target_link_libraries(AddGaussianNoiseToImage ${ITK_LIBRARIES})

add_executable(AddImpulseNoiseToImage AddImpulseNoiseToImage.cxx )
target_link_libraries(AddImpulseNoiseToImage ${ITK_LIBRARIES})

add_executable(AlignImageCenters AlignImageCenters.cxx )
target_link_libraries(AlignImageCenters ${ITK_LIBRARIES})

add_executable(ApproximateImageWithBSplines ApproximateImageWithBSplines.cxx )
target_link_libraries(ApproximateImageWithBSplines ${ITK_LIBRARIES})

add_executable(AttachDicomHeaderToImage AttachDicomHeaderToImage.cxx )
target_link_libraries(AttachDicomHeaderToImage ${ITK_LIBRARIES})

add_executable(BinaryMorphology BinaryMorphology.cxx )
target_link_libraries(BinaryMorphology ${ITK_LIBRARIES})

add_executable(BinaryOperateImages BinaryOperateImages.cxx )
target_link_libraries(BinaryOperateImages ${ITK_LIBRARIES})

add_executable( BSplineBasisFunctions BSplineBasisFunctions.cxx )
target_link_libraries( BSplineBasisFunctions ${ITK_LIBRARIES})

add_executable(CalculateFirstOrderStatisticsFromImage CalculateFirstOrderStatisticsFromImage.cxx )
target_link_libraries(CalculateFirstOrderStatisticsFromImage ${ITK_LIBRARIES})

add_executable(LabelIntensityStatistics LabelIntensityStatistics.cxx )
target_link_libraries(LabelIntensityStatistics ${ITK_LIBRARIES})

#add_executable(CalculateStatisticsFromImagePair CalculateStatisticsFromImagePair.cxx )
#target_link_libraries(CalculateStatisticsFromImagePair ${ITK_LIBRARIES})

add_executable(CalculateHausdorffDistance CalculateHausdorffDistance.cxx )
target_link_libraries(CalculateHausdorffDistance ${ITK_LIBRARIES})

add_executable(CalculateVolumeFromBinaryImage CalculateVolumeFromBinaryImage.cxx )
target_link_libraries(CalculateVolumeFromBinaryImage ${ITK_LIBRARIES})

add_executable(ChangeImageInformation ChangeImageInformation.cxx )
target_link_libraries(ChangeImageInformation ${ITK_LIBRARIES})

add_executable(CompareLandmarks CompareLandmarks.cxx )
target_link_libraries(CompareLandmarks ${ITK_LIBRARIES})

add_executable(ConcatenateAvantsLandmarks ConcatenateAvantsLandmarks.cxx )
target_link_libraries(ConcatenateAvantsLandmarks)

if (USE_VTK_FILE)
add_executable(ConvertAvantsLandmarkFileToVTK ConvertAvantsLandmarkFileToVTK.cxx )
target_link_libraries(ConvertAvantsLandmarkFileToVTK vtkIO)
endif (USE_VTK_FILE)

if (USE_VTK_FILE)
add_executable(ConvertAvantsLandmarkFileToVTKLines ConvertAvantsLandmarkFileToVTKLines.cxx )
target_link_libraries(ConvertAvantsLandmarkFileToVTKLines vtkIO)
endif (USE_VTK_FILE)

add_executable(ConvertDICOMImageSeries ConvertDICOMImageSeries.cxx )
target_link_libraries(ConvertDICOMImageSeries ${ITK_LIBRARIES})

#add_executable(DicomSeriesReadImageWrite DicomSeriesReadImageWrite.cxx )
#target_link_libraries(DicomSeriesReadImageWrite ${ITK_LIBRARIES})

if (USE_VTK_FILE)
add_executable(ConvertScalarImageToVTKSurface ConvertScalarImageToVTKSurface.cxx )
target_link_libraries(ConvertScalarImageToVTKSurface ${ITK_LIBRARIES} vtkIO)
endif (USE_VTK_FILE)

add_executable(ConvertLabeledPointSet ConvertLabeledPointSet.cxx )
target_link_libraries(ConvertLabeledPointSet ${ITK_LIBRARIES})

add_executable(ConvertUnlabeledPointSet ConvertUnlabeledPointSet.cxx )
target_link_libraries(ConvertUnlabeledPointSet ${ITK_LIBRARIES})

add_executable(ConvertDeformationFieldToAvantsLandmarkFiles ConvertDeformationFieldToAvantsLandmarkFiles.cxx )
target_link_libraries(ConvertDeformationFieldToAvantsLandmarkFiles ${ITK_LIBRARIES})

add_executable(ConvertScalarImageToRGB ConvertScalarImageToRGB.cxx )
target_link_libraries(ConvertScalarImageToRGB ${ITK_LIBRARIES})

ADD_EXECUTABLE(ConvertImage ConvertImage.cxx )
TARGET_link_libraries(ConvertImage ${ITK_LIBRARIES})

add_executable(ConvertImageToDicom ConvertImageToDicom.cxx )
target_link_libraries(ConvertImageToDicom ${ITK_LIBRARIES})

add_executable(ConvertImageToDicomRGB ConvertImageToDicomRGB.cxx )
target_link_libraries(ConvertImageToDicomRGB ${ITK_LIBRARIES})

add_executable(GetDicomTagValue GetDicomTagValue.cxx )
target_link_libraries(GetDicomTagValue ${ITK_LIBRARIES})

add_executable(ConvertImageSeries ConvertImageSeries.cxx )
target_link_libraries(ConvertImageSeries ${ITK_LIBRARIES})

add_executable(ConvertImageSeriesToDicom ConvertImageSeriesToDicom.cxx )
target_link_libraries(ConvertImageSeriesToDicom ${ITK_LIBRARIES})

add_executable(ConvolveImage ConvolveImage.cxx )
target_link_libraries(ConvolveImage ${ITK_LIBRARIES})

# add_executable(CreateDeformedGrid CreateDeformedGrid.cxx )
# target_link_libraries(CreateDeformedGrid ${ITK_LIBRARIES})

add_executable(CreateADCImage CreateADCImage.cxx )
target_link_libraries(CreateADCImage ${ITK_LIBRARIES})

add_executable(CreateBoundingBoxMask CreateBoundingBoxMask.cxx )
target_link_libraries(CreateBoundingBoxMask ${ITK_LIBRARIES})

add_executable(CreateDirectionalStrainImages CreateDirectionalStrainImages.cxx )
target_link_libraries(CreateDirectionalStrainImages ${ITK_LIBRARIES})

add_executable(CreateColoredLabeledImage CreateColoredLabeledImage.cxx )
target_link_libraries(CreateColoredLabeledImage ${ITK_LIBRARIES})

add_executable(CreateImageSource CreateImageSource.cxx )
target_link_libraries(CreateImageSource ${ITK_LIBRARIES} )

#add_executable(CreateJacobianDeterminantImage CreateJacobianDeterminantImage.cxx )
#target_link_libraries(CreateJacobianDeterminantImage ${ITK_LIBRARIES} )

add_executable(CreateRidgeMap CreateRidgeMap.cxx )
target_link_libraries(CreateRidgeMap ${ITK_LIBRARIES})

add_executable(CreateImage CreateImage.cxx )
target_link_libraries(CreateImage ${ITK_LIBRARIES})

add_executable(CreatePCAImageDecomposition CreatePCAImageDecomposition.cxx )
target_link_libraries(CreatePCAImageDecomposition ${ITK_LIBRARIES})

add_executable(CreatePrincipalStrainImages CreatePrincipalStrainImages.cxx )
target_link_libraries(CreatePrincipalStrainImages ${ITK_LIBRARIES})

add_executable( CreateBSplineDeformationFieldFromPoints CreateBSplineDeformationFieldFromPoints.cxx )
target_link_libraries( CreateBSplineDeformationFieldFromPoints ${ITK_LIBRARIES} )

add_executable(CreateTopologicalNumberMapFromBinaryImage CreateTopologicalNumberMapFromBinaryImage.cxx )
target_link_libraries(CreateTopologicalNumberMapFromBinaryImage ${ITK_LIBRARIES})

add_executable(DivideLungs DivideLungs.cxx )
target_link_libraries(DivideLungs ${ITK_LIBRARIES})

add_executable( DrawBresenhamLinesFromMaskToPoint DrawBresenhamLinesFromMaskToPoint.cxx )
target_link_libraries(DrawBresenhamLinesFromMaskToPoint ${ITK_LIBRARIES})

add_executable( DynamicProgrammingForSkullLayers DynamicProgrammingForSkullLayers.cxx )
target_link_libraries(DynamicProgrammingForSkullLayers ${ITK_LIBRARIES})

add_executable(ExtractLungs ExtractLungs.cxx )
target_link_libraries(ExtractLungs ${ITK_LIBRARIES})

add_executable(ExtractContours ExtractContours.cxx )
target_link_libraries(ExtractContours ${ITK_LIBRARIES})

add_executable(ExtractMaskImage ExtractMaskImage.cxx )
target_link_libraries(ExtractMaskImage ${ITK_LIBRARIES})

add_executable(ExtractRGBComponent ExtractRGBComponent.cxx )
target_link_libraries(ExtractRGBComponent ${ITK_LIBRARIES})

add_executable( ExtractRegionFromImage ExtractRegionFromImage.cxx )
target_link_libraries( ExtractRegionFromImage ${ITK_LIBRARIES} )

add_executable(ExtractSliceFromImage ExtractSliceFromImage.cxx )
target_link_libraries(ExtractSliceFromImage ${ITK_LIBRARIES})

if(USE_FFTWF)
#  add_executable(ExtractTagLinePoints ExtractTagLinePoints.cxx )
#  target_link_libraries(ExtractTagLinePoints ${ITK_LIBRARIES})
#  add_executable(FourierTransformImage FourierTransformImage.cxx )
#  target_link_libraries(FourierTransformImage ${ITK_LIBRARIES})
endif(USE_FFTWF)

add_executable( FindClosestPointsBetweenTwoObjects FindClosestPointsBetweenTwoObjects.cxx )
target_link_libraries(FindClosestPointsBetweenTwoObjects ${ITK_LIBRARIES})

add_executable( FitBSplineCurveToAvantsLandmarks FitBSplineCurveToAvantsLandmarks.cxx )
target_link_libraries(FitBSplineCurveToAvantsLandmarks ${ITK_LIBRARIES})

add_executable(FlipImage FlipImage.cxx )
target_link_libraries(FlipImage ${ITK_LIBRARIES})

add_executable(FastMarching FastMarching.cxx )
target_link_libraries(FastMarching ${ITK_LIBRARIES})

add_executable(FloodFill FloodFill.cxx )
target_link_libraries(FloodFill ${ITK_LIBRARIES})

add_executable(GetCaminoToPhysicalSpaceAffine GetCaminoToPhysicalSpaceAffine.cxx )
target_link_libraries(GetCaminoToPhysicalSpaceAffine ${ITK_LIBRARIES})

add_executable(GetPhysicalSpaceToParaviewAffine GetPhysicalSpaceToParaviewAffine.cxx )
target_link_libraries(GetPhysicalSpaceToParaviewAffine ${ITK_LIBRARIES})

add_executable(GetConnectedComponents GetConnectedComponents.cxx )
target_link_libraries(GetConnectedComponents ${ITK_LIBRARIES})

add_executable(GenerateAttenuationMask GenerateAttenuationMask.cxx )
target_link_libraries(GenerateAttenuationMask ${ITK_LIBRARIES})

add_executable(GenerateImageFromIndexFile GenerateImageFromIndexFile.cxx )
target_link_libraries(GenerateImageFromIndexFile ${ITK_LIBRARIES})

add_executable(GeneralToBSplineDeformationField GeneralToBSplineDeformationField.cxx )
target_link_libraries(GeneralToBSplineDeformationField ${ITK_LIBRARIES})

add_executable(GenerateCooccurrenceMeasures GenerateCooccurrenceMeasures.cxx )
target_link_libraries(GenerateCooccurrenceMeasures ${ITK_LIBRARIES})

add_executable(GenerateDeformationFieldMagnitudeImage GenerateDeformationFieldMagnitudeImage.cxx )
target_link_libraries(GenerateDeformationFieldMagnitudeImage ${ITK_LIBRARIES})

add_executable(GenerateDistanceImage GenerateDistanceImage.cxx )
target_link_libraries(GenerateDistanceImage ${ITK_LIBRARIES})

add_executable(GenerateFractalImage GenerateFractalImage.cxx )
target_link_libraries(GenerateFractalImage ${ITK_LIBRARIES})

add_executable(GenerateImageFromControlPointLattice GenerateImageFromControlPointLattice.cxx )
target_link_libraries(GenerateImageFromControlPointLattice ${ITK_LIBRARIES})

add_executable(GenerateOEFImages GenerateOEFImages.cxx )
target_link_libraries(GenerateOEFImages ${ITK_LIBRARIES})

#add_executable( GenerateTPSDeformationField GenerateTPSDeformationField.cxx )
#target_link_libraries( GenerateTPSDeformationField ${ITK_LIBRARIES} )

add_executable(GenerateVectorFieldFromControlPointLattice GenerateVectorFieldFromControlPointLattice.cxx )
target_link_libraries(GenerateVectorFieldFromControlPointLattice ${ITK_LIBRARIES})

add_executable(GeneratePercentileAttenuationMask GeneratePercentileAttenuationMask.cxx )
target_link_libraries(GeneratePercentileAttenuationMask ${ITK_LIBRARIES})

add_executable(GenerateRunLengthMeasures GenerateRunLengthMeasures.cxx )
target_link_libraries(GenerateRunLengthMeasures ${ITK_LIBRARIES})

add_executable(GenerateStatisticsFromImage GenerateStatisticsFromImage.cxx )
target_link_libraries(GenerateStatisticsFromImage ${ITK_LIBRARIES})

add_executable(GetImageInformation GetImageInformation.cxx )
target_link_libraries(GetImageInformation ${ITK_LIBRARIES})

add_executable(GetJointSamples GetJointSamples.cxx )
target_link_libraries(GetJointSamples ${ITK_LIBRARIES})

add_executable(GetPhysicalToImageAffineTransform GetPhysicalToImageAffineTransform.cxx )
target_link_libraries(GetPhysicalToImageAffineTransform ${ITK_LIBRARIES})

#add_executable(GradientAnisotropicDiffusionImageFilter GradientAnisotropicDiffusionImageFilter.cxx )
#target_link_libraries(GradientAnisotropicDiffusionImageFilter ${ITK_LIBRARIES})

add_executable(GradientImage GradientImage.cxx )
target_link_libraries(GradientImage ${ITK_LIBRARIES})

add_executable(HessianImage HessianImage.cxx )
target_link_libraries(HessianImage ${ITK_LIBRARIES})

add_executable(GradientMagnitudeImageFilter GradientMagnitudeImageFilter.cxx )
target_link_libraries(GradientMagnitudeImageFilter ${ITK_LIBRARIES})

add_executable(GrayscaleMorphology GrayscaleMorphology.cxx )
target_link_libraries(GrayscaleMorphology ${ITK_LIBRARIES})

add_executable(HessianBasedFeatures HessianBasedFeatures.cxx )
target_link_libraries(HessianBasedFeatures ${ITK_LIBRARIES} )

add_executable(HistogramMatchImages HistogramMatchImages.cxx )
target_link_libraries(HistogramMatchImages ${ITK_LIBRARIES} )

add_executable(HoughTransform HoughTransform.cxx )
target_link_libraries(HoughTransform ${ITK_LIBRARIES})

add_executable(InvertImageIntensity InvertImageIntensity.cxx )
target_link_libraries(InvertImageIntensity ${ITK_LIBRARIES})

#add_executable(JHCT JHCT.cxx)
#target_link_libraries(JHCT ${ITK_LIBRARIES} )

# add_executable(EvaluatePointSetMetric EvaluatePointSetMetric.cxx)
# target_link_libraries(EvaluatePointSetMetric ${ITK_LIBRARIES} )

add_executable(LabelGeometryMeasures LabelGeometryMeasures.cxx )
target_link_libraries(LabelGeometryMeasures ${ITK_LIBRARIES})

#add_executable(MakePDFImageFromPointSet MakePDFImageFromPointSet.cxx )
#target_link_libraries(MakePDFImageFromPointSet ${ITK_LIBRARIES})

add_executable(MakeCheckerBoardImage MakeCheckerBoardImage.cxx )
target_link_libraries(MakeCheckerBoardImage ${ITK_LIBRARIES})

add_executable( MinimalPath MinimalPath.cxx )
target_link_libraries( MinimalPath ${ITK_LIBRARIES} )

add_executable( STAPLEAnalysis STAPLEAnalysis.cxx )
target_link_libraries( STAPLEAnalysis ${ITK_LIBRARIES} )

add_executable(MultipleOperateImages MultipleOperateImages.cxx )
target_link_libraries(MultipleOperateImages ${ITK_LIBRARIES})

add_executable(OtsuThresholdImage OtsuThresholdImage.cxx )
target_link_libraries(OtsuThresholdImage ${ITK_LIBRARIES})

#add_executable(MSQ MSQ.cxx )
#target_link_libraries(MSQ ${ITK_LIBRARIES})

add_executable(OrderDicomSeries OrderDicomSeries.cxx )
target_link_libraries(OrderDicomSeries ${ITK_LIBRARIES})

add_executable(PadImage PadImage.cxx )
target_link_libraries(PadImage ${ITK_LIBRARIES})

add_executable(PasteImageIntoImage PasteImageIntoImage.cxx )
target_link_libraries(PasteImageIntoImage ${ITK_LIBRARIES})

add_executable(PointSetPairStatistics PointSetPairStatistics.cxx )
target_link_libraries(PointSetPairStatistics ${ITK_LIBRARIES})

add_executable(PermuteAxesImage PermuteAxesImage.cxx)
target_link_libraries(PermuteAxesImage ${ITK_LIBRARIES})

add_executable(PropagateLabels PropagateLabels.cxx)
target_link_libraries(PropagateLabels ${ITK_LIBRARIES})

add_executable(RegionalMinMax RegionalMinMax.cxx )
target_link_libraries(RegionalMinMax ${ITK_LIBRARIES})

#add_executable(RegisterPointSetsDMFFD RegisterPointSetsDMFFD.cxx ${SOURCES} )
#target_link_libraries(RegisterPointSetsDMFFD ${ITK_LIBRARIES})

#add_executable(ICP ICP.cxx ${SOURCES} )
#target_link_libraries(ICP ${ITK_LIBRARIES})

#add_executable(RegisterPerfusionImages RegisterPerfusionImages.cxx ${SOURCES})
#target_link_libraries(RegisterPerfusionImages ${ITK_LIBRARIES})

add_executable(ResampleImage ResampleImage.cxx )
target_link_libraries(ResampleImage ${ITK_LIBRARIES})

add_executable(RescaleImageIntensity RescaleImageIntensity.cxx )
target_link_libraries(RescaleImageIntensity ${ITK_LIBRARIES})

add_executable(RigidTransformImage RigidTransformImage.cxx )
target_link_libraries(RigidTransformImage ${ITK_LIBRARIES})

add_executable(SegmentAirways SegmentAirways.cxx )
target_link_libraries(SegmentAirways ${ITK_LIBRARIES})

add_executable(SegmentHeliumLungs SegmentHeliumLungs.cxx )
target_link_libraries(SegmentHeliumLungs ${ITK_LIBRARIES})

add_executable(SeparateLungs SeparateLungs.cxx )
target_link_libraries(SeparateLungs ${ITK_LIBRARIES})

add_executable(ShapeMorphology ShapeMorphology.cxx )
target_link_libraries(ShapeMorphology ${ITK_LIBRARIES})

add_executable(ShapePriorLevelSet ShapePriorLevelSet.cxx )
target_link_libraries(ShapePriorLevelSet ${ITK_LIBRARIES})

add_executable(SimpleRegistration SimpleRegistration.cxx )
target_link_libraries(SimpleRegistration ${ITK_LIBRARIES})

add_executable(SmoothLungs SmoothLungs.cxx )
target_link_libraries(SmoothLungs ${ITK_LIBRARIES})

add_executable(SplitHeliumLungs SplitHeliumLungs.cxx )
target_link_libraries(SplitHeliumLungs ${ITK_LIBRARIES})

add_executable(SSI SSI.cxx )
target_link_libraries(SSI ${ITK_LIBRARIES})

add_executable( StatsForSkullThickness StatsForSkullThickness.cxx )
target_link_libraries( StatsForSkullThickness ${ITK_LIBRARIES})

add_executable(BSplineCalc BSplineCalc.cxx )
target_link_libraries(BSplineCalc ${ITK_LIBRARIES})

add_executable(ThresholdImage ThresholdImage.cxx )
target_link_libraries(ThresholdImage ${ITK_LIBRARIES})

add_executable(ThresholdPlane ThresholdPlane.cxx )
target_link_libraries(ThresholdPlane ${ITK_LIBRARIES})

add_executable(TileImages TileImages.cxx )
target_link_libraries(TileImages ${ITK_LIBRARIES})

add_executable(TruncateImageIntensity TruncateImageIntensity.cxx )
target_link_libraries(TruncateImageIntensity ${ITK_LIBRARIES})

add_executable(UnaryOperateImage UnaryOperateImage.cxx )
target_link_libraries(UnaryOperateImage ${ITK_LIBRARIES})

add_executable(WarpPoints WarpPoints.cxx )
target_link_libraries(WarpPoints ${ITK_LIBRARIES})

add_executable(WarpVTKPolyData WarpVTKPolyData.cxx )
target_link_libraries(WarpVTKPolyData ${ITK_LIBRARIES})

add_executable(WellComposeImage WellComposeImage.cxx )
target_link_libraries(WellComposeImage ${ITK_LIBRARIES})

add_executable(ZeroCrossingImage ZeroCrossingImage.cxx )
target_link_libraries(ZeroCrossingImage ${ITK_LIBRARIES})
