PROJECT( Utilities )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

#  
# Find ITK 
# 
INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake) 
IF (USE_ITK_FILE) 
  INCLUDE(${USE_ITK_FILE}) 
ENDIF(USE_ITK_FILE) 

SET(CMAKE_MODULE_PATH ${ITK_SOURCE_DIR}/CMake)
OPTION(USE_FFTWD "Use double precision FFTW if found" OFF)
OPTION(USE_FFTWF "Use single precision FFTW if found" ON)
FIND_PACKAGE( FFTW )

IF(USE_FFTWF)
   LINK_LIBRARIES(${FFTWF_LIB})
ENDIF(USE_FFTWF)

INCLUDE (${CMAKE_ROOT}/Modules/FindVTK.cmake)
IF (USE_VTK_FILE)
  INCLUDE(${USE_VTK_FILE})
ENDIF (USE_VTK_FILE)
 
LINK_DIRECTORIES(
  ${ITK_LIBRARY_PATH}
)

ADD_EXECUTABLE(AdaptiveHistogramEqualizeImage AdaptiveHistogramEqualizeImage.cxx )
TARGET_LINK_LIBRARIES(AdaptiveHistogramEqualizeImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(AddImage AddImage.cxx )
TARGET_LINK_LIBRARIES(AddImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(AlignImageCenters AlignImageCenters.cxx )
TARGET_LINK_LIBRARIES(AlignImageCenters ITKIO)

ADD_EXECUTABLE(BinaryMorphologicalCloseImage BinaryMorphologicalCloseImage.cxx )
TARGET_LINK_LIBRARIES(BinaryMorphologicalCloseImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(BoundedSpaceDilateImage BoundedSpaceDilateImage.cxx )
TARGET_LINK_LIBRARIES(BoundedSpaceDilateImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(CalculateGradientImage CalculateGradientImage.cxx )
TARGET_LINK_LIBRARIES(CalculateGradientImage ITKNumerics ITKIO 
ITKCommon)

ADD_EXECUTABLE(CalculateHausdorffDistance CalculateHausdorffDistance.cxx )
TARGET_LINK_LIBRARIES(CalculateHausdorffDistance ITKNumerics ITKIO 
ITKCommon)

ADD_EXECUTABLE(CalculateMetricImages CalculateMetricImages.cxx )
TARGET_LINK_LIBRARIES(CalculateMetricImages ITKNumerics ITKIO 
ITKCommon)

ADD_EXECUTABLE(CalculateVolumeFromBinaryImage CalculateVolumeFromBinaryImage.cxx )
TARGET_LINK_LIBRARIES(CalculateVolumeFromBinaryImage ITKNumerics ITKIO 
ITKCommon)

ADD_EXECUTABLE(ContourExtractor ContourExtractor.cxx )
TARGET_LINK_LIBRARIES(ContourExtractor ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(ConvertDeformationFieldToAvantsLandmarkFiles ConvertDeformationFieldToAvantsLandmarkFiles.cxx )
TARGET_LINK_LIBRARIES(ConvertDeformationFieldToAvantsLandmarkFiles ITKNumerics ITKIO ITKCommon vtkIO)

ADD_EXECUTABLE(ConvertDICOMImageSeries ConvertDICOMImageSeries.cxx )
TARGET_LINK_LIBRARIES(ConvertDICOMImageSeries ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(ConvertImage ConvertImage.cxx )
TARGET_LINK_LIBRARIES(ConvertImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(ConvertSegmentedImageToAvantsLandmarkFile ConvertSegmentedImageToAvantsLandmarkFile.cxx )
TARGET_LINK_LIBRARIES(ConvertSegmentedImageToAvantsLandmarkFile ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(CreateDeformedGrid CreateDeformedGrid.cxx )
TARGET_LINK_LIBRARIES(CreateDeformedGrid ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(CreateDifferenceImage CreateDifferenceImage.cxx )
TARGET_LINK_LIBRARIES(CreateDifferenceImage ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(CreateDirectionalStrainImages CreateDirectionalStrainImages.cxx )
TARGET_LINK_LIBRARIES(CreateDirectionalStrainImages ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(CreateGaussianImage CreateGaussianImage.cxx )
TARGET_LINK_LIBRARIES(CreateGaussianImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(CreateGridImage CreateGridImage.cxx )
TARGET_LINK_LIBRARIES(CreateGridImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(CreateImageVolumeFromImageSlices CreateImageVolumeFromImageSlices.cxx )
TARGET_LINK_LIBRARIES(CreateImageVolumeFromImageSlices ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(CreateJacobianDeterminantImage CreateJacobianDeterminantImage.cxx )
TARGET_LINK_LIBRARIES(CreateJacobianDeterminantImage ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(CreatePrincipalStrainImages CreatePrincipalStrainImages.cxx )
TARGET_LINK_LIBRARIES(CreatePrincipalStrainImages ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(CreateRandomImage CreateRandomImage.cxx )
TARGET_LINK_LIBRARIES(CreateRandomImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(CreateStrainTensorComponentImages CreateStrainTensorComponentImages.cxx )
TARGET_LINK_LIBRARIES(CreateStrainTensorComponentImages ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(CreateWarpedGridImage CreateWarpedGridImage.cxx )
TARGET_LINK_LIBRARIES(CreateWarpedGridImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(CreateZeroImage CreateZeroImage.cxx )
TARGET_LINK_LIBRARIES(CreateZeroImage ITKNumerics ITKIO ITKCommon)


ADD_EXECUTABLE(CropImage CropImage.cxx )
TARGET_LINK_LIBRARIES(CropImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(CurvatureAnisotropicDiffusionImageFilter CurvatureAnisotropicDiffusionImageFilter.cxx )
TARGET_LINK_LIBRARIES(CurvatureAnisotropicDiffusionImageFilter ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(DiceMetric DiceMetric.cxx )
TARGET_LINK_LIBRARIES(DiceMetric ITKIO)

ADD_EXECUTABLE(DilateImage DilateImage.cxx )
TARGET_LINK_LIBRARIES(DilateImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(ErodeImage ErodeImage.cxx )
TARGET_LINK_LIBRARIES(ErodeImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(ExtractMaskImage ExtractMaskImage.cxx )
TARGET_LINK_LIBRARIES(ExtractMaskImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(ExtractTagLinePoints ExtractTagLinePoints.cxx )
TARGET_LINK_LIBRARIES(ExtractTagLinePoints ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(FillHoles FillHoles.cxx )
TARGET_LINK_LIBRARIES(FillHoles ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(FlipImage FlipImage.cxx )
TARGET_LINK_LIBRARIES(FlipImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(FourierTransformImage FourierTransformImage.cxx )
TARGET_LINK_LIBRARIES(FourierTransformImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(GaussianSmoothImage GaussianSmoothImage.cxx )
TARGET_LINK_LIBRARIES(GaussianSmoothImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(GetConnectedComponents GetConnectedComponents.cxx )
TARGET_LINK_LIBRARIES(GetConnectedComponents ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(GeneralToBSplineDeformationField GeneralToBSplineDeformationField.cxx )
TARGET_LINK_LIBRARIES(GeneralToBSplineDeformationField ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(GenerateDeformationFieldDerivatives GenerateDeformationFieldDerivatives.cxx )
TARGET_LINK_LIBRARIES(GenerateDeformationFieldDerivatives ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(GenerateDeformationFieldFromAvantsLandmarks GenerateDeformationFieldFromAvantsLandmarks.cxx )
TARGET_LINK_LIBRARIES(GenerateDeformationFieldFromAvantsLandmarks  ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(GenerateDistanceImage GenerateDistanceImage.cxx )
TARGET_LINK_LIBRARIES(GenerateDistanceImage ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(GenerateHistogram GenerateHistogram.cxx )
TARGET_LINK_LIBRARIES(GenerateHistogram ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(GenerateHistogramFromBoxMask GenerateHistogramFromBoxMask.cxx )
TARGET_LINK_LIBRARIES(GenerateHistogramFromBoxMask ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(GenerateInverseBSplineDeformationField GenerateInverseBSplineDeformationField.cxx )
TARGET_LINK_LIBRARIES(GenerateInverseBSplineDeformationField ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(GenerateImageFromControlPointLattice GenerateImageFromControlPointLattice.cxx )
TARGET_LINK_LIBRARIES(GenerateImageFromControlPointLattice ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(GenerateStatisticsFromImage GenerateStatisticsFromImage.cxx )
TARGET_LINK_LIBRARIES(GenerateStatisticsFromImage ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(GradientAnisotropicDiffusionImageFilter GradientAnisotropicDiffusionImageFilter.cxx )
TARGET_LINK_LIBRARIES(GradientAnisotropicDiffusionImageFilter ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(GrayscaleDilateImage GrayscaleDilateImage.cxx )
TARGET_LINK_LIBRARIES(GrayscaleDilateImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(GrayscaleErodeImage GrayscaleErodeImage.cxx )
TARGET_LINK_LIBRARIES(GrayscaleErodeImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(InhomogeneityCorrectImage InhomogeneityCorrectImage.cxx )
TARGET_LINK_LIBRARIES(InhomogeneityCorrectImage ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(JHCT JHCT.cxx )
TARGET_LINK_LIBRARIES(JHCT ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(JHCTMulti JHCTMulti.cxx )
TARGET_LINK_LIBRARIES(JHCTMulti ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(Hough2DTransform Hough2DTransform.cxx )
TARGET_LINK_LIBRARIES(Hough2DTransform ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(LabelAndCloseImage LabelAndCloseImage.cxx )
TARGET_LINK_LIBRARIES(LabelAndCloseImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(LevelSets LevelSets.cxx )
TARGET_LINK_LIBRARIES(LevelSets ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(LiveWire LiveWire.cxx )
TARGET_LINK_LIBRARIES(LiveWire ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(MakeCheckerBoardImage MakeCheckerBoardImage.cxx )
TARGET_LINK_LIBRARIES(MakeCheckerBoardImage ITKNumerics ITKIO ITKCommon)

#ADD_EXECUTABLE(MakePDFImageFromPointSet MakePDFImageFromPointSet.cxx )
#TARGET_LINK_LIBRARIES(MakePDFImageFromPointSet ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(MultiplyImages MultiplyImages.cxx )
TARGET_LINK_LIBRARIES(MultiplyImages ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(OtsuThresholdImage OtsuThresholdImage.cxx )
TARGET_LINK_LIBRARIES(OtsuThresholdImage ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(PermuteAxesImage PermuteAxesImage.cxx )
TARGET_LINK_LIBRARIES(PermuteAxesImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(PadImage PadImage.cxx )
TARGET_LINK_LIBRARIES(PadImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(RegisterImagesFFD RegisterImagesFFD.cxx )
TARGET_LINK_LIBRARIES(RegisterImagesFFD ITKNumerics ITKIO ITKCommon ITKStatistics)

ADD_EXECUTABLE(RegisterMultiplePointSetsFFD RegisterMultiplePointSetsFFD.cxx )
TARGET_LINK_LIBRARIES(RegisterMultiplePointSetsFFD ITKNumerics ITKIO ITKCommon ITKStatistics vtkIO)

ADD_EXECUTABLE(RegisterPointSetsFFD RegisterPointSetsFFD.cxx )
TARGET_LINK_LIBRARIES(RegisterPointSetsFFD ITKNumerics ITKIO ITKCommon ITKStatistics vtkIO)

ADD_EXECUTABLE(RegisterWithMultiResDemons RegisterWithMultiResDemons.cxx )
TARGET_LINK_LIBRARIES(RegisterWithMultiResDemons ITKNumerics ITKIO  ITKCommon)

ADD_EXECUTABLE(ResampleImage ResampleImage.cxx )
TARGET_LINK_LIBRARIES(ResampleImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(RescaleImageIntensity RescaleImageIntensity.cxx )
TARGET_LINK_LIBRARIES(RescaleImageIntensity ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(STAPLEAnalysis STAPLEAnalysis.cxx )
TARGET_LINK_LIBRARIES(STAPLEAnalysis ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(ShapeMorphology ShapeMorphology.cxx )
TARGET_LINK_LIBRARIES(ShapeMorphology ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(SubtractImage SubtractImage.cxx )
TARGET_LINK_LIBRARIES(SubtractImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(ThinTagLines ThinTagLines.cxx )
TARGET_LINK_LIBRARIES(ThinTagLines ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(ThresholdImage ThresholdImage.cxx )
TARGET_LINK_LIBRARIES(ThresholdImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(TileImages TileImages.cxx )
TARGET_LINK_LIBRARIES(TileImages ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(TranslateImage TranslateImage.cxx )
TARGET_LINK_LIBRARIES(TranslateImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(TranslationRegistration TranslationRegistration.cxx )
TARGET_LINK_LIBRARIES(TranslationRegistration ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(WarpAvantsPoints WarpAvantsPoints.cxx )
TARGET_LINK_LIBRARIES(WarpAvantsPoints ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(WarpImage WarpImage.cxx )
TARGET_LINK_LIBRARIES(WarpImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(WellComposeImage WellComposeImage.cxx )
TARGET_LINK_LIBRARIES(WellComposeImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(ZeroCrossingsImage ZeroCrossingsImage.cxx )
TARGET_LINK_LIBRARIES(ZeroCrossingsImage ITKNumerics ITKIO ITKCommon)


