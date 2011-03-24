#!/user/bin/perl -w

$outputDirectory = $ARGV[0];
$uncorrectedHeliumImage = $ARGV[1];

$utilityDirectory = "/home/tustison/Utilities/bin64/";
$antsDirectory = "/home/tustison/ANTS/bin64/";
$shapeIntensityTemplateDirectory = "/home/tustison/Data/HeliumLungStudies/Templates/ANTSShapeIntensityTemplate/Axial/UnbiasedTemplate/";
$shapeIntensityTemplate = $shapeIntensityTemplateDirectory . "OUT3template_biascorrected.nii.gz";
# $shapeIntensityTemplateSegmentation = $shapeIntensityTemplateDirectory . "OUT3template_segmentation.nii.gz";
$pcaShapeTemplateFiles = "/home/tustison/Data/HeliumLungStudies/Templates/PCAShapeTemplate/PCAModes/pca_modes_%03d.nii.gz";



$PID = "run00";

##
# Change image information to have 0x0x0 origin and identity direction
##

$changeInformation = $outputDirectory . "change_image_information_${PID}.nii.gz";

@args = ( $utilityDirectory . "ChangeImageInformation", 3,
  $uncorrectedHeliumImage, $changeInformation, 0, "0x0x0" );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "ChangeImageInformation", 3,
  $changeInformation, $changeInformation, 2 );
system( @args ) == 0 or die "system @args failed: $?";

##
# Get the lung image mask use pca shape template
#   1. affine warp to shape and intensity template
#   2. calculate lung mask
#   3. inverse affine warp lung mask back to original image space
##

$output = $outputDirectory . "registrationToTemplateForPCA";

@args = ( "${antsDirectory}/ANTS", 3,
        "-m", "PR[${shapeIntensityTemplate},${changeInformation},1,4]",
        "-t", "SyN[0.5]",
        "-i", "0",
        "-r", "Gauss[3,0]",
        "-o", "${output}.nii.gz" );
system( @args ) == 0 or die "system @args failed: $?";

$featureImage = $outputDirectory . "feature_image_${PID}.nii.gz";
$backgroundSegmentationImage = $outputDirectory . "background_segmentation_image_${PID}.nii.gz";
$initialLevelSet = $outputDirectory . "initial_levelset_${PID}.nii.gz";
$initialSegmentation = $outputDirectory . "initial_segmentation_${PID}.nii.gz";

@args = ( "${antsDirectory}/WarpImageMultiTransform", 3,
        $changeInformation, $initialSegmentation,
        "-R", $shapeIntensityTemplate,
        "--use-NN",
        "${output}Affine.txt" );
system( @args ) == 0 or die "system @args failed: $?";

@args = ( $utilityDirectory . "GradientMagnitudeImageFilter", 3,
  $initialSegmentation, $featureImage, 1.0 );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "RescaleImageIntensity", 3,
  $featureImage, $featureImage, 0, 1 );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "UnaryOperateImage", 3,
  $featureImage, 'b', 0, $featureImage );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "RescaleImageIntensity", 3,
  $featureImage, $featureImage, 0, 1 );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "SegmentHeliumLungs", 3,
  $initialSegmentation, $backgroundSegmentationImage, 0, 1 );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "ThresholdImage", 3,
  $backgroundSegmentationImage, $backgroundSegmentationImage, 2, 3, 1, 0 );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "BinaryMorphology", 3,
  $backgroundSegmentationImage, $backgroundSegmentationImage, 1, 1, 1 );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "BinaryOperateImages", 3,
  $backgroundSegmentationImage, 'max', $featureImage, $featureImage );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "UnaryOperateImage", 3,
  $featureImage, 's', 0, $featureImage, 0.05, 0.75 );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "RescaleImageIntensity", 3,
  $featureImage, $featureImage, 0, 1 );
system( @args ) == 0 or die "system @args failed: $?";

@args = ( $utilityDirectory . "BinaryMorphology", 3,
  $backgroundSegmentationImage, $initialLevelSet, 1, 3, 1, 1 );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "GenerateDistanceImage", 3,
  $initialLevelSet, $initialLevelSet );
system( @args ) == 0 or die "system @args failed: $?";

@args = ( "${utilityDirectory}/ShapePriorLevelSet", 3,
        $pcaShapeTemplateFiles, 52,
        $featureImage,
        $initialLevelSet,
        $initialSegmentation,
        400, 0.02, 0.5, 1, 20, 1 );
system( @args ) == 0 or die "system @args failed: $?";

@args = ( "${antsDirectory}/WarpImageMultiTransform", 3,
        $initialSegmentation, $initialSegmentation,
        "-R", $changeInformation,
        "--use-NN",
        "-i", "${output}Affine.txt" );
system( @args ) == 0 or die "system @args failed: $?";

#
# Do a slice by slice vesselness filter to remove vessels/airways
#
# $hessianImage = $tmpDirectory . "hessian_mask_${PID}.nii.gz";
#
# @imageInfo = `${utilityDirectory}/GenerateStatisticsFromImage $heliumImage`;
# $begin = index( $imageInfo[1], '[' );
# $end = index( $imageInfo[1], ']' );
# $sizeString = substr( $imageInfo[1], $begin+1, $end - $begin - 1 );
# @sizes = split( /, /, $sizeString );
#
# $smallestSize = $sizes[0];
# $smallestIndex = 0;
# for( $i = 1; $i < @sizes; $i++ )
#   {
#   if( $sizes[$i] < $smallestSize )
#     {
#     $smallestSize = $sizes[$i];
#     $smallestIndex = $i;
#     }
#   }
#
# $smallestSize -= 1;
#
# for( $i = 0; $i <= $smallestSize; $i++ )
#   {
#   $sliceImage = "${tmpDirectory}/slice_${i}.nii.gz";
#   `${utilityDirectory}/ExtractSliceFromImage $heliumImage $sliceImage $smallestIndex $i`;
#   `${utilityDirectory}/HessianBasedFeatures 2 $sliceImage $sliceImage 1 0.5 4 10 0.5 0.5 5 0`;
#   `${utilityDirectory}/OtsuThresholdImage 2 $sliceImage $sliceImage 1 200`;
#   }
# `${utilityDirectory}/ConvertImageSeries $tmpDirectory slice_%d.nii.gz $hessianImage 0 $smallestSize 1`;
# `${utilityDirectory}/ChangeImageInformation 3 $hessianImage $hessianImage 4 $heliumImage`;
# `${utilityDirectory}/ThresholdImage 3 $hessianImage $hessianImage 1 1 1 0`;
# `${utilityDirectory}/BinaryOperateImages 3 $hessianImage x $initialSegmentation $initialSegmentation`;







$prefix[0] = "leftandright";
$mask[0] = $outputDirectory . "/mask_" . $prefix[0] . "_" . $PID . ".nii.gz";

@args = ( "${utilityDirectory}/SplitHeliumLungs", 3,
       $initialSegmentation, $mask[0] );
system( @args ) == 0 or die "system @args failed: $?";

##
# Bias field correction
##

$heliumImage = $outputDirectory . "/he3_bias_corrected_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "InhomogeneityCorrectImage", 3,
           '-i', $changeInformation,
           '-x', $initialSegmentation,
           '-s', 1,
           '-c', '[50x50x50,0.00001]',
           '-b', '[100,3,0,0]',
           '-t', '[0.1,0.1,200]',
           '-o', "[${heliumImage},${outputDirectory}/bias_field_${PID}.nii.gz]" );
system( @args ) == 0 or die "system @args failed: $?";


##
# Ventilation analysis
##

$apocritaImage = $outputDirectory . "apocrita_5class_both_${PID}.nii.gz";

@args = ( $utilityDirectory . "ApocritaSegmentation", 3,
          "-i", "otsu[${heliumImage},5]",
          "-n", 5,
          "-x", $initialSegmentation,
          "-m", "[0.1,2x2x1,0.0,0.0]",
          "-o", $apocritaImage );
system( @args ) == 0 or die "system @args failed: $?";

## Threshold
$apocritaLessVentilatedImage = $outputDirectory . "apocrita_lessventilated_tmp.nii.gz";
$apocritaMoreVentilatedImage = $outputDirectory . "apocrita_moreventilated_tmp.nii.gz";
$apocritaVentilatedImage = $outputDirectory . "apocrita_ventilated.nii.gz";

$command = $utilityDirectory . "ThresholdImage";
@args = ( $command, 3, $apocritaImage, $apocritaLessVentilatedImage, 1, 2, 1, 0 );
system( @args ) == 0 or die "system @args failed: $?";
$command = $utilityDirectory . "ThresholdImage";
@args = ( $command, 3, $apocritaImage, $apocritaMoreVentilatedImage, 3, 5, 2, 0 );
system( @args ) == 0 or die "system @args failed: $?";
$command = $utilityDirectory . "BinaryOperateImages";
@args = ( $command, 3, $apocritaLessVentilatedImage, '+',
  $apocritaMoreVentilatedImage, $apocritaVentilatedImage );
system( @args ) == 0 or die "system @args failed: $?";

$apocritaLeftImage = $outputDirectory  . "apocrita_left_tmp.nii.gz";
$apocritaRightImage = $outputDirectory  . "apocrita_right_tmp.nii.gz";

$command = $utilityDirectory . "ExtractMaskImage";
@args = ( $command, 3, $apocritaVentilatedImage, $mask[0], $apocritaLeftImage, 2, 0 );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $command, 3, $apocritaVentilatedImage, $mask[0], $apocritaRightImage, 3, 0 );
system( @args ) == 0 or die "system @args failed: $?";


$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_2class_ventilation_ratio_both_1_" . $PID . ".txt";
$out = `$command 3 $apocritaVentilatedImage 1 $initialSegmentation 1 > $output`;
 $command = $utilityDirectory . "CalculateFirstOrderStatisticsFromImage";
 $output = $outputDirectory . "/apocrita_2class_ventilation_ratio_both_statistics_1_" . $PID . ".txt";
 $out = `$command 3 $heliumImage $apocritaVentilatedImage 1 > $output`;

$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_2class_ventilation_ratio_both_1_left_" . $PID . ".txt";
$out = `$command 3 $apocritaLeftImage 1 $mask[0] 2 > $output`;
$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_2class_ventilation_ratio_both_1_right_" . $PID . ".txt";
$out = `$command 3 $apocritaRightImage 1 $mask[0] 3 > $output`;

$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_2class_ventilation_ratio_both_2_" . $PID . ".txt";
$out = `$command 3 $apocritaVentilatedImage 2 $initialSegmentation 1 > $output`;
 $command = $utilityDirectory . "CalculateFirstOrderStatisticsFromImage";
 $output = $outputDirectory . "/apocrita_2class_ventilation_ratio_both_statistics_2_" . $PID . ".txt";
 $out = `$command 3 $heliumImage $apocritaVentilatedImage 2 > $output`;

$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_2class_ventilation_ratio_both_2_left_" . $PID . ".txt";
$out = `$command 3 $apocritaLeftImage 2 $mask[0] 2 > $output`;
$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_2class_ventilation_ratio_both_2_right_" . $PID . ".txt";
$out = `$command 3 $apocritaRightImage 2 $mask[0] 3 > $output`;

unlink( $apocritaLeftImage );
unlink( $apocritaRightImage );
unlink( $apocritaLessVentilatedImage );
unlink( $apocritaMoreVentilatedImage );
unlink( $apocritaVentilatedImage );

return;





