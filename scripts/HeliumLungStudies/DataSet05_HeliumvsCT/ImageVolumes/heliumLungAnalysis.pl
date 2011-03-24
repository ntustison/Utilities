#!/user/bin/perl -w

$outputDirectory = $ARGV[0];
$uncorrectedHeliumImage = $ARGV[1];

$utilityDirectory = "/home/tustison/Utilities/bin64/";
$antsDirectory = "/home/tustison/ANTS/bin64/";
$shapeIntensityTemplateDirectory = "/home/tustison/Data/HeliumLungStudies/Templates/ANTSShapeIntensityTemplate/Axial/UnbiasedTemplate/";
$shapeIntensityTemplate = $shapeIntensityTemplateDirectory . "OUT3template_biascorrected.nii.gz";
# $shapeIntensityTemplateSegmentation = $shapeIntensityTemplateDirectory . "OUT3template_segmentation.nii.gz";
$pcaShapeTemplateFiles = "/home/tustison/Data/HeliumLungStudies/Templates/PCAShapeTemplate/PCAModes/pca_modes_%03d.nii.gz";

$numberOfHistogramBins = 200;

$PID = "run03";

##
# Change image information to have 0x0x0 origin and identity direction
##

$changeInformation = $outputDirectory . "change_image_information_${PID}.nii.gz";

@args = ( $utilityDirectory . "ChangeImageInformation", 3,
  $uncorrectedHeliumImage, $changeInformation, 0, "0x0x0" );
#system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "ChangeImageInformation", 3,
  $changeInformation, $changeInformation, 2 );
#system( @args ) == 0 or die "system @args failed: $?";

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
#system( @args ) == 0 or die "system @args failed: $?";

$featureImage = $outputDirectory . "feature_image_${PID}.nii.gz";
$backgroundSegmentationImage = $outputDirectory . "background_segmentation_image_${PID}.nii.gz";
$initialLevelSet = $outputDirectory . "initial_levelset_${PID}.nii.gz";
$initialSegmentation = $outputDirectory . "initial_segmentation_${PID}.nii.gz";

@args = ( "${antsDirectory}/WarpImageMultiTransform", 3,
        $changeInformation, $initialSegmentation,
        "-R", $shapeIntensityTemplate,
        "--use-NN",
        "${output}Affine.txt" );
#system( @args ) == 0 or die "system @args failed: $?";

# @args = ( $utilityDirectory . "GradientMagnitudeImageFilter", 3,
#   $initialSegmentation, $featureImage, 1.0 );
# system( @args ) == 0 or die "system @args failed: $?";
# @args = ( $utilityDirectory . "RescaleImageIntensity", 3,
#   $featureImage, $featureImage, 0, 1 );
# system( @args ) == 0 or die "system @args failed: $?";
# @args = ( $utilityDirectory . "UnaryOperateImage", 3,
#   $featureImage, 'b', 0, $featureImage );
# system( @args ) == 0 or die "system @args failed: $?";
# @args = ( $utilityDirectory . "RescaleImageIntensity", 3,
#   $featureImage, $featureImage, 0, 1 );
# system( @args ) == 0 or die "system @args failed: $?";
# @args = ( $utilityDirectory . "SegmentHeliumLungs", 3,
#   $initialSegmentation, $backgroundSegmentationImage, 0, 1 );
# system( @args ) == 0 or die "system @args failed: $?";
# @args = ( $utilityDirectory . "ThresholdImage", 3,
#   $backgroundSegmentationImage, $backgroundSegmentationImage, 2, 3, 1, 0 );
# system( @args ) == 0 or die "system @args failed: $?";
# @args = ( $utilityDirectory . "BinaryMorphology", 3,
#   $backgroundSegmentationImage, $backgroundSegmentationImage, 1, 1, 1 );
# system( @args ) == 0 or die "system @args failed: $?";
# @args = ( $utilityDirectory . "BinaryOperateImages", 3,
#   $backgroundSegmentationImage, 'max', $featureImage, $featureImage );
# system( @args ) == 0 or die "system @args failed: $?";
# @args = ( $utilityDirectory . "UnaryOperateImage", 3,
#   $featureImage, 's', 0, $featureImage, 0.05, 0.75 );
# system( @args ) == 0 or die "system @args failed: $?";
# @args = ( $utilityDirectory . "RescaleImageIntensity", 3,
#   $featureImage, $featureImage, 0, 1 );
# system( @args ) == 0 or die "system @args failed: $?";
#
# @args = ( $utilityDirectory . "BinaryMorphology", 3,
#   $backgroundSegmentationImage, $initialLevelSet, 1, 3, 1, 1 );
# system( @args ) == 0 or die "system @args failed: $?";
# @args = ( $utilityDirectory . "GenerateDistanceImage", 3,
#   $initialLevelSet, $initialLevelSet );
# system( @args ) == 0 or die "system @args failed: $?";
#
# @args = ( "${utilityDirectory}/ShapePriorLevelSet", 3,
#         $pcaShapeTemplateFiles, 52,
#         $featureImage,
#         $initialLevelSet,
#         $initialSegmentation,
#         400, 0.02, 0.5, 1, 20, 1 );
# system( @args ) == 0 or die "system @args failed: $?";
#
# @args = ( "${antsDirectory}/WarpImageMultiTransform", 3,
#         $initialSegmentation, $initialSegmentation,
#         "-R", $changeInformation,
#         "--use-NN",
#         "-i", "${output}Affine.txt" );
# system( @args ) == 0 or die "system @args failed: $?";

$prefix[0] = "leftandright";
$mask[0] = $outputDirectory . "/mask_" . $prefix[0] . "_" . $PID . ".nii.gz";

@args = ( "${utilityDirectory}/SplitHeliumLungs", 3,
       $initialSegmentation, $mask[0] );
#system( @args ) == 0 or die "system @args failed: $?";

##
# Bias field correction
##

$heliumImage = $outputDirectory . "/he3_bias_corrected_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "InhomogeneityCorrectImage", 3,
           '-i', $changeInformation,
           '-x', $initialSegmentation,
           '-s', 1,
           '-c', '[50,0.00001]',
           '-b', '[4,1,3]',
           '-t', '[0.15,0.1,1,200]',
           '-o', "[${heliumImage},${outputDirectory}/bias_field_${PID}.nii.gz]" );
#system( @args ) == 0 or die "system @args failed: $?";

##
# Create mask images
##


#   $prefix[0] = "leftandright";
#   $output = $outputDirectory . "registrationToTemplate";
#   $mask[0] = $outputDirectory . "/mask_" . $prefix[0] . "_" . $PID . ".nii.gz";
#
#   @args = ( "${antsDirectory}/ANTS", 3,
#           "-m", "PR[${template},${heliumImage},1,4]",
#           "-t", "SyN[0.5]",
#           "-i", "100x50x0",
#           "-r", "Gauss[3,2]",
#           "-o", "${output}.nii.gz" );
#   system( @args ) == 0 or die "system @args failed: $?";
#
#   @args = ( "${antsDirectory}/ComposeMultiTransform", 3,
#           "${output}TotalInverseWarp.nii.gz",
#           "-R", $heliumImage,
#           "-i", "${output}Affine.txt",
#           "${output}InverseWarp.nii.gz" );
#   system( @args ) == 0 or die "system @args failed: $?";
#
#   @args = ( "${antsDirectory}/WarpImageMultiTransform", 3,
#           "${templateSegmentation}", $mask[0],
#           "-R", $heliumImage,
#           "--use-NN",
#           "${output}TotalInverseWarp.nii.gz" );
#   system( @args ) == 0 or die "system @args failed: $?";

##
# Ventilation analysis
##

$apocritaImage = $outputDirectory . "apocrita_2class_both_${PID}.nii.gz";

# @args = ( $utilityDirectory . "ApocritaSegmentation", 3,
#           "-i", "otsu[${heliumImage},2]",
#           "-n", 10,
# #           "-x", $initialSegmentation,
#           "-m", "[0.1,1x1x1,0.1,0.5]",
#           "-o", $apocritaImage );
# system( @args ) == 0 or die "system @args failed: $?";
# @args = ( $utilityDirectory . "BinaryOperateImages", 3, $initialSegmentation,
#   'x', $apocritaImage, $apocritaImage );
# system( @args ) == 0 or die "system @args failed: $?";

$apocritaLeftImage = $outputDirectory  . "apocrita_left_tmp.nii.gz";
$apocritaRightImage = $outputDirectory  . "apocrita_right_tmp.nii.gz";

$command = $utilityDirectory . "ExtractMaskImage";
@args = ( $command, 3, $apocritaImage, $mask[0], $apocritaLeftImage, 2, 0 );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $command, 3, $apocritaImage, $mask[0], $apocritaRightImage, 3, 0 );
system( @args ) == 0 or die "system @args failed: $?";


$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_2class_ventilation_ratio_both_1_" . $PID . ".txt";
$out = `$command 3 $apocritaImage 1 $initialSegmentation 1 > $output`;
 $command = $utilityDirectory . "CalculateFirstOrderStatisticsFromImage";
 $output = $outputDirectory . "/apocrita_2class_ventilation_ratio_both_statistics_1_" . $PID . ".txt";
 $out = `$command 3 $heliumImage $apocritaImage 1 > $output`;

$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_2class_ventilation_ratio_both_1_left_" . $PID . ".txt";
$out = `$command 3 $apocritaLeftImage 1 $mask[0] 2 > $output`;
$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_2class_ventilation_ratio_both_1_right_" . $PID . ".txt";
$out = `$command 3 $apocritaRightImage 1 $mask[0] 3 > $output`;



$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_2class_ventilation_ratio_both_2_" . $PID . ".txt";
$out = `$command 3 $apocritaImage 2 $initialSegmentation 1 > $output`;
 $command = $utilityDirectory . "CalculateFirstOrderStatisticsFromImage";
 $output = $outputDirectory . "/apocrita_2class_ventilation_ratio_both_statistics_2_" . $PID . ".txt";
 $out = `$command 3 $heliumImage $apocritaImage 2 > $output`;

$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_2class_ventilation_ratio_both_2_left_" . $PID . ".txt";
$out = `$command 3 $apocritaLeftImage 2 $mask[0] 2 > $output`;
$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_2class_ventilation_ratio_both_2_right_" . $PID . ".txt";
$out = `$command 3 $apocritaRightImage 2 $mask[0] 3 > $output`;

unlink( $apocritaLeftImage );
unlink( $apocritaRightImage );

##
## also calculate 3 classes
##


$apocritaImage = $outputDirectory . "apocrita_3class_both_${PID}.nii.gz";

 @args = ( $utilityDirectory . "ApocritaSegmentation", 3,
          "-i", "otsu[${heliumImage},3]",
          "-n", 10,
#           "-x", $initialSegmentation,
          "-m", "[0.1,1x1x1,0.1,0.5]",
          "-o", $apocritaImage );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "BinaryOperateImages", 3, $initialSegmentation,
  'x', $apocritaImage, $apocritaImage );
system( @args ) == 0 or die "system @args failed: $?";

$apocritaLeftImage = $outputDirectory  . "apocrita_left_tmp.nii.gz";
$apocritaRightImage = $outputDirectory  . "apocrita_right_tmp.nii.gz";

$command = $utilityDirectory . "ExtractMaskImage";
@args = ( $command, 3, $apocritaImage, $mask[0], $apocritaLeftImage, 2, 0 );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $command, 3, $apocritaImage, $mask[0], $apocritaRightImage, 3, 0 );
system( @args ) == 0 or die "system @args failed: $?";



$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_3class_ventilation_ratio_both_1_" . $PID . ".txt";
$out = `$command 3 $apocritaImage 1 $initialSegmentation 1 > $output`;
 $command = $utilityDirectory . "CalculateFirstOrderStatisticsFromImage";
 $output = $outputDirectory . "/apocrita_3class_ventilation_ratio_both_statistics_1_" . $PID . ".txt";
 $out = `$command 3 $heliumImage $apocritaImage 1 > $output`;

$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_3class_ventilation_ratio_both_1_left_" . $PID . ".txt";
$out = `$command 3 $apocritaLeftImage 1 $mask[0] 2 > $output`;
$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_3class_ventilation_ratio_both_1_right_" . $PID . ".txt";
$out = `$command 3 $apocritaRightImage 1 $mask[0] 3 > $output`;


$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_3class_ventilation_ratio_both_2_" . $PID . ".txt";
$out = `$command 3 $apocritaImage 2 $initialSegmentation 1 > $output`;
 $command = $utilityDirectory . "CalculateFirstOrderStatisticsFromImage";
 $output = $outputDirectory . "/apocrita_3class_ventilation_ratio_both_statistics_2_" . $PID . ".txt";
 $out = `$command 3 $heliumImage $apocritaImage 2 > $output`;

$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_3class_ventilation_ratio_both_2_left_" . $PID . ".txt";
$out = `$command 3 $apocritaLeftImage 2 $mask[0] 2 > $output`;
$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_3class_ventilation_ratio_both_2_right_" . $PID . ".txt";
$out = `$command 3 $apocritaRightImage 2 $mask[0] 3 > $output`;

$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_3class_ventilation_ratio_both_3_" . $PID . ".txt";
$out = `$command 3 $apocritaImage 2 $initialSegmentation 1 > $output`;
 $command = $utilityDirectory . "CalculateFirstOrderStatisticsFromImage";
 $output = $outputDirectory . "/apocrita_3class_ventilation_ratio_both_statistics_3_" . $PID . ".txt";
 $out = `$command 3 $heliumImage $apocritaImage 2 > $output`;

$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_3class_ventilation_ratio_both_3_left_" . $PID . ".txt";
$out = `$command 3 $apocritaLeftImage 3 $mask[0] 2 > $output`;
$command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
$output = $outputDirectory . "/apocrita_3class_ventilation_ratio_both_3_right_" . $PID . ".txt";
$out = `$command 3 $apocritaRightImage 3 $mask[0] 3 > $output`;

unlink( $apocritaLeftImage );
unlink( $apocritaRightImage );

return;






















##
# Create remainder of masks
##

$prefix[1] = "bothlungs";
$mask[1] = $outputDirectory . "/mask_" . $prefix[1] . "_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "ThresholdImage", 3, $mask[0], $mask[1], 2, 3, 1, 0 );
system( @args ) == 0 or die "system @args failed: $?";


$distance = -10;  # mm

$prefix[2] = "inner_core_2";
$mask[2] = $outputDirectory . "/mask_" . $prefix[2] . "_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "GenerateDistanceImage", 3, $mask[0], $mask[2], 2 );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "ThresholdImage", 3, $mask[2], $mask[2], -1000000000, $distance, 2, 0 );
system( @args ) == 0 or die "system @args failed: $?";

$prefix[4] = "inner_core_3";
$mask[4] = $outputDirectory . "/mask_" . $prefix[4] . "_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "GenerateDistanceImage", 3, $mask[0], $mask[4], 3 );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "ThresholdImage", 3, $mask[4], $mask[4], -1000000000, $distance, 3, 0 );
system( @args ) == 0 or die "system @args failed: $?";

$prefix[3] = "outer_rind_2";
$mask[3] = $outputDirectory . "/mask_" . $prefix[3] . "_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "BinaryOperateImages", 3, $mask[0], '-', $mask[2], $mask[3] );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "ThresholdImage", 3, $mask[3], $mask[3], 2, 2, 2, 0 );
system( @args ) == 0 or die "system @args failed: $?";

$prefix[5] = "outer_rind_3";
$mask[5] = $outputDirectory . "/mask_" . $prefix[5] . "_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "BinaryOperateImages", 3, $mask[0], '-', $mask[4], $mask[5] );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "ThresholdImage", 3, $mask[5], $mask[5], 3, 3, 3, 0 );
system( @args ) == 0 or die "system @args failed: $?";

for( $i = 2; $i <= 5; $i++ )
  {
  @args = ( $utilityDirectory . "ThresholdImage", 3, $mask[$i], $mask[$i], 0, 0, 0, 1 );
  system( @args ) == 0 or die "system @args failed: $?";
  }

$prefix[6] = "inner_core_both";
$mask[6] = $outputDirectory . "/mask_" . $prefix[6] . "_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "BinaryOperateImages", 3, $mask[2], '+', $mask[4], $mask[6] );
system( @args ) == 0 or die "system @args failed: $?";

$prefix[7] = "outer_rind_both";
$mask[7] = $outputDirectory . "/mask_" . $prefix[7] . "_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "BinaryOperateImages", 3, $mask[3], '+', $mask[5], $mask[7] );
system( @args ) == 0 or die "system @args failed: $?";

$prefix[8] = "division";
$mask[8] = $outputDirectory . "/mask_" . $prefix[8] . "_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "DivideLungs", $mask[0], $mask[8], 2, 2 );
system( @args ) == 0 or die "system @args failed: $?";

for( $i = 1; $i <= 4; $i++ )
  {
  $prefix[8+$i] = "division_" . $i;
  $mask[8+$i] = $outputDirectory . "/mask_" . $prefix[8+$i] . "_" . $PID . ".nii.gz";
  @args = ( $utilityDirectory . "ThresholdImage", 3, $mask[8], $mask[8+$i], $i, $i, 1, 0 );
  system( @args ) == 0 or die "system @args failed: $?";
  }

for( $i = 2; $i <= 3; $i++ )
  {
  $prefix[11+$i] = "leftright_" . $i;
  $mask[11+$i] = $outputDirectory . "/mask_" . $prefix[11+$i] . "_" . $PID . ".nii.gz";
  @args = ( $utilityDirectory . "ThresholdImage", 3, $mask[0], $mask[11+$i], $i, $i, 1, 0 );
  system( @args ) == 0 or die "system @args failed: $?";
  }

##
# Create fractal images
##

$fractalRadius[0] = "1x1x1";
#  $fractalRadius[1] = "1x1x0";
#  $fractalRadius[2] = "2x2x0";

for ( $i = 0; $i < @fractalRadius; $i++ )
  {
  $fractalImage[$i] = $outputDirectory . "/fractalimage_" . $fractalRadius[$i] . "_" . $PID . ".nii.gz";

  @args = ( $utilityDirectory . "GenerateFractalImage", 3, $heliumImage,
    $fractalImage[$i], $fractalRadius[$i], $mask[1], 1 );
  system( @args ) == 0 or die "system @args failed: $?"
  }

##
# Calculate features
##

# @whichMasks = ( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 );
@whichMasks = ( 1, 6, 7 );

for( $m = 0; $m < @whichMasks; $m++ )
  {

  $i = $whichMasks[$m];

  $command = $utilityDirectory . "CalculateFirstOrderStatisticsFromImage";
  $output = $outputDirectory . "/stats_" . $prefix[$i] . "_" . $PID . ".txt";
  $out = `$command 3 $heliumImage $mask[$i] 1 $numberOfHistogramBins > $output`;

  $command = $utilityDirectory . "GenerateCooccurrenceMeasures";
  $output = $outputDirectory . "/cm_axial_" . $prefix[$i] . "_" . $PID . ".txt";
  $out = `$command 3 $heliumImage 4 1x0x0 1x1x0 0x1x0 -1x-1x0 $mask[$i] 1 > $output`;

  $command = $utilityDirectory . "GenerateCooccurrenceMeasures";
  $output = $outputDirectory . "/cm_" . $prefix[$i] . "_" . $PID . ".txt";
  $out = `$command 3 $heliumImage 13 1x0x0 1x1x0 0x1x0 -1x-1x0 1x1x1 0x1x1 -1x1x1 1x0x1 0x0x1 -1x0x1 1x-1x1 0x-1x1 -1x-1x1 $mask[$i] 1 > $output`;

  $command = $utilityDirectory . "GenerateRunLengthMeasures";
  $output = $outputDirectory . "/rlm_axial_" . $prefix[$i] . "_" . $PID . ".txt";
  $out = `$command 3 $heliumImage 4 1x0x0 1x1x0 0x1x0 -1x-1x0 $mask[$i] 1 > $output`;

  $command = $utilityDirectory . "GenerateRunLengthMeasures";
  $output = $outputDirectory . "/rlm_" . $prefix[$i] . "_" . $PID . ".txt";
  $out = `$command 3 $heliumImage 13 1x0x0 1x1x0 0x1x0 -1x-1x0 1x1x1 0x1x1 -1x1x1 1x0x1 0x0x1 -1x0x1 1x-1x1 0x-1x1 -1x-1x1 $mask[$i] 1 > $output`;

  ##
  # Calculate features from fractal images
  ##

  for( $j = 0; $j < @fractalImage; $j++ )
    {
    $command = $utilityDirectory . "CalculateFirstOrderStatisticsFromImage";
    $output = $outputDirectory . "/stats_fractal_" . $prefix[$i] . "_" . $PID . ".txt";
    $out = `$command 3 $fractalImage[$j] $mask[$i] 1 $numberOfHistogramBins > $output`;

    $command = $utilityDirectory . "GenerateCooccurrenceMeasures";
    $output = $outputDirectory . "/cm_axial_fractal_" . $prefix[$i] . "_" . $PID . ".txt";
    $out = `$command 3 $fractalImage[$j] 4 1x0x0 1x1x0 0x1x0 -1x-1x0 $mask[$i] 1 > $output`;

    $command = $utilityDirectory . "GenerateCooccurrenceMeasures";
    $output = $outputDirectory . "/cm_fractal_" . $prefix[$i] . "_" . $PID . ".txt";
    $out = `$command 3 $fractalImage[$j] 13 1x0x0 1x1x0 0x1x0 -1x-1x0 1x1x1 0x1x1 -1x1x1 1x0x1 0x0x1 -1x0x1 1x-1x1 0x-1x1 -1x-1x1 $mask[$i] 1 > $output`;

    $command = $utilityDirectory . "GenerateRunLengthMeasures";
    $output = $outputDirectory . "/rlm_axial_fractal_" . $prefix[$i] . "_" . $PID . ".txt";
    $out = `$command 3 $fractalImage[$j] 4 1x0x0 1x1x0 0x1x0 -1x-1x0 $mask[$i] 1 > $output`;

    $command = $utilityDirectory . "GenerateRunLengthMeasures";
    $output = $outputDirectory . "/rlm_fractal_" . $prefix[$i] . "_" . $PID . ".txt";
    $out = `$command 3 $fractalImage[$j] 13 1x0x0 1x1x0 0x1x0 -1x-1x0 1x1x1 0x1x1 -1x1x1 1x0x1 0x0x1 -1x0x1 1x-1x1 0x-1x1 -1x-1x1 $mask[$i] 1 > $output`;
    }
  }

##
# Noise statistics
##

$radius = 2;
$noiseMask = $outputDirectory . "/noise_mask_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "BinaryMorphology", 3, $mask[1], $noiseMask, 0, $radius, 0, 1 );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "BinaryOperateImages", 3, $noiseMask, '-', $mask[1], $noiseMask );
system( @args ) == 0 or die "system @args failed: $?";

$command = $utilityDirectory . "CalculateFirstOrderStatisticsFromImage";
$output = $outputDirectory . "/stats_noise_" . $PID . ".txt";
$out = `$command 3 $heliumImage $noiseMask 1 $numberOfHistogramBins > $output`;



