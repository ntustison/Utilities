#!/user/bin/perl -w

$outputDirectory = $ARGV[0];
$heliumImage = $ARGV[1];

$utilityDirectory = "/mnt/data1/tustison/Utilities/bin/";

$numberOfHistogramBins = 200;

$PID = "run01";


## 
# Create mask images
##

$prefix[0] = "leftandright";
$mask[0] = $outputDirectory . "/mask_" . $prefix[0] . "_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "SegmentHeliumLungs", 3, $heliumImage, $mask[0] );
system( @args ) == 0 or die "system @args failed: $?";

##
# Ventilation analysis
##

for( $i = 1; $i < 5; $i++ )
  {  
  $erodedMaskImage = $outputDirectory . "/mask_eroded_" . $PID . ".nii.gz";
  $otsuImage2 = $outputDirectory . "/otsu_" . $i . "_2_" . $PID . ".nii.gz"; 
  $otsuImage3 = $outputDirectory . "/otsu_" . $i . "_3_" . $PID . ".nii.gz"; 

  $otsuImageBoth = $outputDirectory . "/otsu_" . $i . "_both_" . $PID . ".nii.gz";   
  
  @args = ( $utilityDirectory . "BinaryMorphology", 3, $mask[0], $erodedMaskImage, 
    1, 2, 1, 2 );
  system( @args ) == 0 or die "system @args failed: $?";
  @args = ( $utilityDirectory . "BinaryMorphology", 3, $erodedMaskImage, $erodedMaskImage,  
    1, 2, 1, 3 );
  system( @args ) == 0 or die "system @args failed: $?";

  @args = ( $utilityDirectory . "OtsuThresholdImage", 3, $heliumImage, $otsuImage2,  
    $i, $numberOfHistogramBins, $erodedMaskImage, 2 );
  system( @args ) == 0 or die "system @args failed: $?";
  $command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
  $output = $outputDirectory . "/otsu_ventilation_ratio_" . $i . "_2_" . $PID . ".txt";
  $out = `$command 3 $otsuImage2 2 $erodedMaskImage 2 > $output`;

  @args = ( $utilityDirectory . "OtsuThresholdImage", 3, $heliumImage, $otsuImage3,  
    $i, $numberOfHistogramBins, $erodedMaskImage, 3 );
  system( @args ) == 0 or die "system @args failed: $?";
  $command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
  $output = $outputDirectory . "/otsu_ventilation_ratio_" . $i . "_3_" . $PID . ".txt";
  $out = `$command 3 $otsuImage3 3 $erodedMaskImage 3 > $output`;

  @args = ( $utilityDirectory . "ThresholdImage", 3, $erodedMaskImage, $erodedMaskImage,  
    2, 3, 1, 0 );
  system( @args ) == 0 or die "system @args failed: $?";
  @args = ( $utilityDirectory . "OtsuThresholdImage", 3, $heliumImage, $otsuImageBoth,  
    $i, $numberOfHistogramBins, $erodedMaskImage, 1 );
  system( @args ) == 0 or die "system @args failed: $?";
  $command = $utilityDirectory . "CalculateVolumeFromBinaryImage";
  $output = $outputDirectory . "/otsu_ventilation_ratio_" . $i . "_both_" . $PID . ".txt";
  $out = `$command 3 $otsuImageBoth 1 $erodedMaskImage 1 > $output`;

  @args = ( "rm", $erodedMaskImage );
  system( @args ) == 0 or die "system @args failed: $?";
  }

##
# Create remainder of masks
##

$prefix[1] = "bothlungs";
$mask[1] = $outputDirectory . "/mask_" . $prefix[1] . "_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "ThresholdImage", 3, $mask[0], $mask[1], 2, 3, 1, 0 );
system( @args ) == 0 or die "system @args failed: $?";

$radius = 3;

$prefix[2] = "inner_core_2";
$mask[2] = $outputDirectory . "/mask_" . $prefix[2] . "_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "BinaryMorphology", 3, $mask[0], $mask[2], 1, $radius, 1, 2 );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "ThresholdImage", 3, $mask[2], $mask[2], 2, 2, 2, 0 );
system( @args ) == 0 or die "system @args failed: $?";

$prefix[3] = "outer_rind_2";
$mask[3] = $outputDirectory . "/mask_" . $prefix[3] . "_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "SubtractImage", 3, $mask[0], $mask[2], $mask[3] );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "ThresholdImage", 3, $mask[3], $mask[3], 2, 2, 2, 0 );
system( @args ) == 0 or die "system @args failed: $?";

$prefix[4] = "inner_core_3";
$mask[4] = $outputDirectory . "/mask_" . $prefix[4] . "_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "BinaryMorphology", 3, $mask[0], $mask[4], 1, $radius, 1, 3 );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "ThresholdImage", 3, $mask[4], $mask[4], 3, 3, 3, 0 );
system( @args ) == 0 or die "system @args failed: $?";

$prefix[5] = "outer_rind_3";
$mask[5] = $outputDirectory . "/mask_" . $prefix[5] . "_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "SubtractImage", 3, $mask[0], $mask[4], $mask[5] );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "ThresholdImage", 3, $mask[5], $mask[5], 3, 3, 1, 0 );
system( @args ) == 0 or die "system @args failed: $?";

for( $i = 2; $i <= 5; $i++ )
  {
  @args = ( $utilityDirectory . "ThresholdImage", 3, $mask[$i], $mask[$i], 0, 0, 0, 1 );
  system( @args ) == 0 or die "system @args failed: $?";
  }

$prefix[6] = "inner_core_both";
$mask[6] = $outputDirectory . "/mask_" . $prefix[6] . "_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "AddImage", 3, $mask[2], $mask[4], $mask[6] );
system( @args ) == 0 or die "system @args failed: $?";

$prefix[7] = "outer_rind_both";
$mask[7] = $outputDirectory . "/mask_" . $prefix[7] . "_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "AddImage", 3, $mask[3], $mask[5], $mask[7] );
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
$fractalRadius[1] = "1x1x0";
$fractalRadius[2] = "2x2x0";

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

@whichMasks = ( 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14 );

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
@args = ( $utilityDirectory . "SubtractImage", 3, $noiseMask, $mask[1], $noiseMask );
system( @args ) == 0 or die "system @args failed: $?";

$command = $utilityDirectory . "CalculateFirstOrderStatisticsFromImage";
$output = $outputDirectory . "/stats_noise_" . $PID . ".txt";
$out = `$command 3 $heliumImage $noiseMask 1 $numberOfHistogramBins > $output`; 



