#!/user/bin/perl -w

use POSIX;

$outputDirectory = $ARGV[0];
$inputImage = $ARGV[1];

$utilityDirectory = "/home/tustison/Utilities/bin64/";
$antsDirectory = "/home/tustison/ANTS/bin64/";

$correctedImage = $outputDirectory . "/proton_corrected.nii.gz";
$biasField = $outputDirectory . "/proton_biasfield.nii.gz";
if( ! -e $correctedImage )
  {
  @args = ( $utilityDirectory . "InhomogeneityCorrectImage", 3,
             '-i', $inputImage,
             '-s', 1,
             '-c', '[100x50x50,0.00001]',
             '-b', '[100,3]',
             '-t', '[0.15,0.01,200]',
             '-o', "[${correctedImage},${biasField}]" );
  system( @args ) == 0 or die "system @args failed: $?";
  }


$atroposImage = $outputDirectory . "/proton_atropos.nii.gz";

if( ! -e $atroposImage )
  {
  @args = ( "${utilityDirectory}AtroposSegmentation", 3,
    '-i', "KMeans[2]",
    '-a', "${correctedImage}",
    '-c', '[5,0.00001]',
    '-k', 'Gaussian',
    '-m', '[10.0,2x2x1]',
    '-w', 'BoxPlot[0.25,0.75,1.5]',
    '-o', "${atroposImage}" );
  system( @args ) == 0 or die "system @args failed: $?";
  }

$closedMask = $outputDirectory . "/proton_closed.nii.gz";
@args = ( "${utilityDirectory}/BinaryMorphology", 3, $atroposImage, $closedMask,
  2, 3, 0, 2 );
system( @args ) == 0 or die "system @args failed: $?";

$backgroundMask = $outputDirectory . "/proton_background.nii.gz";
@args = ( "${utilityDirectory}/FloodFill", 3, $closedMask, $backgroundMask,
  0.5, 1.5, 0, '0x0x0' );
system( @args ) == 0 or die "system @args failed: $?";

 @args = ( "${utilityDirectory}/BinaryOperateImages", 3, $atroposImage, '-',
   $backgroundMask, $atroposImage );
 system( @args ) == 0 or die "system @args failed: $?";

 @args = ( "${utilityDirectory}/ThresholdImage", 3, $atroposImage, $atroposImage,
   1, 1, 2, 0 );
 system( @args ) == 0 or die "system @args failed: $?";

 @args = ( "${utilityDirectory}/FloodFill", 3, $atroposImage, $atroposImage,
   -0.5, 0.5, 0, '0x0x0' );
 system( @args ) == 0 or die "system @args failed: $?";

 @args = ( "${utilityDirectory}/ThresholdImage", 3, $atroposImage, $atroposImage,
   1, 1, 0, 1 );
 system( @args ) == 0 or die "system @args failed: $?";

 @args = ( "${utilityDirectory}/GetConnectedComponents", 3, $atroposImage, $atroposImage, 0 );
 system( @args ) == 0 or die "system @args failed: $?";

 @args = ( "${utilityDirectory}/ThresholdImage", 3, $atroposImage, $atroposImage, 1, 1, 1, 0 );
 system( @args ) == 0 or die "system @args failed: $?";
