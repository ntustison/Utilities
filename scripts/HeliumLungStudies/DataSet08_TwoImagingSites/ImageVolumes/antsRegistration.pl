#!/usr/bin64/perl

##
# Register point-sets
##

$utilityDirectory = "/home/tustison/Utilities/bin64/";
$antsDirectory = "/home/tustison/ANTS/bin64/";

##
# Perform refinement with point-set registration
##

$PREDIRECTORY = $ARGV[0];
$POSTDIRECTORY = $ARGV[1];
$PRESERIES = $ARGV[2];
$POSTSERIES = $ARGV[3];

#print $PREDIRECTORY . "\n";
#print $POSTDIRECTORY . "\n";
#print $PRESERIES . "\n";
#print $POSTSERIES . "\n";


$fixedImage = $POSTDIRECTORY . "/" . $POSTSERIES . ".nii.gz";
$movingImage = $PREDIRECTORY . "/" . $PRESERIES . ".nii.gz";

$outputDirectory = $POSTDIRECTORY . "/";

$executable = $antsDirectory . "ANTS";

$prefix1 = "ants_";
@args = ( $executable, 3, 
  "--image-metric", "PR[".$fixedImage.",".$movingImage.",1,4]", 
  "--transformation-model", "SyN[0.5]", 
  "--regularization", "Gauss[3,1]", 
  "--number-of-iterations", "200x100x40", 
  "--output-naming", $outputDirectory . $prefix1 . ".nii.gz", 
  "--verbose" );
system( @args ) == 0 or die "system @args failed: $?";

opendir( DIR, $outputDirectory );
@files = grep( /\.nii$/, readdir(DIR) );
closedir( DIR );

foreach $file( @files )
  {
#  @args = ( "/bin/gzip", $outputDirectory . $file );
#  system( @args ) == 0 or die "system @args failed: $?";
  unlink( $outputDirectory . $file );
  }  

$executable = "/home/tustison/ANTS/bin64/WarpImageMultiTransform";
@args = ( $executable, 3, 
  $PREDIRECTORY . "/" . $PRESERIES . ".nii.gz", 
  $POSTDIRECTORY . "/ants_Warped.nii.gz", "-R", 
  $POSTDIRECTORY . "/" . $POSTSERIES . ".nii.gz",
  $POSTDIRECTORY . "/" . "ants_Warp.nii.gz", 
  $POSTDIRECTORY . "/" . "ants_Affine.txt" ); 
system( @args ) == 0 or die "system @args failed: $?";

$executable = "/home/tustison/Utilities/bin64/ConvertScalarImageToRGB";
@args = ( $executable, 3, 
  $POSTDIRECTORY . "/ants_Warped.nii.gz", 
  $POSTDIRECTORY . "/ants_Warped_hot.mha",
  "hot" ); 
system( @args ) == 0 or die "system @args failed: $?";

$executable = "/home/tustison/Utilities/bin64/ConvertScalarImageToRGB";
@args = ( $executable, 3, $POSTDIRECTORY . "/" . $POSTSERIES . ".nii.gz", 
  $POSTDIRECTORY . "/" . $POSTSERIES . "_hot.mha",
  "hot" ); 
system( @args ) == 0 or die "system @args failed: $?";

$executable = "/home/tustison/Utilities/bin64/ConvertScalarImageToRGB";
@args = ( $executable, 3, $PREDIRECTORY . "/" . $PRESERIES . ".nii.gz", 
  $PREDIRECTORY . "/" . $PRESERIES . "_hot.mha",
  "hot" ); 
system( @args ) == 0 or die "system @args failed: $?";


