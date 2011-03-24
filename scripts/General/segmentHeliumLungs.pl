#!/usr/bin/perl

$imageDirectory[0] = "/mnt/data2/PUBLIC/Data/Input/HeliumLungStudies/RegistrationWithCTData/ImageVolumes/Pre Post Exercise/Case 1/Post Ex Ventilation/Series20030609_101358_218000/";

$images[0] = "Series20030609_101358_218000";

$suffix = ".nii.gz";

$utilityDirectory = "/mnt/data1/tustison/Utilities/bin/3D/float/";
$outputDirectory = "/mnt/data1/tustison/Projects/HeliumLungStudies/Registration/segmentation/";

$segmentation = "/mnt/aibs1/suyash/bin/itkSegmentationNMRF"; 

for $i ( 0 .. @images-1 )
  {
  $inputImage = $imageDirectory[$i] . $images[$i] . $suffix;
  $thresholdImage = $outputDirectory . $images[$i] . "_threshold.nii.gz"; 
  $otsuImage = $outputDirectory . $images[$i] . "_otsu.nii.gz"; 
  
  $algorithm = $utilityDirectory . "ThresholdImage";
  @args = ( $algorithm, $inputImage, $thresholdImage, 0, 25, 0, 1 );
  system( @args ) == 0 or die "system @args failed: $?";

  $algorithm = $utilityDirectory . "ShapeMorphology";
  @args = ( $algorithm, $thresholdImage, $thresholdImage, 
    1, 0, 1, 0, 2, 0, 0, 3, 0, 4, 1, 1, 1, 1, 0 );
  system( @args ) == 0 or die "system @args failed: $?";
   
  $algorithm = $utilityDirectory . "OtsuThresholdImage";
  @args = ( $algorithm, $inputImage, $otsuImage, 255, 2 );
  system( @args ) == 0 or die "system @args failed: $?";
  
  $algorithm = $utilityDirectory . "AddImage"; 
  @args = ( $algorithm, $thresholdImage, $otsuImage, $otsuImage );
  system( @args ) == 0 or die "system @args failed: $?";

  $priorityImage1 = $outputDirectory . $images[$i] . "_priority1.nii.gz"; 
  $algorithm = $utilityDirectory . "ThresholdImage";
  @args = ( $algorithm, $otsuImage, $priorityImage1, 1, 2, 1, 0);
  system( @args ) == 0 or die "system @args failed: $?";

  $priorityImage2 = $outputDirectory . $images[$i] . "_priority2.nii.gz"; 
  $algorithm = $utilityDirectory . "ThresholdImage";
  @args = ( $algorithm, $otsuImage, $priorityImage2, 3, 3, 1, 0);
  system( @args ) == 0 or die "system @args failed: $?";

  
  @args = ( $segmentation, "-d", $inputImage, 
                           "-p", $priorityImage1,
                           "-p", $priorityImage2,
                           "-o", $outputDirectory . $images[$i] . "_seg",
                           "-m", $outputDirectory . $images[$i] . "_mem"
                           );
  system( @args ) == 0 or die "system @args failed: $?";
  }
