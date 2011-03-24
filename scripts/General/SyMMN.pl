#!/usr/bin/perl

$utilityDirectory = "/mnt/data1/tustison/Utilities/bin/3D/float/";

$imageDirectory = "/mnt/data1/tustison/Projects/Spiromics/DrewCT/downsampled/";
$fixedImagePrefix = "InspB41f";
$movingImagePrefix = "ExpB41f";

$outputDirectory = "/mnt/data1/tustison/Projects/Spiromics/DrewCT/Results/";


$algorithm = "/mnt/data2/Avants/bin/3D/mstart3dt4_mask";

@args = ( $algorithm, $imageDirectory . $fixedImagePrefix . "_downsampled.nii.gz",
                      $imageDirectory . $movingImagePrefix . "_downsampled.nii.gz",
                      $imageDirectory . $movingImagePrefix . "to" . $fixedImagePrefix . "_affine.nii.gz",
                      $imageDirectory . $movingImagePrefix . "to" . $fixedImagePrefix . "_affine.txt",
                      0,
                      1000,
                      16,
                      9000
                      );
#system( @args ) == 0 or die "system @args failed: $?";


$algorithm = "/mnt/data2/Avants/bin/3D/SyMMN";

$naming = "SyMMN" . $movingImagePrefix . "to" . $fixedImagePrefix;

@args = ( $algorithm, "-f", $imageDirectory . $fixedImagePrefix . "_downsampledx8.mha",
                      "-m", $imageDirectory . $movingImagePrefix . "_downsampledx8.mha",
                      "-c", 5,
                      "-n", 1,
                      "-i", "1",
                      "-o", $outputDirectory . $naming
                      );
system( @args ) == 0 or die "system @args failed: $?";
