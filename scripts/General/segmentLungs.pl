#!/usr/bin/perl

$imageDirectory = "/Users/nick/Desktop/DrewCT/";


$prefix[0] = "ExpB41f";
$prefix[1] = "InspB41f";
$prefix[2] = "ExpB50f";
$prefix[3] = "InspB50f";

$executableDirectory = "/Users/nick/pkg/Utilities/bin/3D/float/";

$extractLungs = $executableDirectory . "ExtractLungs";
$segmentAirways = $executableDirectory . "SegmentAirways";
$separateLungs = $executableDirectory . "SeparateLungs";
$smoothLungs = $executableDirectory . "SmoothLungs";


for $i ( 0 .. @prefix-1 )
  {
  $inputImage = $imageDirectory . $prefix[$i] . "_downsampled.nii.gz";
  $tmpImage = $imageDirectory . $prefix[$i] . "_segmentation.nii.gz";

  print "Processing image " . $inputImage . "\n";
  print "  Extracting lungs ...";
  @args = ( $extractLungs, $inputImage, $tmpImage );
  system( @args ) == 0 or die "system @args failed: $?";
  print " Done. \n";

  print "  Segmenting airways ...";
  @args = ( $segmentAirways, $inputImage, $tmpImage, $tmpImage );
  system( @args ) == 0 or die "system @args failed: $?";
  print " Done. \n";

  print "  Separating lungs ...";
  @args = ( $separateLungs, $inputImage, $tmpImage, $tmpImage );
  system( @args ) == 0 or die "system @args failed: $?";
  print " Done. \n";

  print "  Smoothing lungs ...";
  @args = ( $smoothLungs, $tmpImage, $tmpImage );
  system( @args ) == 0 or die "system @args failed: $?";
  print " Done. \n";
  }
