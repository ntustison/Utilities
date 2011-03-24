#!/usr/bin/perl

$baseDirectory = "/Users/nick/pkg/Projects/HeliumLungStudies/NormalImages/PCALevelSets/PCA_2/";
$directory[0] = "${baseDirectory}/Asthma/";
$directory[1] = "${baseDirectory}/Normal/";
$directory[2] = "${baseDirectory}/Repeatability/";

$outputDirectory = "${baseDirectory}/AxialImages/";

$count = 0;
for( $i = 0; $i < 3; $i++ )
  {
  opendir( DIR, $director[$i] );
  while( defined( $file = readdir( DIR ) ) )
    {
    if( $file =~ m/\.nii\.gz$/ )
      {
      if( $count < 10 )
        {
        $index = "00${count}";
        }
      elsif( $count < 100 )
        {
        $index = "0${count}";
        }
      else
        {
        $index = "${count}";
        }
      $outputFile = "${outputDirectory}/he3_${index}.nii.gz";
      `/bin/cp $file $outputFile`;
      $count++;
      }
    }
  closedir( DIR );
  }
