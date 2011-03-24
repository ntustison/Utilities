#!/usr/bin/perl

#  $baseDirectory = "/Users/nick/pkg/Projects/HeliumLungStudies/NormalImages/PCALevelSets/PCA_2/";
#  $directory[0] = "${baseDirectory}/Asthma/";
#  $directory[1] = "${baseDirectory}/Normal/";
#  $directory[2] = "${baseDirectory}/Repeatability/";
#
#  $outputDirectory = "${baseDirectory}/AxialImages/";
#
#  $count = 0;
#  for( $i = 0; $i < 3; $i++ )
#    {
#    opendir( DIR, $directory[$i] );
#    while( defined( $file = readdir( DIR ) ) )
#      {
#      if( $file =~ m/\.nii\.gz$/ )
#        {
#        if( $count < 10 )
#          {
#          $index = "00${count}";
#          }
#        elsif( $count < 100 )
#          {
#          $index = "0${count}";
#          }
#        else
#          {
#          $index = "${count}";
#          }
#        $outputFile = "${outputDirectory}/he3_${index}.nii.gz";
#         `/bin/cp ${directory[$i]}/${file} $outputFile`;
#  #      print $file . "->" . $outputFile . "\n";
#        $count++;
#        }
#      }
#    closedir( DIR );
#    }
#


$baseDirectory = "/home/tustison/Data/HeliumLungStudies/Templates/PCAShapeTemplate/PCAModes/";
print $baseDirectory . "\n";

opendir( DIR, $baseDirectory );
while( defined( $file = readdir( DIR ) ) )
		{
  print $file . "\n";
		`/home/tustison/Utilities/bin64/ChangeImageInformation 3 ${baseDirectory}/${file} ${baseDirectory}/${file} 0 0x0x0`;
		}
closedir( DIR );
