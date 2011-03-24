#!/user/bin/perl -w

use File::Find;
use File::Basename;

$directory = $ARGV[0];
find( \&wanted, $directory );

sub wanted
  {

  $utilityDirectory = "/home/tustison/Utilities/bin64/";


  ( $filename, $directories, $suffix ) = fileparse( $File::Find::name );

  if( $directories =~ m/HeliumImages/ && $filename =~ m/\.nii\.gz$/
      && $filename !~ m/inspiration/ && $filename !~ m/expiration/ && $filename !~ m/run/ )
    {

    sleep 1;

    if( $directories =~ m/Site1/ )
      {
      ( $base, $path, $suffix ) = fileparse( $directories . $filename, qr/\.nii\.gz/ );
      $imagePrefix = "${path}/${base}";

#       `/usr/bin/perl /home/tustison/Data/HeliumLungStudies/DataSet08_TwoImagingSites/ImageVolumes/splitSite1Images.pl ${directories}/${filename}`;
#
#       system( "${utilityDirectory}/PermuteAxesImage 3 ${imagePrefix}_inspiration.nii.gz ${imagePrefix}_inspiration.nii.gz 0x2x1" );
#       system( "${utilityDirectory}/FlipImage 3 ${imagePrefix}_inspiration.nii.gz ${imagePrefix}_inspiration.nii.gz 0x0x1" );
#       system( "${utilityDirectory}/ChangeImageInformation 3 ${imagePrefix}_inspiration.nii.gz ${imagePrefix}_inspiration.nii.gz 2 0x0x0" );

      @args = ( "qsub", '-p', '-513', $directory . "/ventilationAnalysis.sh",
        $directory, $directories, $imagePrefix . '_inspiration.nii.gz' );
      print "@args\n";
      system( @args ) == 0 or die "system @args failed: $?";
      }
    else
      {
#       `${utilityDirectory}/PermuteAxesImage 3 $filename $filename 0x2x1`;
#       `${utilityDirectory}/FlipImage 3 $filename $filename 0x0x1`;
#       `${utilityDirectory}/ChangeImageInformation 3 $filename $filename 2 0x0x0`;
#
      @args = ( "qsub", '-p', '-513', $directory . "/ventilationAnalysis.sh",
        $directory, $directories, $directories . $filename );
      system( @args ) == 0 or die "system @args failed: $?"
      }
    }
  }
