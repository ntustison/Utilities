#!/user/bin/perl -w

use File::Find;
use File::Basename;

$directory = $ARGV[0];

find( \&wanted, $directory );

sub wanted
  {

  ( $filename, $directories, $suffix ) = fileparse( $File::Find::name );
  @current_vars = split( /\//, $directories );

  if( $directories =~ m/HeliumImages/ &&
    ( $filename =~ m/ndif/ || $filename =~ m/^1.2.840/ )
    && $filename !~ m/.nii.gz/ )
#  if( $filename =~ m/^1.2.840/ && $filename =~ m/nii.gz$/ )
    {
    eval
      {
      $utilitiesDirectory = "/home/tustison/Utilities/bin64/";
      print $directories . $filename . "/" . "\n";
      unlink glob( $directories . $filename . '/*run0*' );
      unlink glob( $directories . $filename . '/registrationToTemplate*' );
#      $changeImageInformation = $utilitiesDirectory . "ChangeImageInformation";
#
#      print "Processing $filename\n";
#
#      @args = ( $changeImageInformation, 3, $filename, $filename, 0, "0x0x0" );
#      system( @args ) == 0 or die "system @args failed: $?";
#      @args = ( $changeImageInformation, 3, $filename, $filename, 1, "1.7188x1.7188x7.5" );
#      system( @args ) == 0 or die "system @args failed: $?";

#      $image = $directories . $filename . "/" . $filename . ".nii.gz";
#      $segmentation = $directories . $filename . "/segmentation.nii.gz";
#
#      @args = ( "snap", "-g", $image, "-s", $segmentation );
#      system( @args ) == 0 or die "system @args failed: $?";
#      print "@args\n";
      }
    }
  }

