#!/user/bin/perl -w

use File::Find;
use File::Basename;

$directory = $ARGV[0];

find( \&wanted, $directory );


sub wanted
  {

  ( $filename, $directories, $suffix ) = fileparse( $File::Find::name );
  @current_vars = split( /\//, $directories ); 

  if( -d && $filename =~ m/^Series/ )
    { 
    eval
      {
      print $directories . $filename . "/" . "\n"; 
      unlink glob( $directories . $filename . '/*run0*' );
     
#      $image = $directories . $filename . "/" . $filename . ".nii.gz"; 
#      $segmentation = $directories . $filename . "/segmentation.nii.gz"; 
#        
#      @args = ( "snap", "-g", $image, "-s", $segmentation );
#      system( @args ) == 0 or die "system @args failed: $?";
#      print "@args\n";
      }
    }
  }

