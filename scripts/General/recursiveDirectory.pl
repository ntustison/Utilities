#!/user/bin/perl -w

use File::Find;
use File::Basename;

$directory = $ARGV[0];
find( \&wanted, $directory );

sub wanted
  {
  if ( -d )
    {
    ( $filename, $directories, $suffix ) = fileparse( $File::Find::name );

    @vars = split( /\//, $directories );

    $newdirectory = "/Users/nick/Desktop/Images/";
    for $i ( 4.. ( @vars-1 ) )
      {
      $newdirectory = $newdirectory . $vars[$i] . "/";
      }
    $newdirectory = $newdirectory . $filename . "/";

    print $newdirectory . "\n";

 if ( !( -e $newdirectory ) )
 {
 @args = ( "mkdir", $newdirectory );
 system( @args ) == 0 or die "system @args failed: $?";
 }

    eval
      {
$algorithm = "/Users/nick/pkg/Utilities/bin/3D/float/ConvertDICOMImageSeries";
@args = ( $algorithm, $File::Find::name, 0, $newdirectory . $filename . ".nii.gz" );
system( @args ) == 0 or die "system @args failed: $?"
      }
    }
  }
