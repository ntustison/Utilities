#!/user/bin/perl -w

use File::Find;
use File::Basename;

$directory = $ARGV[0];
find( \&wanted, $directory );

sub wanted
  {

  ( $filename, $directories, $suffix ) = fileparse( $File::Find::name );
  @current_vars = split( /\//, $directories );

  if( $directories =~ m/Processed/ && $filename =~ m/he3.nii.gz/ )
    {

    @args = ( "qsub", '-p', '-512', $directory . "/ventilationDefectSegmentation.sh",
              $directory, $directories  );
    system( @args ) == 0 or die "system @args failed: $?";
    print "@args\n";
    }
  }
