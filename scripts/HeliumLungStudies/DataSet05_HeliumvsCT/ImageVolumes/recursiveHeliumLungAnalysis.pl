#!/user/bin/perl -w

use File::Find;
use File::Basename;

$directory = $ARGV[0];
find( \&wanted, $directory );

sub wanted
  {

  ( $filename, $directories, $suffix ) = fileparse( $File::Find::name );
  @current_vars = split( /\//, $directories );

  if( -f && $filename =~ m/^$current_vars[@current_vars-1]/
    && $filename =~ m/\.nii\.gz$/ && $filename =~ m/^Series/ )
    {

#    print $filename . "\n";
#    print $directories . "\n";
#    print $directory . "\n";

    eval
      {
      @args = ( "qsub", $directory . "/heliumLungAnalysis.sh",
        $directory, $directories, $directories . $filename );
      system( @args ) == 0 or die "system @args failed: $?"
#      print "@args\n";
      }
    }
  }
