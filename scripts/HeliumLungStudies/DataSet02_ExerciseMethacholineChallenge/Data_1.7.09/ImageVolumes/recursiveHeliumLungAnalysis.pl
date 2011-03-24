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
    && $filename =~ m/\.nii\.gz$/ )
    { 

#    print $filename . "\n";
#    print $directories . "\n";
#     print $directory . "\n";
    
    eval
      {
      @args = ( "qsub", "-q", "x86,mac", $directory . "/heliumLungAnalysis.sh", 
        $directories, $directories . $filename );
#      @args = ( "sh", $directory . "/heliumLungAnalysis.sh", 
#        $directories, $directories . $filename );
      system( @args ) == 0 or die "system @args failed: $?"
#      print "@args\n";
      }
    }
  

#  if ( -d )
#    {
#    ( $filename, $directories, $suffix ) = fileparse( $File::Find::name );
#
#    @base_vars = split( /\//, $directory );
#    
#    @current_vars = split( /\//, $directories ); 
#
#    $newdirectory = "";
#    for $i ( 0..( @base_vars-2 ) )
#      {
#      $newdirectory = $newdirectory . $base_vars[$i] . "/";
#      }  
#    $newdirectory = $newdirectory . "ImageVolumes/";
#
#    for $i ( @base_vars.. ( @current_vars-1 ) )
#      {
#      $newdirectory = $newdirectory . $current_vars[$i] . "/";
#     }
#    $newdirectory = $newdirectory . $filename . "/";
#
#    print $newdirectory . "\n";
#    
#    print $File::Find::name . "\n";
#
#    if ( !( -e $newdirectory ) )
#      {
#      eval
#        { 
#        @args = ( "mkdir", $newdirectory );
#        system( @args ) == 0 or die "system @args failed: $?";
#        }
#      }
#
#    eval
#      {
#      print $File::Find::name . "  " .  "----" . $suffix . "\n";
#      
#      $algorithm = "/mnt/data1/tustison/Utilities/bin/ConvertDICOMImageSeries";
#      @args = ( $algorithm, $File::Find::name, 0, $newdirectory . $filename . ".nii.gz" );
#      system( @args ) == 0 or die "system @args failed: $?"
#      }
#    }
  }
