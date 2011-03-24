#!/user/bin/perl -w

use File::Find;
use File::Basename;

$directory = $ARGV[0];
find( \&wanted, $directory );

sub wanted
  {

  ( $filename, $directories, $suffix ) = fileparse( $File::Find::name );
  @current_vars = split( /\//, $directories ); 

  if( -d && $filename =~ m/^Post$/ )
    { 
    print $directories . "\n"; 
    opendir( preDIR, $directories . "Pre/" );
    @preseries = grep( /^Series/, readdir( preDIR ) );
    closedir( preDIR );

    opendir( postDIR, $directories . "Post/" );
    @postseries = grep( /^Series/, readdir( postDIR ) );
    closedir( postDIR );
        
    foreach $postserie ( @postseries )
      {
      foreach $preserie ( @preseries )
        {
        eval
          {
          $preDirectory = $directories . "Pre/" . $preserie; 
          $postDirectory = $directories . "Post/" . $postserie; 
           
          @args = ( "qsub", "-q", "x86", $directory . "/warpImages.sh", 
            $preDirectory, $postDirectory, $preserie, $postserie );
          system( @args ) == 0 or die "system @args failed: $?"
          }
        }
      }
    }
  }
