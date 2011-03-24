#!/user/bin/perl -w

use File::Find;
use File::Basename;

$directory = $ARGV[0];
find( \&wanted, $directory );

sub wanted
  {
  ( $filename, $directories, $suffix ) = fileparse( $File::Find::name );
  @current_vars = split( /\//, $directories ); 

  if( -d && $filename =~ m/^Series/ && $current_vars[@current_vars-1] =~ m/^Post$/ )
    {
    print $directories  . "\n";

    @args = ( "mv", $directories . $filename . "psr_ControlPointLattice__0_3xvec.nii.gz", $directories . $filename . "/psr__0_3xvec.nii.gz" ); 
    system( @args ) == 0;
    @args = ( "mv", $directories . $filename . "psr_ControlPointLattice__0_3yvec.nii.gz", $directories . $filename . "/psr__0_3yvec.nii.gz" ); 
    system( @args ) == 0;
    @args = ( "mv", $directories . $filename . "psr_ControlPointLattice__0_3zvec.nii.gz", $directories . $filename . "/psr__0_3zvec.nii.gz" ); 
    system( @args ) == 0;
    @args = ( "mv", $directories . $filename . "psr_ControlPointLattice__1_3xvec.nii.gz", $directories . $filename . "/psr__1_3xvec.nii.gz" ); 
    system( @args ) == 0;
    @args = ( "mv", $directories . $filename . "psr_ControlPointLattice__1_3yvec.nii.gz", $directories . $filename . "/psr__1_3yvec.nii.gz" ); 
    system( @args ) == 0;
    @args = ( "mv", $directories . $filename . "psr_ControlPointLattice__1_3zvec.nii.gz", $directories . $filename . "/psr__1_3zvec.nii.gz" ); 
    system( @args ) == 0;
    @args = ( "mv", $directories . $filename . "psr_ControlPointLattice_3xvec.nii.gz", $directories . $filename . "/psr_3xvec.nii.gz" ); 
    system( @args ) == 0;
    @args = ( "mv", $directories . $filename . "psr_ControlPointLattice_3yvec.nii.gz", $directories . $filename . "/psr_3yvec.nii.gz" ); 
    system( @args ) == 0;
    @args = ( "mv", $directories . $filename . "psr_ControlPointLattice_3zvec.nii.gz", $directories . $filename . "/psr_3zvec.nii.gz" ); 
    system( @args ) == 0;
    @args = ( "rm", $directories . $filename . "psr_WarpedPoints_3.vtk" );
    system( @args ) == 0;
    }

     
#    opendir( DIR, $directories . $filename );
#    @files = grep( /^Series/, readdir( DIR ) );
#    closedir( DIR );
#    
#    print @files . "\n";
#
#    foreach $file ( @files )
#      {
#      if( $file =~ m/^Series/ && ( $file =~ m/\.nii\.gz$/ || $file =~ m/\.vtk$/ ) )
#        {
#        print "   " . $file . "\n";
#        }
#      }  
    
#      {
#      foreach $preserie ( @preseries )
#        {
#        eval
#          {
#          $preDirectory = $directories . "Pre/" . $preserie; 
#          $postDirectory = $directories . "Post/" . $postserie; 
#           
#          @args = ( "qsub", $directory . "/pointSetRegistration.sh", 
#            $preDirectory, $postDirectory );
#          system( @args ) == 0 or die "system @args failed: $?"
#          }
#        }
#      }
  }
