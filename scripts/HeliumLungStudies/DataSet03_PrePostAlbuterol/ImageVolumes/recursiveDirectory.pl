#!/user/bin/perl -w

use File::Find;
use File::Basename;

$directory = $ARGV[0];
find( \&wanted, $directory );

$index = 0;


sub wanted
  {
  ( $filename, $directories, $suffix ) = fileparse( $File::Find::name );
  @current_vars = split( /\//, $directories );

  if( -d && $filename =~ m/^Series/ && $filename !~ m/\.nii\.gz$/ )
    {
    print $filename . "\n";

    $outputDirectory = "/home/tustison/Data/HeliumLungStudies/ExerciseMethacholineChallenge/Data_2.12.09/ImageVolumes/TestImages/";

    @args = ( "/bin/cp", "${directories}/${filename}/he3_bias_corrected_run03.nii.gz", "${outputDirectory}/he3_${index}.nii.gz" );
    system( @args ) == 0;
    @args = ( "/bin/cp", "${directories}/${filename}/initial_segmentation_run03.nii.gz", "${outputDirectory}/initial_segmentation_${index}.nii.gz" );
    system( @args ) == 0;

    $index++;
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
