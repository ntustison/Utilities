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

    my $class = "";
    if( $directories =~ m/Asthma/ )
      {
      $class = "asthma";
      }
    else
      {
      $class = "normal";
      }

    my $time = "";
    if( $directories =~ m/Post/ )
      {
      $time = "post";
      }
    else
      {
      $time = "pre";
      }

    @args = ( "/bin/cp", "${directories}/${filename}/he3_bias_corrected_run03.nii.gz", "${outputDirectory}/he3_${class}_${time}_${index}.nii.gz" );
    system( @args ) == 0;
    @args = ( "/bin/cp", "${directories}/${filename}/apocrita_both_run03.nii.gz", "${outputDirectory}/apocrita_${class}_${time}_${index}.nii.gz" );
    system( @args ) == 0;

    @args = ( "/bin/cp", "${directories}/${filename}/${filename}.nii.gz", "${outputDirectory}/he3_uncorrected_${index}.nii.gz" );
#    system( @args ) == 0;
    @args = ( "/bin/cp", "${directories}/${filename}/bias_field_run03.nii.gz", "${outputDirectory}/bias_field_${index}.nii.gz" );
#    system( @args ) == 0;

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
