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
#     `/home/tustison/Utilities/bin64/CreateZeroImage 3 ${directories}/${filename} ${directories}/ones_tmp.nii.gz 1`;
#
#     @fileContent = `/home/tustison/Utilities/bin64/CalculateFirstOrderStatisticsFromImage 3 ${directories}/${filename} ${directories}/ones_tmp.nii.gz 1 200 ${directories}/hist_tmp.txt`;
#     @statistics = split( / /, $fileContent[1] );
#     print $statistics[9] . "\n";

#     if( $statistics[9] > 200 )
#       {
#       print "${directories}/${filename}\n";
#       @args = ( "qsub", $directory . "/registerToTemplate.sh",
#         $directory, $directories, $directories . $filename );
#       system( @args ) == 0 or die "system @args failed: $?";
#       @args = ( "qsub", '-p', '-512', $directory . "/heliumLungAnalysis.sh",
#         $directory, $directories, $directories . $filename );
#       system( @args ) == 0 or die "system @args failed: $?";
      @args = ( "qsub", '-p', '-512', $directory . "/ventilationDefectSegmentation.sh",
        $directory, $directories, $filename );
#       @args = ( 'sh', $directory . "/ventilationDefectSegmentation.sh",
#         $directory, $directories, $filename );
      system( @args ) == 0 or die "system @args failed: $?";
#        exit( 0 );
#        print "@args\n";
#       }
#     unlink( "${directories}/ones_tmp.nii.gz" );
#     unlink( "${directories}/hist_tmp.txt" );
#      system( "/bin/rm ${directories}/*_run0*" );
#      system( "/bin/rm ${directories}/*latex*" );
#      system( "/bin/rm ${directories}/*ants*" );
#      system( "/bin/rm ${directories}/*registration*" );
    }
  }
