#!/user/bin/perl -w

$directory = $ARGV[0];

print "Converting point set files\n";

$convertPoints = "/Users/nick/pkg/Utilities/bin/3D/float/ConvertAvantsLandmarkFileToVTK";

opendir( DIR, $directory );
while( defined( $file = readdir( DIR ) ) )
  {
  next unless $file =~ /Points/ and $file =~ /txt/;

		$basefilename = $file;
		$basefilename =~ s!^(?:.*/)?(.+?)(?:\.[^.]*)?$!$1!;

		@args = ( $convertPoints, $directory . "/" . $file, $directory . "/" . $file, 
            $directory . "/" . $basefilename, 2 );
		system( @args ) == 0 or die "system @args failed: $?";
  print "  converting " . $basefilename . "\n";  

		@args = ( "rm", $directory . "/" . $file );
		system( @args ) == 0 or die "system @args failed: $?";
  }
closedir( DIR );

print "Converting gradients\n";

opendir( DIR, $directory );
while( defined( $file = readdir( DIR ) ) )
  {
  next unless $file =~ /Gradient/ and $file =~ /Moving/;

		$basefilename = $file;
		$basefilename =~ s!^(?:.*/)?(.+?)(?:\.[^.]*)?$!$1!;

  $basefilename =~ s/Moving//g;

  print "  converting " . $basefilename . "\n";  
		@args = ( $convertPoints, 
            $directory . "/" . $basefilename . "Fixed.txt", 
            $directory . "/" . $basefilename . "Moving.txt", 
            $directory . "/" . $basefilename, 2 );
		system( @args ) == 0 or die "system @args failed: $?";

		@args = ( "rm", $directory . "/" . $basefilename . "Fixed.txt" );
		system( @args ) == 0 or die "system @args failed: $?";
		@args = ( "rm", $directory . "/" . $basefilename . "Moving.txt" );
		system( @args ) == 0 or die "system @args failed: $?";
		@args = ( "rm", $directory . "/" . $basefilename . "Fixed.vtk" );
		system( @args ) == 0 or die "system @args failed: $?";
		@args = ( "mv", $directory . "/" . $basefilename . "Moving.vtk",  
                  $directory . "/" . $basefilename . ".vtk" );
		system( @args ) == 0 or die "system @args failed: $?";
  }
closedir( DIR );

print "Converting probability images\n";
$convertImage = "/Users/nick/pkg/Utilities/bin/3D/float/ConvertImage";

opendir( DIR, $directory );
while( defined( $file = readdir( DIR ) ) )
  {
  next unless $file =~ /ProbabilityImage/;

		$basefilename = $file;
		$basefilename =~ s!^(?:.*/)?(.+?)(?:\.[^.]*)?$!$1!;

		@args = ( $convertImage, $directory . "/" . $file,  
            $directory . "/" . $basefilename . ".mha" );
		system( @args ) == 0 or die "system @args failed: $?";
  print "  converting " . $basefilename . "\n";  

		@args = ( "rm", $directory . "/" . $file );
		system( @args ) == 0 or die "system @args failed: $?";
  }
closedir( DIR );

print "Converting deformed grids\n";
$convertLines = "/Users/nick/pkg/Utilities/bin/3D/float/ConvertAvantsLandmarkFileToVTKLines";

opendir( DIR, $directory );
while( defined( $file = readdir( DIR ) ) )
  {
  next unless $file =~ /Grid/ and $file =~ /txt/;

		$basefilename = $file;
		$basefilename =~ s!^(?:.*/)?(.+?)(?:\.[^.]*)?$!$1!;

		@args = ( $convertLines, $directory . "/" . $file, 
    $directory . "/" . $basefilename . ".vtk" );
		system( @args ) == 0 or die "system @args failed: $?";
  print "  converting " . $basefilename . "\n";  

		@args = ( "rm", $directory . "/" . $file );
		system( @args ) == 0 or die "system @args failed: $?";
  }
closedir( DIR );


