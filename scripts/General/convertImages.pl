#!/user/bin/perl -w

$directory = $ARGV[0];

opendir( DIR, $directory );
while( defined( $file = readdir( DIR ) ) )
  {
  next unless $file =~ /tif/;

		$basefilename = $file;
		$basefilename =~ s!^(?:.*/)?(.+?)(?:\.[^.]*)?$!$1!;

  print "  converting " . $basefilename . "\n";
		$out = `convert -colorspace RGB -units PixelsPerInch ${directory}/${file} -verbose -resample 600 -depth 8 -background white -flatten +matte ${directory}/${file}`;

# 		@args = ( "rm", $directory . "/" . $file );
# 		system( @args ) == 0 or die "system @args failed: $?";
  }
closedir( DIR );

