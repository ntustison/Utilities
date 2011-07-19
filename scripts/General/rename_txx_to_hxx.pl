#!/user/bin/perl -w

use File::Find;
use File::Basename;

$directory = $ARGV[0];
find( \&wanted, $directory );

sub wanted
  {
  if ( -d )
    {
    my ( $filename, $directories, $suffix ) = fileparse( $File::Find::name );

    my $localDirectory = "${directories}/${filename}";
    opendir( DIR, $localDirectory ) || die "can't opendir: $!";

    while( my $file = readdir( DIR ) )
      {
      next if( $file !~ /\.txx$/ && $file !~ /\.cxx$/ && $file !~ /\.h$/ &&
        $file !~ /^CMakeLists.txt$/ && $file !~ /\.hxx$/ );

#       print $file . "\n";

      my $filePath = "${localDirectory}/${file}";
      open( DAT, $filePath ) || die( "Could not open file!" );
      my @contents=<DAT>;
      close( DAT );

      open( FILE, ">${filePath}" );

      for( my $i = 0; $i < @contents; $i++ )
        {
        $contents[$i] =~ s/txx/hxx/g;
#         $contents[$i] =~ s/^\s//;
        print FILE "${contents[$i]}";
        }
      close( FILE );

      if( $filePath =~ /\.txx$/ )
        {
        my $newFilePath = $filePath;
        $newFilePath =~ s/\.txx/\.hxx/;
        `mv $filePath $newFilePath`;
        }
      }
    closedir( DIR );
    }
  }
