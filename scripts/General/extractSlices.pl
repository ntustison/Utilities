#!/user/bin/perl -w

$image = $ARGV[0];
$outputPrefix = $ARGV[1];
$outputExtension = $ARGV[2];
$direction = $ARGV[3];
$start = $ARGV[4];
$end = $ARGV[5];

for( $i = $start; $i <= $end; $i++ )
  {
  if( $i < 10 )
    {
    $outputImagePrefix = "${outputPrefix}00${i}";
    }
  elsif( $i < 100 )
    {
    $outputImagePrefix = "${outputPrefix}0${i}";
    }
  elsif( $i < 1000 )
    {
    $outputImagePrefix = "${outputPrefix}${i}";
    }

  print "Writing ${outputImagePrefix}\n";

  `/Users/nick/pkg/Utilities/bin/ExtractSliceFromImage $image ${outputImagePrefix}tmp.nii.gz $direction $i`;


  if( $outputExtension =~ m/jpg/ || $outputExtension =~ m/tiff/ )
    {
    `ConvertImage 2 ${outputImagePrefix}tmp.nii.gz ${outputImagePrefix}.${outputExtension} 1`;
    }
  else
    {
    `ConvertImage 2 ${outputImagePrefix}tmp.nii.gz ${outputImagePrefix}.${outputExtension} 0`;
    }
  unlink( "${outputImagePrefix}tmp.nii.gz" )
  }
