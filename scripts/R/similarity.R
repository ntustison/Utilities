#@alpha = ( 1.0, 1.5, 2.0 );
#@pointSetSigma = ( 5 );
#@noiseSigma = ( 50, 100 );
#@noisePercentage = ( 200, 300 );

alpha_idx = '*';
pssig_idx = 0;
nssig_idx = 1;
nsper_idx = 1;

patternJHCT <- paste( "translationJHCT_", alpha_idx, "_", pssig_idx, "_",
  nssig_idx, "_", nsper_idx, ".csv", sep = "", collapse = ' ' );
patternICP <- paste( "translationICP_", 0, "_", pssig_idx, "_",
  nssig_idx, "_", nsper_idx, ".csv", sep = "", collapse = ' ' );

jhctFiles <- list.files( path = ".", pattern = glob2rx( patternJHCT ) );

icp <- read.csv( patternICP, header = FALSE );

plot( icp[,1], ( icp[,2] - min( icp[,2] ) ) / ( max( icp[,2] ) - min( icp[,2] ) ), type = "l" );
for( i in seq( 1, length( jhctFiles ) ) )
  {
  jhct <- read.csv( jhctFiles[i] );
  lines( jhct[,1], ( jhct[,2] - min( jhct[,2] ) ) / ( max( jhct[,2] ) - min( jhct[,2] ) ) );
  }

