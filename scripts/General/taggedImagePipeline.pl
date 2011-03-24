#!/usr/bin/perl

$imageDirectory = "/Users/nick/pkg/Projects/HeliumLungStudies/He3-833/TrinityFiles/";
$outputDirectory = $imageDirectory . "PipelineResults/";

$images[0] = "He3_833_ins_resampled.nii.gz";
$images[1] = "He3_833_exp_resampled.nii.gz";

$dimension = 3;
$keepImageFiles = 1;
$thinTagLines = 1;

#####################################################################
#####################################################################

if ( $dimension == 2 )
  {
  $utilityDirectory = "/Users/nick/pkg/Utilities/bin/2D/float/";
  }
else
  {
  $utilityDirectory = "/Users/nick/pkg/Utilities/bin/3D/float/";
  }

for $i ( 0..@images-1 )
  {
  print "Processing " . $images[$i] . "\n";

		##
		# Anisotropic smoothing
		##
  print "Anisotropic smoothing\n";
		$executable = $utilityDirectory . "GradientAnisotropicDiffusionImageFilter";
  $anis_input = $imageDirectory . $images[$i];
  $anis_output = $outputDirectory . "anis_" . $i . "_" . $images[$i];		
  $numberOfIterations = 5;
  $timeStep = 0.125;
  if ( $dimension == 3 )
    {
    $timeStep = 0.0625;
    } 
  $conductance = 7;
		@args = ( $executable, $anis_input, $anis_output, $numberOfIterations, $timeStep, $conductance );
		system( @args ) == 0 or die "system @args failed: $?";

		##
		# Grayscale dilation
		##
  print "Grayscale dilation\n";
		$executable = $utilityDirectory . "GrayscaleDilateImage";
  $gdil_input = $anis_output;
  $gdil_output = $outputDirectory . "gdil_" . $i . "_" . $images[$i];		
  $radius = 1;
		@args = ( $executable, $gdil_input, $gdil_output, $radius, 0 );
		system( @args ) == 0 or die "system @args failed: $?";

		##
		# Otsu thresholding
		##
  print "Otsu thresholding\n";
		$executable = $utilityDirectory . "OtsuThresholdImage";
  $otsu_input = $gdil_output;
  $otsu_output = $outputDirectory . "otsu_" . $i . "_" . $images[$i];		
  $numberOfBins = 100;
  $numberOfThresholds = 1;
		@args = ( $executable, $otsu_input, $otsu_output, $numberOfBins, $numberOfThresholds );
		system( @args ) == 0 or die "system @args failed: $?";

		##
		# Label and close
		##
  print "Label and close\n";
		$executable = $utilityDirectory . "LabelAndCloseImage";
  $landc_input = $otsu_output;
  $landc_output[$i] = $outputDirectory . "lac_" . $i . "_" . $images[$i];		
		@args = ( $executable, $landc_input, $landc_output[$i] );
		system( @args ) == 0 or die "system @args failed: $?";

  ##
  # Inhomogeneity correction
  ##
  print "Inhomogeneity correction\n";
  $executable = $utilityDirectory . "InhomogeneityCorrectImage";
  $ic_input = $gdil_output;
  $ic_output = $outputDirectory . "ic_" . $i . "_" . $images[$i];
  $ic_mask = $otsu_output;
  @args = ( $executable, $ic_input, $ic_output, $ic_mask );
		system( @args ) == 0 or die "system @args failed: $?";

  ##
  # Extract tag planes
  ##
  print "Extract tag planes\n";
  $executable = $utilityDirectory . "ExtractTagLinePoints";
  $tag_input = $ic_output;
  $tag_mask = $landc_output[$i];
  $tag_output = $outputDirectory . "tag_" . $i . "_";
  $spacing = 28;
  $percentage = 0.85;   # use automated thresholding
  $numberOfAngleSteps = 3;
  $numberOfTagSpacingSteps = 3;
  $angleOffset = 5;
  $tagSpacingFactor = 0.3;
  if ( $i == 0 )
    {
    $numberOfAngleSteps = 3;
    $numberOfTagSpacingSteps = 3;
    $tagSpacingFactor = 0.1;
    $angleOffset = 2;
    } 
  @args = ( $executable, $tag_input, $tag_mask, $tag_output, $spacing, $percentage,
            $numberOfAngleSteps, $numberOfTagSpacingSteps, $angleOffset, $tagSpacingFactor );
		system( @args ) == 0 or die "system @args failed: $?";

  for $d ( 0..$dimension-1 )
    {
    $tagLabelImage[$d] = $tag_output . "GaborLabel." . $d . ".nii.gz";
    $tagGaborResponseImage[$d] = $tag_output . "GaborMaximalResponse." . $d . ".nii.gz";
    $tagPoints[$d] = $tag_output . "GaborPoints." . $d . ".txt";
    }  

  if ( $thinTagLines )
    {
				##
				# Thin the point-sets
				##
				print "Thin tag planes\n";
				for $d ( 0..$dimension-1 )
						{
  				$executable = $utilityDirectory . "ThinTagLines";
						@args = ( $executable, $tagLabelImage[$d], $tagLabelImage[$d], $d );
  	  	system( @args ) == 0 or die "system @args failed: $?";

  				$executable = $utilityDirectory . "ConvertSegmentedImageToAvantsLandmarkFile";
						@args = ( $executable, $tagLabelImage[$d], $tagPoints[$d] );
  	  	system( @args ) == 0 or die "system @args failed: $?";
      }   
    }  

  ##
  # Convert point-set files
  ##
  $executable = "/Users/nick/pkg/Utilities/bin/3D/float/ConvertAvantsLandmarkFileToVTK";
  for $d ( 0..$dimension-1 )
    {
  		@args = ( $executable, $tagPoints[$d], $tagPoints[$d], 
                $tag_output . "GaborPoints." . $d, 2  );
	  	system( @args ) == 0 or die "system @args failed: $?";
    }    

  ##
  # Remove temporary files
  ##
  if ( ! $keepImageFiles )
    {
    $executable = "rm";
     
				@args = ( $executable, $anis_output, $gdil_output, $otsu_output, $ic_output );
				system( @args ) == 0 or die "system @args failed: $?";

				for $d ( 0..$dimension-1 )
						{
      $executable = "rm";
						@args = ( $executable, $tagLabelImage[$d], $tagGaborResponseImage[$d], $tagPoints[$d] );
						system( @args ) == 0 or die "system @args failed: $?";
						}  
    } 
  }


##
# Register point-sets
##

$executable = $utilityDirectory . "RegisterPointSetsFFD";
$prefix = "registered";
$order = 3;

for $i ( 0..0 )
#for $i ( 1..@images-1 )
  {
  for $d ( 0..$dimension-1 )
    {
				
				$movingPointSet = $outputDirectory . "tag_" . $i . "_GaborPoints." . $d . ".vtk";
				$fixedPointSet = $outputDirectory . "tag_0_GaborPoints." . $d . ".vtk";
				
    $directionality = "0";
    if ( $d == 0 )
      {
      $directionality = "1";
      } 
    for $f ( 1..$dimension-1 )
      {
      if ( $d == $f )
        {
        $directionality = $directionality . "x1"
        }
      else
        {
        $directionality = $directionality . "x0"
        }     
      }     

				$naming = $prefix . "_" . $i . "_" . $directionality;
				
				@args = ( $executable, "-f", $fixedPointSet, 
																										"-m", $movingPointSet,
																										"-a", 2.0,                                # alpha [1.0, 2.0]      
																										"-g", 0,                                  # generate mean shape 
																										"-c", 6,                                  # regularization sigma
																										"-d", 1,                                  # initialization sigma
																										"-r", 0.9,
																										"-h", 1,                                  # use input as samples
																										"-k", 4,    
																										"-n", 3, 
																										"-i", "10x5x2", 
																										"-B", $order,
																										"-R", "3x3x3",
																										"-l", 20,  
																										"-o", $outputDirectory . $naming,
																										"-p", 0,
																										"-D", $directionality, 
																										"-I", $imageDirectory . $images[$i]
																										);
				system( @args ) == 0 or die "system @args failed: $?";
    }

		##
		# Combine control point results
		##
		
		for $d ( 0..$dimension-1 )
				{
				$directionality = "0";
				if ( $d == 0 )
						{
						$directionality = "1";
						} 
				for $f ( 1..$dimension-1 )
						{
						if ( $d == $f )
								{
								$directionality = $directionality . "x1"
								}
						else
								{
								$directionality = $directionality . "x0"
								}     
						}     
				if ( $d == 0 )
						{
						$extension = "xvec.nii";
						} 
				elsif ( $d == 1 )
						{
						$extension = "yvec.nii";
						} 
				elsif ( $d == 2 )
						{
						$extension = "zvec.nii";
						} 

				$naming = $prefix . "_" . $i . "_" . $directionality;
			
				$executable = "cp";
				$source = $outputDirectory . $naming . "ControlPointLattice_1_" . $order . $extension;
				$target = $outputDirectory . $prefix . "_" . $i . "_" . "ControlPointLattice_1_" . $order . $extension;
				@args = ( $executable, $source, $target );
				system( @args ) == 0 or die "system @args failed: $?";
		
				##
				# Create image from control points
				##
				$executable = $utilityDirectory . "GenerateImageFromControlPointLattice";
				$controlPointLattice = $target;
				$output = $outputDirectory . $prefix . "_" . $i . "_" . "DeformationField" . $extension . ".gz";  
				@args = ( $executable, $controlPointLattice, $output, $imageDirectory . $images[0], $order );
				system( @args ) == 0 or die "system @args failed: $?";
		
				##
				# Remove extraneous results
				##
				opendir( DIR, $outputDirectory );
				while( defined( $file = readdir( DIR ) ) )
						{
						next unless $file =~ /$prefix/ and $file =~ /$directionality/ and $file =~ /ControlPointLattice/;
				
						@args = ( "rm", $outputDirectory . $file );
						system( @args ) == 0 or die "system @args failed: $?";
						}
				closedir( DIR );
		
				##
				# Convert warped points to vtk and rm txt file
				##
				$executable = "/Users/nick/pkg/Utilities/bin/3D/float/ConvertAvantsLandmarkFileToVTK";
				$filePrefix = $outputDirectory . $prefix . "_" . $i . "_" . $directionality  . "WarpedPoints_1_" . $order;
				@args = ( $executable, $filePrefix . ".txt", $filePrefix . ".txt", $filePrefix, 2 );
				system( @args ) == 0 or die "system @args failed: $?";
		
				$executable = "rm";
				$filePrefix = $outputDirectory . $prefix . "_" . $i . "_" . $directionality  . "WarpedPoints_1_" . $order;
				@args = ( $executable, $filePrefix . ".txt" );
				system( @args ) == 0 or die "system @args failed: $?";
		
				$filePrefix = $outputDirectory . $prefix . "_" . $i . "_" . $directionality  . "WarpedPoints_0_" . $order;
				@args = ( $executable, $filePrefix . ".txt" );
				system( @args ) == 0 or die "system @args failed: $?";
				}
		
		##
		# Generate inverse deformation field
		##
		print "Generating inverse deformation field.\n";
		$executable = $utilityDirectory . "GenerateInverseBSplineDeformationField";
		$ifield_input = $outputDirectory . $prefix . "_" . $i . "_" . "DeformationField.nii.gz";
		$ifield_output = $outputDirectory . $prefix . "_" . $i . "_" . "InverseDeformationField.nii.gz";
		$nlevels = 6;
		$mask = "";  
		@args = ( $executable, $ifield_input, $ifield_output, $order, $nlevels );
		system( @args ) == 0 or die "system @args failed: $?";
		
		##
		# Generate principal strains from inverse deformation field
		##
		print "Creating principal strain images.\n";
		$executable = $utilityDirectory . "CreatePrincipalStrainImages";
		$input = $ifield_output;
		$outputPrefix = $outputDirectory . $prefix . "_" . $i . "_" . "PrincipalStrain";
		$mask = $landc_output[0];
		@args = ( $executable, $input, $outputPrefix, $mask );
		system( @args ) == 0 or die "system @args failed: $?";

		## 
		# Create vector fields
		##
		$executable = $utilityDirectory . "ConvertDeformationFieldToAvantsLandmarkFiles";
		$input_field = $outputDirectory . $prefix . "_" . $i . "_" . "DeformationField.nii.gz";
  $output_prefix = $outputDirectory . $prefix . "_" . $i . "_" . "DeformationField";
  $mask = $landc_output[$i];
		@args = ( $executable, $input_field, $output_prefix, 1, $mask );
		system( @args ) == 0 or die "system @args failed: $?";

		$executable = "/Users/nick/pkg/Utilities/bin/3D/float/ConvertAvantsLandmarkFileToVTK";
  $moving_file = $output_prefix . "Moving.txt";
  $fixed_file = $output_prefix . "Fixed.txt";
		@args = ( $executable, $fixed_file, $moving_file, $output_prefix, 2 );
		system( @args ) == 0 or die "system @args failed: $?";

		@args = ( "mv", $output_prefix . "Moving.vtk", $output_prefix . ".vtk" );
		system( @args ) == 0 or die "system @args failed: $?";
		@args = ( "rm", $output_prefix . "Fixed.vtk" );
		system( @args ) == 0 or die "system @args failed: $?";
		@args = ( "rm", $output_prefix . "Moving.txt" );
		system( @args ) == 0 or die "system @args failed: $?";
		@args = ( "rm", $output_prefix . "Fixed.txt" );
		system( @args ) == 0 or die "system @args failed: $?";

		$executable = $utilityDirectory . "ConvertDeformationFieldToAvantsLandmarkFiles";
		$input_field = $outputDirectory . $prefix . "_" . $i . "_" . "InverseDeformationField.nii.gz";
  $output_prefix = $outputDirectory . $prefix . "_" . $i . "_" . "InverseDeformationField";
  $mask = $landc_output[0];
		@args = ( $executable, $input_field, $output_prefix, 1, $mask );
		system( @args ) == 0 or die "system @args failed: $?";

		$executable = "/Users/nick/pkg/Utilities/bin/3D/float/ConvertAvantsLandmarkFileToVTK";
  $moving_file = $output_prefix . "Moving.txt";
  $fixed_file = $output_prefix . "Fixed.txt";
		@args = ( $executable, $fixed_file, $moving_file, $output_prefix, 2 );
		system( @args ) == 0 or die "system @args failed: $?";

		@args = ( "mv", $output_prefix . "Moving.vtk", $output_prefix . ".vtk" );
		system( @args ) == 0 or die "system @args failed: $?";
		@args = ( "rm", $output_prefix . "Fixed.vtk" );
		system( @args ) == 0 or die "system @args failed: $?";
		@args = ( "rm", $output_prefix . "Moving.txt" );
		system( @args ) == 0 or die "system @args failed: $?";
		@args = ( "rm", $output_prefix . "Fixed.txt" );
		system( @args ) == 0 or die "system @args failed: $?";

  if ( $dimension == 2 )
    {
    $start = 2;
    }
  else
    {
    $start = 1;
    }
  for $d ( $start..3 )
    {
  		$executable = $utilityDirectory . "ConvertDeformationFieldToAvantsLandmarkFiles";
  		$input_field = $outputDirectory . $prefix . "_" . $i . "_" . "PrincipalStrain" . $d . ".nii.gz";
    $output_prefix = $outputDirectory . $prefix . "_" . $i . "_" . "PrincipalStrain" . $d;
    $mask = $landc_output[0];
		  @args = ( $executable, $input_field, $output_prefix, 1, $mask );
		  system( @args ) == 0 or die "system @args failed: $?";
		
				$executable = "/Users/nick/pkg/Utilities/bin/3D/float/ConvertAvantsLandmarkFileToVTK";
				$moving_file = $output_prefix . "Moving.txt";
				$fixed_file = $output_prefix . "Fixed.txt";
				@args = ( $executable, $fixed_file, $moving_file, $output_prefix, 2 );
				system( @args ) == 0 or die "system @args failed: $?";
		
				@args = ( "mv", $output_prefix . "Moving.vtk", $output_prefix . ".vtk" );
				system( @args ) == 0 or die "system @args failed: $?";
				@args = ( "rm", $output_prefix . "Fixed.vtk" );
				system( @args ) == 0 or die "system @args failed: $?";
				@args = ( "rm", $output_prefix . "Moving.txt" );
				system( @args ) == 0 or die "system @args failed: $?";
				@args = ( "rm", $output_prefix . "Fixed.txt" );
				system( @args ) == 0 or die "system @args failed: $?";
    }  
		
		## 
		# Create deformed grids
		##
		print "Create deformed grids\n";
		$executable = $utilityDirectory . "CreateDeformedGrid";
		$input = $outputDirectory . $prefix . "_" . $i . "_" . "DeformationField.nii.gz";
		$output = $outputDirectory . $prefix . "_" . $i . "_" . "DeformedGrid.txt";
		@args = ( $executable, $input, $output, 12, 12 );
		system( @args ) == 0 or die "system @args failed: $?";
		
		$input = $outputDirectory . $prefix . "_" . $i . "_" . "InverseDeformationField.nii.gz";
		$output = $outputDirectory . $prefix . "_" . $i . "_" . "InverseDeformedGrid.txt";
		@args = ( $executable, $input, $output, 12, 12 );
		system( @args ) == 0 or die "system @args failed: $?";
		
		$executable = "/Users/nick/pkg/Utilities/bin/3D/float/ConvertAvantsLandmarkFileToVTKLines";
		$input = $outputDirectory . $prefix . "_" . $i . "_" . "DeformedGrid.txt";
		$output = $outputDirectory . $prefix . "_" . $i . "_" . "DeformedGrid.vtk";
		@args = ( $executable, $input, $output );
		system( @args ) == 0 or die "system @args failed: $?";
		
		$input = $outputDirectory . $prefix . "_" . $i . "_" . "InverseDeformedGrid.txt";
		$output = $outputDirectory . $prefix . "_" . $i . "_" . "InverseDeformedGrid.vtk";
		@args = ( $executable, $input, $output );
		system( @args ) == 0 or die "system @args failed: $?";
		
		$executable = "rm";
		$input = $outputDirectory . $prefix . "_" . $i . "_" . "DeformedGrid.txt";
		@args = ( $executable, $input );
		system( @args ) == 0 or die "system @args failed: $?";
		
		$input = $outputDirectory . $prefix . "_" . $i . "_" . "InverseDeformedGrid.txt";
		@args = ( $executable, $input );
		system( @args ) == 0 or die "system @args failed: $?";
		
  }
