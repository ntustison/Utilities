#!/usr/bin/perl

$baseDirectory = "/Users/nick/pkg/Projects/HeliumLungStudies/NormalImages/PCALevelSets/PCA_2/";
$directory = "${baseDirectory}/PCAModes/";

for( $i = 1; $i <= 52; $i++ )
  {
  print "Calculating mode ${i} images.\n";
  for( $j = 1; $j <= 3; $j++ )
    {
    if( $i < 10 )
      {
      $index = "00${i}";
      }
    else
      {
      $index = "0${i}";
      }

    `UnaryOperateImage 3 ${directory}/pca_modes_${index}.nii.gz x $j temp.nii.gz`;
    `BinaryOperateImages 3 ${directory}/pca_modes_000.nii.gz + temp.nii.gz ${directory}/pca_sample_modes_${index}_plus${j}.nii.gz`;
    `ThresholdImage 3 ${directory}/pca_sample_modes_${index}_plus${j}.nii.gz ${directory}/pca_sample_modes_${index}_plus${j}.nii.gz -10000000 0 1 0`;

    `UnaryOperateImage 3 ${directory}/pca_modes_${index}.nii.gz x -${j} temp.nii.gz`;
    `BinaryOperateImages 3 ${directory}/pca_modes_000.nii.gz + temp.nii.gz ${directory}/pca_sample_modes_${index}_minus${j}.nii.gz`;
    `ThresholdImage 3 ${directory}/pca_sample_modes_${index}_minus${j}.nii.gz ${directory}/pca_sample_modes_${index}_minus${j}.nii.gz -10000000 0 1 0`;

    }
  }
