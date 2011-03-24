# request Bourne shell as shell for job
#$ -S /bin/sh

PREDIRECTORY=$1
POSTDIRECTORY=$2

PREPREFIX=$3
POSTPREFIX=$4


UTILITIESDIR=/mnt/data1/tustison/Utilities/bin/
ANTSDIR=/mnt/data1/tustison/PICSL/ANTS/bin/
TMPDIR=/mnt/data1/tustison/tmp

##
# Generate the vector fields from the control points

 ${UTILITIESDIR}/GenerateVectorFieldFromControlPointLattice 3 \
  ${POSTDIRECTORY}/psr_3.nii.gz \
  ${TMPDIR}/${POSTPREFIX}_Warp.nii.gz \
  ${POSTDIRECTORY}/${POSTPREFIX}.nii.gz 3

##
# Warp the expiration image
##

 ${ANTSDIR}/WarpImageMultiTransform 3 \
  ${PREDIRECTORY}/${PREPREFIX}.nii.gz \
  ${POSTDIRECTORY}/${PREPREFIX}_Warped.nii.gz \
  -R ${POSTDIRECTORY}/${POSTPREFIX}.nii.gz \
  ${TMPDIR}/${POSTPREFIX}_Warp.nii.gz

##
# Convert to RGB
##

 ${UTILITIESDIR}/ConvertScalarImageToRGB 3 \
  ${PREDIRECTORY}/${PREPREFIX}.nii.gz \
  ${PREDIRECTORY}/${PREPREFIX}_hot.mha \
  hot

 ${UTILITIESDIR}/ConvertScalarImageToRGB 3 \
  ${POSTDIRECTORY}/${PREPREFIX}_Warped.nii.gz \
  ${POSTDIRECTORY}/${PREPREFIX}_Warped_hot.mha \
  hot

##
# Create warped grid image
##

 ${ANTSDIR}/CreateWarpedGridImage 3 \
  ${TMPDIR}/${POSTPREFIX}_Warp.nii.gz \
  ${POSTDIRECTORY}/${PREPREFIX}_WarpedGrid.nii.gz \
  1x1x0


rm ${TMPDIR}/${POSTPREFIX}_Warp?vec.nii.gz


