#!/bin/bash

VERSION="0.0"

## need to change to put everything in ANTs (ImageMath?)
UTILPATH=/Users/ntustison/Pkg/Utilities/bin

if [[ ! -s ${ANTSPATH}/N4BiasFieldCorrection ]]; then
  echo we cant find the N4 program -- does not seem to exist.  please \(re\)define \$ANTSPATH in your environment.
  exit
fi
if [[ ! -s ${ANTSPATH}/Atropos ]]; then
  echo we cant find the Atropos program -- does not seem to exist.  please \(re\)define \$ANTSPATH in your environment.
  exit
fi
if [[ ! -s ${ANTSPATH}/ImageMath ]]; then
  echo we cant find the Atropos program -- does not seem to exist.  please \(re\)define \$ANTSPATH in your environment.
  exit
fi

function Usage {
    cat <<USAGE

`basename $0` creates feature images for tumor segmentation using random decision forests.
   See <reference goes here>

Usage:

`basename $0` -d imageDimension
              -a inputImage(s)
              -r neighgorhoodRadius
              -s smoothSigma
              -t symmetricTemplate(s)
              -o outputPrefix

Example:

  bash $0 -d 3 -a t1.nii.gz -t symmetricTemplate.nii.gz -o output

Required arguments:

     -d:  image dimension                       2 or 3 (for 2- or 3-dimensional image)
     -a:  input image                           Anatomical image(s) comprised of flair, t2, t1, or t1 contrast.
                                                It is assumed that the following preprocessing steps have
                                                already been performed:
                                                   * brain extraction
                                                   * bias correction
                                                   * intensity normalized (i.e. histogram matching and intensity truncation)
     -t:  symmetric anatomical templates        Symmetric templates.  Need to be specified in the same order as
                                                the input anatomical images.
     -o:  output prefix                         The following images are created:
                                                  * ${OUTPUT_PREFIX}N4Corrected.${OUTPUT_SUFFIX}
                                                  * ${OUTPUT_PREFIX}Segmentation.${OUTPUT_SUFFIX}
                                                  * ${OUTPUT_PREFIX}SegmentationPosteriors.${OUTPUT_SUFFIX}

Optional arguments:

     -r   radius                                Neighborhood voxel radius for local statistics image (default = 2)
     -s   smoothingSigma                        Standard deviation of gaussian smoothing kernel (in physical space) for
                                                the following images:
                                                  * symmetric template difference
                                                  * contralateral difference

USAGE
    exit 1
}

echoParameters() {
    cat <<PARAMETERS

    Using createFeatureImages with the following arguments:
      image dimension         = ${DIMENSION}
      anatomical image        = ${ANATOMICAL_IMAGES[@]}
      symmetric template      = ${SYMMETRIC_TEMPLAGE}
      radius                  = ${RADIUS}
      smoothing sigma         = ${SMOOTHING_SIGMA}
      output prefix           = ${OUTPUT_PREFIX}

PARAMETERS
}


#    local  myresult='some value'
#    echo "$myresult"

# Echos a command to both stdout and stderr, then runs it
function logCmd() {
  cmd="$*"
  echo "BEGIN >>>>>>>>>>>>>>>>>>>>"
  echo $cmd
  logCmdOutput=$( $cmd | tee /dev/tty )
  echo "END   <<<<<<<<<<<<<<<<<<<<"
  echo
  echo
}

################################################################################
#
# Main routine
#
################################################################################

HOSTNAME=`hostname`
DATE=`date`

CURRENT_DIR=`pwd`/
OUTPUT_DIR=${CURRENT_DIR}/tmp$RANDOM/
OUTPUT_PREFIX=${OUTPUT_DIR}/tmp
OUTPUT_SUFFIX="nii.gz"

DIMENSION=3

ANATOMICAL_IMAGES=()
SYMMETRIC_TEMPLATE=()

RADIUS=2
SMOOTHING_SIGMA=0

################################################################################
#
# Programs and their parameters
#
################################################################################

if [[ $# -lt 3 ]] ; then
  Usage >&2
  exit 1
else
  while getopts "a:d:h:o:r:s:t:" OPT
    do
      case $OPT in
          a) #anatomical t1 image
       ANATOMICAL_IMAGES[${#ANATOMICAL_IMAGES[@]}]=$OPTARG
       ;;
          d) #dimensions
       DIMENSION=$OPTARG
       if [[ ${DIMENSION} -gt 4 || ${DIMENSION} -lt 2 ]];
         then
           echo " Error:  ImageDimension must be 2, 3, or 4 "
           exit 1
         fi
       ;;
          h) #help
       Usage >&2
       exit 0
       ;;
          o) #output prefix
       OUTPUT_PREFIX=$OPTARG
       ;;
          r)
       RADIUS=$OPTARG
       ;;
          s)
       SMOOTHING_SIGMA=$OPTARG
       ;;
          t)
       SYMMETRIC_TEMPLATES[${#SYMMETRIC_TEMPLATES[@]}]=$OPTARG
       ;;
          *) # getopts issues an error message
       echo "ERROR:  unrecognized option -$OPT $OPTARG"
       exit 1
       ;;
      esac
  done
fi

################################################################################
#
# Preliminaries:
#  1. Check existence of inputs
#  2. Figure out output directory and mkdir if necessary
#
################################################################################

for (( i = 0; i < ${#ANATOMICAL_IMAGES[@]}; i++ ))
  do
  if [[ ! -f ${ANATOMICAL_IMAGES[$i]} ]];
    then
      echo "The specified image \"${ANATOMICAL_IMAGES[$i]}\" does not exist."
      exit 1
    fi
  done

for (( i = 0; i < ${#SYMMETRIC_TEMPLATES[@]}; i++ ))
  do
  if [[ ! -f ${SYMMETRIC_TEMPLATES[$i]} ]];
    then
      echo "The specified image \"${SYMMETRIC_TEMPLATES[$i]}\" does not exist."
      exit 1
    fi
  done

if [[ ! ${ANATOMICAL_IMAGES[$i]} -neq ${SYMMETRIC_TEMPLATES[$i]} ]]
  then
      echo "The number of symmetric templates does not match the number of anatomical images."
      exit 1
  fi

OUTPUT_DIR=${OUTPUT_PREFIX%\/*}
if [[ ! -e $OUTPUT_PREFIX ]];
  then
    echo "The output directory \"$OUTPUT_DIR\" does not exist. Making it."
    mkdir -p $OUTPUT_DIR
  fi

echoParameters >&2

echo "---------------------  Running `basename $0` on $HOSTNAME  ---------------------"

time_start=`date +%s`

################################################################################
#
# Calculate statistics images
#
################################################################################

STATS=${UTILPATH}/CalculateStatisticsImage

for (( i = 0; i < ${#ANATOMICAL_IMAGES[@]}; i++ ))
  do

    # mean image
    OUTPUT_IMAGE=${OUTPUT_PREFIX}IMAGE${i}_MEAN_RADIUS_${RADIUS}.${OUTPUT_SUFFIX}
    if [[ ! -f ${OUTPUT_IMAGE} ]];
      then
        logCmd $STATS ${DIMENSION} ${ANATOMICAL_IMAGES[$i]} $OUTPUT_IMAGE 0 ${RADIUS}
      fi

    # standard deviation image
    OUTPUT_IMAGE=${OUTPUT_PREFIX}IMAGE${i}_SIGMA_RADIUS_${RADIUS}.${OUTPUT_SUFFIX}
    if [[ ! -f ${OUTPUT_IMAGE} ]];
      then
        logCmd $STATS ${DIMENSION} ${ANATOMICAL_IMAGES[$i]} $OUTPUT_IMAGE 4 ${RADIUS}
      fi

    # skewness image
    OUTPUT_IMAGE=${OUTPUT_PREFIX}IMAGE${i}_SKEWNESS_RADIUS_${RADIUS}.${OUTPUT_SUFFIX}
    if [[ ! -f ${OUTPUT_IMAGE} ]];
      then
        logCmd $STATS ${DIMENSION} ${ANATOMICAL_IMAGES[$i]} $OUTPUT_IMAGE 5 ${RADIUS}
      fi

  done

################################################################################
#
# Calculate normalized distance image (inside mask values between [0,1])
#
################################################################################

DISTANCE=${UTILPATH}/GenerateDistanceImage

OUTPUT_IMAGE=${OUTPUT_PREFIX}NORMALIZED_DISTANCE.${OUTPUT_SUFFIX}
TMP_IMAGE=${OUTPUT_PREFIX}TMP.${OUTPUT_SUFFIX}
if [[ -f ${OUTPUT_IMAGE} ]];
  then
    logCmd ${ANTSPATH}ThresholdImage 3 ${ANATOMICAL_IMAGES[0]} $TMP_IMAGE 0 0 0 1
    logCmd $DISTANCE ${DIMENSION} $TMP_IMAGE $OUTPUT_IMAGE 0
    logCmd ${ANTSPATH}ImageMath 3 $OUTPUT_IMAGE m $TMP_IMAGE $OUTPUT_IMAGE
    logCmd ${ANTSPATH}ImageMath 3 $OUTPUT_IMAGE Normalize $OUTPUT_IMAGE
    logCmd rm $TMP_IMAGE
  fi

################################################################################
#
# Compute transform between anatomical image and symmetric template
#
################################################################################

OUTPUT_REGISTRATION_PREFIX=${OUTPUT_PREFIX}ANTs_REGISTRATION

REG_BASE="${ANTSPATH}/antsRegistration -d ${DIMENSION} -w [0.025,0.975] -o ${OUTPUT_REGISTRATION_PREFIX}"
REG_LEV0="-t Rigid[0.2] -m MI[${ANATOMICAL_IMAGES[0]},${SYMMETRIC_TEMPLATES[0]},1,32] -s 2x1x0 -f 4x2x1 -c [500x250x100,1e-8,15]"
REG_LEV1="-t Affine[0.2] -m MI[${ANATOMICAL_IMAGES[0]},${SYMMETRIC_TEMPLATES[0]},1,32] -s 2x1x0 -f 4x2x1 -c [500x250x100,1e-8,15]"
REG_LEV2="-t SyN[0.1,3,0.2] -m CC[${ANATOMICAL_IMAGES[0]},${SYMMETRIC_TEMPLATES[0]},1,4] -s 3x2x1x0 -f 8x4x2x1 -c [100x50x40x0,1e-8,15]"

logCmd $REG_BASE $REG_LEV0 $REG_LEV1 $REG_LEV2

AFFINE=${OUTPUT_REGISTRATION_PREFIX}0GenericAffine.mat
WARP=${OUTPUT_REGISTRATION_PREFIX}1Warp.nii.gz
INVERSE_WARP=${OUTPUT_REGISTRATION_PREFIX}1InverseWarp.nii.gz

# log jacobian image
OUTPUT_IMAGE=${OUTPUT_PREFIX}LOG_JACOBIAN.${OUTPUT_SUFFIX}
if [[ ! -f ${OUTPUT_IMAGE} ]];
  then
    logCmd ${ANTSPATH}/ANTSJacobian ${DIMENSION} $WARP $OUTPUT_REGISTRATION_PREFIX
    logCmd mv ${OUTPUT_REGISTRATION_PREFIX}logjacobian.nii.gz $OUTPUT_IMAGE
  fi


for (( i = 0; i < ${#ANATOMICAL_IMAGES[@]}; i++ ))
  do

    # symmetric template difference images
    OUTPUT_IMAGE=${OUTPUT_PREFIX}IMAGE${i}_SYMMETRIC_TEMPLATE_DIFFERENCE.${OUTPUT_SUFFIX}
    if [[ ! -f ${OUTPUT_IMAGE} ]];
      then
        logCmd ${ANTSPATH}/antsApplyTransforms -d ${DIMENSION} -n BSpline -r ${ANATOMICAL_IMAGES[$i]} -i ${SYMMETRIC_TEMPLATES[$i]} -o ${OUTPUT_IMAGE} -t $WARP -t $AFFINE
        logCmd ${ANTSPATH}/ImageMath ${DIMENSION} $OUTPUT_IMAGE - ${ANATOMICAL_IMAGES[$i]} $OUTPUT_IMAGE
        logCmd ${ANTSPATH}/SmoothImage ${DIMENSION} $OUTPUT_IMAGE $SMOOTHING_SIGMA $OUTPUT_IMAGE 1
      fi

    # contralateral difference images
    OUTPUT_IMAGE=${OUTPUT_PREFIX}CONTRALATERAL_DIFFERENCE.${OUTPUT_SUFFIX}
    if [[ ! -f ${OUTPUT_IMAGE} ]];
      then
        logCmd ${ANTSPATH}/antsApplyTransforms -d ${DIMENSION} -n BSpline -i ${ANATOMICAL_IMAGES[$i]} -r ${SYMMETRIC_TEMPLATES[$i]} -o ${OUTPUT_IMAGE} -t [$AFFINE,1] -t $INVERSE_WARP
        logCmd ${ANTSPATH}/PermuteFlipImageOrientationAxes -d ${DIMENSION} $OUTPUT_IMAGE 0 1 2 1 0 0
        logCmd ${UTILPATH}/ChangeImageInformation ${DIMENSION} $OUTPUT_IMAGE $OUTPUT_IMAGE 4 ${SYMMETRIC_TEMPLATES[$i]}
        logCmd ${ANTSPATH}/antsApplyTransforms -d ${DIMENSION} -n BSpline -i $OUTPUT_IMAGE -r ${ANATOMICAL_IMAGES[$i]} -o ${OUTPUT_IMAGE} -t $WARP -t $AFFINE
        logCmd ${ANTSPATH}/ImageMath ${DIMENSION} $OUTPUT_IMAGE - ${ANATOMICAL_IMAGES[$i]} $OUTPUT_IMAGE
        logCmd ${ANTSPATH}/SmoothImage 3 $OUTPUT_IMAGE $SMOOTHING_SIGMA $OUTPUT_IMAGE 1
      fi

  done

################################################################################
#
# Construct GMM for each class
#
################################################################################






################################################################################
#
# Construct skeleton for each connected component and assign each pixel the
# value of the number of pixels comprising the skeleton for that particular
# component.
#
################################################################################




################################################################################
#
# End of main routine
#
################################################################################

time_end=`date +%s`
time_elapsed=$((time_end - time_start))

echo
echo "--------------------------------------------------------------------------------------"
echo " Done with creating feature images"
echo " Script executed in $time_elapsed seconds"
echo " $(( time_elapsed / 3600 ))h $(( time_elapsed %3600 / 60 ))m $(( time_elapsed % 60 ))s"
echo "--------------------------------------------------------------------------------------"

