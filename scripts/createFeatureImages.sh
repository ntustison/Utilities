#!/bin/bash

VERSION="0.0"

## need to change to put everything in ANTs (ImageMath?)

if [[ ! -d "$UTILPATH" ]];
  then
    echo We can\'t find the Utilities path -- does not seem to exist.  Please \(re\)define \$UTILPATH in your environment.
    exit 1
  fi
if [[ ! -d "$ANTSPATH" ]];
  then
    echo We can\'t find the ANTs path -- does not seem to exist.  Please \(re\)define \$ANTSPATH in your environment.
    exit 1
  fi

function Usage {
    cat <<USAGE

`basename $0` creates feature images for tumor segmentation using random decision forests.
   See <reference goes here>

Usage:

`basename $0` -d imageDimension
              -a inputImage(s)
              -t symmetricTemplate(s)
              -c clusterCenters
              -n imageNames
              -r neighgorhoodRadius
              -s smoothSigma
              -f differencePair
              -x maskImage
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
     -c:  cluster centers                       Array describing the intensity centers of the intensity normalized images.
                                                Need one for each input image.  Should be of the form:  e.g. 0.14x0.57x0.37x0.83x0.95
                                                (for 7 classes: csf, gm, wm, edema, tumor, core, enhancement).  Note that either the cluster
                                                centers are specified (for testing) or the truth labels (for training) but not both.
     -g:  truth labels                          Truth labels.  Note that either the cluster centers are specified (for testing)
                                                or the truth labels (for training) but not both.
     -t:  symmetric anatomical templates        Symmetric templates.  Need to be specified in the same order as
                                                the input anatomical images.
     -x:  mask image                            Mask image defining the region of interest.
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
     -p:  Brain segmentation priors             Tissue *probability* priors. Specified using c-style formatting, e.g.
                                                -p labelsPriors%02d.nii.gz.
     -f:  difference pair                       pair of indices \"i.e. 3x1\" to create difference image \"image[3] - image[1]\"
     -l:  tumor core label                      used to create distance feature map (default = 5)
     -n   imageNames                            used in the naming of the images (otherwise, labeled IMAGE0, IMAGE1, etc)

USAGE
    exit 1
}

echoParameters() {
    cat <<PARAMETERS

    Using createFeatureImages with the following arguments:
      image dimension         = ${DIMENSION}
      anatomical image        = ${ANATOMICAL_IMAGES[@]}
      symmetric templates     = ${SYMMETRIC_TEMPLATES[@]}
      cluster centers         = ${CLUSTER_CENTERS[@]}
      image names             = ${IMAGE_NAMES[@]}
      radius                  = ${RADIUS}
      smoothing sigma         = ${SMOOTHING_SIGMA}
      priors                  = ${SEGMENTATION_PRIOR}
      difference pairs        = ${DIFFERENCE_PAIRS[@]}
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
CLUSTER_CENTERS=()
IMAGE_NAMES=()
DIFFERENCE_PAIRS=()

SEGMENTATION_PRIOR=""

MASK_IMAGE=""

RADIUS=2
SMOOTHING_SIGMA=0
CORE_LABEL=5

################################################################################
#
# Programs and their parameters
#
################################################################################

if [[ $# -lt 3 ]] ; then
  Usage >&2
  exit 1
else
  while getopts "a:c:d:f:g:h:l:n:o:p:r:s:t:x:" OPT
    do
      case $OPT in
          a) #anatomical image
       ANATOMICAL_IMAGES[${#ANATOMICAL_IMAGES[@]}]=$OPTARG
       ;;
          c) # cluster centers
       CLUSTER_CENTERS[${#CLUSTER_CENTERS[@]}]=$OPTARG
       ;;
          d) #dimensions
       DIMENSION=$OPTARG
       if [[ ${DIMENSION} -gt 4 || ${DIMENSION} -lt 2 ]];
         then
           echo " Error:  ImageDimension must be 2, 3, or 4 "
           exit 1
         fi
       ;;
          f)
       DIFFERENCE_PAIRS[${#DIFFERENCE_PAIRS[@]}]=$OPTARG
       ;;
          g)
       TRUTH_LABELS=$OPTARG
       ;;
          h) #help
       Usage >&2
       exit 0
       ;;
          l)
       CORE_LABEL=$OPTARG
       ;;
          n)
       IMAGE_NAMES[${#IMAGE_NAMES[@]}]=$OPTARG
       ;;
          o) #output prefix
       OUTPUT_PREFIX=$OPTARG
       ;;
          p) #brain segmentation label prior image
       SEGMENTATION_PRIOR=$OPTARG
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
          x)
       MASK_IMAGE=$OPTARG
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
    if [[ -z ${IMAGE_NAMES[$i]} ]];
      then
        IMAGE_NAMES[$i]=IMAGE${i}
      fi
  done

for (( i = 0; i < ${#SYMMETRIC_TEMPLATES[@]}; i++ ))
  do
    if [[ ! -f ${SYMMETRIC_TEMPLATES[$i]} ]];
      then
        echo "The specified template image \"${SYMMETRIC_TEMPLATES[$i]}\" does not exist."
        exit 1
      fi
  done

if [[ ! -f ${MASK_IMAGE} ]];
  then
    echo "The specified mask image \"${MASK_IMAGE}\" does not exist."
    exit 1
  fi

if [[ ${#ANATOMICAL_IMAGES[@]} -ne ${#SYMMETRIC_TEMPLATES[@]} ]];
  then
      echo "The number of symmetric templates does not match the number of anatomical images."
      exit 1
  fi

if [[ ${#ANATOMICAL_IMAGES[@]} -ne ${#CLUSTER_CENTERS[@]} ]]
  then
      echo "The number of cluster center arrays does not match the number of anatomical images."
      exit 1
  fi

CLUSTERS=( `echo ${CLUSTER_CENTERS[0]} | tr 'x' ' '` )
NUMBER_OF_LABELS=${#CLUSTERS[@]}
for (( i = 1; i < ${#CLUSTER_CENTERS[@]}; i++ ))
  do
    CLUSTERS=( `echo ${CLUSTER_CENTERS[i]} | tr 'x' ' '` )
    if [[ ${#CLUSTERS[@]} -ne NUMBER_OF_LABELS ]];
      then
        echo "The number of labels is not equal across the cluster center arrays."
        exit 1
      fi
  done

OUTPUT_DIR=${OUTPUT_PREFIX%\/*}
if [[ ! -d $OUTPUT_DIR ]];
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
    OUTPUT_IMAGE=${OUTPUT_PREFIX}${IMAGE_NAMES[$i]}_MEAN_RADIUS_${RADIUS}.${OUTPUT_SUFFIX}
    if [[ ! -f ${OUTPUT_IMAGE} ]];
      then
        logCmd $STATS ${DIMENSION} ${ANATOMICAL_IMAGES[$i]} $OUTPUT_IMAGE 0 ${RADIUS}
      fi

    # standard deviation image
    OUTPUT_IMAGE=${OUTPUT_PREFIX}${IMAGE_NAMES[$i]}_SIGMA_RADIUS_${RADIUS}.${OUTPUT_SUFFIX}
    if [[ ! -f ${OUTPUT_IMAGE} ]];
      then
        logCmd $STATS ${DIMENSION} ${ANATOMICAL_IMAGES[$i]} $OUTPUT_IMAGE 4 ${RADIUS}
      fi

    # skewness image
    OUTPUT_IMAGE=${OUTPUT_PREFIX}${IMAGE_NAMES[$i]}_SKEWNESS_RADIUS_${RADIUS}.${OUTPUT_SUFFIX}
    if [[ ! -f ${OUTPUT_IMAGE} ]];
      then
        logCmd $STATS ${DIMENSION} ${ANATOMICAL_IMAGES[$i]} $OUTPUT_IMAGE 5 ${RADIUS}
      fi

    # entropy image
    OUTPUT_IMAGE=${OUTPUT_PREFIX}${IMAGE_NAMES[$i]}_ENTROPY_RADIUS_${RADIUS}.${OUTPUT_SUFFIX}
    if [[ ! -f ${OUTPUT_IMAGE} ]];
      then
        logCmd $STATS ${DIMENSION} ${ANATOMICAL_IMAGES[$i]} $OUTPUT_IMAGE 7 ${RADIUS}
      fi

  done

################################################################################
#
# Calculate normalized distance image (inside mask values between [0,1])
#
################################################################################

DISTANCE=${UTILPATH}/GenerateDistanceImage

OUTPUT_IMAGE=${OUTPUT_PREFIX}NORMALIZED_DISTANCE.${OUTPUT_SUFFIX}
if [[ ! -f ${OUTPUT_IMAGE} ]];
  then
    logCmd $DISTANCE ${DIMENSION} $MASK_IMAGE $OUTPUT_IMAGE 0
    logCmd ${ANTSPATH}ImageMath 3 $OUTPUT_IMAGE m $MASK_IMAGE $OUTPUT_IMAGE
    logCmd ${ANTSPATH}ImageMath 3 $OUTPUT_IMAGE Normalize $OUTPUT_IMAGE
  fi

################################################################################
#
# Compute transform between anatomical image and symmetric template
#
################################################################################

OUTPUT_REGISTRATION_PREFIX=${OUTPUT_PREFIX}ANTs_REGISTRATION

REG_BASE="${ANTSPATH}/antsRegistration -d ${DIMENSION} -w [0.01,0.995] -o ${OUTPUT_REGISTRATION_PREFIX}"
REG_LEV0="-r [${ANATOMICAL_IMAGES[0]},${SYMMETRIC_TEMPLATES[0]},1]"
REG_LEV1="-t Rigid[0.1] -m MI[${ANATOMICAL_IMAGES[0]},${SYMMETRIC_TEMPLATES[0]},1,32,Regular,0.25] -s 2x1x0 -f 4x2x1 -c [500x250x100,1e-8,15]"
REG_LEV2="-t Affine[0.1] -m MI[${ANATOMICAL_IMAGES[0]},${SYMMETRIC_TEMPLATES[0]},1,32,Regular,0.25] -s 2x1x0 -f 4x2x1 -c [500x250x100,1e-8,15]"
REG_LEV3="-t BSplineSyN[0.1,10x10x11,0x0x0] -m CC[${ANATOMICAL_IMAGES[0]},${SYMMETRIC_TEMPLATES[0]},1,4] -s 2x1x0 -f 4x2x1 -c [70x50x10,1e-8,15]"

AFFINE=${OUTPUT_REGISTRATION_PREFIX}0GenericAffine.mat
WARP=${OUTPUT_REGISTRATION_PREFIX}1Warp.nii.gz
INVERSE_WARP=${OUTPUT_REGISTRATION_PREFIX}1InverseWarp.nii.gz

if [[ ! -f $WARP ]];
  then
    logCmd $REG_BASE $REG_LEV0 $REG_LEV1 $REG_LEV2 $REG_LEV3
  fi

# log jacobian image
OUTPUT_IMAGE=${OUTPUT_PREFIX}LOG_JACOBIAN.${OUTPUT_SUFFIX}
if [[ ! -f ${OUTPUT_IMAGE} ]];
  then
    logCmd ${ANTSPATH}/ANTSJacobian ${DIMENSION} $WARP $OUTPUT_REGISTRATION_PREFIX 1
    logCmd mv ${OUTPUT_REGISTRATION_PREFIX}logjacobian.nii.gz $OUTPUT_IMAGE
  fi


for (( i = 0; i < ${#ANATOMICAL_IMAGES[@]}; i++ ))
  do

    # symmetric template difference images
    OUTPUT_IMAGE=${OUTPUT_PREFIX}${IMAGE_NAMES[$i]}_SYMMETRIC_TEMPLATE_DIFFERENCE.${OUTPUT_SUFFIX}
    if [[ ! -f ${OUTPUT_IMAGE} ]];
      then
        logCmd ${ANTSPATH}/antsApplyTransforms -d ${DIMENSION} -n BSpline -r ${ANATOMICAL_IMAGES[$i]} -i ${SYMMETRIC_TEMPLATES[$i]} -o ${OUTPUT_IMAGE} -t $WARP -t $AFFINE
        logCmd ${ANTSPATH}/ImageMath ${DIMENSION} $OUTPUT_IMAGE - ${ANATOMICAL_IMAGES[$i]} $OUTPUT_IMAGE
        logCmd ${ANTSPATH}/SmoothImage ${DIMENSION} $OUTPUT_IMAGE $SMOOTHING_SIGMA $OUTPUT_IMAGE 1
      fi

    # contralateral difference images
    OUTPUT_IMAGE=${OUTPUT_PREFIX}${IMAGE_NAMES[$i]}_CONTRALATERAL_DIFFERENCE.${OUTPUT_SUFFIX}
    if [[ ! -f ${OUTPUT_IMAGE} ]];
      then
        logCmd ${ANTSPATH}/antsApplyTransforms -d ${DIMENSION} -n BSpline -i ${ANATOMICAL_IMAGES[$i]} -r ${SYMMETRIC_TEMPLATES[$i]} -o ${OUTPUT_IMAGE} -t [$AFFINE,1] -t $INVERSE_WARP
        logCmd ${UTILPATH}/FlipImage ${DIMENSION} $OUTPUT_IMAGE $OUTPUT_IMAGE 1x0x0
        logCmd ${UTILPATH}/ChangeImageInformation ${DIMENSION} $OUTPUT_IMAGE $OUTPUT_IMAGE 4 ${SYMMETRIC_TEMPLATES[$i]}
        logCmd ${ANTSPATH}/antsApplyTransforms -d ${DIMENSION} -n BSpline -i $OUTPUT_IMAGE -r ${ANATOMICAL_IMAGES[$i]} -o ${OUTPUT_IMAGE} -t $WARP -t $AFFINE
        logCmd ${ANTSPATH}/ImageMath ${DIMENSION} $OUTPUT_IMAGE - ${ANATOMICAL_IMAGES[$i]} $OUTPUT_IMAGE
        logCmd ${ANTSPATH}/SmoothImage 3 $OUTPUT_IMAGE $SMOOTHING_SIGMA $OUTPUT_IMAGE 1
      fi

  done

################################################################################
#
# Construct GMM or MAP-MRF probability images for each anatomical image
# Also create geometric feature images for each atropos labeled output
#
################################################################################

for (( i = 0; i < ${#ANATOMICAL_IMAGES[@]}; i++ ))
  do

    if [[ -n ${SEGMENTATION_PRIOR} ]];
      then

        OUTPUT_ATROPOS_IMAGE=${OUTPUT_PREFIX}${IMAGE_NAMES[$i]}_ATROPOS_MAP_MRF.${OUTPUT_SUFFIX}
        OUTPUT_ATROPOS_IMAGE_POSTERIORS=${OUTPUT_PREFIX}${IMAGE_NAMES[$i]}_ATROPOS_MAP_MRF_POSTERIORS%d.${OUTPUT_SUFFIX}

        if [[ ! -f ${OUTPUT_ATROPOS_IMAGE} ]];
          then

            bash ${ANTSPATH}/antsAtroposN4.sh \
              -d ${DIMENSION} \
              -a ${ANATOMICAL_IMAGES[$i]} \
              -x ${MASK_IMAGE} \
              -m 1 \
              -n 5 \
              -c ${NUMBER_OF_LABELS} \
              -l 3 \
              -l 2 \
              -p ${SEGMENTATION_PRIOR} \
              -w 0.0 \
              -o ${OUTPUT_PREFIX}${IMAGE_NAMES[$i]} \
              -k 0 \
              -s ${OUTPUT_SUFFIX}

            f=${OUTPUT_PREFIX}${IMAGE_NAMES[$i]}Segmentation.${OUTPUT_SUFFIX}
            newfile=${f/Segmentation/_ATROPOS_MAP_MRF};
            logCmd mv $f $newfile

            for f in `ls ${OUTPUT_PREFIX}${IMAGE_NAMES[$i]}SegmentationPosteriors*.nii.gz`;
              do
                newfile=${f/SegmentationPosteriors/_ATROPOS_MAP_MRF_POSTERIORS};
                logCmd mv $f $newfile
              done
          fi

        OUTPUT_ATROPOS_FEATURES_PREFIX=${OUTPUT_PREFIX}${IMAGE_NAMES[$i]}_ATROPOS_MAP_MRF_
        OUTPUT_ATROPOS_DISTANCE_IMAGE=${OUTPUT_ATROPOS_FEATURES_PREFIX}LABEL${CORE_LABEL}_DISTANCE.${OUTPUT_SUFFIX}

        if [[ ! -f ${OUTPUT_ATROPOS_FEATURES_PREFIX}ECCENTRICITY.nii.gz ]];
          then
            logCmd ${UTILPATH}/GetConnectedComponentsFeatureImages ${DIMENSION} ${OUTPUT_ATROPOS_IMAGE} ${OUTPUT_ATROPOS_FEATURES_PREFIX}
            logCmd ${ANTSPATH}/ThresholdImage ${DIMENSION} ${OUTPUT_ATROPOS_IMAGE} ${OUTPUT_ATROPOS_DISTANCE_IMAGE} ${CORE_LABEL} ${CORE_LABEL} 1 0
            logCmd ${UTILPATH}/GenerateDistanceImage ${DIMENSION} ${OUTPUT_ATROPOS_DISTANCE_IMAGE} ${OUTPUT_ATROPOS_DISTANCE_IMAGE} 1
          fi

      else

        OUTPUT_ATROPOS_IMAGE=${OUTPUT_PREFIX}${IMAGE_NAMES[$i]}_ATROPOS_GMM.${OUTPUT_SUFFIX}
        OUTPUT_ATROPOS_IMAGE_POSTERIORS=${OUTPUT_PREFIX}${IMAGE_NAMES[$i]}_ATROPOS_GMM_POSTERIORS%d.${OUTPUT_SUFFIX}

        SEG_BASE="${ANTSPATH}/Atropos -d ${DIMENSION} -a ${ANATOMICAL_IMAGES[$i]}"
        SEG_0="-i KMeans[${NUMBER_OF_LABELS},${CLUSTER_CENTERS[i]}] -p Socrates[1] -x $MASK_IMAGE"
        SEG_1="-c [0,0] -k Gaussian -m [0.1,1x1x1]"
        SEG_2="-o [${OUTPUT_ATROPOS_IMAGE},${OUTPUT_ATROPOS_IMAGE_POSTERIORS}]"

        if [[ ! -f ${OUTPUT_ATROPOS_IMAGE} ]];
          then
            logCmd $SEG_BASE $SEG_0 $SEG_1 $SEG_2
          fi

        OUTPUT_ATROPOS_FEATURES_PREFIX=${OUTPUT_PREFIX}${IMAGE_NAMES[$i]}_ATROPOS_GMM_
        OUTPUT_ATROPOS_DISTANCE_IMAGE=${OUTPUT_ATROPOS_FEATURES_PREFIX}LABEL${CORE_LABEL}_DISTANCE.${OUTPUT_SUFFIX}

        if [[ ! -f ${OUTPUT_ATROPOS_FEATURES_PREFIX}ECCENTRICITY.nii.gz ]];
          then
            logCmd ${UTILPATH}/GetConnectedComponentsFeatureImages ${DIMENSION} ${OUTPUT_ATROPOS_IMAGE} ${OUTPUT_ATROPOS_FEATURES_PREFIX}
            logCmd ${ANTSPATH}/ThresholdImage ${DIMENSION} ${OUTPUT_ATROPOS_IMAGE} ${OUTPUT_ATROPOS_DISTANCE_IMAGE} ${CORE_LABEL} ${CORE_LABEL} 1 0
            logCmd ${UTILPATH}/GenerateDistanceImage ${DIMENSION} ${OUTPUT_ATROPOS_DISTANCE_IMAGE} ${OUTPUT_ATROPOS_DISTANCE_IMAGE} 1
          fi

      fi
  done

################################################################################
#
# Create difference pairs images
#
################################################################################

for (( i = 0; i < ${#DIFFERENCE_PAIRS[@]}; i++ ))
  do
    DIFFERENCE_PAIR_ARRAY=( ${DIFFERENCE_PAIRS//x/ } )

    if [[ ${#DIFFERENCE_PAIR_ARRAY[@]} -eq 2 ]];
      then
        MINUEND_INDEX=${DIFFERENCE_PAIR_ARRAY[0]}
        SUBTRAHEND_INDEX=${DIFFERENCE_PAIR_ARRAY[1]}

        OUTPUT_DIFFERENCE_PAIR_IMAGE=${OUTPUT_PREFIX}${IMAGE_NAMES[${MINUEND_INDEX}]}_${IMAGE_NAMES[${SUBTRAHEND_INDEX}]}_DIFFERENCE.${OUTPUT_PREFIX}

        logCmd ${ANTSPATH}/ImageMath ${DIMENSION} ${OUTPUT_DIFFERENCE_PAIR_IMAGE} - ${ANATOMICAL_IMAGES[${MINUEND_INDEX}]} ${ANATOMICAL_IMAGES[${SUBTRAHEND_INDEX}]}
      fi
  done

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
