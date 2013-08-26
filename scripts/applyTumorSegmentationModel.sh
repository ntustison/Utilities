#!/bin/bash

VERSION="0.0"

## need to change so that SCRIPTSPATH is where createFeatureImages.sh
## and the other scripts are located.

SCRIPTSPATH=/Users/ntustison/Documents/Academic/SubmittedPapers/BRATS2013/Scripts/

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

`basename $0`
              -d imageDimension
              -a inputImagesPrefix
              -m randomForestModel
              -n modalityString
              -t symmetricTemplate(s)
              -c clusterCenters
              -r neighgorhoodRadius
              -s smoothSigma
              -x maskImage
              -f differencePair
              -o outputPrefix

Example:

  bash $0
          -d 3
          -m /home/njt4n/share/Data/Tumor/BRATS-1/Model/highGradeModel.RData
          -x /home/njt4n/share/Data/Tumor/BRATS-1/Images/BRATS_HG0001_MASK.nii.gz
          -a /home/njt4n/share/Data/Tumor/BRATS-1/Images/BRATS_HG0001_FLAIR.nii.gz
          -n FLAIR
          -c clusterCenterFLAIR
          -t /home/njt4n/share/Data/Kirby/FLAIR_template.nii.gz
          -a /home/njt4n/share/Data/Tumor/BRATS-1/Images/BRATS_HG0001_T1.nii.gz
          -n T1
          -c clusterCenterT1
          -t /home/njt4n/share/Data/Kirby/T1_template.nii.gz
          -a /home/njt4n/share/Data/Tumor/BRATS-1/Images/BRATS_HG0001_T1C.nii.gz
          -n T1C
          -c clusterCenterT1C
          -t /home/njt4n/share/Data/Kirby/T1C_template.nii.gz
          -a /home/njt4n/share/Data/Tumor/BRATS-1/Images/BRATS_HG0001_T2.nii.gz
          -n T2
          -c clusterCenterT2
          -t /home/njt4n/share/Data/Kirby/T2_template.nii.gz
          -o /home/njt4n/share/Data/Tumor/BRATS-1/BRATS_HG0001/BRATS_HG0001

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
                                                (for 5 classes: csf, gm, wm, edema, and tumor)
     -b:  number of clusters                    If -c is specified, this option is not needed.
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

    Using applyTumorSegmentationModel with the following arguments:
      image dimension         = ${DIMENSION}
      anatomical image        = ${ANATOMICAL_IMAGES[@]}
      symmetric templates     = ${SYMMETRIC_TEMPLATES[@]}
      cluster centers         = ${CLUSTER_CENTERS[@]}
      image names             = ${IMAGE_NAMES[@]}
      radii                   = ${RADII[@]}
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

MODEL=""
SEGMENTATION_PRIOR=""

MASK_IMAGE=""

RADII=()
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
  while getopts "a:b:c:d:f:h:l:m:n:o:p:r:s:t:x:" OPT
    do
      case $OPT in
          a) #anatomical image
       ANATOMICAL_IMAGES[${#ANATOMICAL_IMAGES[@]}]=$OPTARG
       ;;
          b)
       NUMBER_OF_LABELS=$OPTARG
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
          h) #help
       Usage >&2
       exit 0
       ;;
          m)
       MODEL=$OPTARG
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
       RADII[${#RADII[@]}]=$OPTARG
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

if [[ ! -f ${MODEL} ]];
  then
    echo "The specified model \"${MODEL}\" does not exist."
    exit 1
  fi

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

if [[ ${#CLUSTER_CENTERS[@]} -gt 0 ]];
  then
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
  fi

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
# Create feature images
#
################################################################################

COMMAND_LINE="-d ${DIMENSION} -x ${MASK_IMAGE} -o ${OUTPUT_PREFIX}"
for (( i = 0; i < ${#ANATOMICAL_IMAGES[@]}; i++ ))
  do
    COMMAND_LINE="${COMMAND_LINE} -a ${ANATOMICAL_IMAGES[$i]}"
    COMMAND_LINE="${COMMAND_LINE} -t ${SYMMETRIC_TEMPLATES[$i]}"
    if [[ ${#CLUSTER_CENTERS[@]} -gt 0 ]];
      then
        COMMAND_LINE="${COMMAND_LINE} -c ${CLUSTER_CENTERS[$i]}"
      fi
    COMMAND_LINE="${COMMAND_LINE} -n ${IMAGE_NAMES[$i]}"
  done

if [[ ${#CLUSTER_CENTERS[@]} -eq 0 ]];
  then
    COMMAND_LINE="${COMMAND_LINE} -b ${NUMBER_OF_LABELS}"
  fi



for (( i = 0; i < ${#RADII[@]}; i++ ))
  do
    COMMAND_LINE="${COMMAND_LINE} -r ${RADII[$i]}"
  done

for (( i = 0; i < ${#DIFFERENCE_PAIRS[@]}; i++ ))
  do
    COMMAND_LINE="${COMMAND_LINE} -f ${DIFFERENCE_PAIRS[$i]}"
  done

COMMAND_LINE="${COMMAND_LINE} -s ${SMOOTHING_SIGMA} -l ${CORE_LABEL}"

if [[ ! -z "${SEGMENTATION_PRIOR}" ]];
  then
    COMMAND_LINE="${COMMAND_LINE} -p ${SEGMENTATION_PRIOR}"
  fi

sh ${SCRIPTSPATH}/createFeatureImages.sh ${COMMAND_LINE}

################################################################################
#
# Create csv file
#
################################################################################
CSV_FILE=${OUTPUT_PREFIX}FeatureImageList.csv

logCmd Rscript ${SCRIPTSPATH}/createCSVFileFromModel.R ${MODEL} ${MASK_IMAGE} ${OUTPUT_PREFIX} ${CSV_FILE}

################################################################################
#
# Apply the model
#
################################################################################

logCmd Rscript ${SCRIPTSPATH}/applyModel.R ${DIMENSION} ${MODEL} ${CSV_FILE} ${OUTPUT_PREFIX}RF_POSTERIORS 1

################################################################################
#
# Need to refine the model
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
echo " Done with applying model"
echo " Script executed in $time_elapsed seconds"
echo " $(( time_elapsed / 3600 ))h $(( time_elapsed %3600 / 60 ))m $(( time_elapsed % 60 ))s"
echo "--------------------------------------------------------------------------------------"

