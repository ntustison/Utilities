#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=18:00:00
#PBS -m ea
#PBS -V
#PBS -M ntustison@gmail.com

# variables that change between calls
T1=${HOME}/share/Data/dp04_001_x_t1.nii.gz
REGISTRATION_TEMPLATE=""
OUTPUT_DATA_DIR=${HOME}/share/Data/TestOutput2/

# ants brain processing script

ABP=${HOME}/Pkg/scripts/abp.sh

# Define variable for local storage on compute nodes associated with the job
LS="/jobtmp/pbstmp.$PBS_JOBID"

# Copy executables and data from your home directory to local storage on the
# master compute node

cd $LS
/bin/cp -R ${HOME}/Pkg/ANTS/bin/ .
/bin/cp $ABP .

export ANTSPATH="${LS}/bin/"


if [[ ! -d "$OUTPUT_DATA_DIR" ]];
  then
    mkdir -p $OUTPUT_DATA_DIR
  fi

T1_DATA_DIR=${LS}/T1Data/
EXTRACTION_DATA_DIR=${LS}/ExtractionData/
SEGMENTATION_DATA_DIR=${LS}/SegmentationData/
LOCAL_OUTPUT_DATA_DIR=${LS}/Output/

mkdir $T1_DATA_DIR
mkdir $EXTRACTION_DATA_DIR
mkdir $SEGMENTATION_DATA_DIR
mkdir $LOCAL_OUTPUT_DATA_DIR

/bin/cp $T1 $T1_DATA_DIR

FILES=( `ls ${T1_DATA_DIR}/*.nii.gz` )
T1_IMAGE=${FILES[0]}

HOME_DATA_DIR=${HOME}/share/Data/Public/

/bin/cp ${HOME_DATA_DIR}/LONI/ANTS_templates/Template1to20/T_template.nii.gz $EXTRACTION_DATA_DIR
/bin/cp ${HOME_DATA_DIR}/LONI/ANTS_templates/Template1to20/T_templateProbabilityMask.nii.gz $EXTRACTION_DATA_DIR
/bin/cp ${HOME_DATA_DIR}/NIREP/ANTS_templates/Template1to8/T_template.nii.gz $SEGMENTATION_DATA_DIR
/bin/cp ${HOME_DATA_DIR}/NIREP/ANTS_templates/Template1to8/3TissueProbabilityImages/smooth*nii.gz $SEGMENTATION_DATA_DIR

# Run program


sh ./abp.sh -d 3 \
          -a ${T1_IMAGE} \
          -o ${LOCAL_OUTPUT_DATA_DIR}/apb_ \
          -w 3 \
          -g 2 \
          -k 0 \
          -e ${EXTRACTION_DATA_DIR}/T_template.nii.gz \
          -m ${EXTRACTION_DATA_DIR}/T_templateProbabilityMask.nii.gz \
          -l ${SEGMENTATION_DATA_DIR}/T_template.nii.gz \
          -p ${SEGMENTATION_DATA_DIR}/smoothprior\%d.nii.gz \
          -t $REGISTRATION_TEMPLATE

# cp data to it's final location

DIRS=( `ls /gluster/pbstemp/*${PBS_JOBID}*` )
TMP_OUTPUT_DIR=${DIRS[0]}

echo "Outputting to $TMP_OUTPUT_DIR"

if [[ -d "$TMP_OUTPUT_DIR" ]];
  then
    /bin/cp ${TMP_OUTPUT_DIR}/Output/apb_* $OUTPUT_DATA_DIR
  else
    echo "$TMP_OUTPUT_DIR does not exist."
  fi
