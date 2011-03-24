# !/bin/bash 
# 
#
NUMPARAMS=$#
#
#
if [ $NUMPARAMS -lt 2  ]
    then 
    echo " USAGE ::  "
    echo "  bash   $0   ImageToExtractBrainFrom  (PipeDream)ParamFile OPTIONAL-OUTPREFIX OPTIONAL-GroundTruthSegmentationImage " 
    echo " will output the file ::   \${OUTPREFIX}brain.nii.gz " 
    echo " make sure you know whether you are using the forward or inverse transform to propagate template information --- here we assume inverse. " 
    echo " be sure to set ANTSPATH environment variable -- Default is: $ANTSPATH " 
    echo
    echo the param file should define the following variable, at minimum, 
    echo  \$TEMPLATEBRAINMASK
    echo you should also define 
    echo \$TEMPLATE_SPHERE_INIT  -- initialization for g-0 segmentation 
    echo \$TEMPLATE_RELIABLE_BRAIN_EXTRACTION -- where you think the template based initialization is totally reliable 
    echo  e.g. in cerebellum / brainstem 
    echo if you do not define those, then we will estimate them from \$TEMPLATEBRAINMASK
    exit
fi
#
# initialization, here, is unbiased 
DIM=3
INPUT_FILE=$1
OUTPUTNAME=${INPUT_FILE%.*.*}
PIPEDREAMPARAMFILE=$2
echo $2

if [ -s $2 ] ; then 
    source $2
else
    echo " No Param File !!  Exiting!! "
    exit 1
fi
echo OK 
if [ $NUMPARAMS -gt 2 ] ; then
  OUTPUTNAME=${3}
fi  
if [ $NUMPARAMS -gt 3 ] ; then
  GROUNDTRUTH=${4}
  echo groundtruth $4 
fi  
if [ ! -s ${OUTPUTNAME}Affine.txt ] ; then 
  echo you need to compute ${OUTPUTNAME}Affine.txt and ${OUTPUTNAME}Warp.nii.gz !!! 
  exit 1 
fi 
# echo Mask is $TEMPLATEBRAINMASK

# for individual to template 
WARPCOMMANDA=" ${ANTSPATH}WarpImageMultiTransform $DIM  "
WARPCOMMANDB="  -R $TEMPLATEBRAINMASK  ${OUTPUTNAME}Warp.nii.gz ${OUTPUTNAME}Affine.txt  "

# for template to individual 
INVWARPCOMMANDA="  ${ANTSPATH}WarpImageMultiTransform $DIM "
INVWARPCOMMANDB=" -R $INPUT_FILE  -i ${OUTPUTNAME}Affine.txt   ${OUTPUTNAME}InverseWarp.nii.gz   "

# Code for extracting a brain from a brain image, using SyN and Atropos.
#
# Usage:   "sh extract_brain.sh INPUT_FILE OUT_PATH fix_headers"
# Example: "sh extract_brain.sh ./subjects/S01.nii.gz ./output/S01 0"
#
# NOTE: Make sure to set the ANTSPATH environment variable
#       either locally or on the cluster:
#         export ANTSPATH="/Users/arno/ANTS-build/"
#         export PATH="${PATH}:${ANTSPATH}"
#
#       Set fix_headers=1 below only if subject and template
#       headers do not correspond.
#
# Steps of the algorithm:
# 0. Run N3 bias field correction.
# 1. Deform a template image to a subject image that we want to extract,
#    and apply the transform to template segmentation class images.
# 2. Create a mask about the boundary of the deformed template.
# 3. Run Atropos with the masked, transformed segmented images as priors.
# 4. Extract the largest components of class 1, 2, and 3 to extract the brain.
# 5. Mask the resulting brain with a transformed template brain mask.
#
# Note: Atropos should be used in a data-specific manner:
# Build a local template, model relevant brain tissues, and apply to new data.
# Ideally, one would localize errors using training data and refine the template.
#
# (c) Brian Avants, arno klein (2010) arno@mindboggle.info
#
# Arguments:

###################################
# EDIT INFORMATION BELOW WITH CARE:
###################################
#export ANTSPATH="/Users/arno/ANTS-build/"
#export ANTSPATH="/data/BI/Toolbox/software/ants_mac_svn460/bin/"
#export PATH="${PATH}:${ANTSPATH}"

####################
# SET STAGES TO RUN:
####################
run_n3=1
run_register=0
run_transform=1
run_postpostprocess=1

FORMAT='.nii.gz' # format of the template and output files (NIfTI: .nii.gz)

####################
# parameters begin 
####################
  ITERS='30x90x20'
  ITERS='30x0x0'
 ## max distance of gm surface from wm surface 
  GMSURFDISTFROMWM=5
 ## dilation of initial template mask 
  BM_TEM_INIT_DILATION=2
 ## fast marching propagation - non critical 
  FMP=" 5000 2 " # preserve topology
  FMP=" 1 0 " # do not preserve topo, but smooth fit 
 ## atropos segmentation 
  ATROPOS_PARAMS1=" -c [1,0.001]  -m [0.,1] -i Kmeans[3] "
  ATROPOS_PARAMS2=" -c [1,0.001]  -m [0.,1] -i Kmeans[4] "
 ## basically dont add any csf to the gm mask ... erode it away 
 ## if this is zero, then you will add csf back to the mask
  CSFERODE=10
 ## amount to erode before getting largest components 
 ## this isolates non-brain regions 
  EROSION=2
  let DILATION=$EROSION+2
  MORPHCLOSE=5
# parameters end

####################
# Template files in: T1, segmented, and partial-volum/smoothed individual classes:
####################
BM_IN=${TEMPLATEBRAINMASK}
BM_SP_IN=${TEMPLATE_SPHERE_INIT}  # initialization for g-0 segmentation 
BM_TOPO_IN=${TEMPLATE_RELIABLE_BRAIN_EXTRACTION} # reliable brain extraction zone
if [[ ! -s $BM_SP_IN ]] ; then 
echo  create topo images from $BM_IN  
  if [ ${#BM_SP_IN} -lt 1 ] ; then 
    BM_SP_IN=${OUTPUTNAME}sphere.nii.gz
    BM_TOPO_IN=$BM_SP_IN 
  fi 
  ${ANTSPATH}ImageMath $DIM $BM_SP_IN MD $BM_IN 10 
  ${ANTSPATH}ImageMath $DIM $BM_SP_IN ME $BM_SP_IN 10 
  ${ANTSPATH}ImageMath $DIM $BM_SP_IN ME $BM_SP_IN 15 
  BM_TOPO_IN=$BM_SP_IN 
fi 
if [ ! -s $BM_TOPO_IN ] ; then 
  BM_TOPO_IN=$BM_SP_IN 
fi 

for x in $BM_IN $BM_SP_IN $BM_TOPO_IN ; do 
  if [ ! -s $x ] ; then 
    echo $x does not exist 
    exit 1
    else 
      echo $x 
  fi 
done 

 
BM_TEM=${OUTPUTNAME}brainmasktemplate${FORMAT}
BM_OUT=${OUTPUTNAME}brainmaskinit${FORMAT}
BM_TOPO_OUT=${OUTPUTNAME}brainmasktopoinit${FORMAT}
BM_TOPO_FIT=${OUTPUTNAME}brainmasktopofit${FORMAT}
BM_SP_OUT=${OUTPUTNAME}sphere${FORMAT}
BM_SEG_OUT=${OUTPUTNAME}seg${FORMAT}
BM_OUT_MASK=${OUTPUTNAME}topomask${FORMAT}
OUT_HEAD_N3=${OUTPUTNAME}head_n3${FORMAT}  
GM=${OUTPUTNAME}GM${FORMAT}        # Segmented grey matter   (post-postprocessing)
WM=${OUTPUTNAME}WM${FORMAT}        # Segmented white matter  (post-postprocessing)
CSF=${OUTPUTNAME}CSF${FORMAT}      # Segmented CSF           (post-postprocessing)
OUT_BRAIN=${OUTPUTNAME}brain${FORMAT}      # Segmented CSF           (post-postprocessing)

########
# Begin:
########


##################################
# 0. Run N3 bias field correction:
##################################
if  [[ ! -s $OUT_HEAD_N3  ]] ; then
  ${ANTSPATH}ImageMath $DIM $OUT_HEAD_N3  TruncateImageIntensity $INPUT_FILE 0.025 0.975 64
  ${ANTSPATH}N3BiasFieldCorrection $DIM  $OUT_HEAD_N3 $OUT_HEAD_N3 4
fi

########################################################################
# 1. Deform a template image to a subject image that we want to extract,
#    and apply the transform to template segmentation class images:
########################################################################


# 1b. Apply the registration transform to the segmented class images of the template:
#  we have two masks --- the initial brain mask and the topologically correct mask.
if [ $run_transform -eq 1 ]; then
  c_start=$INVWARPCOMMANDA
  c_end=$INVWARPCOMMANDB
  $c_start $BM_IN $BM_OUT $c_end  
  ${ANTSPATH}ThresholdImage $DIM $BM_OUT $BM_TEM 0.5 1 
  ${ANTSPATH}ImageMath $DIM $BM_OUT MD $BM_TEM  $BM_TEM_INIT_DILATION
  $c_start $BM_TOPO_IN $BM_TOPO_OUT $c_end  
  $c_start $BM_SP_IN $BM_SP_OUT $c_end 
  ${ANTSPATH}ImageMath $DIM $BM_SP_OUT GetLargestComponent $BM_SP_OUT
  ${ANTSPATH}SmoothImage $DIM $BM_SP_OUT 0.5 $BM_SP_OUT 
  ${ANTSPATH}ThresholdImage $DIM $BM_SP_OUT $BM_SP_OUT 0.25 1 
  ${ANTSPATH}ImageMath $DIM $BM_SP_OUT GetLargestComponent $BM_SP_OUT

fi

##############################
# 6. Run Atropos segmentation:
##############################
if [ $run_postpostprocess -eq 1 ]; then
  # 3-tissue segmentation -- can be replaced with N-tissue segmentation   
    if [ ! -s $BM_SEG_OUT ] ; then 
	${ANTSPATH}Atropos $DIM -a $OUT_HEAD_N3 -x $BM_OUT  -o $BM_SEG_OUT $ATROPOS_PARAMS1
	cp $BM_SEG_OUT ${OUTPUTNAME}temp.nii.gz 
	${ANTSPATH}ThresholdImage $DIM $BM_SEG_OUT ${OUTPUTNAME}temp1.nii.gz    3 3
	${ANTSPATH}MultiplyImages $DIM $BM_TOPO_OUT ${OUTPUTNAME}temp1.nii.gz ${OUTPUTNAME}temp1.nii.gz 
	${ANTSPATH}Atropos $DIM -a $OUT_HEAD_N3 -x $BM_OUT  -o $BM_SEG_OUT $ATROPOS_PARAMS2
	${ANTSPATH}ThresholdImage $DIM $BM_SEG_OUT ${OUTPUTNAME}temp2.nii.gz    3 3
	${ANTSPATH}MultiplyImages $DIM $BM_TOPO_OUT ${OUTPUTNAME}temp2.nii.gz ${OUTPUTNAME}temp2.nii.gz 
	brainvol1=` ${ANTSPATH}ImageMath $DIM   ${OUTPUTNAME}temp1.nii.gz  GetLargestComponent ${OUTPUTNAME}temp1.nii.gz | grep max |  cut -d ':' -f 2  `
	brainvol2=` ${ANTSPATH}ImageMath $DIM   ${OUTPUTNAME}temp2.nii.gz  GetLargestComponent ${OUTPUTNAME}temp2.nii.gz | grep max |  cut -d ':' -f 2  `
	rm ${OUTPUTNAME}temp1.nii.gz ${OUTPUTNAME}temp2.nii.gz 
	if [ $brainvol1 -gt $brainvol2 ] ; then 
	    cp  ${OUTPUTNAME}temp.nii.gz $BM_SEG_OUT
	fi 
    else
	echo already have seg
    fi
    
  ${ANTSPATH}ThresholdImage $DIM $BM_SEG_OUT $WM  3 3
  ${ANTSPATH}ThresholdImage $DIM $BM_SEG_OUT $GM  2 2
  ${ANTSPATH}ThresholdImage $DIM $BM_SEG_OUT $CSF 1 1

  # Get the largest connected component of WM and of GM and make sure each has no holes
  ${ANTSPATH}ImageMath $DIM $WM GetLargestComponent $WM
  ${ANTSPATH}ImageMath $DIM $GM GetLargestComponent $GM
  ${ANTSPATH}ImageMath $DIM $WM FillHoles $WM
  ${ANTSPATH}ImageMath $DIM ${OUTPUTNAME}temp.nii.gz D $WM 
  ${ANTSPATH}ThresholdImage $DIM ${OUTPUTNAME}temp.nii.gz ${OUTPUTNAME}temp.nii.gz 0 $GMSURFDISTFROMWM # prior for surf distance from wm 
  ${ANTSPATH}ImageMath $DIM $GM FillHoles $GM
  ${ANTSPATH}MultiplyImages $DIM $GM ${OUTPUTNAME}temp.nii.gz $GM 

 # Recreate the segmentation with these "repaired" images
  ${ANTSPATH}MultiplyImages $DIM $WM 3 $WM
  ${ANTSPATH}ImageMath $DIM ${OUTPUTNAME}temp.nii.gz ME $CSF $CSFERODE
  ${ANTSPATH}ImageMath $DIM $GM addtozero $GM ${OUTPUTNAME}temp.nii.gz
  rm -f ${OUTPUTNAME}temp.nii.gz 
  ${ANTSPATH}MultiplyImages $DIM $GM 2 $GM
  ${ANTSPATH}ImageMath $DIM $BM_SEG_OUT addtozero $WM $GM
  ${ANTSPATH}ImageMath $DIM $BM_SEG_OUT addtozero $BM_SEG_OUT $CSF 


  ###############################################
  # 7. Mask the brain with the segmentation mask:
  ###############################################
  ${ANTSPATH}ThresholdImage $DIM $BM_SEG_OUT $BM_OUT_MASK 2 3
  ${ANTSPATH}ImageMath $DIM $BM_OUT_MASK ME $BM_OUT_MASK $EROSION
  ${ANTSPATH}ImageMath $DIM $BM_OUT_MASK GetLargestComponent $BM_OUT_MASK
  ${ANTSPATH}ImageMath $DIM $BM_OUT_MASK MD $BM_OUT_MASK $DILATION
  ${ANTSPATH}ImageMath $DIM $BM_OUT_MASK FillHoles $BM_OUT_MASK
# add "reliable" parts of the brain ... 
  ${ANTSPATH}ImageMath $DIM   $BM_OUT_MASK  addtozero  $BM_OUT_MASK $BM_TOPO_OUT
  ${ANTSPATH}ImageMath $DIM $BM_OUT_MASK MD $BM_OUT_MASK $MORPHCLOSE
  ${ANTSPATH}ImageMath $DIM $BM_OUT_MASK ME $BM_OUT_MASK $MORPHCLOSE
# make the topology correct - NO!
# instead propagate toward template-based solution 
  ${ANTSPATH}SmoothImage $DIM  $BM_TEM  0.5 ${OUTPUTNAME}temp2.nii.gz 
  ${ANTSPATH}MultiplyImages $DIM $BM_TEM  ${OUTPUTNAME}temp2.nii.gz  ${OUTPUTNAME}temp2.nii.gz 
  ${ANTSPATH}ImageMath $DIM  $BM_OUT_MASK FastMarchingSegmentation  ${OUTPUTNAME}temp2.nii.gz  $BM_OUT_MASK   $FMP 
# ${ANTSPATH}CheckTopology  $BM_TOPO_FIT 

  ${ANTSPATH}MultiplyImages $DIM $BM_OUT_MASK $OUT_HEAD_N3 $OUT_BRAIN

fi

# evaluate 
# if [ -s $GROUNDTRUTH ] ; then 
#  ${ANTSPATH}ThresholdImage $DIM $BM_OUT_MASK  ${OUTPUTNAME}temp.nii.gz  1 99999999
#  ${ANTSPATH}ThresholdImage $DIM $GROUNDTRUTH ${OUTPUTNAME}temp2.nii.gz  1 99999999
#  ${ANTSPATH}LabelOverlapMeasures $DIM ${OUTPUTNAME}temp.nii.gz  ${OUTPUTNAME}temp2.nii.gz  > ${OUTPUTNAME}eval.txt 
#  ${ANTSPATH}LabelOverlapMeasures $DIM $BM_TEM ${OUTPUTNAME}temp2.nii.gz  >> ${OUTPUTNAME}eval.txt 
# fi 

# clean-up  --- but not $BM_OUT_MASK 
rm $BM_OUT $BM_TOPO_OUT $BM_TOPO_FIT $BM_SP_OUT $OUT_HEAD_N3 $GM $WM $CSF ${OUTPUTNAME}temp.nii.gz  ${OUTPUTNAME}temp2.nii.gz 


exit 0

