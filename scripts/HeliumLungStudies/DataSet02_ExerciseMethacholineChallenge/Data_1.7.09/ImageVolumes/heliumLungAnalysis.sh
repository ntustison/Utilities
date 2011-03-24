# request Bourne shell as shell for job
#$ -S /bin/sh

OUTPUTDIR="$1"
HELIUMIMAGE="$2"

PID="021209"

UTILITYDIR="/mnt/data1/tustison/Utilities/bin/"


TMPDIR=${OUTPUTDIR}

SEGMENTATION=${OUTPUTDIR}/segmentation_${PID}.nii.gz

${UTILITYDIR}/SegmentHeliumLungs 3 $HELIUMIMAGE $SEGMENTATION

##
# Convert segmentation image to point set for later registration
##
${UTILITYDIR}/ConvertLabeledPointSet 3 \
  ${SEGMENTATION} ${OUTPUTDIR}/points.vtk 1 1


echo "Both lungs analysis"

PREFIX=both_

MASK=${OUTPUTDIR}/segmentation_both_${PID}.nii.gz

${UTILITYDIR}/ThresholdImage 3 ${SEGMENTATION} ${MASK} 2 3 1 0
${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${HELIUMIMAGE} \
  ${MASK} 1 100 > ${OUTPUTDIR}/${PREFIX}statistics_${PID}.txt
${UTILITYDIR}/GenerateCooccurrenceMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 1 > ${OUTPUTDIR}/${PREFIX}cooccurrence_${PID}.txt
${UTILITYDIR}/GenerateRunLengthMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 1 > ${OUTPUTDIR}/${PREFIX}runlength_${PID}.txt
${UTILITYDIR}/CalculateVolumeFromBinaryImage 3 \
  ${MASK} 1 > ${OUTPUTDIR}/${PREFIX}volume_${PID}.txt

##
# Calculate stochastic fractal dimension image
##

FRACTALIMAGE=${TMPDIR}/fractalimage_${PID}.nii.gz

${UTILITYDIR}/GenerateFractalImage 3 ${HELIUMIMAGE} ${FRACTALIMAGE} 1 \
  ${MASK} 1
${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${FRACTALIMAGE} \
  ${MASK} 1 100 > ${OUTPUTDIR}/${PREFIX}statistics_${PID}.txt


##
# whole lung analysis
##

echo "Whole lung analysis"

MASK=${SEGMENTATION}

PREFIX=whole_2_

${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${HELIUMIMAGE} \
  ${MASK} 2 100 > ${OUTPUTDIR}/${PREFIX}statistics_${PID}.txt
${UTILITYDIR}/GenerateCooccurrenceMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 2 > ${OUTPUTDIR}/${PREFIX}cooccurrence_${PID}.txt
${UTILITYDIR}/GenerateRunLengthMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 2 > ${OUTPUTDIR}/${PREFIX}runlength_${PID}.txt
${UTILITYDIR}/CalculateVolumeFromBinaryImage 3 \
  ${MASK} 2 > ${OUTPUTDIR}/${PREFIX}volume_${PID}.txt
${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${FRACTALIMAGE} \
  ${MASK} 2 100 > ${OUTPUTDIR}/${PREFIX}fractal_statistics_${PID}.txt

PREFIX=whole_3_

${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${HELIUMIMAGE} \
  ${MASK} 3 100 > ${OUTPUTDIR}/${PREFIX}statistics_${PID}.txt
${UTILITYDIR}/GenerateCooccurrenceMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 3 > ${OUTPUTDIR}/${PREFIX}cooccurrence_${PID}.txt
${UTILITYDIR}/GenerateRunLengthMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 3 > ${OUTPUTDIR}${PREFIX}runlength_${PID}.txt
${UTILITYDIR}/CalculateVolumeFromBinaryImage 3 \
  ${MASK} 3 > ${OUTPUTDIR}/${PREFIX}volume_${PID}.txt
${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${FRACTALIMAGE} \
  ${MASK} 3 100 > ${OUTPUTDIR}/${PREFIX}fractal_statistics_${PID}.txt

##
# Inner rind / outer core analysis
##

echo "Inner rind / outer core analysis"


RADIUS=5

## inner core 2

PREFIX=inner_core_2

MASK=${TMPDIR}/${PREFIX}_${PID}.nii.gz

${UTILITYDIR}/BinaryMorphology 3 ${SEGMENTATION} ${MASK} \
  1 ${RADIUS} 1 2
${UTILITYDIR}/ThresholdImage 3 ${MASK} ${MASK} 2 2 2 0

${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${HELIUMIMAGE} \
  ${MASK} 2 100 > ${OUTPUTDIR}/${PREFIX}_statistics_${PID}.txt
${UTILITYDIR}/GenerateCooccurrenceMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 2 > ${OUTPUTDIR}/${PREFIX}_cooccurrence_${PID}.txt
${UTILITYDIR}/GenerateRunLengthMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 2 > ${OUTPUTDIR}/${PREFIX}_runlength_${PID}.txt
${UTILITYDIR}/CalculateVolumeFromBinaryImage 3 \
  ${MASK} 2 > ${OUTPUTDIR}/${PREFIX}_volume_${PID}.txt
${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${FRACTALIMAGE} \
  ${MASK} 2 100 > ${OUTPUTDIR}/${PREFIX}_fractal_statistics_${PID}.txt


## outer rind 2


PREFIX=outer_rind_2

MASK=${TMPDIR}/${PREFIX}_${PID}.nii.gz

${UTILITYDIR}/SubtractImage 3 ${SEGMENTATION} ${TMPDIR}/inner_core_2_${PID}.nii.gz ${MASK}
${UTILITYDIR}/ThresholdImage 3 ${MASK} ${MASK} 2 2 2 0

MASK=${TMPDIR}/${PREFIX}_${PID}.nii.gz

${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${HELIUMIMAGE} \
  ${MASK} 2 100 > ${OUTPUTDIR}/${PREFIX}_statistics_${PID}.txt
${UTILITYDIR}/GenerateCooccurrenceMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 2 > ${OUTPUTDIR}/${PREFIX}_cooccurrence_${PID}.txt
${UTILITYDIR}/GenerateRunLengthMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 2 > ${OUTPUTDIR}/${PREFIX}_runlength_${PID}.txt
${UTILITYDIR}/CalculateVolumeFromBinaryImage 3 \
  ${MASK} 2 > ${OUTPUTDIR}/${PREFIX}_volume_${PID}.txt
${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${FRACTALIMAGE} \
  ${MASK} 2 100 > ${OUTPUTDIR}/${PREFIX}_fractal_statistics_${PID}.txt

## inner core 3

PREFIX=inner_core_3

MASK=${TMPDIR}/${PREFIX}_${PID}.nii.gz

${UTILITYDIR}/BinaryMorphology 3 ${SEGMENTATION} ${MASK} \
  1 ${RADIUS} 1 3
${UTILITYDIR}/ThresholdImage 3 ${MASK} ${MASK} 3 3 3 0

${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${HELIUMIMAGE} \
  ${MASK} 3 100 > ${OUTPUTDIR}/${PREFIX}_statistics_${PID}.txt
${UTILITYDIR}/GenerateCooccurrenceMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 3 > ${OUTPUTDIR}/${PREFIX}_cooccurrence_${PID}.txt
${UTILITYDIR}/GenerateRunLengthMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 3 > ${OUTPUTDIR}/${PREFIX}_runlength_${PID}.txt
${UTILITYDIR}/CalculateVolumeFromBinaryImage 3 \
  ${MASK} 3 > ${OUTPUTDIR}/${PREFIX}_volume_${PID}.txt
${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${FRACTALIMAGE} \
  ${MASK} 3 100 > ${OUTPUTDIR}/${PREFIX}_fractal_statistics_${PID}.txt


## outer rind 3

PREFIX=outer_rind_3

MASK=${TMPDIR}/${PREFIX}_${PID}.nii.gz

${UTILITYDIR}/SubtractImage 3 ${SEGMENTATION} ${TMPDIR}/inner_core_3_${PID}.nii.gz \
  ${MASK}
${UTILITYDIR}/ThresholdImage 3 ${MASK} ${MASK} 3 3 3 0

MASK=${TMPDIR}/${PREFIX}_${PID}.nii.gz

${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${HELIUMIMAGE} \
  ${MASK} 3 100 > ${OUTPUTDIR}/${PREFIX}_statistics_${PID}.txt
${UTILITYDIR}/GenerateCooccurrenceMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 3 > ${OUTPUTDIR}/${PREFIX}_cooccurrence_${PID}.txt
${UTILITYDIR}/GenerateRunLengthMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 3 > ${OUTPUTDIR}/${PREFIX}_runlength_${PID}.txt
${UTILITYDIR}/CalculateVolumeFromBinaryImage 3 \
  ${MASK} 3 > ${OUTPUTDIR}/${PREFIX}_volume_${PID}.txt
${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${FRACTALIMAGE} \
  ${MASK} 3 100 > ${OUTPUTDIR}/${PREFIX}_fractal_statistics_${PID}.txt

${UTILITYDIR}/ThresholdImage 3 ${TMPDIR}/${PREFIX}.nii.gz ${TMPDIR}/${PREFIX}.nii.gz \
  3 3 3 0

## inner rind both

PREFIX=inner_core_both
MASK=${TMPDIR}/${PREFIX}_${PID}.nii.gz

${UTILITYDIR}/AddImage 3 ${TMPDIR}/inner_core_2_${PID}.nii.gz ${TMPDIR}/inner_core_3_${PID}.nii.gz \
  ${MASK}
${UTILITYDIR}/ThresholdImage 3 ${MASK} ${MASK} 2 3 1 0

${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${HELIUMIMAGE} \
  ${MASK} 1 100 > ${OUTPUTDIR}/${PREFIX}_statistics_${PID}.txt
${UTILITYDIR}/GenerateCooccurrenceMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 1 > ${OUTPUTDIR}/${PREFIX}_cooccurrence_${PID}.txt
${UTILITYDIR}/GenerateRunLengthMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 1 > ${OUTPUTDIR}/${PREFIX}_runlength_${PID}.txt
${UTILITYDIR}/CalculateVolumeFromBinaryImage 3 \
  ${MASK} 1 > ${OUTPUTDIR}/${PREFIX}_volume_${PID}.txt
${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${FRACTALIMAGE} \
  ${MASK} 1 100 > ${OUTPUTDIR}/${PREFIX}_fractal_statistics_${PID}.txt

## outer rind both

PREFIX=outer_rind_both
MASK=${TMPDIR}/${PREFIX}_${PID}.nii.gz

${UTILITYDIR}/AddImage 3 ${TMPDIR}/outer_rind_2_${PID}.nii.gz ${TMPDIR}/outer_rind_3_${PID}.nii.gz \
  ${MASK}
${UTILITYDIR}/ThresholdImage 3 ${MASK} ${MASK} 2 3 1 0

${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${HELIUMIMAGE} \
  ${MASK} 1 100 > ${OUTPUTDIR}/${PREFIX}_statistics_${PID}.txt
${UTILITYDIR}/GenerateCooccurrenceMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 1 > ${OUTPUTDIR}/${PREFIX}_cooccurrence_${PID}.txt
${UTILITYDIR}/GenerateRunLengthMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 1 > ${OUTPUTDIR}/${PREFIX}_runlength_${PID}.txt
${UTILITYDIR}/CalculateVolumeFromBinaryImage 3 \
  ${MASK} 1 > ${OUTPUTDIR}/${PREFIX}_volume_${PID}.txt
${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${FRACTALIMAGE} \
  ${MASK} 1 100 > ${OUTPUTDIR}/${PREFIX}_fractal_statistics_${PID}.txt

##
# Lung division analysis
##

echo "lung division analysis"


MASK=${TMPDIR}/division_${PID}.nii.gz

${UTILITYDIR}/DivideLungs ${SEGMENTATION} ${MASK} 2 2

## region 1

PREFIX=region_1

${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${HELIUMIMAGE} \
  ${MASK} 1 100 > ${OUTPUTDIR}/${PREFIX}_statistics_${PID}.txt
${UTILITYDIR}/GenerateCooccurrenceMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 1 > ${OUTPUTDIR}/${PREFIX}_cooccurrence_${PID}.txt
${UTILITYDIR}/GenerateRunLengthMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 1 > ${OUTPUTDIR}/${PREFIX}_runlength_${PID}.txt
${UTILITYDIR}/CalculateVolumeFromBinaryImage 3 \
  ${MASK} 1 > ${OUTPUTDIR}/${PREFIX}_volume_${PID}.txt
${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${FRACTALIMAGE} \
  ${MASK} 1 100 > ${OUTPUTDIR}/${PREFIX}_fractal_statistics_${PID}.txt

## region 2

PREFIX=region_2

${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${HELIUMIMAGE} \
  ${MASK} 2 100 > ${OUTPUTDIR}/${PREFIX}_statistics_${PID}.txt
${UTILITYDIR}/GenerateCooccurrenceMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 2 > ${OUTPUTDIR}/${PREFIX}_cooccurrence_${PID}.txt
${UTILITYDIR}/GenerateRunLengthMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 2 > ${OUTPUTDIR}/${PREFIX}_runlength_${PID}.txt
${UTILITYDIR}/CalculateVolumeFromBinaryImage 3 \
  ${MASK} 2 > ${OUTPUTDIR}/${PREFIX}_volume_${PID}.txt
${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${FRACTALIMAGE} \
  ${MASK} 2 100 > ${OUTPUTDIR}/${PREFIX}_fractal_statistics_${PID}.txt

## region 3

PREFIX=region_3

${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${HELIUMIMAGE} \
  ${MASK} 3 100 > ${OUTPUTDIR}/${PREFIX}_statistics_${PID}.txt
${UTILITYDIR}/GenerateCooccurrenceMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 3 > ${OUTPUTDIR}/${PREFIX}_cooccurrence_${PID}.txt
${UTILITYDIR}/GenerateRunLengthMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 3 > ${OUTPUTDIR}/${PREFIX}_runlength_${PID}.txt
${UTILITYDIR}/CalculateVolumeFromBinaryImage 3 \
  ${MASK} 3 > ${OUTPUTDIR}/${PREFIX}_volume_${PID}.txt
${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${FRACTALIMAGE} \
  ${MASK} 3 100 > ${OUTPUTDIR}/${PREFIX}_fractal_statistics_${PID}.txt

## region 4

PREFIX=region_4

${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${HELIUMIMAGE} \
  ${MASK} 4 100 > ${OUTPUTDIR}/${PREFIX}_statistics_${PID}.txt
${UTILITYDIR}/GenerateCooccurrenceMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 4 > ${OUTPUTDIR}/${PREFIX}_cooccurrence_${PID}.txt
${UTILITYDIR}/GenerateRunLengthMeasures 3 ${HELIUMIMAGE} \
  ${MASK} 4 > ${OUTPUTDIR}/${PREFIX}_runlength_${PID}.txt
${UTILITYDIR}/CalculateVolumeFromBinaryImage 3 \
  ${MASK} 4 > ${OUTPUTDIR}/${PREFIX}_volume_${PID}.txt
${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${FRACTALIMAGE} \
  ${MASK} 4 100 > ${OUTPUTDIR}/${PREFIX}_fractal_statistics_${PID}.txt

##
# Ventilation analysis
##

${UTILITYDIR}/BinaryMorphology 3 ${SEGMENTATION} ${OUTPUTDIR}/segmentation_eroded_${PID}.nii.gz \
  1 2 1 2
${UTILITYDIR}/BinaryMorphology 3 ${OUTPUTDIR}/segmentation_eroded_${PID}.nii.gz \
  ${OUTPUTDIR}/segmentation_eroded_${PID}.nii.gz 1 2 1 3

${UTILITYDIR}/OtsuThresholdImage 3 ${HELIUMIMAGE} ${OUTPUTDIR}/segmentation_otsu_2_${PID}.nii.gz \
  1 200 ${OUTPUTDIR}/segmentation_eroded_${PID}.nii.gz 2
${UTILITYDIR}/OtsuThresholdImage 3 ${HELIUMIMAGE} ${OUTPUTDIR}/segmentation_otsu_3_${PID}.nii.gz \
  1 200 ${OUTPUTDIR}/segmentation_eroded_${PID}.nii.gz 3

${UTILITYDIR}/CalculateVolumeFromBinaryImage 3 \
  ${OUTPUTDIR}/segmentation_otsu_2_${PID}.nii.gz 2 \
  ${OUTPUTDIR}/segmentation_eroded_${PID}.nii.gz 2 \
  > ${OUTPUTDIR}/otsu_wellventilated_volume_2_${PID}.txt
${UTILITYDIR}/CalculateVolumeFromBinaryImage 3 \
  ${OUTPUTDIR}/segmentation_otsu_3_${PID}.nii.gz 2 \
  ${OUTPUTDIR}/segmentation_eroded_${PID}.nii.gz 3 \
  > ${OUTPUTDIR}/otsu_wellventilated_volume_3_${PID}.txt

${UTILITYDIR}/BinaryMorphology 3 ${OUTPUTDIR}/segmentation_both_${PID}.nii.gz \
  ${OUTPUTDIR}/segmentation_eroded_${PID}.nii.gz 1 2 1 1
${UTILITYDIR}/OtsuThresholdImage 3 ${HELIUMIMAGE} ${OUTPUTDIR}/segmentation_otsu_both_${PID}.nii.gz \
  1 200 ${OUTPUTDIR}/segmentation_eroded_${PID}.nii.gz 1
${UTILITYDIR}/CalculateVolumeFromBinaryImage 3 \
  ${OUTPUTDIR}/segmentation_otsu_both_${PID}.nii.gz 2 \
  ${OUTPUTDIR}/segmentation_eroded_${PID}.nii.gz 1 \
  > ${OUTPUTDIR}/otsu_wellventilated_volume_both_${PID}.txt


##
# Noise statistics
##

RADIUS=5
PREFIX=noise
MASK=${TMPDIR}/${PREFIX}_mask_${PID}.nii.gz

${UTILITYDIR}/BinaryMorphology 3 ${TMPDIR}/segmentation_both_${PID}.nii.gz \
  ${MASK} 0 ${RADIUS} 0 1
${UTILITYDIR}/SubtractImage 3  ${MASK} ${TMPDIR}/segmentation_both_${PID}.nii.gz ${MASK}
  
${UTILITYDIR}/CalculateFirstOrderStatisticsFromImage 3 ${HELIUMIMAGE} \
  ${MASK} 1 100 > ${OUTPUTDIR}/${PREFIX}_statistics_${PID}.txt



