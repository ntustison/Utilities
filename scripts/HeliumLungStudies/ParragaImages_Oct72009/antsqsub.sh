#!/bin/bash
#$ -S /bin/bash

DIR=/home/tustison/Data/HeliumLungStudies/ParragaImages_Oct72009/
ANTS=/home/tustison/ANTS/bin64/ANTS

$ANTS 3 -m PR[${DIR}C3proton_corrected_inverted.nii.gz,${DIR}C3helium_corrected.nii.gz,1,4] \
        -t SyN[0.5] \
        -i 50x50x10 \
        -o ${DIR}C3.nii.gz \
        -r Gauss[3.0,1.0]
