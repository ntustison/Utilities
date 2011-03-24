# request Bourne shell as shell for job
#$ -S /bin/sh

ANTS=/home/tustison/ANTS/bin64/ANTS
fixedImage=/home/tustison/Data/HeliumLungStudies/SevereAsthma/ImageVolumes/VIDAMasks/SR6_resampled.nii.gz
movingImage=/home/tustison/Data/HeliumLungStudies/SevereAsthma/ImageVolumes/HeliumProtonAndVentilation/Series20070914_100413_656000/Series20070914_100413_656000_flipped.nii.gz

output=/home/tustison/Data/HeliumLungStudies/SevereAsthma/ImageVolumes/RegistrationResults/antsPR4_he3_

$ANTS 3 -m PR[${fixedImage},${movingImage},1,4] \
        -t SyN[0.5] \
        -i 100x100x50x10 \
        -r Gauss[3,1] \
        -o $output 
