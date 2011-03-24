# request Bourne shell as shell for job
#$ -S /bin/sh

ANTS=/home/tustison/ANTS/bin64/ANTS
fixedImage=$1
movingImage=$2
output=$3

$ANTS 3 -m PR[${fixedImage},${movingImage},1,4] \
        -t SyN[0.5] \
        -i 100x100x50x10 \
        -r Gauss[3,1] \
        -o $output 
