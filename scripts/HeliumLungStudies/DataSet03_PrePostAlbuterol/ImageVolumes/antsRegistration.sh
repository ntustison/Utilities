# request Bourne shell as shell for job
#$ -S /bin/sh

PREDIRECTORY=$1
POSTDIRECTORY=$2
PRESERIES=$3
POSTSERIES=$4

sleep 1

perl /home/tustison/Data/HeliumLungStudies/PrePostAlbuterol/ImageVolumes/antsRegistration.pl \
  $PREDIRECTORY $POSTDIRECTORY $PRESERIES $POSTSERIES
