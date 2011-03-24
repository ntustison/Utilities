# request Bourne shell as shell for job
#$ -S /bin/sh

PREDIRECTORY=$1
POSTDIRECTORY=$2

perl /mnt/data2/PUBLIC/Data/Input/HeliumLungStudies/ExerciseMethacholineChallenge/ImageVolumes/pointSetRegistration.pl \
  $PREDIRECTORY $POSTDIRECTORY
