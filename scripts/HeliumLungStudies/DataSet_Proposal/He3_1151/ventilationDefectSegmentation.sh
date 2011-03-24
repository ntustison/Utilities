# request Bourne shell as shell for job
#$ -S /bin/sh

SCRIPTDIR="$1"
OUTPUTDIR="$2"

perl ${SCRIPTDIR}/ventilationDefectSegmentation.pl $OUTPUTDIR
