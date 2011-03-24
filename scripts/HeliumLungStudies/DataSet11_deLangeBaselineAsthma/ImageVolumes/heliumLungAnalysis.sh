# request Bourne shell as shell for job
#$ -S /bin/sh

SCRIPTDIR="$1"
OUTPUTDIR="$2"
HELIUMIMAGE="$3"

perl ${SCRIPTDIR}/heliumLungAnalysis.pl $OUTPUTDIR $HELIUMIMAGE
