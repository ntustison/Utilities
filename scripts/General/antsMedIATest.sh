# request Bourne shell as shell for job#$ -S /bin/sh

ANTS=/home/tustison/ANTS/bin/ANTS

sleep 1

echo $ANTS $1 --output $2 --iterations $3 --transformation $4 --regularization $5 --metric $6
