# request Bourne shell as shell for job
#$ -S /bin/sh


for (( i = $1; i <= $2; i++ ))
  do
    qdel $i
  done 
