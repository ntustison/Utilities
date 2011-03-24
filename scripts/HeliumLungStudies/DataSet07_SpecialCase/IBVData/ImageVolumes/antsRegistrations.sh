# request Bourne shell as shell for job
#$ -S /bin/sh

# qsub antsPR.sh \
#  /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/CT_03_23_09/CT_592.nii.gz \
#  /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/He3_1125_03_23_09/he3_aligned.nii.gz \
#  /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/RegistrationResults/he3ToCT_2009

# qsub antsPR.sh \
#  /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/CT_03_23_09/CT_592.nii.gz \
#  /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/He3_1125_03_23_09/proton_aligned.nii.gz \
#  /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/RegistrationResults/protonToCT_2009

# qsub antsPR.sh \
#  /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/CT_08_19_08/CT_707.nii.gz \
#  /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/He3_1108_09_09_2008/he3_aligned.nii.gz \
#  /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/RegistrationResults/he3ToCT_2008

# qsub antsPR.sh \
#  /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/CT_08_19_08/CT_707.nii.gz \
#  /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/He3_1108_09_09_2008/proton_aligned.nii.gz \
#  /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/RegistrationResults/protonToCT_2008

# qsub antsPR.sh \
#   /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/CT_08_19_08/CT_707.nii.gz \
#   /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/CT_03_23_09/CT_592.nii.gz \
#   /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/RegistrationResults/CT592ToCT707
 
 qsub antsPR.sh \
   /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/He3_1108_09_09_2008/proton.nii.gz \
   /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/He3_1125_03_23_09/proton.nii.gz \
   /home/tustison/Data/HeliumLungStudies/IBVData/ImageVolumes/RegistrationResults/proton2009ToProton2008
