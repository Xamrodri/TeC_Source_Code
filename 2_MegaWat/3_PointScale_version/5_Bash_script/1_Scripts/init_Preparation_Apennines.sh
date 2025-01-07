#!/bin/bash

#SBATCH --job-name=Max_2024	# Job name  (max. 12345678 characters)
#SBATCH -o /nfs/scistore18/pelligrp/mrodrigu/1_TC/2_MegaWat/3_PointScale_version/5_Bash_script/2_Script_results/Results_Example_%j.out 			# Standard output
#SBATCH -e /nfs/scistore18/pelligrp/mrodrigu/1_TC/2_MegaWat/3_PointScale_version/5_Bash_script/2_Script_results/Results_Example_%j.err 			# Standard error output

### NODE (CENTOS) ###
#SBATCH --mail-user=mrodrigu@ista.ac.at 	# My email
#SBATCH --mail-type=SUBMIT,END,FAIL 	 	# Notification circumstances
#SBATCH --mem=65G 			# RAM request
#SBATCH --time=01-00:00 			# Run time DD-HH:MM
#SBATCH -c 30			# Number of cores
#SBATCH --constraint=matlab 
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV


module load matlab/R2024b
srun matlab -nodisplay -nodesktop -nosplash -r "clear all; delete(gcp('nocreate')); s=2; fnum=1; run('/nfs/scistore18/pelligrp/mrodrigu/1_TC/2_MegaWat/4_Preparation_files/Apennine_preparation_Achille_Method_PointScale_Forcing/Prepare_forcing_MultiPoint.m');exit;" | tail -n +11