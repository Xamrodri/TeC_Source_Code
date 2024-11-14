#!/bin/bash

#SBATCH --job-name=Kyzylsu_2021_2022	# Job name  (max. 12345678 characters)
#SBATCH -o /nfs/scistore18/pelligrp/ajoubert/Jobs/logs/TC_distributed/TC_Dist_KYZ_%j.out 			# Standard output
#SBATCH -e /nfs/scistore18/pelligrp/ajoubert/Jobs/logs/TC_distributed/TC_Dist_KYZ_%j.err 			# Standard error output

### NODE (CENTOS) ###
#SBATCH --mail-user=achille.jouberton@ist.ac.at 	# My email, put yours instead !
#SBATCH --mail-type=BEGIN,END,FAIL 		# Notification circumstances
#SBATCH --mem=65G 				# RAM request
#SBATCH --time=10-00:00 			# Run time DD-HH:MM
#SBATCH -c 30			# Number of cores
#SBATCH --constraint=matlab 
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

module load matlab/R2023b
srun matlab -nodisplay -nodesktop -nosplash -r "clear all; delete(gcp('nocreate')); s=2; fnum=1;sim_comment = 'Test_studycase_HPC'; date_start = '01-Oct-2021 00:00:00'; date_end = '30-Sep-2022 23:00:00'; run('/nfs/scistore18/pelligrp/ajoubert/TeC_Source_Code/Case_study/Kyzylsu_distributed/Launcher_Kyzylsu_Distributed.m');exit;" | tail -n +11