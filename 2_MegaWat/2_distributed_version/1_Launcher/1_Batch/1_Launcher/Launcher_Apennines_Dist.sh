#!/bin/bash

#SBATCH --job-name=Distributed_model	# Job name  (max. 12345678 characters)
#SBATCH -o /nfs/scistore18/pelligrp/mrodrigu/1_TC/3_Model_Source/2_MegaWat/2_distributed_version/1_Launcher/1_Batch/2_Outputs/HPC_output_%j.out 			# Standard output
#SBATCH -e /nfs/scistore18/pelligrp/mrodrigu/1_TC/3_Model_Source/2_MegaWat/2_distributed_version/1_Launcher/1_Batch/2_Outputs/HPC_output_%j.err 			# Standard error output

### NODE (CENTOS) ###
#SBATCH --mail-user=mrodrigu@ista.ac.at 	# My email
#SBATCH --mail-type=SUBMIT,END,FAIL 	 	# Notification circumstances
#SBATCH --mem=260G 			                # RAM request
#SBATCH --time=01-00:00 			        # Run time DD-HH:MM
#SBATCH -c 160			                    # Number of cores
#SBATCH --constraint=matlab 
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV


module load matlab/R2024b
srun matlab -nodisplay -nodesktop -nosplash -r "clear all; delete(gcp('nocreate')); run('/nfs/scistore18/pelligrp/mrodrigu/1_TC/3_Model_Source/2_MegaWat/2_distributed_version/1_Launcher/Launcher_Apennines_Distributed.m');exit;" | tail -n +11

