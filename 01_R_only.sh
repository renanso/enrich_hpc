#!/bin/bash
#SBATCH --job-name=enrich	# Job name
#SBATCH --nodes=1	        # N of nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem="512G"			# Memory per node; by default using M as unit
#SBATCH --time=48:00:00     # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=ALL        #export user’s explicit environment variables to compute node
#SBATCH --output=%x_%j.out		# Standard output log
#SBATCH --error=%x_%j.err		# Standard error log
#SBATCH --partition=normal

#running the enrich pipeline on R

eval "$(conda shell.bash hook)"
conda activate r_env

#R CMD BATCH 03_gc_filter_v1.R
R CMD BATCH 03_enrich_v6.3.R
