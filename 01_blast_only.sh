#!/bin/bash
#SBATCH --job-name=blast	# Job name
#SBATCH --nodes=1	        # N of nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem="256G"			# Memory per node; by default using M as unit
#SBATCH --time=48:00:00     # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=ALL        #export user’s explicit environment variables to compute node
#SBATCH --output=%x_%j.out		# Standard output log
#SBATCH --error=%x_%j.err		# Standard error log
#SBATCH --partition=plant


ulimit -S 200000000

#1-running blast on the candidates

# blast needs to be in the PATH:  export PATH=$PATH:$HOME/ncbi-blast-2.10.1+/bin
# there is an installation on the cluster

#export PATH=/cluster/software/blast-2.9.0/bin
export PATH=/cluster/home/rsouza/apps/blast/ncbi-blast-2.15.0+/bin
# A local database for your genome needs to be created with:

#makeblastdb -in genome.fasta -parse_seqids -dbtype nucl -out genome

# database needs to be in the PATH: export BLASTDB=$HOME/blastdb
export BLASTDB=/cluster/home/rsouza/apps/blast/blastdb/sintegrifolium

blastn -db sintegrifolium -query candidates2.fasta -out blast_results.txt -outfmt 6
