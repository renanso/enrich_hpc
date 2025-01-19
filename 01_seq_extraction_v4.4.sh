#!/bin/bash
#SBATCH --job-name=enrich	# Job name
#SBATCH --nodes=1	        # N of nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem="256G"			# Memory per node; by default using M as unit
#SBATCH --time=48:00:00     # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=ALL        #export user’s explicit environment variables to compute node
#SBATCH --output=%x_%j.out		# Standard output log
#SBATCH --error=%x_%j.err		# Standard error log
#SBATCH --partition=normal

# Load module
ml cluster/bedtools/2.28.0

#directories and inputs
#reference genome
REF=/cluster/home/rsouza/silphium/Silphium_integrifolium_var_Bad_Astra_V1_release/Silphium_integrifolium_var_Bad_Astra/sequences/Silphium_integrifolium_var_Bad_Astra.mainGenome.fasta
#Input VCF
VCF=/cluster/home/rsouza/enrichment/snpeff_files/vcfs/final_srt_merge_smiss001_dp5_maf001_77m_ann.vcf.gz
#Annotations
#ANN=/cluster/home/rsouza/enrichment/snpeff_files/vcfs/sintegrifolium.ann.tsv

# 1- extract the information from the annotated vcf

#selecting targets

mkdir tmp

#ignoring intergenic candidates
#cat ${VCF} | grep -v 'intergenic' | grep -E 'HIGH|MODERATE|LOW|MODIFIER' | sed 's/|/\t/g' | awk -v OFS='\t' '{ print $14, $1, $2, $2-250, $2+250, $3, "REF="$4,"ALT="$5, $9, $10}' | sed 's/CHR_START-//g' |  sed 's/-CHR_END//g' | sed 's/\.t[1-9]//g' > tmp/targets.txt

#only high impact
#cat ${VCF} | grep 'HIGH' | sed 's/|/\t/g' | awk -v OFS='\t' '{ print $14, $1, $2, $2-150, $2+150, $3, "REF="$4,"ALT="$5, $9, $10}' | sed 's/CHR_START-//g' |  sed 's/-CHR_END//g' | sed 's/\.t[1-9]//g' > tmp/targets.txt

# If SNP ID is not present for some reason

#cat ${VCF} | grep 'HIGH' | sed 's/|/\t/g' | awk -v OFS='\t' '{ print $14, $1, $2, $2-150, $2+150, $1"_"$2, "REF="$4,"ALT="$5, $9, $10}' | sed 's/CHR_START-//g' |  sed 's/-CHR_END//g' | sed 's/\.t[1-9]//g' > tmp/targets.txt

cat ${VCF} | grep 'intergenic' | sed 's/|/\t/g' | awk -v OFS='\t' '{ print $14, $1, $2, $2-150, $2+150, $1"_"$2, "REF="$4,"ALT="$5, $9, $10}' | sed 's/CHR_START-//g' |  sed 's/-CHR_END//g' | sed 's/\.t[1-9]//g' > tmp/targets.txt

#high, moderate and low impact
#cat ${VCF} | grep -E 'HIGH|MODERATE|LOW' | sed 's/|/\t/g' | awk -v OFS='\t' '{ print $14, $1, $2, $2-150, $2+150, $3, "REF="$4,"ALT="$5, $9, $10}' | sed 's/CHR_START-//g' |  sed 's/-CHR_END//g' | sed 's/\.t[1-9]//g' > tmp/targets.txt

#sort list and dropping targets with more than 1 allele
cat tmp/targets.txt | sort -nk 1 | grep -v ',' > tmp/targets_srt.txt

#simplify the gene file
#cat ${ANN} | sed 's/\.t[1-9]//g'  | awk -v OFS='\t' '{ print $1, $6, $7, $8, $9, $10}' > tmp/gene.ann2.tsv

#simplify the gene list
#cat tmp/targets_srt.txt | cut -f1 > tmp/gene_list.txt

# filter the annotation file to get the meaninful annotations. Removing 'consensus disorder'
#cat tmp/gene.ann2.tsv |  grep -v 'consensus' > tmp/gene.ann3.tsv

# Get the genes from annotation [The -f option instructs grep to read the patterns to look for from file2.txt]
#grep -f tmp/gene_list.txt tmp/gene.ann3.tsv > tmp/gene_ann.txt

#sort
#cat tmp/gene_ann.txt | sort -nk 1 > tmp/gene_ann_srt.txt

#3- join snp data with annotation file
#join tmp/targets_srt.txt tmp/gene_ann_srt.txt > tmp/gene_targets.txt

# remove missing annotations
#cat tmp/gene_targets.txt | grep -v ' \- ' >  tmp/gene_targets2.txt

# keep only one variant per gene
#awk '!a[$1]++' tmp/gene_targets2.txt > tmp/gene_targets3.txt

#4- prepare the final bed file
#cat tmp/gene_targets3.txt  | awk -v OFS='\t' '{ print $2, $4, $5, $6"|"$7"|"$8"|"$10"|"$1"|"$11"-"$12"-"$13"-"$14"-"$15}' > final_sites.bed

cat tmp/targets_srt.txt  | awk -v OFS='\t' '{ print $2, $4, $5, $6"|"$7"|"$8"|"$10"|"$1}' > final_sites.bed

#5- extract the sequences 150 up and down with bedtools (this step can be also done on R)

bedtools getfasta -fi ${REF} -bed final_sites.bed -name > final_seq_300bp.fasta

#6- delete intermediate files

#rm -r tmp

#7- run sliding window with 120bp and 10 bp window
## using fasta_windows_v1.1.sh from https://github.com/kdillmcfarland/sliding_windows

source /cluster/projects/khufu/qtl_seq_II/khufu_II/utilities/load_modules.sh 
fasta2kmer final_seq_300bp.fasta 120 2 > candidates.fasta

#bash 02_fasta_windows_v1.1.sh  final_seq_250bp.fasta 120 10

## this will result in a new file: windows_final_seq_250bp.fasta

#8-run blast on the result
 
# blast needs to be in the PATH:  export PATH=$PATH:$HOME/ncbi-blast-2.10.1+/bin
#export PATH=/cluster/home/rsouza/apps/blast/ncbi-blast-2.15.0+/bin

# database needs to be in the PATH: export BLASTDB=$HOME/blastdb
#export BLASTDB=/cluster/home/rsouza/apps/blast/blastdb/sintegrifolium

#blastn -db sintegrifolium -query candidates.fasta -out blast_results.txt -outfmt 6

#9-running the enrich pipeline on R

#eval "$(conda shell.bash hook)"
#conda activate r_env

#R CMD BATCH 03_enrich_v6.R

