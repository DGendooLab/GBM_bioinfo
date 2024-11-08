#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=50G
#SBATCH --time=2-00
#SBATCH --qos=bbdefault
#SBATCH --account=cazierj-msc-bioinf

set -e

module purge; module load bluebear
module load bear-apps/2021b/live; module load Salmon/1.9.0-GCC-11.2.0


transcriptome="/rds/projects/g/gendood-preclinomics/Marcello_Thesis/Project_1/REF_GENOME/GRCh38.p14/Transcripts/gencode.v45.transcripts.fa"  # Path to your transcriptome fasta file
index_dir="/rds/projects/g/gendood-preclinomics/Marcello_Thesis/Project_1/REF_GENOME/GRCh38.p14/"          # Directory to store the index
index_name="Salmon_index" 

salmon index -t $transcriptome -i $index_dir/$index_name --gencode
