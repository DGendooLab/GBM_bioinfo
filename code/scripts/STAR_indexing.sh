#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=50G
#SBATCH --time=2-00
#SBATCH --qos=bbdefault
#SBATCH --account=cazierj-msc-bioinf

set -e

module purge; module load bluebear

module load bear-apps/2022b/live
module load STAR/2.7.11a-GCC-12.2.0

STAR --runThreadN 30 --runMode genomeGenerate --genomeDir ./data/REF_GENOME/GRCh38.p14/indices --genomeFastaFiles ./data/REF_GENOME/GRCh38.p14/GRCh38.p14.genome.fa --sjdbGTFfile ./data/REF_GENOME/GRCh38.p14/gencode.v45.basic.annotation.gtf  --sjdbOverhang 50
