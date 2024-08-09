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


fastqdir="./data/Trimmed"
ref_ind="./data/REF_GENOME/GRCh38.p14/indices"
bam_dir="./results/DEG/temp/bams"
gtf_ref="./data/REF_GENOME/GRCh38.p14/gencode.v45.basic.annotation.gtf"

mkdir -p ${bam_dir}

#Improvement should be made so this list is automatically generated.
declare -a unique_fl=("BT594" "BT972" "BT241" "BT935" "MBT168" "MBT357" "MBT373")

#Extract unique file names 
full_names=()

for file in "${fastqdir}"/*.gz; do
	echo ${file}
	a=$(echo ${file} | cut -d'_' -f 1)
	b=$(echo ${file} | cut -d'_' -f 2)
	names+=("${a}_${b}")
	#gzip -d $file

done


#Carries out alignment and produces read counts table

for name in "${unique_fl[@]}"; do
	echo ${name}
	STAR --runThreadN 40\
	--genomeDir "${ref_ind}" \
	--genomeLoad NoSharedMemory\
	--readFilesIn "${fastqdir}/${name}.R1.fastq.gz" "${fastqdir}/${name}.R2.fastq.gz" \
	--readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode GeneCounts\
	--sjdbGTFfile ${gtf_ref}\
	--outFileNamePrefix "${bam_dir}/${name}/"\
	 --limitBAMsortRAM 1904111303 \
	--quantTranscriptomeBan IndelSoftclipSingleend

	echo "${name} aligned!!"
done
