#!/bin/bash

#SBATCH --ntasks=20
#SBATCH --time=6-00
#SBATCH --qos=bbdefault
#SBATCH --account=cazierj-msc-bioinf

set -e


module purge; module load bluebear
module load bear-apps/2019b/live
module load bear-apps/2019a/live
module load VarScan/2.4.4-Java-1.8
module load SAMtools/1.9-GCC-8.2.0-2.31.1

# Defines variables. see README.md for file storage.
ref_genome="./REF_GENOME/GRCh38.p14/GRCh38.p14.genome.fa"

SNP_file_name="./results/var_call/rm_S_BT972_SNVS.vcf"
INDEL_file_name="./results/var_call/rm_S_BT972_INDELS.vcf"
Output_name="./results/var_call/rm_S_BT972_merged.vcf"



# Lists all files excluding outlier
bam_files=$(ls ./results/DEG/temp/bams/*.bam | grep -v -E 'S_BT972') # Excludes outliers

# Runs variant calling determining SNVs and INDELs
samtools mpileup  -B -f ${ref_genome} ${bam_files} | java -jar  $EBROOTVARSCAN/VarScan.v2.4.4.jar mpileup2snp  /Pileup_wout_outliers/rm_S_BT972_samples.mpileup --vcf-sample-list all_samples.txt --p-value 0.01 --output-vcf 1 > "${SNP_file_name}"
samtools mpileup  -B -f ${ref_genome} ${bam_files} | java -jar  $EBROOTVARSCAN/VarScan.v2.4.4.jar mpileup2indel /Pileup_wout_outliers/rm_S_BT972_samples.mpileup --vcf-sample-list all_samples.txt --p-value 0.01 --output-vcf 1 > "${INDEL_file_name}"

echo "run successful, analysis completed. GoodNight Node :)"

# Zips and concatenates SNV and INDEL vcfs
bgzip ${SNP_file_name}
bgzip ${INDEL_file_name}

bcftools index "${SNP_file_name}.gz"
bcftools index "${INDEL_file_name}.gz"


bcftools concat --allow-overlaps "${SNP_file_name}.gz"  "${INDEL_file_name}.gz" > "${Output_name}"