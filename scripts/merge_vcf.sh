SNP_file_name=$1
INDEL_file_name=$2
Output_name=$3

bgzip ${SNP_file_name}
bgzip ${INDEL_file_name}

bcftools index "${SNP_file_name}.gz"
bcftools index "${INDEL_file_name}.gz"


bcftools concat --allow-overlaps "${SNP_file_name}.gz"  "${INDEL_file_name}.gz" > "${Output_name}.vcf"%
