#!/bin/sh

#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=50G
#SBATCH --time=2-00
#SBATCH --qos=bbdefault
#SBATCH --account=cazierj-msc-bioinf


set -e

module purge; module load bluebear
module load bear-apps/2021b/live
module load Salmon/1.9.0-GCC-11.2.0

out_dir="./results/DEG/temp/salmon_star_3.0/"
mkdir -p ${out_dir}


declare -a unique_fl=("BT594" "BT972" "BT241" "BT935" "MBT168" "MBT357" "MBT373")
index_dir="./data/REF_GENOME/GRCh38.p14/Salmon_index"

for name in "${unique_fl[@]}"; do

	salmon quant -i ${index_dir}  --threads 40  -l A -1 ./data/Trimmed/${name}.R1.fastq.gz -2 ./data/Trimmed/${name}.R2.fastq.gz -o ${out_dir}/${name}

done 


