#!/bin/bash
#SBATCH --job-name sc-atac-seq_mbc_bamslice
#SBATCH --cpus-per-task 8
#SBATCH -c 8
#SBATCH --mem 32g
#SBATCH --partition allnodes
#SBATCH --time 1-00:00:00
#SBATCH --output=log_slurm/slurm_mbc_bamslice.out
#SBATCH --error=log_slurm/slurm_mbc_bamslice.err



cancer_type="male-bc"


source /home/hkim77/share/compsci/python/anaconda3/etc/profile.d/conda.sh
conda activate archr




dir_meme=/home/hkim77/share/compsci/meme/5.4.1
export PATH=$dir_meme/bin:$dir_meme/libexec/meme-5.4.1:$PATH


dir_run=$PWD
dir_bam_subset_archr=$dir_run/output_p2g_male-bc/bam_subset_archr





##### bamslice




dir_cellranger_dna=/home/hkim77/share/compsci/cellranger-dna/cellranger-dna-1.1.0
dir_count=../count_${cancer_type}

for filename in $dir_bam_subset_archr/csv/*_config.csv; do
	name=$(basename $filename)
	name=$(echo $name | cut -d'_' -f 2)
	echo $name
	$dir_cellranger_dna/cellranger-dna bamslice --id=BAMs_${name} \
	     --csv=${filename} \
	     --bam=$dir_count/${name}/outs/possorted_bam.bam
done 
mkdir -p $dir_bam_subset_archr/bam
mv -f BAMs_* $dir_bam_subset_archr/bam

# output
# ./bam_subset_archr/bam/BAMs_446B7L/outs/subsets/
#    446B7L_2-Epithelial_cells.bam
#    446B7L_2-Epithelial_cells.bam.bai







