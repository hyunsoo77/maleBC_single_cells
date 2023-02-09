#!/bin/bash
#SBATCH --job-name sc-atac-seq_mbc_fimo
#SBATCH --cpus-per-task 32
#SBATCH -c 32
#SBATCH --mem 128g
#SBATCH --partition allnodes
#SBATCH --time 1-00:00:00
#SBATCH --output=log_slurm/slurm_mbc_fimo.out
#SBATCH --error=log_slurm/slurm_mbc_fimo.err




dir_meme=/home/hkim77/share/compsci/meme/5.4.1
export PATH=$dir_meme/bin:$dir_meme/libexec/meme-5.4.1:$PATH


dir_work=./output_p2g_male-bc/bam_subset_archr
dir_bam=$dir_work/bam
dir_bed=$dir_work/bed
dir_fa=$dir_work/fa
dir_fimo=$dir_work/fimo
dir_vcf=$dir_work/vcf

mkdir -p $dir_fa
mkdir -p $dir_fimo
mkdir -p $dir_vcf

sample="4CC61L"



bcftools mpileup -d 1000 -Ou -f /datastore/nextgenout5/share/labs/francolab/Data/refdata-cellranger-atac-GRCh38-1.2.0/fasta/genome.fa $dir_bam/BAMs_$sample/outs/subsets/$sample.bam | bcftools call -mv -Oz -o $dir_vcf/${sample}_variants.vcf.gz



bcftools index $dir_vcf/${sample}_variants.vcf.gz
bcftools view --types snps $dir_vcf/${sample}_variants.vcf.gz > $dir_vcf/${sample}_variants_new.vcf
bgzip -c $dir_vcf/${sample}_variants_new.vcf > $dir_vcf/${sample}_variants_new.vcf.gz
tabix -p vcf $dir_vcf/${sample}_variants_new.vcf.gz




# apply cluster variants to output cluster reference genome
cat /datastore/nextgenout5/share/labs/francolab/Data/refdata-cellranger-atac-GRCh38-1.2.0/fasta/genome.fa | bcftools consensus $dir_vcf/${sample}_variants_new.vcf.gz > $dir_fa/${sample}_consensus.fa



echo "get sequence of enhancers and promoter"

# get sequence of enhancers and promoter
echo -e "\t${sample}_ANXA2"
bedtools getfasta -fi $dir_fa/${sample}_consensus.fa -fo $dir_fa/${sample}_ANXA2_enhancer1.fa -bed $dir_bed/ANXA2_enhancer1.bed
bedtools getfasta -fi $dir_fa/${sample}_consensus.fa -fo $dir_fa/${sample}_ANXA2_promoter.fa -bed $dir_bed/ANXA2_promoter.bed

echo -e "\t${sample}_PRDX4"
bedtools getfasta -fi $dir_fa/${sample}_consensus.fa -fo $dir_fa/${sample}_PRDX4_enhancer1.fa -bed $dir_bed/PRDX4_enhancer1.bed
bedtools getfasta -fi $dir_fa/${sample}_consensus.fa -fo $dir_fa/${sample}_PRDX4_promoter.fa -bed $dir_bed/PRDX4_promoter.bed






echo "run FIMO motif scanning"

filename_motif=/home/hkim77/share/bio/protein/tf_motif_database/jaspar/jaspar2020/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme

echo -e "\t${sample}_ANXA2"
fimo --bgfile motif-file --oc $dir_fimo/${sample}_ANXA2_enhancer1_fimo $filename_motif $dir_fa/${sample}_ANXA2_enhancer1.fa
fimo --bgfile motif-file --oc $dir_fimo/${sample}_ANXA2_promoter_fimo $filename_motif $dir_fa/${sample}_ANXA2_promoter.fa

echo -e "\t${sample}_PRDX4"
fimo --bgfile motif-file --oc $dir_fimo/${sample}_PRDX4_enhancer1_fimo $filename_motif $dir_fa/${sample}_PRDX4_enhancer1.fa
fimo --bgfile motif-file --oc $dir_fimo/${sample}_PRDX4_promoter_fimo $filename_motif $dir_fa/${sample}_PRDX4_promoter.fa



