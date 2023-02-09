# maleBC_single_cells

This repository is for male breast cancer scRNA-seq and scATAC-seq data processing and figure generation.

The general flow of data processing and figure generation is (1) scRNA-seq data processing from 10x Genomics scRNA-seq FASTQ files, (2) scATAC-seq data processing from 10x Genomics scATAC-seq FASTQ files, (3) figure generation.



### Installing

```
git clone https://github.com/hyunsoo77/maleBC_single_cells.git
```



### scRNA-seq data processing

Step 1: align sequences in scRNA-seq FASTQ files to GRCh38 reference transcriptome by 10x Genomics [cellranger count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct) to obtain two filtered_feature_bc_matrix.h5 files for two samples.

Step 2: Make the following directoy structure with copy or link.

```
../count_male-bc
├── Patient1
│   ├── outs
│   │   └── filtered_feature_bc_matrix.h5
└── Patient2
    └── outs
        └── filtered_feature_bc_matrix.h5
```


Step 3: Make Seurat object for each sample with the following command:

```
./make_sc-rna-seq_seurat_obj.R --dir_count ../count_male-bc --dir_output ./output_male-bc --dir_seurat_obj ./output_male-bc/rds_male-bc --type_qc arguments --min_ncount_rna 5000 --min_nfeature_rna 2000 --th_percent.mt 25 --max_dimstouse 30 --seurat_resolution 0.8 --method_to_update_cell_types epithelial_cell_types --method_to_identify_subtypes none --type_infercnv_argset vignettes --method_to_determine_th_cna_value_corr fixed --th_cna_value 0.05 --th_cna_corr 0.35 male-bc Patient1
```

The above example is only for Patient1, you can make another Seurat object for Patient2 by changing the last argument. The contents of the output directory of "./output_male-bc" follows:

```
output_male-bc/
├── infercnv
│   ├── male-bc_Patient1_cnv_postdoublet
│   └── male-bc_Patient2_cnv_postdoublet
├── output
│   └── log
├── rds_male-bc
│   ├── male-bc_Patient1_sc-rna-seq_sample_seurat_obj.rds
│   ├── male-bc_Patient2_sc-rna-seq_sample_seurat_obj.rds
│   └── wilcox_degs
├── tsv
│   ├── infercnv_input_barcode_group_male-bc_Patient1.tsv
│   └── infercnv_input_barcode_group_male-bc_Patient2.tsv
└── xlsx
    ├── male-bc_Patient1_sc-rna-seq_pipeline_summary.xlsx
    └── male-bc_Patient2_sc-rna-seq_pipeline_summary.xlsx
```

Step 4: Merge Seurat objects for multiple samples to make merged Seurat object by the following command:

```
./make_sc-rna-seq_merged_seurat_obj.R --dir_output ./output_male-bc --dir_seurat_obj ./output_male-bc/rds_male-bc --type_parsing_rds_filename unc-male-bc --method_integration none --max_dimstouse 30 --seurat_resolution 0.2 --harmony_theta 0 male-bc
```

The output file is located under ./output_male-bc/rds_male-bc that was defined by an argument of --dir_seurat_obj.


```
output_male-bc/
│   ...
├── rds_male-bc
│   ├── male-bc_Patient1_sc-rna-seq_sample_seurat_obj.rds
│   ├── male-bc_Patient2_sc-rna-seq_sample_seurat_obj.rds
│   ├── male-bc_sc-rna-seq_merged_seurat_obj.rds
│   └── wilcox_degs
...
```















### scATAC-seq data processing


Step 1: align sequences in scATAC-seq FASTQ files to GRCh38 reference genome by 10x Genomics [cellranger-atac count](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count) to obtain two fragments.tsv.gz files for two samples.


Step 2: Make the following directoy structure with copy or link.

```
../count_male-bc
├── Patient1
│   ├── outs
│   │   ├── fragments.tsv.gz
│   │   └── fragments.tsv.gz.tbi
└── Patient2
    └── outs
        ├── fragments.tsv.gz
        └── fragments.tsv.gz.tbi
```

Step 3: Make an ArchRProject object for two samples with the following command:

```
./make_sc-atac-seq_archr_obj.R --n_cores 60 --dir_count ../count_male-bc --dir_output ./output_male-bc --dir_seurat_obj ./dir_seurat_obj_male-bc --keep_most_common_cell_type_for_each_cluster --max_dimstouse 30 --seurat_resolution 0.2 --min_tss 0 --min_frags 1000 --colorlim 1.5 --umap_n_neighbors 30 --umap_min_dist 0.3 --umap_metric cosine --umap_metric_doubletscores cosine --harmony_theta 0 --harmony_lambda 20 male-bc
```

This will take care of all samples under ../count_male-bc. The contents of the output directory of "./output_male-bc" follows:


```
output_male-bc
├── archr_output
│   ├── Annotations
│   ├── ArrowFiles
│   ├── Embeddings
│   ├── LSI_ATAC
│   ├── Peak2GeneLinks
│   ├── PeakCalls
│   │   ├── InsertionBeds
│   │   └── ReplicateCalls
│   ├── Plots
│   └── RNAIntegration
│       └── GeneIntegrationMatrix_ArchR
├── log
├── pdf
├── png
│   ├── male-bc_Patient1_qc.png
│   ├── male-bc_Patient2_qc.png
│   └── ...
├── qc
│   ├── Patient1
│   │   ├── Patient1-Doublet-Summary.pdf
│   │   ├── Patient1-Fragment_Size_Distribution.pdf
│   │   └── Patient1-TSS_by_Unique_Frags.pdf
│   └── Patient2
│       ├── Patient2-Doublet-Summary.pdf
│       ├── Patient2-Fragment_Size_Distribution.pdf
│       └── Patient2-TSS_by_Unique_Frags.pdf
├── rds
│   └── male-bc_archrproj_obj_final.rds
├── tmp
└── xlsx
    └── male-bc_sc-atac-seq_pipeline_summary.xlsx
```


Let's review the content of QC files (male-bc_Patient1_qc.png and male-bc_Patient2_qc.png).


<p align="center">
<img src="https://github.com/hyunsoo77/maleBC_single_cells/blob/main/png/male-bc_Patient1_qc.png" width="300" height="300">
<img src="https://github.com/hyunsoo77/maleBC_single_cells/blob/main/png/male-bc_Patient2_qc.png" width="300" height="300">
</p>



Step 4: Add peak-to-gene links to the final ArchRProject object and perform motif enrichment analysis by the following command:



```
./analyze_cancer_specific_p2g.R --n_cores 4 --dir_output ./output_p2g_male-bc --dir_seurat_obj ./dir_seurat_obj_male-bc --dir_archr_output ./output_male-bc/archr_output --subset_archrproject_force_to_update --max_dimstouse 30 --harmony_theta 0 --seed_kmeans 4 --exclude_cluster_epi_unassigned --exclude_cluster_highly_overlapped_with_normal_peaks --exclude_cluster_low_nfrags male-bc
```


This will subset the final ArchRProject in order to remove clustes of unassigned epithelail cells by InferCNV and low mean of nFrags. The peak-to-gene links were obtained by addPeak2GeneLinks() (see document [ArchR book/peak2genelinkage-with-archr](https://www.archrproject.com/bookdown/peak2genelinkage-with-archr.html)). The contents of the output directory of "./output_p2g_male-bc" follows:



```
output_p2g_male-bc
├── archr_output
│   ├── Annotations
│   ├── ArrowFiles
│   ├── Embeddings
│   ├── LSI_ATAC
│   ├── Peak2GeneLinks
│   │   ├── seATAC-Group-KNN.rds
│   │   └── seRNA-Group-KNN.rds
│   ├── PeakCalls
│   │   ├── ...
│   │   ├── InsertionBeds
│   │   ├── ReplicateCalls
│   │   ├── X0.Epi..Tumor-reproduciblePeaks.gr.rds
│   │   ├── X2.Epi..Tumor-reproduciblePeaks.gr.rds
│   │   └── ...
│   ├── Plots
│   ├── RNAIntegration
│   └── Save-ArchR-Project.rds
├── log
├── pdf
│   ├── scatterplot_cisbp_motif_up.normal-vs-cancer_cancer-specific_enhancer.pdf
│   ├── venn_chippeakanno_overlaps_of_enhancer_peaks.pdf
│   └── ...
├── png
├── rds
│   ├── all_p2g_observed.rds
│   ├── cancer_enriched_enhancer_p2g_table.rds
│   ├── cancer_specific_enhancer_p2g_table.rds
│   ├── cancer_specific_p2g_table_degs.rds
│   ├── find_overlaps_of_enhancer_peaks_output_overlappingpeaks_obj.rds
│   ├── male-bc_archrproj_obj_p2gs.rds
│   ├── markerpeaks.for_comparison.normal-vs-cancer_cancer-specific_enhancer.rds
│   ├── p2g.df.sub.plot_enhancer.rds
│   ├── proj.archr.for_comparison.normal-vs-cancer_cancer-specific_enhancer.rds
│   ├── proj.archr.for_comparison.normal-vs-cancer_non-cancer-specific_enhancer.rds
│   ├── venn_chippeakanno_overlaps_of_enhancer_peaks.rds
│   ├── ...
├── tmp
├── tsv
│   ├── cancer_specific_enhancer_p2g_table.tsv
│   ├── peakannoenrichment_cisbp_motif_up.normal-vs-cancer_non-cancer-specific_enhancer.tsv
│   └── ...
└── xlsx
    └── male-bc_sc-atac-seq_cancer_specific_p2g_summary.xlsx

```


The final ArchProject object (male-bc_archrproj_obj_p2gs.rds) is located under output_p2g_male-bc/rds. The peak2gene data.frame for all enhancers related with epithelial cells was stored at cancer_enriched_enhancer_p2g_table.rds. The peak2gene data.frame for cancer specific enhancers was stored at cancer_specific_enhancer_p2g_table.rds. The cancer specific enhancers were defined by subtracting cancer enriched enhancers by peak2gene links overlapped with enhancer peaks in normal epithelial cells (i.e. human mamary epithelial cells (HMEC) H3K27ac peaks in our case).







### Jupyter notebook

Figures were generated by Jupyter notebook scripts. In order to install Jupyter notebook/lab, see [jupyter.org](https://jupyter.org/). You need to change dir_rna and/or dir_atac to locate the merged Seurat object or final ArchRProject object you generated. The output files include PDF files that will be located at the directory of "pdf".


```
./
├── figure1.ipynb
├── figure2_dge.ipynb
├── figure3_01_piedonut.ipynb
├── figure3_02_venn.ipynb
├── figure3_03_heatmap_peak2gene.ipynb
├── figure3_04_normal-vs-cancer_enhancers.ipynb
├── figure4_01_boxplot_for_browser_track.ipynb
├── figure4_02_volcanoplot_dge.cse.ipynb
├── figure4_03_heatmap_dge.cse.ipynb
├── figure4_04_browser_track.ipynb
├── figure_s1_01_sc-rna-seq_qc.ipynb
├── figure_s1_02_sc-atac-seq_qc.ipynb
├── figure_s1_03_sc-atac-seq_qc.ipynb
├── log
├── pdf
│   ├── ...
│   ├── featureplot_male-bc_er+bc-epi_er+bc_vs_male-bc.pdf
│   ├── heatmap_male-bc_er+bc-epi_er+bc_vs_male-bc.pdf
│   ├── heatmap_male-bc_er+bc-epi_er+bc_vs_male-bc_zscore.pdf
│   ├── heatmap_peak2gene_legend.pdf
│   ├── piedonut_peak_call_summary.pdf
│   ├── umap_male-bc_cluster_labels_atac.pdf
│   ├── umap_male-bc_cluster_labels_rna.pdf
│   ├── umap_male-bc_cluster_types_atac.pdf
│   ├── umap_male-bc_cluster_types_rna.pdf
│   ├── ...
│   └── volcanoplot_male-bc_sc-rna-seq_female-bc-vs-male-bc.enhancer_overlap.pdf
├── r
├── rds
├── reference
├── tsv
│   ├── df_dge.cse.tsv
│   └── male-bc_er+bc-epi_er+bc_vs_male-bc.tsv
├── txt
└── xlsx
    ├── ...
    └── male-bc_er+bc-epi_er+bc_vs_male-bc.xlsx
```



The scRNA-seq pipeline and scATAC-seq pipeline are actively developed. Other single cell data analysis projects will use the current version with different parameters or upgraded version of these pipelines.  





