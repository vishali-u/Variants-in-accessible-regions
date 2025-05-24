# Variants in accessible regions

Prioritizing psychiatric GWAS variants in neurodevelopmental cell types using chromatin accessibility

### chromatin_accessibility
* seurat_tutorial
  * Example analysis of ATAC-seq data using Signac and Seurat: https://stuartlab.org/signac/articles/pbmc_vignette.html
* trevino_2021
  * Paper: https://doi.org/10.1016/j.cell.2021.07.039
  * Data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162170 (GEO: GSE162170)
  * Included single-ATAC-seq data and multiome ATAC-seq data
* ziffra_2021
  * Paper: https://doi.org/10.1038/s41586-021-03209-8
  * Data: supplementary table 2

### gwas
* Paper: https://doi.org/10.1038/s41586-022-04434-5
* Data: supplementary table 11
* Note: SNP locations were converted to reference genome hg38 using https://genome.ucsc.edu/cgi-bin/hgLiftOver

### gwas_dar
* GWAS and ATAC-seq data intersection
* intersect.R: searches for overlapping SNPs and ATAC-seq peaks
* create_scripts.sh: generates shell scripts needed to submit a batch job (saves shell scripts in shell_scripts)
* intersected_files: .rds files generated from batch job
* combine_files.R: combines the .rds files and saves results in intersection.rds

