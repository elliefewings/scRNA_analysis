# scRNA_analysis - s02_cluster_identity
## R script to perform Seurat based clustering
### Requirements
 * R >= version 3.4
 * scRNAseq data processed by Cellranger count
 * .Rdata file produced by [stage1 of pipeline](https://github.com/elliefewings/scRNA_analysis/tree/master/s01_qc_processing) 

    If you don't have an appropriate Cellranger count output (outs/filtered_feature_bc_matrix), you can create this from fastqs using [this cellranger wrapper](https://github.com/elliefewings/cellranger_wrapper)

### Basic process:
  1. Find clusters and generate UMAP
  2. Find marker genes for clusters
  3. Create heatmap of cluster markers
  4. Apply cell identities to clusters if supplied
  5. Save data
  6. Create PDF and interactive data report
  
## Usage
```
$ s02_cluster_identity/s02_cluster_identity.R -h
Usage: s02_cluster_identity/s02_cluster_identity.R [options]


Options:
        -i INPUT, --input=INPUT
                Path to Rdata output of s01_qc_processing.R (s01_qc_processing.Rdata) [required]

        -o OUTPUT, --output=OUTPUT
                Path to desired output directory [default = directory of input]

        -n NPC, --npc=NPC
                Number of principle components for clustering [default = see elbow plot in s01 report]

        -k RES, --res=RES
                Resolution for clustering (see clustree output in s01 report)[default = 0.5]

        -m MARKERS, --markers=MARKERS
                Path to text file containing marker genes, see README for example

        -h, --help
                Show this help message and exit

```

## Output
By default, data is saved to an `pipeline_output` directory within the input directory (see above for how to set this option).

This directory will include:


    1. s02_cluster_identity.Rdata - Rdata containing Seurat object, clusters and plots
    
    2. sample_name.clusteridentity.report.pdf - Static shiny report showing heatmap of markers and clusters
    
    3. sample_name.clusteridentity.report.html - Link to interactive shiny report showing showing heatmap of markers and clusters
    
    
Shiny interactive report can be accessed [here](https://saezlab.shinyapps.io/cluster_identity_report) or via the link in the output directory.

Simply upload the s02 Rdata file using the "Browse" button or drag and drop and click "Generate report" 
