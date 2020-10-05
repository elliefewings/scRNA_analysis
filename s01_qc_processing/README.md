# scRNA_analysis - s01_qc_processing
## R script to perform Seurat based processing and qc of cellranger output

### Requirements
 * R >= version 3.4
 * scRNAseq data processed by Cellranger count

    If you don't have an appropriate Cellranger count output (outs/filtered_feature_bc_matrix), you can create this from fastqs using [this cellranger wrapper](https://github.com/elliefewings/scRNA_raw_data_toolkit)

### Basic process:
  1. Creation of quality control plots
  2. Filtering of cells and features based on default or input criteria
  3. Demultiplex if hashtag option applied
  4. Data normalisation
  5. Find variable features
  6. Scale data
  7. Run PCA and find appropriate principle components
  8. Save data
  9. Create PDF and interactive data report
  
## Usage
```
 $ ./s01_qc_processing.R [options]

 Options:
  -i INPUT, --input=INPUT
  Path to sample directory containing results of cellranger count [required]
  
  -o OUTPUT, --output=OUTPUT
  Path to desired output directory [default = input]
  
  -c MINCELLS, --mincells=MINCELLS
  Feature filter: Minimum number of cells expressing a feature for it to be included [default = 3]
  
  -n MINFEATURES, --minfeatures=MINFEATURES
  Cell filter: Minimum number of features a cell should express [default = 200]
  
  -x MAXFEATURES, --maxfeatures=MAXFEATURES
  Cell filter: Maximum number of features a cell should express [default = 2500]
  
  -m MAXPERCENTMT, --maxpercentmt=MAXPERCENTMT
  Cell filter: Maximum percentage of mitochondrial features a cell should express [default = 5]
  
  -r HASHTAG, --hashtag=HASHTAG
  Path to UMI hashtag data from CITE-seq
  
  -h, --help
  Show this help message and exit
```

## Output
By default, data is saved to an `pipeline_output` directory within the input directory (see above for how to set this option).

This directory will include:


    1. s01_qc_processing.Rdata - Rdata containing Seurat object and generated QC figures
    
    2. sample_name.qcprocessing.report.pdf - Static shiny report showing raw data metrics and principle component selection
    
    3. sample_name.qcprocessing.report.html - Link to interactive shiny report showing raw and filtered data metrics and interactive princple component exploration
    
    
Shiny interactive report can be accessed [here](https://saezlab.shinyapps.io/qc_processing_report/) via the link in the output directory.

Simply upload s01 Rdata file using the "Browse" button or drag and drop and click "Generate report" 
