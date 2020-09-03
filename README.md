# scRNA_analysis

## R script to perform Seurat based analysis of cellranger output
### Basic process:
  1. Creation of quality control plots
  2. Filtering of cells and features based on default or input criteria
  3. Data normalisation
  4. Find variable features
  5. Scale data
  6. Run PCA and find appropriate principle components
  7. Save data
  8. Create PDF and interactive data report
  
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
  
  -h, --help
  Show this help message and exit
```

## Output
By default, data is saved to an `pipeline_output` directory within the input directory (see above for how to set this option).

This directory will include:


    1. s01_qc_processing.Rdata - Rdata containing Seurat object and generated QC figures
    
    2. sample_name.qcprocessing.report.pdf - Static shiny report showing raw data metrics and principle component selection
    
    3. sample_name.qcprocessing.report.html - Link to interactive shiny report showing raw and filtered data metrics and interactive princple component exploration
    
    
Shiny interactive report can be accessed [here](http://saezlab.shinyapps.io/qc_processing_report) or via the link in the output directory.

Simply upload Rdata file using the "Browse" button or drag and drop and click "Generate report" 
