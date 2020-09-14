#!/usr/bin/env Rscript
# coding: utf-8
# Copyright (C) 2020 Eleanor Fewings
#
# Contact: eleanor.fewings@bioquant.uni-heidelberg.de
#
# ====================
# GNU-GLPv3:
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# A full copy of the GNU General Public License can be found on
# http://www.gnu.org/licenses/.
#
# ============================================================
# DESCRIPTION:
# Processes scRNAseq cellranger output and generates QC report
#   
#   Seurat based process to load cellranger output, generate
#   QC plots, filter, normalise data, and find variable
#   features. 
# ============================================================
# Usage: (from terminal)
#
# $ ./s01_qc_processing.R [options]
#
# Options:
#  -i INPUT, --input=INPUT
#  Path to sample directory containing results of cellranger count [required]
#  
#  -o OUTPUT, --output=OUTPUT
#  Path to desired output directory [default = input]
#  
#  -c MINCELLS, --mincells=MINCELLS
#  Feature filter: Minimum number of cells expressing a feature for it to be included [default = 3]
#  
#  -n MINFEATURES, --minfeatures=MINFEATURES
#  Cell filter: Minimum number of features a cell should express [default = 200]
#  
#  -x MAXFEATURES, --maxfeatures=MAXFEATURES
#  Cell filter: Maximum number of features a cell should express [default = 2500]
#  
#  -m MAXPERCENTMT, --maxpercentmt=MAXPERCENTMT
#  Cell filter: Maximum percentage of mitochondrial features a cell should express [default = 5]
#  
#  -r HASHTAG, --hashtag=HASHTAG
#  Path to UMI hashtag data from CITE-seq
#
#  -h, --help
#  Show this help message and exit
#
# ============================================================

#############
## Startup ##
#############

## Load libraries
libs <- c("Seurat", "dplyr", "GetoptLong", "optparse", "magrittr", "stringr", "ggplot2", "webshot", "shiny", "gridExtra", "RColorBrewer", "clustree")

for (i in libs) {
  if (! suppressPackageStartupMessages(suppressWarnings(require(i, character.only = TRUE, quietly = TRUE)))) { 
    install.packages(i)
    if (! suppressPackageStartupMessages(suppressWarnings(require(i, character.only = TRUE, quietly = TRUE)))) {
      stop(paste("Unable to install package: ", i, ". Please install manually and restart.", sep=""))
      }
    }
}

if (is.null(webshot:::find_phantom())) {
  webshot::install_phantomjs()
}

## Find script directory
initial.options <- commandArgs(trailingOnly = FALSE)
script.dir <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))

## Source functions
source(paste(script.dir, "/../source/source.R", sep=""))

## Get options
option_list <- list(
  make_option(c("--input", "-i"), action="store", default=NULL, type='character',
              help="Path to sample directory containing results of cellranger count [required]"),
  make_option(c("--output", "-o"), action="store", default=NULL, type='character',
              help="Path to desired output directory [default = input]"),
  make_option(c("--mincells", "-c"), action="store", default=3, type='integer',
              help="Feature filter: Minimum number of cells expressing a feature for it to be included [default = 3]"),
  make_option(c("--minfeatures", "-n"), action="store", default=200, type='integer',
              help="Cell filter: Minimum number of features a cell should express [default = 200]"),
  make_option(c("--maxfeatures", "-x"), action="store", default=2500, type='integer',
              help="Cell filter: Maximum number of features a cell should express [default = 2500]"),
  make_option(c("--maxpercentmt", "-m"), action="store", default=5, type='numeric',
              help="Cell filter: Maximum percentage of mitochondrial features a cell should express [default = 5]"),
  make_option(c("--hashtag", "-r"), action="store", default=NULL, type='character',
              help="Path to UMI hashtag data from CITE-seq")
)

opt <- parse_args(OptionParser(option_list=option_list))

## Check for input and output options
if (is.null(opt$input)) {
  message("ERROR: Input missing, please specify input directory with --input, -i flags.")
  stop(parse_args(OptionParser(option_list=option_list), args = c("--help")))
}

if (is.null(opt$output)) {
  opt$output <- paste(opt$input, "/pipeline_output", sep="")
}

## Find filtered feature bc matrix directory
indir <- paste(opt$input, "/outs/filtered_feature_bc_matrix", sep="") %>% str_replace_all("/outs/filtered_feature_bc_matrix/outs/filtered_feature_bc_matrix", "/outs/filtered_feature_bc_matrix") %>% str_replace_all("/outs/filtered_feature_bc_matrix//outs/filtered_feature_bc_matrix", "/outs/filtered_feature_bc_matrix") 

if (!file.exists(indir)) {
  message(paste("ERROR: Input directory doesn't exist or doesn't contain cellranger output: ", indir, sep=""))
  stop(parse_args(OptionParser(option_list=option_list), args = c("--help")))
}

## Check other options
# Minimum number of cells per feature
if (!is.numeric(opt$mincells) | opt$mincells < 0 ) {
  message("ERROR: Option --mincells/-c is not numeric or negative. Please supply a positive integer")
  stop(parse_args(OptionParser(option_list=option_list), args = c("--help")))
}

if (opt$mincells > 10) {
  warning("WARNING: Option --mincells/-c filter is HIGH, consider changing to a number less than 10")
}

if (opt$mincells < 1) {
  warning("WARNING: Option --mincells/-c filter is LOW, consider changing to a number greater than 0")
}

# Minimum number of features per cell

if (!is.numeric(opt$minfeatures)) {
  message("ERROR: Option --minfeatures/-n is not numeric")
  stop(parse_args(OptionParser(option_list=option_list), args = c("--help")))
}

if (opt$minfeatures < 100) {
  warning("WARNING: Option --minfeatures/-n filter is LOW, consider changing to a number greater than 100")
}

if (opt$minfeatures > 1000) {
  warning("WARNING: Option --minfeatures/-n filter is HIGH, consider changing to a number less than 1000")
}

# Maximum number of features per cell

if (!is.numeric(opt$maxfeatures)) {
  message("ERROR: Option --maxfeatures/-x is not numeric")
  stop(parse_args(OptionParser(option_list=option_list), args = c("--help")))
}

if (opt$maxfeatures < 1000) {
  warning("WARNING: Option --maxfeatures/-x filter is LOW, consider changing to a number greater than 1000")
}

if (opt$maxfeatures > 10000) {
  warning("WARNING: Option --maxfeatures/-x filter is HIGH, consider changing to a number less than 10000")
}

if (opt$maxfeatures < opt$minfeatures) {
  message("ERROR: Option --maxfeatures/-x is less than option --minfeatures/-n")
  stop(parse_args(OptionParser(option_list=option_list), args = c("--help")))
}

# Maximum percentage of features from mitochondrial genes

if (!is.numeric(opt$maxpercentmt)) {
  message("ERROR: Option --maxpercentmt/-m is not numeric")
  stop(parse_args(OptionParser(option_list=option_list), args = c("--help")))
}

if (opt$maxpercentmt > 10) {
  warning("WARNING: Option --maxfeatures/-x filter is HIGH, consider changing to a number less than 10")
}

if (opt$maxpercentmt < 2) {
  warning("WARNING: Option --maxfeatures/-x filter is LOW, consider changing to a number greater than 1")
}

# If hashtag mode implemented, check if file exists

if (!is.null(opt$hashtag)) {
  hashdir <- paste(opt$hashtag, "/umi_count", sep="") %>% str_replace_all("/umi_count/umi_count", "/umi_count") %>% str_replace_all("/umi_count//umi_count", "/umi_count") 
  if (!file.exists(hashdir)){
    message(paste("ERROR: Hashtag directory doesn't exist or doesn't contain CITE-seq output: ", hashdir, sep=""))
    stop(parse_args(OptionParser(option_list=option_list), args = c("--help"))) 
  }
}

###############
## Load data ##
###############

# Infer sample name from directory name
sample <- indir %>% str_replace_all("/outs/filtered_feature_bc_matrix", "") %>% basename()

# Read 10x data
input_data <- Read10X(data.dir = indir)

# Find if hashtag mode is implemented and load data including hashtags
if (!is.null(opt$hashtag)) {
  
  umis <- Read10X(data.dir = hashdir, gene.column = 1)
  
  #Remove barcodes with few counts and unmapped
  umis <- umis[rowSums(as.data.frame(umis)) > 100 & row.names(umis) != "unmapped",]
  
  # Add -1 to header
  colnames(umis) <- paste(colnames(umis), "-1", sep="") %>% str_replace_all("-1-1", "-1")
  
  # Find intersection
  joint.bcs <- intersect(colnames(input_data), colnames(umis))
  
  # Subset RNA and HTO counts by the joint cell barcodes
  input_data <- input_data[,joint.bcs]
  umis <- as.matrix(umis[,joint.bcs])

  # Create Seurat Object
  data <- CreateSeuratObject(counts = input_data, 
                             project = sample, 
                             min.cells = opt$mincells)
  
  # Add HTO data
  data[["HTO"]] <- CreateAssayObject(counts = umis)
} else {

  data <- CreateSeuratObject(counts = input_data, 
                           project = sample, 
                           min.cells = opt$mincells)
}
##############
## QC Plots ##
##############

# Find out how mitochrondrial genes are labelled (by selecting MT- with case insensitivity)
mt.patt <- rownames(data)[grepl("^MT-", rownames(data), ignore.case = TRUE)] %>% gsub("-.*", "-", .) %>% unique()

# Create percent.mt metric
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = mt.patt)

# Create QC plots
qc1 <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, cols="grey")

qc2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt", cols="steelblue4")

qc3 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols="steelblue4")

# Subset meta data
data.meta <- data@meta.data

# Create summary table of meta data
data.meta.summ <- summarise(data.meta, ncells = length(orig.ident),
                   med_nCount_RNA = median(nCount_RNA),
                   min_nCount_RNA = min(nCount_RNA),
                   max_nCount_RNA = max(nCount_RNA),
                   med_nFeature_RNA = median(nFeature_RNA),
                   min_nFeature_RNA = min(nFeature_RNA),
                   max_nFeature_RNA = max(nFeature_RNA),
                   med_percent.mt = median(percent.mt),
                   min_percent.mt = min(percent.mt),
                   max_percent.mt = max(percent.mt)) %>% t() %>% as.data.frame()

colnames(data.meta.summ) <- "pre-filtering"

###############
## Filtering ##
###############

data <- data[,(data$nFeature_RNA >= opt$minfeatures) & (data$nFeature_RNA <= opt$maxfeatures) & (data$percent.mt < opt$maxpercentmt)]

# Subset meta data
data.meta <- data@meta.data

# Create summary table of meta data
data.meta.summ$'post-filtering' <- summarise(data.meta, ncells = length(orig.ident),
                            med_nCount_RNA = median(nCount_RNA),
                            min_nCount_RNA = min(nCount_RNA),
                            max_nCount_RNA = max(nCount_RNA),
                            med_nFeature_RNA = median(nFeature_RNA),
                            min_nFeature_RNA = min(nFeature_RNA),
                            max_nFeature_RNA = max(nFeature_RNA),
                            med_percent.mt = median(percent.mt),
                            min_percent.mt = min(percent.mt),
                            max_percent.mt = max(percent.mt)) %>% t()


# Create post filter QC plots
qc1.f <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, cols="lightskyblue3")

qc2.f <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt", cols="steelblue4")

qc3.f <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols="steelblue4")

#################
## Demultiplex ##
#################
if (!is.null(opt$hashtag)) {
  
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  data <- NormalizeData(data, assay = "HTO", normalization.method = "CLR")
  
  # Demultiplex
  data <- HTODemux(data, assay="HTO", positive.quantile = 0.99)
  
  db.count <- table(data$HTO_classification.global) %>% as.data.frame()
  
  db.count$Var1 <- factor(db.count$Var1, levels=c("Doublet", "Singlet", "Negative"))
  
  # Plot number of doublets and singlets
  doublet <- ggplot(db.count, aes(x=Var1, y=Freq, fill=Var1)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("darkblue","lightblue", "grey")) +
    xlab("") +
    ylab("Cell Count") +
    theme(legend.position="none")
  
  # Plot expression amongst different barcodes
  Idents(data) <- 'HTO_maxID'
  
  ridge <- RidgePlot(data, assay = 'HTO', features = rownames(umis), ncol = 1, cols=brewer.pal(n=nrow(umis),"PuBuGn"))
  
  # Extract singlets
  Idents(data) <- 'HTO_classification.global'
  
  data <- subset(data, idents = 'Singlet')
  
  # Add to data table after filtering doublets
  data.meta <- data@meta.data
  data.meta.summ$Singlets <- summarise(data.meta, ncells = length(orig.ident),
                                               med_nCount_RNA = median(nCount_RNA),
                                               min_nCount_RNA = min(nCount_RNA),
                                               max_nCount_RNA = max(nCount_RNA),
                                               med_nFeature_RNA = median(nFeature_RNA),
                                               min_nFeature_RNA = min(nFeature_RNA),
                                               max_nFeature_RNA = max(nFeature_RNA),
                                               med_percent.mt = median(percent.mt),
                                               min_percent.mt = min(percent.mt),
                                               max_percent.mt = max(percent.mt)) %>% t()
  
  }

###################
## Normalisation ##
###################

# Normalise with default log transform and scale
data <- NormalizeData(data, verbose = FALSE)

#######################
## Variable Features ##
#######################

# Find variable features (can be adjusted if requested)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Scale data on variable features (can be changed if requested)
data <- ScaleData(data, verbose = FALSE)

# Run PCA
data <- RunPCA(data, features = VariableFeatures(object = data), verbose = FALSE)

npcs <- get_npcs(seurat_object = data, create_plot = TRUE)

pca <- DimPlot(data, reduction = "pca")

########################
## Initial Clustering ##
########################

data <- invisible(FindNeighbors(data, reduction="pca", dims=1:npcs$npcs, verbose=FALSE))

data <- invisible(FindClusters(data, resolution = seq(from=0.1, to=1.5, by=0.1), verbose=FALSE))

clust <- suppressWarnings(clustree(data)) + theme(
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(size=5)))

#######################
## Clean up and Save ##
#######################

# Remove old data
suppressWarnings(rm(i, data.meta, input_data, libs, initial.options, hashdir, joint.bcs, option_list))

# Create output directory
dir.create(opt$output, showWarnings = FALSE)

# Save image
rdata <- paste(opt$output, "/s01_qc_processing.Rdata", sep="")

save.image(rdata)


######################
## Run Shiny Report ##
######################

# Source app functions
source(paste(script.dir, "/qc_processing_report.pdf/app.R", sep=""))

# Create app
app <- shinyApp(ui = ui, server = server)

# Create PDF screenshot of app
appshot(app,  paste(opt$output, "/", sample, ".qcprocessing.report.pdf", sep=""),  envvars = c(rdata = rdata), delay=10, port = getOption("shiny.port"), vwidth = 1500)

# Create HTML with link to shiny app io
out.html <- paste(opt$output, "/", sample, ".qcprocessing.report.html", sep="")

fileConn <- file(out.html)
writeLines(c("<html>",
             "<head>",
             '<meta http-equiv="refresh" content="0; url=http://saezlab.shinyapps.io/qc_processing_report" />',
             "</head>",
             "</html>"), fileConn)
close(fileConn)
