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
# $ ./s02_cluster_identity.R [options]
#
#Options:
#  -i INPUT, --input=INPUT
#Path to Rdata output of s01_qc_processing.R (s01_qc_processing.Rdata) [required]
#
#-o OUTPUT, --output=OUTPUT
#Path to desired output directory [default = directory of input]
#
#-n NPC, --npc=NPC
#Number of principle components for clustering [default = see elbow plot in s01 report]
#
#-k RES, --res=RES
#Resolution for clustering (see clustree output in s01 report)[default = 0.5]
#
#-m MARKERS, --markers=MARKERS
#Path to text file containing marker genes, see README for example
#
#-h, --help
#Show this help message and exit
#
# ============================================================

#############
## Startup ##
#############

## Load libraries
libs <- c("Seurat", "dplyr", "GetoptLong", "optparse", "magrittr", "stringr", "ggplot2", "webshot", "shiny", "clustree", "tidyr")

for (i in libs) {
  if (! suppressPackageStartupMessages(suppressWarnings(require(i, character.only = TRUE, quietly = TRUE)))) { 
    stop(paste("Unable to find package: ", i, ". Please install and restart.", sep=""))
  }
}


## Find script directory
initial.options <- commandArgs(trailingOnly = FALSE)
script.dir <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))

#For testing
#script.dir <- "C:/Users/ellie/OneDrive/Saez/Pipeline/github/scRNA_analysis/s02_cluster_identity/"

## Get options
option_list <- list(
  make_option(c("--input", "-i"), action="store", default="", type='character',
              help="Path to Rdata output of s01_qc_processing.R (s01_qc_processing.Rdata) [required]"),
  make_option(c("--output", "-o"), action="store", default=NULL, type='character',
              help="Path to desired output directory [default = directory of input]"),
  make_option(c("--npc", "-n"), action="store", default="", type='integer',
              help="Number of principle components for clustering [default = see elbow plot in s01 report]"),
  make_option(c("--res", "-k"), action="store", default=0.5, type='numeric',
              help="Resolution for clustering (see clustree output in s01 report)[default = 0.5]"),
  make_option(c("--markers", "-m"), action="store", default="", type='character',
              help="Path to text file containing marker genes, see README for example"),
  make_option(c("--species", "-s"), action="store", default="mouse", type='character',
              help="Species for assigning cell type identity ('mouse' or 'human')[default = 'mouse']")
)

opt2 <- parse_args(OptionParser(option_list=option_list))

## Check for input and output options
if (!grepl(".Rdata", opt2$input) | !file.exists(opt2$input)) {
  message("ERROR: Input missing or not an .Rdata file, please specify input with --input, -i flags.")
  stop(parse_args(OptionParser(option_list=option_list), args = c("--help")))
}

if (is.null(opt2$output)) {
  opt2$output <- dirname(opt2$input)
}

# Generate warning messing if npcs or resolution is outside expected values

if ( !is.na(opt2$npc) & opt2$npc < 2 ) {
  warning("WARNING: Option --npc/-n is LOW, consider changing to a number greater than 2")
}

if ( !is.na(opt2$npc) & opt2$npc >= 20 ) {
  warning("WARNING: Option --npc/-n is HIGH, consider changing to a number less than 20")
}

if ( opt2$res < 0.1 ) {
  warning("WARNING: Option --res/-k is LOW, consider changing to a number greater than 0.1")
}

if ( opt2$res >= 20 ) {
  warning("WARNING: Option --res/-k is HIGH, consider changing to a number less than 20")
}

# Check if marker file is supplied and if it exists

if (opt2$markers != "" & !file.exists(opt2$markers)) {
  message("ERROR: Markers file does not exist.")
  stop(parse_args(OptionParser(option_list=option_list), args = c("--help")))
}

###############
## Load data ##
###############
#opt2$input <- "C:/Users/ellie/OneDrive/Saez/Pipeline/github/data/CK114/pipeline_output/s01_qc_processing.Rdata"
#opt2$output <- dirname(opt2$input)

# Load Rdata from input
load(opt2$input)

# Reload script directory
initial.options <- commandArgs(trailingOnly = FALSE)
script.dir <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))
#script.dir <- "C:/Users/ellie/OneDrive/Saez/Pipeline/github/scRNA_analysis/s02_cluster_identity/"

# Clean-up previous script data
rm(data.meta.summ, i, indir, libs, pca, qc1, qc1.f, qc2, qc2.f, qc3, qc3.f)


## Source functions
source(paste(script.dir, "/../source/source.R", sep=""))

######################
## Cluster and UMAP ##
######################

# Set NPC from elbow plot if not set in params
if (is.na(opt2$npc)) {
  opt2$npc <- npcs$npcs
}

# Find clusters based on set resolution
data <- FindClusters(data, resolution = opt2$res, verbose = FALSE)

# Run UMAP
data <- RunUMAP(data, reduction = "pca", dims = 1:opt2$npc, verbose = FALSE)

# Plot clusters
pca <- DimPlot(data, reduction = "umap")

##################
## Find Markers ##
##################

# Find markers for all clusters
all.markers <- suppressMessages(FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE))

# Find number of genes needed per cluster to achieve 100 genes total for heatmap
ngenes <- round(100/nlevels(all.markers$cluster))

# Gather 100 genes total from top of each cluster
top <- all.markers %>% group_by(cluster) %>% top_n(ngenes, avg_log2FC)

# Plot heatmap
heat <- suppressWarnings(DoHeatmap(object = data, features = top$gene, label = TRUE))

##########################
## Find Cell Identities ##
##########################

if (opt2$markers != "") { 

  # Read in marker genes
  marker <- read.table(opt2$markers, sep="\t", header = TRUE)
  
  # Assign cluster number to cell identities
  ids <- assign.identity(data, marker)
  
  # If found, apply cell type to cluster
  ids$label <- ifelse(is.na(ids$label1), levels(ids$cluster), paste(ids$cluster, ids$label1,  sep="-"))
  
  # Set new ids
  new.ids <- ids$label
  
  names(new.ids) <- levels(data)
  
  data <- RenameIdents(data, new.ids)
  
  # Replot with new cell identities
  pca2 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 1, label.size = 6) 
  }

if (opt2$markers == "") { 

  # Use Panglao database if no markers
  pldb <- read.table(gzfile(paste(dirname(script.dir), "/source/PanglaoDB_markers_27_Mar_2020.tsv.gz", sep="")), sep = "\t", header=TRUE)

  #Filter database by species
  sp <- ifelse(opt2$species == "mouse", "Mm", "Hs")
  pldb <- pldb %>% filter(grepl(sp, pldb$species))
  
  #Reformat db
  ids <- assign.pldb(data, pldb)
  
  # If found, apply cell type to cluster
  ids$label <- ifelse(ids$label1 == "NA(1)", levels(ids$cluster), paste(ids$cluster, ids$label1,  sep="-"))
  
  # Set new ids
  new.ids <- ids$label
  
  names(new.ids) <- levels(data)
  
  data <- RenameIdents(data, new.ids)
  
  # Replot with new cell identities
  pca2 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 1, label.size = 6) 
  
  }

#######################
## Clean up and Save ##
#######################

# Remove old data
suppressWarnings(rm(marker, top, new.ids, ngenes, option_list, rdata))

# Create output directory
dir.create(opt2$output, showWarnings = FALSE)

# Save image
rdata <- paste(opt2$output, "/s02_cluster_identity.Rdata", sep="")

save.image(rdata)


######################
## Run Shiny Report ##
######################

# Source app functions
source(paste(script.dir, "/cluster_identity_report.pdf/app.R", sep=""))

# Create app
app <- shinyApp(ui = ui, server = server)

# Create PDF screenshot of app
appshot(app,  paste(opt2$output, "/", sample, ".clusteridentity.report.pdf", sep=""),  envvars = c(rdata = rdata), delay=10, port = getOption("shiny.port"), vwidth = 1500)

# Create HTML with link to shiny app io
out.html <- paste(opt2$output, "/", sample, ".clusteridentity.report.html", sep="")

fileConn <- file(out.html)
writeLines(c("<html>",
             "<head>",
             '<meta http-equiv="refresh" content="0; url=https://saezlab.shinyapps.io/cluster_identity_report" />',
             "</head>",
             "</html>"), fileConn)
close(fileConn)
