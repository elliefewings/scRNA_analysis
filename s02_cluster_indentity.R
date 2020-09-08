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
libs <- c("Seurat", "dplyr", "GetoptLong", "optparse", "magrittr", "stringr", "ggplot2", "webshot", "shiny", "clustree")

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

#For testing
script.dir <- "C:/Users/ellie/OneDrive/Saez/Pipeline/github/scRNA_analysis/"

## Source functions
source(paste(script.dir, "/source/source.R", sep=""))

## Get options
option_list <- list(
  make_option(c("--input", "-i"), action="store", default="", type='character',
              help="Path to Rdata output of s01_qc_processing.R (s01_qc_processing.Rdata) [required]"),
  make_option(c("--output", "-o"), action="store", default=NULL, type='character',
              help="Path to desired output directory [default = directory of input]"),
  make_option(c("--npc", "-n"), action="store", default=NULL, type='character',
              help="Number of principle components for clustering [default = decided by elbow plot in previous stage of pipeline]"),
  make_option(c("--hashtag", "-r"), action="store", default=NULL, type='character',
              help="Path to UMI hashtag data from CITE-seq")
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


###############
## Load data ##
###############

# Load Rdata from input
load(opt2$input)

# Reload script directory
initial.options <- commandArgs(trailingOnly = FALSE)
script.dir <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))

# Clean-up previous script data
rm(data.meta.summ, i, indir, libs, pca, qc1, qc1.f, qc2, qc2.f, qc3, qc3.f)

###################
## Find Clusters ##
###################

data <- FindNeighbors(data, reduction="pca", dims=1:npcs$npcs)

data <- FindClusters(data, resolution = seq(from=0.1, to=1.5, by=0.1))

clustree(data)

#######################
## Clean up and Save ##
#######################

# Remove old data
suppressWarnings(rm(i, data.meta, input_data, libs, initial.options, hashdir, joint.bcs))

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
