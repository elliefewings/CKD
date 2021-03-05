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
# Generates QC report from seurat object created in create_seurat.R
#   
# ============================================================
# Usage: (from terminal)
#
# $ ./create_qc_report.R [options]
#
# Options:
#  -i INPUT, --input=INPUT
#  Path to normalised mtx (e.g. ./data/01a_Report/sample.mtx) [required]
#  -o OUTPUT, --output=OUTPUT
#  Path to desired output [default = /path/to/input/sample_QC.pdf]
#  -h, --help
#  Show this help message and exit
#
# ============================================================

#############
## Startup ##
#############

# Source functions
source("./src/seurat_fx.R")

## Load libraries
libs <- c("webshot", "shiny", "gridExtra", "stringr", "optparse", "ggplot2", "ggraph", "Matrix", "Seurat", "dplyr", "clustree")

for (i in libs) {
  suppressPackageStartupMessages(suppressWarnings(require(i, character.only = TRUE, quietly = TRUE)))
}


if (is.null(webshot:::find_phantom())) {
  webshot::install_phantomjs()
}

## Find script directory
initial.options <- commandArgs(trailingOnly = FALSE)
script.dir <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))

## Get options
option_list <- list(
  make_option(c("--input", "-i"), action="store", default="./data/01a_Report/GSM4191941.mtx", type='character',
              help="Path to normalised mtx (e.g. ./data/01a_Report/sample.mtx) [required]"),
  make_option(c("--output", "-o"), action="store", default="./data/01a_Report/GSM4191941_QC.pdf", type='character',
              help="Path to desired output [default = /path/to/input/sample_QC.pdf]")
)

opt2 <- parse_args(OptionParser(option_list=option_list))

#Test input
#opt2$input <- "C:/Users/ellie/OneDrive/Saez/Pipeline/github/CKD/data/testmatrix.mtx"

## Check for input and output options
if (is.null(opt2$input)) {
  message("ERROR: Input missing, please specify input mtx with --input, -i flags.")
  stop(parse_args(OptionParser(option_list=option_list), args = c("--help")))
} else {
  data.in <- read.table(opt2$input)
}

# Extract sample name from input
sample <- basename(opt2$input) %>% str_replace_all(".mtx", "")

if (is.null(opt2$output)) {
  opt2$output <- str_replace_all(opt2$input, ".mtx", "_QC.pdf")
}

############################
## Generate quality plots ##
############################

# Convert matrix to seurat object
data <- CreateSeuratObject(counts = as.matrix(data.in), 
                           project = sample)

# Find out how mitochrondrial genes are labelled (by selecting MT- with case insensitivity)
mt.patt <- rownames(data)[grepl("^MT-", rownames(data), ignore.case = TRUE)] %>% gsub("-.*", "-", .) %>% unique()

# Create percent.mt metric
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = mt.patt)

# Create QC plots
qc1 <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, cols="grey")

qc1.1 <- qc1[[1]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
qc1.2 <- qc1[[2]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
qc1.3 <- qc1[[3]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank())


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

colnames(data.meta.summ) <- "post-filtering"

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

######################
## Run Shiny Report ##
######################

# Source app functions
source("./src/s01_app.R")

# Create app
app <- shinyApp(ui = ui, server = server)

# Create PDF screenshot of app
appshot(app,  opt2$output, delay=10, port = getOption("shiny.port"), vwidth = 1500)

