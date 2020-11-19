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
#  Path to report RDS created in 'create_seurat.R' (e.g. ./data/01a_Report/sample_R.rds) [required]
#  -o OUTPUT, --output=OUTPUT
#  Path to desired output [default = /path/to/input/sample_QC.pdf]
#  -h, --help
#  Show this help message and exit
#
# ============================================================

#############
## Startup ##
#############

## Load libraries
libs <- c("webshot", "shiny", "gridExtra", "stringr", "optparse", "ggplot2", "ggraph")

for (i in libs) {
  if (! suppressPackageStartupMessages(suppressWarnings(require(i, character.only = TRUE, quietly = TRUE)))) { 
    install.packages(i, repos = "https://ftp.fau.de/cran/")
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

## Get options
option_list <- list(
  make_option(c("--input", "-i"), action="store", default="./data/01a_Report/GSM4191941_R.rds", type='character',
              help="Path to report RDS created in 'create_seurat.R' (e.g. ./data/01a_Report/sample_R.rds) [required]"),
  make_option(c("--output", "-o"), action="store", default="./data/01a_Report/GSM4191941_QC.pdf", type='character',
              help="Path to desired output [default = /path/to/input/sample_QC.pdf]")
)

opt2 <- parse_args(OptionParser(option_list=option_list))

## Check for input and output options
if (is.null(opt2$input)) {
  message("ERROR: Input missing, please specify input directory with --input, -i flags.")
  stop(parse_args(OptionParser(option_list=option_list), args = c("--help")))
} else {
  report <- readRDS(opt2$input)
}

# Extract sample name from input
sample <- basename(opt2$input) %>% str_replace_all("_R.rds", "")

if (is.null(opt2$output)) {
  opt2$output <- str_replace_all(opt2$input, "_R.rds", "_QC.pdf")
}

# Unlist report
for(objs in names(report)) {
  assign(objs, report[[objs]])
}


######################
## Run Shiny Report ##
######################

# Source app functions
source("./src/s01_app.R")

# Create app
app <- shinyApp(ui = ui, server = server)

# Create PDF screenshot of app
appshot(app,  opt2$output, delay=10, port = getOption("shiny.port"), vwidth = 1500)

