# Copyright (C) 2020 Javier Perales-Patón 
#
# Contact: javier.perales@bioquant.uni-heidelberg.de
#
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
# Create Seurat Object following standard workflow with SCTransform normalization
# ===============================================================================
#
#  It creates a Seurat Object and few txt files. These could follow downstream analysis of
#  QC test and cell clustering/identity annotation.
#  The input is the YAML file and the cellranger's matrices output of UMI counts and barcodes.
#  It is assumed that cellranger is run in another in-house facility, so genomics data
#  (seq data) is not needed to be copied/shared externally here.
#  to be copied/shared externally.
#  The output data is exported as follow:
#      - Seurat Object: a RDS file in the OUTDIR.
#	- nPCs : plain-text file with Number of Principal Components selected.
#	- nCells: plain-text file with Number of cells before and after filtering.

#--- Env
set.seed(1234)

#--- Libraries
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(require(yaml))
source("./src/seurat_fx.R")

#--- Get input parameters
option_list = list(
	make_option(c("--MATRIX"), action="store", 
	  default="./data/00_scRNAseq/GSM4191941/matrix.mtx.gz", type='character',
	  help="Input: CellRanger's output matrix with UMI counts per cell."),
	make_option(c("--FEATURES"), action="store", 
	  default="./data/00_scRNAseq/GSM4191941/features.tsv.gz", type='character',
	  help="Input: CellRanger's output genes table."),
	make_option(c("--BARCODES"), action="store", 
	  default="./data/00_scRNAseq/GSM4191941/barcodes.tsv.gz", type='character',
	  help="Input: CellRanger's output cell barcodes."),
	make_option(c("--SEURATOBJ"), action="store", 
	  default="./data/01_Seurat/GSM4191941/S.rds", type='character',
	  help="Output file to store the Seurat Object as RDS."),
	make_option(c("--NPCS"), action="store", 
	  default="./data/01_Seurat/GSM4191941/nPCs.txt", type='character',
	  help="Output plain-text file with selected N principal components."),
	make_option(c("--NCELLS"), action="store", 
	  default="./data/01_Seurat/GSM4191941/nCells.txt", type='character',
	  help="Output plain-text file with recorded N cells before and after filtering."),
	make_option(c("--YAML"), action="store", 
	  default="./index/samples/GSM4191941.yaml", type='character',
	  help="YAML file with settings for the sample.")
)

# Parse the parameters
opt = parse_args(OptionParser(option_list=option_list))

# Cat the input parameters
cat("[INFO] Input parameters\n", file=stdout())
for(user_input in names(opt)) {
	  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
    assign(user_input,opt[[user_input]])
}

# Prepare INPUT directory for Seurat
# That is, the directory containing the three cellranger output files
in1 <- dirname(MATRIX)
in2 <- dirname(FEATURES)
in3 <- dirname(BARCODES)
stopifnot(in1==in2 && in1 == in3)
INDIR <- in1
rm(in1, in2, in3)

# Prepare OUTPUT directory
out1 <- dirname(SEURATOBJ)
out2 <- dirname(NCELLS)
out3 <- dirname(NPCS)
stopifnot(out1==out2 && out1 == out3)
OUTDIR <- out1
rm(out1, out2, out3)

# Sanity check
if(!dir.exists(INDIR)) {
	stop("[ERROR] Input directory does _NOT_ exist.\n")
}
if(!dir.exists(OUTDIR)) {
	cat("[INFO] Creating output directory.\n", file=stdout())
	dir.create(OUTDIR)
}
if(!file.exists(YAML)) {
	stop("[ERROR] Input YAML file does _NOT_ exist.\n")
}

#--- Read index
# Read index file
idx <- yaml::read_yaml(YAML)
# Reformat for easy access to info
idx$seurat_cutoff <- unlist(lapply(idx$seurat_cutoff, function(z) setNames(z[[1]], names(z))))
idx$metadata <- unlist(lapply(idx$metadata, function(z) setNames(z[[1]], names(z))))

#--- Creating Seurat
# Getting ID of sample
# It is a standard practice to name the prefix of YAML as the sample ID.
ID <- strsplit(basename(YAML), split="\\.")[[1]][0]
cat(paste0("Processing ",ID, "\n"), file=stdout())
# Find folder of 10x cellranger count
names(INDIR) <- ID

dat <- Read10X(data.dir=INDIR)
S <- CreateSeuratObject(counts = dat, project = ID, 
			min.cells=5, 
			min.features=idx$seurat_cutoff["MIN_nFeatures"])
S[["percent.mt"]] <- PercentageFeatureSet(S, pattern = "^MT-")

#--- Filter-out
# Before
before_nCells <- ncol(S)
# apply Filter
S <- subset(S, subset= nFeature_RNA > idx$seurat_cutoff["MIN_nFeatures"] & 
	    		nFeature_RNA < idx$seurat_cutoff["MAX_nFeatures"] &
	    percent.mt < idx$seurat_cutoff["MAX_percMT"])
# After
after_nCells <- ncol(S)

#--- Normalization using scTransform
#NOTE: this single command replaces NormalizeData, ScaleData, and FindVariableFeatures.
# Transformed data will be available in the SCT assay, which is set as the default after running sctransform
# During normalization, we can also remove confounding sources of variation, for example, 
# mitochondrial mapping percentage
S <- SCTransform(S, vars.to.regress = "percent.mt", verbose = FALSE)

#--- Add metadata
for(feat in names(idx$metadata)) {
	    S[[feat]] <- idx$metadata[[feat]]
}

#--- Finding Variable genes
S <- FindVariableFeatures(S, selection.method = "vst", nfeatures = 2000)

#--- Run PCA, UMAP
S <- RunPCA(S, npcs=50, features=VariableFeatures(S))
calcNPCs <- get_npcs(S)
nPCs <- calcNPCs$npcs

S <- RunUMAP(S, dims= 1:nPCs)

#--- Doublet detection (using DoubletFinder)
# No ground-truth: pK identification
sweep.res.list <- paramSweep_v3(S, PCs = 1:nPCs, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
# Do I have to select here pK based on something??
# > Optimal pK for any scRNA-seq data can be manually discerned as maxima in BCmvn distributions
# NOTE: pK is encoded as a factor. But actually is a number. If you transform
# a factor level to a number directly you will get level position. So as.char is needed.
pK_sel <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))

# Doublet proportion estimation
## Assuming 0.75% doublet formation rate per each 1000 cells in 10x run.
## src: https://assets.ctfassets.net/an68im79xiti/1eX2FPdpeCgnCJtw4fj9Hx/7cb84edaa9eca04b607f9193162994de/CG000204_ChromiumNextGEMSingleCell3_v3.1_Rev_D.pdf
## NOTE: Just search by 'multiplet' in the PDF above.
nExp_poi <- (before_nCells/1e3 * 0.0075)  
S <- doubletFinder_v3(S, PCs = 1:nPCs, pN = 0.25, 
			       pK = pK_sel, 
			       nExp = nExp_poi, reuse.pANN = FALSE, 
			       sct = TRUE)
# Rename cols
DF_idx <- grep("^(pANN_|DF\\.)", colnames(S@meta.data))
colnames(S@meta.data)[DF_idx] <- gsub("^(pANN|DF)(_|\\.)?.*","\\1",colnames(S@meta.data)[DF_idx])

#--- Save object
saveRDS(S, SEURATOBJ)
cat(c(before_nCells, after_nCells), sep="\n", 
    file=NCELLS)
cat(nPCs, sep="\n", 
    file=NPCS)

#--- Show sessionInfo
sessionInfo()
