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
#  It creates a Seurat Object and an image with a bundle of objects to follow-up QC.
#  The input are indexes that point-out to the output of cellranger. It is assumed that
#  cellranger is processed in an in-house external facility, so genomics data is not needed
#  to be copied/shared externally.
#  The output data is exported as follow:
#      - Seurat Object: a RDS file in the OUTDIR.
#      - QC_data : a RData file in the OUTDIR.

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
	make_option(c("--INDIR"), action="store", 
	  default="./data/00_scRNAseq/", type='character',
	  help="Input directory that stores CellRanger's output for each library in an individual directory named with sampleID."),
	make_option(c("--OUTDIR"), action="store", 
	  default="./data/01_Seurat/", type='character',
	  help="Output directory that stores the resulting RDS/RData file for Seurat Object/QC image."),
	make_option(c("--ID"), action="store", 
	  default="GSM4191941", type='character',
	  help="sampleID unique of a library pointing-out to INDIR and'./index/samples/$ID.yaml' file.")
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

# Build path to YAML
YAML <- paste0("./index/samples/",ID,".yaml")
# TAG for the sample subdirectory of OUTDIR
INPUT_TAG <- "_Seurat"

# Sanity check
if(!dir.exists(INDIR)) {
	stop("[ERROR] Input directory does _NOT_ exist.\n")
}
if(!dir.exists(OUTDIR)) {
	cat("[INFO] Creating output directory.\n")
	dir.create(OUTDIR)
}
if(!file.exists(YAML)) {
	stop("[ERROR] Input YAML file does _NOT_ exist.\n")
}

#--- Read index
yaml_fl <- YAML
# Find index file
idx <- yaml::read_yaml(yaml_fl)
# Reformat for easy access to info
idx$seurat_cutoff <- unlist(lapply(idx$seurat_cutoff, function(z) setNames(z[[1]], names(z))))
idx$metadata <- unlist(lapply(idx$metadata, function(z) setNames(z[[1]], names(z))))

#--- Creating Seurat
cat(paste0("Processing ",ID, "\n"), file=stdout())
# Find folder of 10x cellranger count
input_dir <- paste0(INDIR, ID)
names(input_dir) <- ID

dat <- Read10X(data.dir=input_dir)
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
nExp_poi <- (before_nCells/1e3 * 0.0075)  
S <- doubletFinder_v3(S, PCs = 1:nPCs, pN = 0.25, 
			       pK = pK_sel, 
			       nExp = nExp_poi, reuse.pANN = FALSE, 
			       sct = TRUE)
# Rename cols
DF_idx <- grep("^(pANN_|DF\\.)", colnames(S@meta.data))
colnames(S@meta.data)[DF_idx] <- gsub("^(pANN|DF)(_|\\.)?.*","\\1",colnames(S@meta.data)[DF_idx])

##### DEVELOPING STAGE
#NOTE: add plots for QC

#--- Save object
saveRDS(S, paste0(OUTDIR,"/S.rds"))
save(list=c("before_nCells", "after_nCells", "nPCs"))
cat(c(before_nCells, after_nCells), sep="\n", 
    file=paste0(OUTDIR,"nCells_filtering.txt"))
cat(nPCs, sep="\n", 
    file=paste0(OUTDIR,"nPCs.txt"))

#--- Show sessionInfo
sessionInfo()
