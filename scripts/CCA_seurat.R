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
# Data integration of Seurat Objects with SCTransform normalization
# =================================================================
#
#  It performs data integration of multiple samples from Seurat Objects with
#  SCTransform normalization. The input is a YAML pointing-out samples to the 
#  indicated directory (INDIR). The output is path to a directory to store
#  another Seurat object with the integrated data.

#--- Env
set.seed(1234)

#--- Libraries
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(yaml))
source("./src/seurat_fx.R")
# src: https://stackoverflow.com/questions/40536067/how-to-adjust-future-global-maxsize-in-r
# 6 globals that need to be exported for the future expression (‘FUN()’) is 2.54 GiB
# This exceeds the maximum allowed size of 500.00 MiB
options(future.globals.maxSize= (1024*3)*1024^2)

#--- Get input parameters
option_list = list(
	make_option(c("--INDIR"), action="store", 
	  default="./data/01_Seurat/", type='character',
	  help="Input directory that stores Seurat Objects."),
	make_option(c("--OUTDIR"), action="store", 
	  default="./data/02_SCTintegration/", type='character',
	  help="Output directory that stores the resulting RDS for Seurat Object."),
	make_option(c("--ID"), action="store", 
	  default="LDvsTN", type='character',
	  help="group/contrasti ID pointing-out to INDIR and './index/*/$ID.yaml' file.")
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
YAML <- list.files("./index", pattern=paste0(ID,".yaml"), recursive=TRUE, full.names=TRUE)

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
stopifnot(file.exists(yaml_fl))
# Find index file
idx <- yaml::read_yaml(yaml_fl)
if(grepl("contrast", yaml_fl)) {
	# Get groups
	GRx_idx <- setNames(vector("list", length(idx$groups)),
			    idx$groups)
	for(gr in idx$groups) {
		yaml_fl <- paste0("./index/groups/",gr,".yaml")
		stopifnot(file.exists(yaml_fl))
		GRx_idx[[gr]] <- yaml::read_yaml(yaml_fl)$samples
	}
	# Get for samples
	sampleFLs <- vector("character", length=length(unlist(GRx_idx)))
	names(sampleFLs) <- unlist(GRx_idx)
	for(sampleID in unlist(GRx_idx)) {
		S_rds <- paste0(INDIR,"/",sampleID,INPUT_TAG,"/S.rds")
		# 	stopifnot(file.exists(S_rds))
		sampleFLs[sampleID] <- ifelse(file.exists(S_rds), S_rds, NA)
	}
} else if (grepl("group", yaml_fl)) {
 #NOTE: to be dev
}

# Sanity check in detail for individual samples for the contrast
if(any(is.na(sampleFLs))) {
	cat("[ERROR] : Certain samples that belong to any of the groups do not have Seurat objects:\n", file=stderr())
	cat(paste0("[ERROR] : DO exist:", paste0(na.omit(sampleFLs), collapse=", "), "\n"), file=stderr())
	cat(paste0("[ERROR] : do NOT exist:", paste0(names(sampleFLs)[is.na(sampleFLs)], collapse=", "), "\n"), file=stderr())
	stop("Please, generate the samples that do NOT exist and revise indexes for groups and contrasts.\n")
}

#--- Load samples
# Create list of seurat objects
SL <- setNames(vector("list", length=length(sampleFLs)), 
	       names(sampleFLs))
# For every sample, add it to the list
for(sname in names(sampleFLs)) {
	sname_fl <- sampleFLs[sname]
	cat(paste0("[INFO] loading SeuratObject for ",sname,"\n"),file=stdout())
	S <- readRDS(sname_fl)
	cat(paste0("[INFO] Add ",sname," to list of SeuratObjects\n"),file=stdout())
	SL[[sname]] <- S
	rm(S)
}

#--- Integration
features <- SelectIntegrationFeatures(object.list = SL, nfeatures = 3000)
SL <- PrepSCTIntegration(object.list = SL, anchor.features = features, verbose = FALSE)

anchors <- FindIntegrationAnchors(object.list = SL, normalization.method = "SCT", anchor.features=features, verbose=FALSE)
S <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
rm(SL)

#--- High dimensionality reduction
DefaultAssay(S) <- "integrated"

S <- ScaleData(S, verbose=FALSE)
S <- RunPCA(S, npcs=50, verbose= FALSE)
calcNPCs <- get_npcs(S)
(nPCs <- calcNPCs$npcs)
ElbowPlot(S, ndims=50) + geom_vline(xintercept = nPCs)

S <- RunUMAP(S, reduction = "pca", dims = 1:nPCs)

#--- Save object
saveRDS(S, paste0(OUTDIR,"/S.rds"))
cat(nPCs, sep="\n", 
    file=paste0(OUTDIR,"nPCs.txt"))

#--- Show sessionInfo
sessionInfo()
