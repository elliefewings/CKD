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
# CCA Data integration of Seurat Objects with SCTransform normalization
# =====================================================================
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
	make_option(c("--INPUT"), action="store", 
	  default="data/01_Seurat/GSM4191952_S.rds,data/01_Seurat/GSM4191953_S.rds,data/01_Seurat/GSM4191954_S.rds", 
	  type='character',
	  help="Comma separated string with a list of seurat objects as RDS files."),
	make_option(c("--SEURATOBJ"), action="store", 
	  default="data/02_SCTintegration/test_S.RDS",
	  type='character',
	  help="Output directory that stores the resulting RDS for Seurat Object."),
	make_option(c("--MIN_res"), action="store", 
	  default=0.1, type='numeric',
	  help="Input: Minimum resolution to explore in the cell clustering."),
	make_option(c("--MAX_res"), action="store", 
	  default=2.0, type='numeric',
	  help="Input: Maximum resolution to explore in the cell clustering."),
	make_option(c("--BY_res"), action="store", 
	  default=0.1, type='numeric',
	  help="Input: Increased resolution by each step during exploration of the cell clustering."),
	make_option(c("--ACTUAL_res"), action="store", 
	  default=1.0, type='numeric',
	  help="Input: Chosen resolution to define idents in the cell clustering."),
	make_option(c("--NPCS"), action="store", 
	  default="./data/02_CCAintegration/test_nPCs.txt", type='character',
	  help="Output plain-text file with selected N principal components."),
	make_option(c("--CLUST"), action="store", 
	  default="./data/02_CCAintegration/test_clust.txt", type='character',
	  help="Output plain-text file with cell clustering outcome among resolutions."),
	make_option(c("--YAML"), action="store", 
	  default="index/contrasts/LDvsTN.yaml", type='character',
	  help="YAML file with additional metadata for the data integration, if any.")
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

# Seurat objects
sampleFLs <- strsplit(INPUT, split=",")[[1]]
names(sampleFLs) <- gsub("_S.rds","", basename(sampleFLs))

# Ranging of resolutions to explore
EXPLORE_res <- seq(from=MIN_res, to=MAX_res, by=BY_res)

# Sanity check
for(fl in sampleFLs) {
	if(!file.exists(fl)) stop(paste0("[ERROR] fl '",fl,"' does not exist\n"))
}
if(!file.exists(YAML)) {
	stop("[ERROR] Input YAML file does _NOT_ exist.\n")
}
if(!ACTUAL_res %in% EXPLORE_res) {
	cat("[WARN] INPUT ACTUAL_res is not present in the ranging values of resolution exploration.\n", file=stdout())
	cat("Adding ACTUAL_res to the ranging to explore.\n", file=stdout())

	EXPLORE_res <- c(EXPLORE_res, ACTUAL_res)
}


#--- Read index
yaml_fl <- YAML
stopifnot(file.exists(yaml_fl))
#NOTE: DEV : somehow link to metadata from contrasts->by group
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
} else if (grepl("group", yaml_fl)) {
 #NOTE: to be dev
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

#--- Cell Clustering
# Shared-Nearest Neighbour with Graph partitioning, Louvain algorithm (default)
cat(paste0("[INFO] Cell clustering of sample '", Project(S),
	   "' on assay '", DefaultAssay(S), "'\n"), file=stdout())
S <- FindNeighbors(S, dims = 1:nPCs)
S <- FindClusters(S, resolution = EXPLORE_res)

# Rename resolutions as ending with 2 digits for a perfect match
S <- ColRenameSNN(S)

# Actual colname based on input
ACTUAL_col <- paste0(DefaultAssay(S),
		     "_snn_res.", 
		     formatC(ACTUAL_res, digits=1, format="f"))

# Sanity check
if(!ACTUAL_col %in% colnames(S@meta.data)) {
	cat("[WARN]: It was not found the outcome of the actual chosen resolution among metadata.\n", file=stdout())
	cat(paste0("Actual resolution column name: ",ACTUAL_col,".\n"), file=stdout())
	cat(paste0("Possible columns with no match: ",
		   paste(colnames(S@meta.data), collapse="\n"),
		   ".\n"), 
	    file=stdout())
	stop("[ERROR] Not possible to match resolutions.\n")
}

# Set actual chosen resolution
Idents(S) <- S@meta.data[, ACTUAL_col]
S@meta.data$seurat_clusters <- S@meta.data[, ACTUAL_col]

#--- Save object
saveRDS(S, SEURATOBJ)
cat(nPCs, sep="\n", 
    file=NPCS)
write.table(S@meta.data[,c("seurat_clusters", 
			   grep(paste0("^", DefaultAssay(S), "_snn_res"),
				colnames(S@meta.data), value=TRUE)
			   )
			],
            file=CLUST,
            sep=",", col.names = NA, row.names=TRUE, quote=TRUE)

#--- Show sessionInfo
sessionInfo()
