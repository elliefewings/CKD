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
# Data integration using scanorama for a list of Seurat objects
# =============================================================
#
#  scanorama is a method implemented in python. Herein scanorama is run over 
#  a list of seurat objects in R powered by reticulate.
#  NOTE: Scability might be an issue with increasing number of datasets since 
#  R-python interoperability.

#--- Env
set.seed(1234)

#--- Libraries
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(yaml))
suppressPackageStartupMessages(require(reticulate))

#--- Get input parameters
option_list = list(
	make_option(c("--INPUT"), action="store", 
	  default="data/01_Seurat/GSM4191952_S.rds,data/01_Seurat/GSM4191953_S.rds,data/01_Seurat/GSM4191954_S.rds", 
	  type='character',
	  help="Comma separated string with a list of seurat objects as RDS files."),
	make_option(c("--SEURATOBJ"), action="store", 
	  default="data/02_scanorama/test_S.RDS",
	  type='character',
	  help="Output directory that stores the resulting RDS for Seurat Object."),
	make_option(c("--NPCS"), action="store", 
	  default="./data/02_scanorama/test_nPCs.txt", type='character',
	  help="Output plain-text file with selected N principal components."),
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

# Sanity check
for(fl in sampleFLs) {
	if(!file.exists(fl)) stop(paste0("[ERROR] fl '",fl,"' does not exist\n"))
}
if(!file.exists(YAML)) {
	stop("[ERROR] Input YAML file does _NOT_ exist.\n")
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


#--- Prepare input for scanorama
#NOTE: input must be lists, but names of elements must be not set. Otherwise reticulate export it incorrectly.
assayList <- lapply(SL, function(S) t(as.matrix(GetAssayData(S, "data"))))
names(assayList) <- NULL

geneList <- lapply(SL, function(S) rownames(S))
names(geneList) <- NULL

#--- Run scanorama pipeline via reticulate
# Ack: @gdagstn, github. 
# src: https://github.com/brianhie/scanorama/issues/38#issuecomment-551446738
scanorama <- import("scanorama")

integrated.data <- scanorama$integrate(assayList, geneList)
corrected.data <- scanorama$correct(assayList, geneList, return_dense=TRUE)
integrated.corrected.data <- scanorama$correct(assayList, geneList, return_dimred=TRUE, return_dense=TRUE)

#--- Get back to Seurat
# From the latest object, we extract:
# - corrected counts: element 2 as cells x genes. Thus transposition is required.
# - dimensional reduction embeddings: element 1
# - common genes: element 3, used to name dimmension to the integrated, batch-corrected matrix.
intdata <- lapply(integrated.corrected.data[[2]], t)
panorama <- do.call(cbind, intdata)
rownames(panorama) <- as.character(integrated.corrected.data[[3]])
colnames(panorama) <- unlist(sapply(assayList, rownames))

intdimred <- do.call(rbind, integrated.corrected.data[[1]])
colnames(intdimred) <- paste0("PC_", 1:100)

#Add standard deviations in order to draw Elbow Plots in Seurat
stdevs <- apply(intdimred, MARGIN = 2, FUN = sd)

# Create Seurat object with corrected counts
S <- CreateSeuratObject(counts = panorama, assay = "pano",  project = "CKD")

# Normalization.
#NOTE: SCTransformation crashes out somehow, so we use just the log and scaling
S <- NormalizeData(S)

# Variable feature selection could be skipped since PCA embeddings were calculated
# by scanorama.

#NOTE: make sure that the rownames of your metadata slot 
# are the same as the colnames of your integrated expression matrix
stopifnot(all(unlist(lapply(SL, colnames)) == colnames(S)))

# Adding metadata from all previous objects 
tmp <- lapply(SL, function(x) x@meta.data)
names(tmp) <- NULL
S@meta.data <- do.call(rbind, tmp)
rm(tmp)
stopifnot(all(rownames(S@meta.data)==colnames(S)))
# rownames(S@meta.data) <- colnames(S)

# Adding PCA embeddings
rownames(intdimred) <- colnames(S)
S[["pca"]] <- CreateDimReducObject(embeddings = intdimred, stdev = stdevs, key = "PC_", assay = "pano")

#--- Save object
saveRDS(S, SEURATOBJ)
cat(nPCs, sep="\n", 
    file=NPCS)

#--- Show sessionInfo
sessionInfo()
