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

#E# May want to add remotes::install_github('chris-mcginnis-ucsf/DoubletFinder') somewhere to documentation
suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(require(yaml))
suppressPackageStartupMessages(require(magrittr))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(clustree))
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
	make_option(c("--MAX_percMT"), action="store", 
	  default=1, type='numeric',
	  help="Input: Maximum percentage of mitochondrial expression genes to remain after filtering-out"),
	make_option(c("--MIN_nFeatures"), action="store", 
	  default=1, type='numeric',
	  help="Input: Minimum number of Features per cell to remain after filtering-out."),
	make_option(c("--MAX_nFeatures"), action="store", 
	  default=1, type='numeric',
	  help="Input: Maximum number of Features per cell to remain after filtering-out."),
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
	make_option(c("--SEURATOBJ"), action="store", 
	  default="./data/01_Seurat/GSM4191941_S.rds", type='character',
	  help="Output file to store the Seurat Object as RDS."),
	make_option(c("--REPORTOBJ"), action="store", 
	  default="./data/01a_Report/GSM4191941_R.rds", type='character',
	  help="Output file to store the Report Object as RDS."),
	make_option(c("--NPCS"), action="store", 
	  default="./data/01_Seurat/GSM4191941_nPCs.txt", type='character',
	  help="Output plain-text file with selected N principal components."),
	make_option(c("--NCELLS"), action="store", 
	  default="./data/01_Seurat/GSM4191941_nCells.txt", type='character',
	  help="Output plain-text file with recorded N cells before and after filtering."),
	make_option(c("--CLUST"), action="store", 
	  default="./data/01_Seurat/GSM4191941_clust.txt", type='character',
	  help="Output plain-text file with cell clustering outcome among resolutions."),
	make_option(c("--YAML"), action="store", 
	  default="./index/samples/GSM4191941.yaml", type='character',
	  help="YAML file with settings for the sample.")
)

# Parse the parameters
opt = parse_args(OptionParser(option_list=option_list))

# Create report object
report <- c()

#Add parameters to report
report$opt <- opt

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
out1a <- dirname(REPORTOBJ)
out2 <- dirname(NCELLS)
out3 <- dirname(NPCS)
stopifnot(out1==out2 && out1 == out3)
OUTDIR <- out1
REPORTDIR <- out1a
rm(out1, out2, out3)

# Ranging of resolutions to explore
EXPLORE_res <- seq(from=MIN_res, to=MAX_res, by=BY_res)

# Sanity check
if(!dir.exists(INDIR)) {
	stop("[ERROR] Input directory does _NOT_ exist.\n")
}
if(!dir.exists(OUTDIR)) {
	cat("[INFO] Creating output directory.\n", file=stdout())
	dir.create(OUTDIR)
}
if(!dir.exists(REPORTDIR)) {
  cat("[INFO] Creating report directory.\n", file=stdout())
  dir.create(REPORTDIR)
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
# Read index file
idx <- yaml::read_yaml(YAML)
# Reformat for easy access to info
idx$seurat_cutoff <- unlist(lapply(idx$seurat_cutoff, function(z) setNames(z[[1]], names(z))))
idx$metadata <- unlist(lapply(idx$metadata, function(z) setNames(z[[1]], names(z))))

#--- Creating Seurat
# Getting ID of sample
# It is a standard practice to name the prefix of YAML as the sample ID.
ID <- strsplit(basename(YAML), split="\\.")[[1]][1]
cat(paste0("Processing ",ID, "\n"), file=stdout())
# Find folder of 10x cellranger count
names(INDIR) <- ID

dat <- Read10X(data.dir=INDIR)
S <- CreateSeuratObject(counts = dat, project = ID, 
			min.cells=5, 
			min.features=idx$seurat_cutoff["MIN_nFeatures"])
S[["percent.mt"]] <- PercentageFeatureSet(S, pattern = "^MT-")

#--- Create summary table and figures of meta data pre-filtering
# Subset meta data
report$data.meta.summ <- summarise(S@meta.data, ncells = length(orig.ident),
                            med_nCount_RNA = median(nCount_RNA),
                            min_nCount_RNA = min(nCount_RNA),
                            max_nCount_RNA = max(nCount_RNA),
                            med_nFeature_RNA = median(nFeature_RNA),
                            min_nFeature_RNA = min(nFeature_RNA),
                            max_nFeature_RNA = max(nFeature_RNA),
                            med_percent.mt = median(percent.mt),
                            min_percent.mt = min(percent.mt),
                            max_percent.mt = max(percent.mt)) %>% t() %>% as.data.frame()

colnames(report$data.meta.summ) <- "pre-filtering"

# Create QC plots
report$qc1 <- VlnPlot(S, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, cols="grey")

report$qc1.1 <- report$qc1[[1]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank()) + geom_hline(yintercept=MAX_nFeatures, color="red") + geom_hline(yintercept=MIN_nFeatures, color="blue")
report$qc1.2 <- report$qc1[[2]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
report$qc1.3 <- report$qc1[[3]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank()) + geom_hline(yintercept=MAX_percMT, color="red")

report$qc2 <- FeatureScatter(S, feature1 = "nCount_RNA", feature2 = "percent.mt", cols="steelblue4")

report$qc3 <- FeatureScatter(S, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols="steelblue4")


#--- Filter-out
# Before
before_nCells <- ncol(S)
# apply Filter
S <- subset(S, subset= nFeature_RNA > idx$seurat_cutoff["MIN_nFeatures"] & 
	    		nFeature_RNA < idx$seurat_cutoff["MAX_nFeatures"] &
	    percent.mt < idx$seurat_cutoff["MAX_percMT"])
# After
after_nCells <- ncol(S)

#--- Create summary table of meta data post-filtering
report$data.meta.summ$'post-filtering' <- summarise(S@meta.data, ncells = length(orig.ident),
                                             med_nCount_RNA = median(nCount_RNA),
                                             min_nCount_RNA = min(nCount_RNA),
                                             max_nCount_RNA = max(nCount_RNA),
                                             med_nFeature_RNA = median(nFeature_RNA),
                                             min_nFeature_RNA = min(nFeature_RNA),
                                             max_nFeature_RNA = max(nFeature_RNA),
                                             med_percent.mt = median(percent.mt),
                                             min_percent.mt = min(percent.mt),
                                             max_percent.mt = max(percent.mt)) %>% t()

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

# Get npc elbow plot for report
report$npcs <- calcNPCs
# Get pca plot for report
report$pca <- DimPlot(S, reduction = "pca")

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

#--- Cell Clustering
# Shared-Nearest Neighbour with Graph partitioning, Louvain algorithm (default)
cat(paste0("[INFO] Cell clustering of sample '", Project(S),
	   "' on assay '", DefaultAssay(S), "'\n"), file=stdout())
S <- FindNeighbors(S, dims = 1:nPCs)
S <- FindClusters(S, resolution = EXPLORE_res)

report$clust <- suppressWarnings(clustree(S)) + theme(
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(size=5)))

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
cat(c(before_nCells, after_nCells), sep="\n", file=NCELLS)
cat(nPCs, sep="\n", file=NPCS)
write.table(S@meta.data[,c("seurat_clusters", 
			   grep(paste0("^", DefaultAssay(S), "_snn_res"), 
				colnames(S@meta.data), value=TRUE)
			   )
			],
            file=CLUST,
            sep=",", col.names = NA, row.names=TRUE, quote=TRUE)

#--- Save object for report
saveRDS(report, REPORTOBJ)

#--- Show sessionInfo
sessionInfo()
