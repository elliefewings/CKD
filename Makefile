
# Samples
YAML = $(shell find ./index/samples -mindepth 1 -name "*.yaml")
ID = $(basename $(notdir $(YAML)))
SEURAT = $(ID:%=data/01_Seurat/%_Seurat)
	
# Contrasts
YAML2 = $(shell find ./index/contrasts -mindepth 1 -name "*.yaml")
ID2 = $(basename $(notdir $(YAML2)))
SCT = $(ID2:%=data/02_SCTintegration/%_SCT)
PAN = $(ID2:%=data/02_scanorama/%_PAN)

all: seurat sct logs

clean:
	test -d logs && rm -f logs
	test -s data/01_Seurat && rm -rf data/01_Seurat
	test -s data/02_SCTintegration && rm -rf data/02_SCTintegration
	test -s data/02_scanorama && rm -rf data/02_scanorama

#--- Single-Sample Seurat runs
seurat: $(SEURAT)

data/01_Seurat/%_Seurat: scripts/create_seurat.R data/01_Seurat
	qsub -q short -l nodes=1:ppn=16,mem=32gb -v SRP=$<,OUTDIR=$@,SAMPLEID=$(@:data/01_Seurat/%_Seurat=%),INDIR=data/00_scRNAseq/ -e logs/$(@:data/01_Seurat/%=%).err -o logs/$(@:data/01_Seurat/%=%).txt ./src/PBS-RscriptIdx.sh 

data/01_Seurat: logs
	test -d $@ || mkdir $@

#--- Multi-Sample Seurat integration based on contrasts
#--- CCA integration with SCTransformation normalization powered by Seurat
sct: $(SCT)

data/02_SCTintegration/%_SCT: scripts/SCTintegration_seurat.R data/02_SCTintegration 
	qsub -q short -l nodes=1:ppn=16,mem=32gb -v SRP=$<,OUTDIR=$@,SAMPLEID=$(@:data/02_SCTintegration/%_SCT=%),INDIR=data/01_Seurat/ -e logs/$(@:data/02_SCTintegration/%=%).err -o logs/$(@:data/02_SCTintegration/%=%).txt ./src/PBS-RscriptIdx.sh 

data/02_SCTintegration: logs
	test -d $@ || mkdir $@


#--- Scanorama integration with SCTransformation normalization
pan: $(PAN)

data/02_scanorama/%_PAN: scripts/scanorama_seurat.R data/02_scanorama
	qsub -q short -l nodes=1:ppn=16,mem=32gb -v SRP=$<,OUTDIR=$@,SAMPLEID=$(@:data/02_scanorama/%_SCT=%),INDIR=data/01_Seurat/ -e logs/$(@:data/02_scanorama/%=%).err -o logs/$(@:data/02_scanorama/%=%).txt ./src/PBS-RscriptIdx.sh 

data/02_scanorama: logs
	test -d $@ || mkdir $@

#--- Just the logs
logs:
	mkdir -p $@

