
# Samples
YAML = $(shell find ../index/samples -mindepth 1 -name "*.yaml")
ID = $(basename $(notdir $(YAML)))
SEURAT = $(ID:%=01_Seurat/%_Seurat)
	
# Contrasts
YAML2 = $(shell find ../index/contrasts -mindepth 1 -name "*.yaml")
ID2 = $(basename $(notdir $(YAML2)))
SCT = $(ID2:%=02_SCTintegration/%_SCT)

all: seurat sct logs

clean:
	rm -rf logs
	rm -rf 01_Seurat

#--- Single-Sample Seurat runs
seurat: $(SEURAT)

01_Seurat/%_Seurat: 01_Seurat/01_Seurat.Rmd  01_Seurat 
	qsub -q short -l nodes=1:ppn=16,mem=32gb -v RMD=$<,PREFIX=$@,SAMPLEID=$(@:01_Seurat/%_Seurat=%) -e logs/$(@:01_Seurat/%=%).err -o logs/$(@:01_Seurat/%=%).txt ../src/PBSrenderIdx.sh 

01_Seurat/01_Seurat.Rmd: 01_Seurat.Rmd 01_Seurat
	test -s $@ || ln -s $(realpath 01_Seurat.Rmd) $@

01_Seurat: logs
	test -d $@ || mkdir $@

#--- Multi-Sample Seurat integration based on contrasts
sct: $(SCT)

02_SCTintegration/%_SCT: 02_SCTintegration/02_SCTintegration.Rmd  02_SCTintegration 
	qsub -q short -l nodes=1:ppn=16,mem=32gb -v RMD=$<,PREFIX=$@,SAMPLEID=$(@:02_SCTintegration/%_SCT=%) -e logs/$(@:02_SCTintegration/%=%).err -o logs/$(@:02_SCTintegration/%=%).txt ../src/PBSrenderIdx.sh 

02_SCTintegration/02_SCTintegration.Rmd: 02_SCTintegration.Rmd 02_SCTintegration
	test -s $@ || ln -s $(realpath 02_SCTintegration.Rmd) $@

02_SCTintegration: logs
	test -d $@ || mkdir $@

#--- Just the logs
logs:
	mkdir -p $@

