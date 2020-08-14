
YAML = $(shell find ../index/samples -mindepth 1 -name "*.yaml")
ID = $(basename $(notdir $(YAML)))
SEURAT = $(ID:%=01_Seurat/%_Seurat)

all: seurat logs

clean:
	rm -rf logs
	rm -rf 01_Seurat

seurat: $(SEURAT)

01_Seurat/%_Seurat: 01_Seurat/01_Seurat.Rmd  01_Seurat 
	qsub -q short -l nodes=1:ppn=16,mem=32gb -v RMD=$<,PREFIX=$@,SAMPLEID=$(@:01_Seurat/%_Seurat=%) -e logs/$(@:01_Seurat/%=%).err -o logs/$(@:01_Seurat/%=%).txt PBSrender.sh 

01_Seurat/01_Seurat.Rmd: 01_Seurat.Rmd 01_Seurat
	test -s $@ || ln -s $(realpath 01_Seurat.Rmd) $@

01_Seurat: logs
	test -d $@ || mkdir $@

logs:
	mkdir -p $@

