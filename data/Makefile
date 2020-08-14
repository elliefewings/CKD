
YAML = $(shell find ../index/samples -mindepth 1 -name "*.yaml")
ID = $(basename $(notdir $(YAML)))
SEURAT = $(ID:%=%_Seurat)

all: seurat logs

clean:
	rm -rf 01_Seurat

seurat: $(SEURAT)

%_Seurat: 01_Seurat.Rmd  01_Seurat 
#	${CONDA_PREFIX}/bin/R -e "rmarkdown::render('$<', output_file='$@.md', params=list(ID='$(@:%_Seurat=%)',OUTDIR='$@'))"
#	mv $@ $(filter-out $<,$^)
#	mv $@.md $(filter-out $<,$^)
#	mv $@.html $(filter-out $<,$^)
	qsub -q short -l nodes=1:ppn=16,mem=32gb -v RMD=$<,PREFIX=$@,SAMPLEID=$(@:%_Seurat=%),DEST=$(filter-out $<,$^) -e logs/$@.err -o logs/$@.txt PBSrender.sh 

01_Seurat: logs
	mkdir -p $@

logs:
	mkdir -p $@

