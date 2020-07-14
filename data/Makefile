
YAML = $(shell find ../index/samples -mindepth 1 -name "*.yaml")
ID = $(basename $(notdir $(YAML)))
SEURAT = $(ID:%=%_Seurat)

all: seurat

clean:
	rm -rf 01_Seurat

seurat: $(SEURAT)

%_Seurat: 01_Seurat.Rmd  01_Seurat
	${CONDA_PREFIX}/bin/R -e "rmarkdown::render('$<', output_file='$@.md', params=list(ID='$(@:%_Seurat=%)',OUTDIR='$@'))"
	mv $@ $(filter-out $<,$^)
	mv $@.md $(filter-out $<,$^)
	mv $@.html $(filter-out $<,$^)

01_Seurat:
	mkdir $@

