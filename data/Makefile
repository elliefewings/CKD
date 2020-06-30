
YAML = $(shell find ../index -mindepth 2 -name "*.yaml")
ID = $(basename $(notdir $(YAML)))
SEURAT = $(ID:%=%_Seurat)

all: seurat

clean:
#	rm -f $(SEURAT:%=/01_Seurat/%.md) $(SEURAT:%=/01_Seurat/%.html) 
	rm -rf 01_Seurat

seurat: $(SEURAT)

%_Seurat: 01_Seurat.Rmd  01_Seurat
	${CONDA_PREFIX}/bin/R -e "rmarkdown::render('$<', output_dir='$(filter-out $<,$^)', output_file='01_Seurat/$@.md', params=list(ID='$(@:%_01_Seurat=%)',OUTDIR='$(filter-out $<,$^)'))"

01_Seurat:
	mkdir $@

