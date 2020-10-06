
import os
import os.path
import re
import glob

DATASETS = glob.glob("index/samples/*.yaml")
IDS = [os.path.splitext(os.path.basename(k))[0] for k in DATASETS]

def getCellRanger(wildcards):
	fls = [ "data/00_scRNAseq/"+wildcards.id+"/matrix.mtx.gz",
			"data/00_scRNAseq/"+wildcards.id+"/features.tsv.gz",
			"data/00_scRNAseq/"+wildcards.id+"/barcodes.tsv.gz"]
	return fls

rule seurat:
	input:
		expand("data/01_Seurat/{id}_{fl}", id=IDS, fl=["S.rds", "nPCs.txt", "nCells.txt"])

rule run_seurat:
	input:
		CR=getCellRanger,
	#	"data/00_scRNAseq/{id}/matrix.mtx.gz"
	#	"data/00_scRNAseq/{id}/features.tsv.gz"
	#	"data/00_scRNAseq/{id}/barcodes.tsv.gz"
		YAML = lambda wildcards: "index/samples/"+wildcards.id+".yaml"
	output:
		"data/01_Seurat/{id}_S.rds",
		"data/01_Seurat/{id}_nPCs.txt",
		"data/01_Seurat/{id}_nCells.txt"
	
	params:
		conda_env = "envs/01_Seurat"

	shell:
         	'set +eu '
 		" && (test -d data/01_Seurat/ || mkdir data/01_Seurat/) "
		' && . $(conda info --base)/etc/profile.d/conda.sh '
		' && conda activate {params.conda_env} '
		"&& $CONDA_PREFIX/bin/Rscript scripts/create_seurat.R --MATRIX {input.CR[0]} --FEATURES {input.CR[1]} --BARCODES {input.CR[2]} --YAML {input.YAML} --SEURATOBJ {output[0]} --NPCS {output[1]} --NCELLS {output[2]}"
# 		"(test -d data/01_Seurat/ || mkdir data/01_Seurat/) "
# 		"&& touch {output[0]} "
# 		"&& touch {output[1]} "
# 		"&& touch {output[2]}"
