
import os
import os.path
import re
import glob
import yaml

DATASETS = glob.glob("index/samples/*.yaml")
IDS = [os.path.splitext(os.path.basename(k))[0] for k in DATASETS]

def getCellRanger(wildcards):
	fls = [ "data/00_scRNAseq/"+wildcards.id+"/matrix.mtx.gz",
			"data/00_scRNAseq/"+wildcards.id+"/features.tsv.gz",
			"data/00_scRNAseq/"+wildcards.id+"/barcodes.tsv.gz"]
	return fls

def getSeuratParams(wildcards):
	yaml_fl = open("index/samples/"+wildcards.id+".yaml")
	parsed_yaml = yaml.load(yaml_fl, Loader=yaml.FullLoader)
	seurat_cutoff = parsed_yaml["seurat_cutoff"]

	if ['MAX_percMT', 'MIN_nFeatures', 'MAX_nFeatures'] == list(seurat_cutoff.keys()):
		seurat_cutoff = list(seurat_cutoff.values())
	else:
		sys.exit("ERROR: cutoff format is not correct")

	return seurat_cutoff

rule seurat:
	input:
		expand("data/01_Seurat/{id}_{fl}", id=IDS, fl=["S.rds", "nPCs.txt", "nCells.txt"])

rule run_seurat:
	input:
		script = "scripts/create_seurat.R"
		CR=getCellRanger,
		YAML = lambda wildcards: "index/samples/"+wildcards.id+".yaml"
	output:
		"data/01_Seurat/{id}_S.rds",
		"data/01_Seurat/{id}_nPCs.txt",
		"data/01_Seurat/{id}_nCells.txt"
	
	params:
		# Activate existing conda env
		conda_env = "envs/01_Seurat",
		# Get seurat params from sample YAML file
		seurat_cutoff = getSeuratParams

	shell:
         	'set +eu '
 		" && (test -d data/01_Seurat/ || mkdir data/01_Seurat/) "
		' && . $(conda info --base)/etc/profile.d/conda.sh '
		' && conda activate {params.conda_env} '
		"&& $CONDA_PREFIX/bin/Rscript {input.script}"
		" --MATRIX {input.CR[0]} --FEATURES {input.CR[1]} --BARCODES {input.CR[2]}"
		" --MAX_percMT {params.seurat_cutoff[0]} --MIN_nFeatures {params.seurat_cutoff[1]} --MAX_nFeatures {params.seurat_cutoff[2]}"
		" --YAML {input.YAML}"
		" --SEURATOBJ {output[0]} --NPCS {output[1]} --NCELLS {output[2]}"
