
import os
import os.path
import re
import glob
import yaml

DATASETS = glob.glob("index/samples/*.yaml")
SIDS = [os.path.splitext(os.path.basename(k))[0] for k in DATASETS]

CONTRASTS = glob.glob("index/contrasts/*.yaml")
CIDS = [os.path.splitext(os.path.basename(k))[0] for k in CONTRASTS]

### All
rule all:
	input:
		expand("data/01_Seurat/{id}_{fl}", id=SIDS, fl=["S.rds", "nPCs.txt", "nCells.txt"]),
		expand("data/02_CCAintegration/{id}_{fl}", id=CIDS, fl=["S.rds", "nPCs.txt"])


### Single-Sample Seurat
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

def getSeuratClust(wildcards):
	YAML = glob.glob("index/*/"+wildcards.id+".yaml")[0]
	yaml_fl = open(YAML)
	parsed_yaml = yaml.load(yaml_fl, Loader=yaml.FullLoader)
	seurat_clust = parsed_yaml["seurat_clust"]

	if ['MIN_res', 'MAX_res', 'BY_res', 'ACTUAL_res'] == list(seurat_clust.keys()):
		seurat_clust = list(seurat_clust.values())
	else:
		sys.exit("ERROR: clustering format is not correct")

	return seurat_clust

rule seurat:
	input:
		expand("data/01_Seurat/{id}_{fl}", id=SIDS, fl=["S.rds", "nPCs.txt", "nCells.txt"])

rule run_seurat:
	input:
		script = "scripts/create_seurat.R",
		CR = getCellRanger,
		YAML = lambda wildcards: "index/samples/"+wildcards.id+".yaml"
	output:
		"data/01_Seurat/{id}_S.rds",
		"data/01_Seurat/{id}_nPCs.txt",
		"data/01_Seurat/{id}_nCells.txt",
		"data/01_Seurat/{id}_clust.csv"
	
	params:
		# Activate existing conda env
		conda_env = "envs/01_Seurat",
		# Get seurat params for cell filtering from sample YAML file
		seurat_cutoff = getSeuratParams,
		# Get seurat params for cell clustering from sample YAML file
		seurat_clust = getSeuratClust

	message:
		"Building single-sample Seurat for {wildcards.id}"
	shell:
         	'set +eu '
 		" && (test -d data/01_Seurat/ || mkdir data/01_Seurat/) "
		' && . $(conda info --base)/etc/profile.d/conda.sh '
		' && conda activate {params.conda_env} '
		"&& $CONDA_PREFIX/bin/Rscript {input.script}"
		# Inputs
		" --MATRIX {input.CR[0]} --FEATURES {input.CR[1]} --BARCODES {input.CR[2]}"
		" --MAX_percMT {params.seurat_cutoff[0]} --MIN_nFeatures {params.seurat_cutoff[1]} --MAX_nFeatures {params.seurat_cutoff[2]}"
		" --MIN_res {params.seurat_clust[0]} --MAX_res {params.seurat_clust[1]} --BY_res {params.seurat_clust[2]} --ACTUAL_res {params.seurat_clust[3]}"
		" --YAML {input.YAML}"
		# Outputs
		" --SEURATOBJ {output[0]} --NPCS {output[1]} --NCELLS {output[2]} --CLUST {output[3]}"

### Data Integration
def getIDSfromGROUP(yaml_fl2):
	parsed_yaml2 = yaml.load(yaml_fl2, Loader=yaml.FullLoader)
	gr_ids = parsed_yaml2.values()
	gr_ids = list(gr_ids)
	if len(gr_ids) == 1:
		gr_ids = gr_ids[0]
	return(gr_ids)

def getSeuratObjects(wildcards):
	yaml_flname = "index/contrasts/"+wildcards.id+".yaml"
	yaml_fl = open(yaml_flname)
	parsed_yaml = yaml.load(yaml_fl, Loader=yaml.FullLoader)

	IDS_tmp = list()
	if "contrasts" in yaml_flname:
		grs = parsed_yaml["groups"]
		if len(grs) == 1:
			grs = grs[0]

		for gr in grs:
			yaml_fl2 = open("index/groups/"+gr+".yaml")
			gr_ids = getIDSfromGROUP(yaml_fl2)
			IDS_tmp.extend(gr_ids)
	
	if "groups" in yaml_fl:
		gr_ids = getIDSfromGROUP(yaml_fl)
		IDS_tmp = gr_ids

	SL = ["data/01_Seurat/"+k+"_S.rds" for k in IDS_tmp]
	return(SL)

rule cca:
	input:
		expand("data/02_CCAintegration/{id}_{fl}", id=CIDS, fl=["S.rds", "nPCs.txt"])

rule run_cca:
	input:
		script = "scripts/CCA_seurat.R",
		SL = getSeuratObjects,
		YAML = lambda wildcards: "index/contrasts/"+wildcards.id+".yaml"
	output:
		"data/02_CCAintegration/{id}_S.rds",
		"data/02_CCAintegration/{id}_nPCs.txt"
	params:
		# Activate existing conda env
		conda_env = "envs/01_Seurat",
		# Get seurat params for cell clustering from sample YAML file
		seurat_clust = getSeuratClust
	message: "Data integration for contrast {wildcards.id} using CCA SCT seurat"
	run:
		seurats = ",".join(input.SL)
		shell('set +eu '\
		" && (test -d data/02_CCAintegration/ || mkdir data/02_CCAintegration/)"\
		' && . $(conda info --base)/etc/profile.d/conda.sh '\
		' && conda activate {params.conda_env} '\
		"&& $CONDA_PREFIX/bin/Rscript {input.script} --INPUT {seurats} --SEURATOBJ {output[0]} --NPCS {output[1]}" 
		" --MIN_res {params.seurat_clust[0]} --MAX_res {params.seurat_clust[1]} --BY_res {params.seurat_clust[2]} --ACTUAL_res {params.seurat_clust[3]}"
		" --YAML {input.YAML}")

rule pano:
	input:
		expand("data/02_scanorama/{id}_{fl}", id=CIDS, fl=["S.rds", "nPCs.txt"])

rule run_pano:
	input:
		script = "scripts/scanorama_seurat.R",
		SL = getSeuratObjects,
		YAML = lambda wildcards: "index/contrasts/"+wildcards.id+".yaml"
	output:
		"data/02_scanorama/{id}_S.rds",
		"data/02_scanorama/{id}_nPCs.txt"
	params:
		# Activate existing conda env
		conda_env = "envs/01_Seurat",
	message: "Data integration for contrast {wildcards.id} using SCANORAMA"
	run:
		seurats = ",".join(input.SL)
		shell('set +eu '\
		" && (test -d data/02_scanorama/ || mkdir data/02_scanorama/)"\
		' && . $(conda info --base)/etc/profile.d/conda.sh '
		' && conda activate {params.conda_env} '
		"&& $CONDA_PREFIX/bin/Rscript {input.script} --INPUT {seurats} --SEURATOBJ {output[0]} --NPCS {output[1]}" 
		" --MIN_res {params.seurat_clust[0]} --MAX_res {params.seurat_clust[1]} --BY_res {params.seurat_clust[2]} --ACTUAL_res {params.seurat_clust[3]}"
		" --YAML {input.YAML}")
