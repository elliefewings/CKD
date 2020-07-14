#!/bin/bash

########################################################
######## RAW DATA ######################################
########################################################

## Processed data of both (sc)RNA-seq
if [ ! -e ./raw ];then mkdir ./raw;fi
url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140989/suppl/GSE140989_RAW.tar"
dest_fl="./raw/$(basename $url)"

if [ ! -e $dest_fl ];then wget -P ./raw $url;fi
tar -xf $dest_fl -C ./raw

### Organize single-cell data and bulk RNA-seq
## RNA-seq
if [ ! -e ./raw/RNAseq ];then mkdir ./raw/RNAseq;fi
# The HTseq-count files are txt.gz
find ./raw -name *.txt.gz -exec mv -t ./raw/RNAseq {} +

## Single-Cell data
# 10x files are genes, barcodes, matrix
if [ ! -e ./raw/scRNAseq ];then mkdir ./raw/scRNAseq;fi

# Some files are old-fashion Cellranger output that crashes in Seurat, where
# files are called genes instead of features. So let's rename them
for fl in `find ./raw -name *genes.tsv.gz`;do
	fl2=`echo $fl | sed -e "s/genes/features/g"`;
	mv $fl $fl2
done

# Organize by sample
for fl in `find ./raw -name *matrix.mtx.gz`;do 
	sname=`basename $fl | cut -d"_" -f 1`; 
	sname_dir="./raw/scRNAseq/$sname";
	echo -e "Processing $sname";
	if [ ! -e $sname_dir ];then mkdir $sname_dir;fi
	# Move matrix to sample folder
	find ./raw -name $sname*"matrix.mtx.gz" -exec mv -t $sname_dir {} +
	# Move features.tsv to sample folder
	find ./raw -name $sname*"features.tsv.gz" -exec mv -t $sname_dir {} +
	# Move barcodes to sample folder
	find ./raw -name $sname*"barcodes.tsv.gz" -exec mv -t $sname_dir {} +
done

# Because Seurat just expect files such as "barcodes.tsv.gz" without prefix,
# we create soft links to them
for fl in `find ./raw -name *barcodes.tsv.gz`;do ln -s "./$(basename $fl)" "$(dirname $fl)/barcodes.tsv.gz";done
for fl in `find ./raw -name *features.tsv.gz`;do ln -s "./$(basename $fl)" "$(dirname $fl)/features.tsv.gz";done
for fl in `find ./raw -name *matrix.mtx.gz`;do ln -s "./$(basename $fl)" "$(dirname $fl)/matrix.mtx.gz";done


########################################################
######## PROCESSED DATA ################################
########################################################
if [ ! -e ./data ];then mkdir ./data;fi
## Single-Cell data
url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140989/suppl/GSE140989_PREMIERE_CombinedDataMatrixFile_Normalized.txt.gz"
dest_fl="./data/$(basename $url)"

if [ ! -e $dest_fl ];then wget -P ./data $url;fi


