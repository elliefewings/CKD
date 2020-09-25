# The transcriptome regulation of Chronic Kidney Disease at single-cell resolution
> Maintainers:  
 Javier Perales-Patón - javier.perales@bioquant.uni-heidelberg.de - ORCID: 0000-0003-0780-6683  

> A joint project of [Saez Lab](http://www.saezlab.org), 
[Kramann Lab](http://www.kramannlab.com) and 
[Stegle lab](https://www.dkfz.de/en/bioinformatik-genomik-systemgenetik/). Project funded by _XXX_.

## Abstract 
Chronic kidney disease (CKD) is a heterogeneous group of diseses that gradually impair kidney function 
over a long period of time. There are different diseases and conditions causing CKD such as diabetes, 
high blood preasure, glomerulonephritis, polycystic kidney, kidney infections, or prolongued 
obstruction of urinary tract. We have profiled the transcriptome and epigenome of individual cells from 
a large cohort of CKD patients to characterize the molecular factors driving the disease. Patient 
clinical data will be used to find associations with disease outcome. A complete characterization of 
transcriptome regulation in CKD will be done.

## Structure

| Folder  | Description                                                                                    |
| :---    | :----                                                                                          |
| raw     | Raw seq data. It is kept in the local for preprocessing - never pushed.                        |
| data    | Processed data for data analysis, i.e. count matrices, objects, gene sets collections, etc.    |
| scripts | Main scripts                                                                                   |
| ext     | External (published) datasets that has been re-analysed incl DOI or source.                    |
| envs    | Docker container / conda environments and install instructions                                 |
| img     | Main images for the repository/project, i.e. graphical abstract, workflow, etc.                |
| src     | Useful auxiliary functions to be loaded in different scripts                                   |
| utils   | External scripts, packages and github submodules used in the project                           |
| doc     | Material and methods section, essential documentation, access the cluster, conda environments. |
| index   | YAML files for individual samples describing parameters of analysis                            |

```bash
# Clone repository
git clone https://github.com/saezlab/CKD.git

# Setup repo structure from root if any folder does not exist already
arr=(raw data scripts ext envs img src utils doc index);
for i in "${arr[@]}";do if [ ! -e $i ];then mkdir $i; fi; done
```

## Environment
Please visit [`env`](./env) for further details.

1. Data preprocessing : Raw data is preprocessed in dedicated HPC at RWTH Aachen University.
2. Data analysis: `R-base 4.0` and `python 3.7` are used as base environments with different packages 
wrapped in conda environments. The processed data is analysed in BinAC bw-cluster (Tuebingen University).

- - -
<img src="./img/unihd_logo.png" height="100" align="left">
