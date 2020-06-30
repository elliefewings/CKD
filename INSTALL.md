## Environment for CKD study

### 1 . Install miniconda3
> NOTE: only done once (1st time) in your local machine

`miniconda3` is the framework to create and manage `envs`.
Before starting to work with virtual environments (`envs`), you need to install `miniconda3`.

Download the installer of `miniconda3`
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
Create a directory where you are going to install your SOFTWARE
inside your working folder in the cluster

```
mkdir /home/$USER/SOFTWARE/
```

Run the installer and define the installation path. See below.
```
bash Miniconda3-latest-Linux-x86_64.sh

# Press <ENTER> to continue an read Miniconda License
# Go down and enter 'yes' to accept Miniconda License
# Importantly! Miniconda3 will be installed into the default location,
# that is your $HOME. However, there is no space in your $HOME at BioQuant.
# Thus, specify another location. For instance:
#	/home/<username>/SOFTWARE/miniconda3
# run conda init when asked during the installation
```
If you didnt run `conda init` during the installation, this is a good time to do so. This will make changes into your `~/.bashrc` file needed to setup the conda environment when needed.

Remove the installer
```
rm Miniconda3-latest-Linux-x86_64.sh
```

### 2. Setup of environment

```bash
conda create --prefix ./envs/CKD --file ./envs/CKD_spec.txt
# Alternative: Only when same machine and environment.yml avail. Otherwise omit
## if [ ! -e ./envs/ ];then mkdir ./envs;fi
## conda env create --prefix ./envs/CKD -f environment.yml
```

### 3. External devtools::packages
There are a couple of packages that are not available in conda channels.
These must be installed manually in an interactive R session from the
conda environment using `devtools`.

```bash
# Terminal
conda activate ./envs/CKD
$CONDA_PREFIX/bin/R
```

```r
# R session
# Setup for devtools installation
getOption("unzip")
Sys.getenv("TAR")
options(unzip = "/opt/conda/bin/unzip")
Sys.setenv(TAR = "/bin/tar")
# Actual required packages
devtools::install_github("mahmoudibrahim/genesorteR") 
```
