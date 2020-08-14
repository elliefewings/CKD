# Ensures that the Linux environment for the job is the same as the one we're working in:
#PBS -V

# Change directory to the project
if [ ! -z ${PBS_O_WORKDIR+x} ];then
  cd ${PBS_O_WORKDIR};
fi

source ~/.bashrc
conda activate ../envs/scrna
export RSTUDIO_PANDOC="$CONDA_PREFIX/bin/pandoc"

${CONDA_PREFIX}/bin/R -e "rmarkdown::render('${RMD}', output_file='${PREFIX}.md', params=list(ID='$SAMPLEID',OUTDIR='${PREFIX}'))"
