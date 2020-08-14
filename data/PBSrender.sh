# Ensures that the Linux environment for the job is the same as the one we're working in:
#PBS -V

# Change directory to the project
if [ ! -z ${PBS_O_WORKDIR+x} ];then
  cd ${PBS_O_WORKDIR};
fi

source ~/.bashrc
conda activate ../envs/scrna
export RSTUDIO_PANDOC="$CONDA_PREFIX/bin/pandoc"

# Original in local
# ${CONDA_PREFIX}/bin/R -e "rmarkdown::render('$<', output_file='$@.md', params=list(ID='$(@:%_Seurat=%)',OUTDIR='$@'))"
# mv $@ $(filter-out $<,$^)
# mv $@.md $(filter-out $<,$^)
# mv $@.html $(filter-out $<,$^)

${CONDA_PREFIX}/bin/R -e "rmarkdown::render('${RMD}', output_file='${PREFIX}.md', params=list(ID='$SAMPLEID',OUTDIR='${PREFIX}'))"
mv ${PREFIX} ${DEST}
mv ${PREFIX}.md ${DEST}
mv ${PREFIX}.html ${DEST}
