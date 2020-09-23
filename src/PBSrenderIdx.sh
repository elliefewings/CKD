# Ensures that the Linux environment for the job is the same as the one we're working in:
#PBS -V

# Change directory to the project
if [ ! -z ${PBS_O_WORKDIR+x} ];then
  cd ${PBS_O_WORKDIR};
fi

source ~/.bashrc
conda activate ../envs/01_Seurat
export RSTUDIO_PANDOC="$CONDA_PREFIX/bin/pandoc"

${CONDA_PREFIX}/bin/R -e "rmarkdown::render('${RMD}', output_file='${PREFIX}.md', params=list(ID='$SAMPLEID',OUTDIR='${PREFIX}'))"

# This should be done by the rmarkdown header: github_document(html_preview=TRUE, keep_html=TRUE). 
# However it is not possible due to this current issue with pandoc, https://github.com/rstudio/rmarkdown/issues/1268
# In brief, html cannot be rendered if files are in a network drive.
# Then, we render the HTML manually only if the md exists
if [ -e "${PREFIX}.md" ];then
$CONDA_PREFIX/bin/pandoc +RTS -K512m -RTS "${PREFIX}.md" --to html4 --from gfm --output "${PREFIX}.html" --standalone --self-contained --highlight-style pygments --template $CONDA_PREFIX/lib/R/library/rmarkdown/rmarkdown/templates/github_document/resources/preview.html --variable "github-markdown-css:$CONDA_PREFIX/lib/R/library/rmarkdown/rmarkdown/templates/github_document/resources/github.css" --email-obfuscation none --metadata pagetitle=PREVIEW 
fi

