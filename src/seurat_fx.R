 
# Find inflection points for selection of principle components
get_npcs <- function(seurat_object, create_plot = T){
  library(ggplot2)
  
  std_data = Stdev(object = seurat_object, reduction = "pca")
  ndims = length(x = std_data)
  elbow_data = data.frame(dims = 1:ndims, stdev = std_data[1:ndims])
  
  reference = elbow_data[,2]
  difference_vect = c(elbow_data[2:nrow(elbow_data),2],0)
  difference = which((reference - difference_vect) < 0.05)
  
  difference_consec = c(difference[2:length(difference)],0) - difference
  names(difference_consec) = difference
  
  npcs = as.numeric(names(which(difference_consec ==1)[1]))
  
  if(create_plot){
    
    plt = ggplot(elbow_data, aes(x = dims, y = stdev)) +
      geom_point() + 
      geom_vline(xintercept=npcs) +
      geom_text(aes(npcs+10, max(stdev), label=paste("Inferred NPCs:", npcs), vjust=0), size=6) +
      xlab("Number of principle components") +
      ylab("Standard deviation of principle components") +
      theme(text = element_text(size=17), axis.text=element_text(size=12))
    
  }
  
  out <- list(npcs=npcs, plot=plt)
  return(out)
  
}

#' ColRenameSNN rename columns related to SNN cell clustering in metadata Seurat Object
#' S Seurat Object

ColRenameSNN <- function(S) {

	assay_tag <- DefaultAssay(S)
	snn_idx <- grep(paste0("^", assay_tag, "_snn_res"),
			colnames(S@meta.data), value=FALSE)

	if(length(snn_idx)==0) {
		cat(paste0("[INFO] METADATA cols:", paste(colnames(S@meta.data), collapse=", "),"\n"), file=stdout())
		stop("[ERROR] : Not a single column related to SNN in metadata\n")
	}

	snn_cols <- colnames(S@meta.data)[snn_idx]
	cat(paste0("[INFO] SNN cols:", paste(snn_cols, collapse=", "),"\n"), file=stdout())
	snn_tagRES <- formatC(as.numeric(gsub(paste0(assay_tag, "_snn_res\\."), "", snn_cols)), 
			      digits=1, format="f")
	snn_rename <- paste0(assay_tag, "_snn_res.",snn_tagRES)

	colnames(S@meta.data)[snn_idx] <- snn_rename

	return(S)	
}
