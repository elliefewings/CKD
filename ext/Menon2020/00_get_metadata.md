Get metadata for Menon 2020 dataset
================
Javier Perales-Patón - <javier.perales@bioquant.uni-heidelberg.de>

Two objectives are addressed here:

  - Getting clinical metadata for each sample. For both technologies,
    bulk rnaseq and single-cell rnaseq.
  - Reverse Engineer what are the library chemistry of single-cell
    samples. The original manuscript described two (v2 and v3). However
    original authors do not provide this information in the metadata. We
    could get which samples are which libraries because they used
    different thresholds in the cutoff for mitochondrial genes. (see
    later).

## Load libraries

``` r
library(GEOquery)
```

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Setting options('download.file.method.GEOquery'='auto')

    ## Setting options('GEOquery.inmemory.gpl'=FALSE)

``` r
library(Biobase)
```

## Get metadata from GEO

For convenience, we use GEOquery to parse the `Series Matrix` file from
the GEO accession number.

``` r
eset_list <- GEOquery::getGEO("GSE140989")
```

    ## Found 2 file(s)

    ## GSE140989-GPL11154_series_matrix.txt.gz

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## File stored at:

    ## /tmp/RtmpjZ9F2A/GPL11154.soft

    ## GSE140989-GPL20301_series_matrix.txt.gz

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## File stored at:

    ## /tmp/RtmpjZ9F2A/GPL20301.soft

``` r
length(eset_list)
```

    ## [1] 2

There are two expression-sets. Actually one for bulk and another for
single-cell data. This is due to the fact that each bunch of samples
were profiled in different instruments (2000 vs 4000, Illumina HiSeq).

``` r
lapply(eset_list, ncol)
```

    ## $`GSE140989-GPL11154_series_matrix.txt.gz`
    ## Samples 
    ##      68 
    ## 
    ## $`GSE140989-GPL20301_series_matrix.txt.gz`
    ## Samples 
    ##      24

``` r
lapply(eset_list, function(eset) pData(eset)$title)
```

    ## $`GSE140989-GPL11154_series_matrix.txt.gz`
    ##  [1] "FSG1 [bulk RNA-seq]"  "FSG2 [bulk RNA-seq]"  "FSG3 [bulk RNA-seq]" 
    ##  [4] "FSG4 [bulk RNA-seq]"  "FSG5 [bulk RNA-seq]"  "FSG6 [bulk RNA-seq]" 
    ##  [7] "FSG7 [bulk RNA-seq]"  "FSG8 [bulk RNA-seq]"  "FSG9 [bulk RNA-seq]" 
    ## [10] "FSG10 [bulk RNA-seq]" "FSG11 [bulk RNA-seq]" "FSG12 [bulk RNA-seq]"
    ## [13] "FSG13 [bulk RNA-seq]" "FSG14 [bulk RNA-seq]" "FSG15 [bulk RNA-seq]"
    ## [16] "FSG16 [bulk RNA-seq]" "FSG17 [bulk RNA-seq]" "FSG18 [bulk RNA-seq]"
    ## [19] "FSG19 [bulk RNA-seq]" "FSG20 [bulk RNA-seq]" "FSG21 [bulk RNA-seq]"
    ## [22] "FSG22 [bulk RNA-seq]" "FSG23 [bulk RNA-seq]" "FSG24 [bulk RNA-seq]"
    ## [25] "FSG25 [bulk RNA-seq]" "FSG26 [bulk RNA-seq]" "FSG27 [bulk RNA-seq]"
    ## [28] "FSG28 [bulk RNA-seq]" "FSG29 [bulk RNA-seq]" "FSG30 [bulk RNA-seq]"
    ## [31] "FSG31 [bulk RNA-seq]" "FSG32 [bulk RNA-seq]" "FSG33 [bulk RNA-seq]"
    ## [34] "FSG34 [bulk RNA-seq]" "FSG35 [bulk RNA-seq]" "FSG36 [bulk RNA-seq]"
    ## [37] "FSG37 [bulk RNA-seq]" "FSG38 [bulk RNA-seq]" "FSG39 [bulk RNA-seq]"
    ## [40] "FSG40 [bulk RNA-seq]" "FSG41 [bulk RNA-seq]" "FSG42 [bulk RNA-seq]"
    ## [43] "FSG43 [bulk RNA-seq]" "FSG44 [bulk RNA-seq]" "FSG45 [bulk RNA-seq]"
    ## [46] "FSG47 [bulk RNA-seq]" "FSG48 [bulk RNA-seq]" "FSG49 [bulk RNA-seq]"
    ## [49] "FSG50 [bulk RNA-seq]" "FSG51 [bulk RNA-seq]" "FSG52 [bulk RNA-seq]"
    ## [52] "FSG54 [bulk RNA-seq]" "FSG56 [bulk RNA-seq]" "FSG57 [bulk RNA-seq]"
    ## [55] "FSG58 [bulk RNA-seq]" "FSG59 [bulk RNA-seq]" "FSG60 [bulk RNA-seq]"
    ## [58] "FSG61 [bulk RNA-seq]" "FSG62 [bulk RNA-seq]" "FSG63 [bulk RNA-seq]"
    ## [61] "FSG64 [bulk RNA-seq]" "FSG65 [bulk RNA-seq]" "FSG66 [bulk RNA-seq]"
    ## [64] "FSG67 [bulk RNA-seq]" "FSG69 [bulk RNA-seq]" "FSG70 [bulk RNA-seq]"
    ## [67] "FSG71 [bulk RNA-seq]" "FSG74 [bulk RNA-seq]"
    ## 
    ## $`GSE140989-GPL20301_series_matrix.txt.gz`
    ##  [1] "PREMIERE-17-1606-2-0 [scRNA-seq]"        
    ##  [2] "PREMIERE-17-1606-2-1 [scRNA-seq]"        
    ##  [3] "PREMIERE-18-139 [scRNA-seq]"             
    ##  [4] "PREMIERE-18142-5-1 [scRNA-seq]"          
    ##  [5] "PREMIERE-18142-5-2 [scRNA-seq]"          
    ##  [6] "PREMIERE-18-162 [scRNA-seq]"             
    ##  [7] "PREMIERE-2017-06-28-TN-1 [scRNA-seq]"    
    ##  [8] "PREMIERE-2017-06-28-TN-2 [scRNA-seq]"    
    ##  [9] "PREMIERE-Hkid-PH-2017-08-25 [scRNA-seq]" 
    ## [10] "PREMIERE-Hkid-YY-2017-08-25 [scRNA-seq]" 
    ## [11] "PREMIERE-KPMP_Pilot_18342 [scRNA-seq]"   
    ## [12] "PREMIERE-SamplePRE027-1 [scRNA-seq]"     
    ## [13] "PREMIERE-SamplePRE038-1 [scRNA-seq]"     
    ## [14] "PREMIERE-SampleTransPRE19025 [scRNA-seq]"
    ## [15] "PREMIERE-TN109-18-242-Lib [scRNA-seq]"   
    ## [16] "PREMIERE-TN116-17-116-Lib-4a [scRNA-seq]"
    ## [17] "PREMIERE-TN117-17-116-Lib-4b [scRNA-seq]"
    ## [18] "PREMIERE-TN120-18-696-4 [scRNA-seq]"     
    ## [19] "PREMIERE-TN121-18-696-10 [scRNA-seq]"    
    ## [20] "PREMIERE-Trans-SV-004 [scRNA-seq]"       
    ## [21] "PREMIERE-Trans-SV-005 [scRNA-seq]"       
    ## [22] "PREMIERE-Trans-SV-006 [scRNA-seq]"       
    ## [23] "PREMIERE-Trans-SV-007 [scRNA-seq]"       
    ## [24] "PREMIERE-Trans-SV-008 [scRNA-seq]"

So actually we could split it

``` r
bulk_mdat <- pData(eset_list[["GSE140989-GPL11154_series_matrix.txt.gz"]])
sc_mdat <- pData(eset_list[["GSE140989-GPL20301_series_matrix.txt.gz"]])
```

## singleCell data

### Load normalized data

We are going to use the normalized data to find library chemistry. As it
is described in the paper method section (see below), different cutoff
of mitochondrial gene expression was used depending on library. So we
expect to see a clear signal.

``` r
sc_mdat$data_processing[1]
```

    ## [1] "Further data analyses were performed on the CellRanger output data files using Seurat 3 R package. As a quality control step, cells with less than 500 genes per cell were filtered out. The percentage of mitochondrial gene read content has emerged as a quality measure for the cell viability in single-cell studies. High percentage of mitochondrial reads indicates that the cells are less viable, as much of their cytoplasmic mRNAs were lost due to the breakage of plasma membrane during dissociation process. For the tumor-nephrectomy data generated by 10Xgenomics version2, a threshold of < 20% mitochondrial reads per cell was used, whereas, a cutoff of < 50% mitochondrial reads per cell was applied for the data from 10X genomics version 3. It has been observed by us and other groups that for certain tissues, including kidney, 10X genomics version 3 chemistry results in higher percentage of mitochondrial reads irrespective of the viability status of the cells. This issue was circumvented by deeper sequencing of the samples to obtain more genes per cell. After the quality control step, the data from the three sample types (tumor-nephrectomy, surveillance biopsy and pre-perfusion biopsy) were merged and further analyzed. In order to exclude cells that were doublets or multiplets from the analysis, we used a cutoff of > 500 and < 5000 genes per cell."

``` r
gz_con <- gzfile("data/GSE140989_PREMIERE_CombinedDataMatrixFile_Normalized.txt.gz", 'rt')
normDat <- read.table(gz_con, sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
close(gz_con) 

# Check how it looks the Cell Ids, which actually contain the sample of origin
head(colnames(normDat))
```

    ## [1] "18-139_AGAATAGTCACTTACT" "18-139_GACTACATCCATGCTC"
    ## [3] "18-139_GCAATCAGTCGGCATC" "18-139_TCACGAACACCGAATT"
    ## [5] "18-139_TGCGTGGGTATGAATG" "18-162_AAGCCGCCAAAGGAAG"

``` r
# Get the sample of origin
cellOrigin <- sapply(colnames(normDat), function(x) strsplit(x, split="_")[[1]][1])

sampleIDs_suppl <- unique(sort(cellOrigin))
```

If we extract the sample IDs reported inthe metadata by the author,
there is almost perfect match. We just have to tweeck a little bit one
sample and they match perfectly: just to remove a prefix, and by
alphanumerical order they match.

``` r
sc_mdat$sampleID <- gsub("^PREMIERE-","", gsub(" \\[scRNA-seq\\]", "",sc_mdat$title))

sidx <- order(gsub("KPMP_Pilot_18342","18-342", sc_mdat$sampleID))
sID_dic <- data.frame(sampleID=sc_mdat$sampleID[sidx], sampleID_suppl=sampleIDs_suppl )
print(sID_dic)
```

    ##               sampleID      sampleID_suppl
    ## 1          17-1606-2-0         17-1606-2-0
    ## 2          17-1606-2-1         17-1606-2-1
    ## 3               18-139              18-139
    ## 4               18-162              18-162
    ## 5     KPMP_Pilot_18342              18-342
    ## 6            18142-5-1           18142-5-1
    ## 7            18142-5-2           18142-5-2
    ## 8      2017-06-28-TN-1     2017-06-28-TN-1
    ## 9      2017-06-28-TN-2     2017-06-28-TN-2
    ## 10  Hkid-PH-2017-08-25  Hkid-PH-2017-08-25
    ## 11  Hkid-YY-2017-08-25  Hkid-YY-2017-08-25
    ## 12      SamplePRE027-1      SamplePRE027.1
    ## 13      SamplePRE038-1      SamplePRE038.1
    ## 14 SampleTransPRE19025 SampleTransPRE19025
    ## 15    TN109-18-242-Lib    TN109-18-242-Lib
    ## 16 TN116-17-116-Lib-4a TN116-17-116-Lib-4a
    ## 17 TN117-17-116-Lib-4b TN117-17-116-Lib-4b
    ## 18      TN120-18-696-4      TN120-18-696-4
    ## 19     TN121-18-696-10     TN121-18-696-10
    ## 20        Trans-SV-004          Trans004-1
    ## 21        Trans-SV-005          Trans005-1
    ## 22        Trans-SV-006          Trans006-1
    ## 23        Trans-SV-007          Trans007-2
    ## 24        Trans-SV-008          Trans008-2

This does not make sense at all. It looks the percentage of
mitochondrial genes is much lower than the expected one for the data
processing (expected 20 or 50% vs. observed 0-10%)

``` r
# Get mitochondrial genes
mito <- grep("^MT-", rownames(normDat), value=TRUE)
print(mito)
```

    ##  [1] "MT-ATP6" "MT-ATP8" "MT-CO1"  "MT-CO2"  "MT-CO3"  "MT-CYB"  "MT-ND1" 
    ##  [8] "MT-ND2"  "MT-ND3"  "MT-ND4"  "MT-ND4L" "MT-ND5"  "MT-ND6"

``` r
# Calculate the total expression of mitochondrial genes per cell
nMito <- colSums(normDat[mito,])/colSums(normDat) * 100

# Regroup by sample
sMito_bySample <- split(nMito, cellOrigin)

# Visualize approach
boxplot(sMito_bySample)
```

![](00_get_metadata_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
unlist(lapply(sMito_bySample, median))
```

    ##         17-1606-2-0         17-1606-2-1              18-139              18-162 
    ##            2.434621            2.842842            2.900768            2.009844 
    ##              18-342           18142-5-1           18142-5-2     2017-06-28-TN-1 
    ##            2.542562            2.086907            3.426722            1.586766 
    ##     2017-06-28-TN-2  Hkid-PH-2017-08-25  Hkid-YY-2017-08-25      SamplePRE027.1 
    ##            3.440920            2.590429            2.314063            5.226543 
    ##      SamplePRE038.1 SampleTransPRE19025    TN109-18-242-Lib TN116-17-116-Lib-4a 
    ##            5.175588            5.296576            3.715420            2.528922 
    ## TN117-17-116-Lib-4b      TN120-18-696-4     TN121-18-696-10          Trans004-1 
    ##            2.497692            2.849025            2.340055            3.627157 
    ##          Trans005-1          Trans006-1          Trans007-2          Trans008-2 
    ##            4.287728            3.743239            3.967328            4.117901

Actually supplementary Table 2 provides really a link of biopsy type to
sample ID. Source: <https://insight.jci.org/articles/view/133267/sd/1>,
suppl tab 2.

``` r
# Curated by hand... suppl tab 2
biopsy_dic <- data.frame(sampleID=c("SamplePRE027-1", "SamplePRE038-1", "SampleTransPRE19025",
              "Trans-SV-004", "Trans-SV-005", "Trans-SV-006", "Trans-SV-007", "Trans-SV-008",
              "17-1606-2-0", "17-1606-2-1", "18-139", "18-162", "KPMP_Pilot_18342", 
              "18142-5-1", "18142-5-2", "2017-06-28-TN-1", "2017-06-28-TN-2",
              "Hkid-PH-2017-08-25", "Hkid-YY-2017-08-25", "TN109-18-242-Lib", 
              "TN116-17-116-Lib-4a", "TN117-17-116-Lib-4b",
              "TN120-18-696-4", "TN121-18-696-10"),
            biopsy=c(rep("preperfusion", 3), rep("surveillance", 5), rep("tumornephrectomy", 16))
       )
stopifnot(all(biopsy_dic$sampleID %in% sID_dic$sampleID))
sID_dic <- merge(sID_dic, biopsy_dic, by.x="sampleID", by.y="sampleID", all.x=TRUE, all.y=FALSE)
```

We try again to figure out the impact of biopsy type on the
mitochondrial expression

``` r
sID_suppl2biopsy <- setNames(sID_dic$biopsy,sID_dic$sampleID_suppl) 
cols <- as.integer(factor(sID_suppl2biopsy[names(sMito_bySample)]))+1

boxplot(sMito_bySample, col=cols)
```

![](00_get_metadata_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

> Wow\! actually there is a signal there regarding type of biopsy and
> mitochondrial expression.

Fortunately we have now a link of type of sample to type of biopsy that
represents. So we could arrange tables and save everything. The most
important thing is that:

  - Extraction protocol: "10X Single Cell 3’ GEX – version 2 chemistry
    for all tumor-nephrectomy samples and the enhanced 10X Single Cell
    3’ GEX – version3 for the more recently processed transplant
    surveillance biopsies and pre-perfusion biopsies.
  - mitochondrial cutoff: tumor-nephrectomy data generated by
    10Xgenomics version2, a threshold of \< 20% mitochondrial reads per
    cell was used, whereas, a cutoff of \< 50% mitochondrial reads per
    cell was applied for the data from 10X genomics version 3

<!-- end list -->

``` r
sc_mdat2 <- merge(sc_mdat, sID_dic, by.x="sampleID", by.y="sampleID", all.x=TRUE, all.y=FALSE)
sc_mdat2$chemistry10X <- ifelse(sc_mdat2$biopsy=="tumornephrectomy", "v2","v3")
sc_mdat2$processing.mtcutoff <- ifelse(sc_mdat2$biopsy=="tumornephrectomy", 20, 50)
sc_mdat2$processing.nFeature.min <- 500
sc_mdat2$processing.nFeature.max <- 5000

write.table(sc_mdat2, file="./data/scRNAseq_metadata.tsv", sep="\t", 
        quote = FALSE, row.names = FALSE, col.names=TRUE)
```

## RNAseq data

``` r
# We remove data_processing since this information is more related to the other dataset (main dataset of the study)
# tha is the singlecell data (cell ranger, cutoffs for cell viability, etc).
write.table(bulk_mdat[,colnames(bulk_mdat)!="data_processing"], 
        file="./data/RNAseq_metadata.tsv", sep="\t", 
        quote = FALSE, row.names = FALSE, col.names=TRUE)
```

## Session

``` r
sessionInfo()
```

    ## R version 4.0.0 (2020-04-24)
    ## Platform: x86_64-conda_cos6-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /home/jperales/miniconda3/envs/CKD0/lib/libopenblasp-r0.3.9.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] GEOquery_2.56.0     Biobase_2.48.0      BiocGenerics_0.34.0
    ## [4] rmarkdown_2.1       nvimcom_0.9-82     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.4.6     xml2_1.3.2       knitr_1.28       magrittr_1.5    
    ##  [5] hms_0.5.3        tidyselect_1.1.0 R6_2.4.1         rlang_0.4.6     
    ##  [9] stringr_1.4.0    dplyr_0.8.5      tools_4.0.0      xfun_0.14       
    ## [13] htmltools_0.4.0  ellipsis_0.3.1   yaml_2.2.1       digest_0.6.25   
    ## [17] assertthat_0.2.1 tibble_3.0.1     lifecycle_0.2.0  crayon_1.3.4    
    ## [21] tidyr_1.1.0      readr_1.3.1      purrr_0.3.4      vctrs_0.3.0     
    ## [25] curl_4.3         glue_1.4.1       evaluate_0.14    limma_3.44.1    
    ## [29] stringi_1.4.6    compiler_4.0.0   pillar_1.4.4     pkgconfig_2.0.3
