---
title: "Expression analysis for FBXO11 knockdown RNA-seq data"
author: "Jennifer Grants"
date: "1/9/2019"
output:
  html_document:
    keep_md: true
    toc: true
---





## Load saved data   
`expr` is mRNA-seq TPM values, `libs` is metadata

```r
load("./FBXO11/expr.RData")
load("./FBXO11/libs.RData")
```


```r
libs$treatment <- ifelse(test = grepl(pattern = "shSCR", x = libs$specimen_subset_external_id), yes = "shSCR", no = 
                           ifelse(test = grepl(pattern = "shFBXO11#9", x = libs$specimen_subset_external_id), yes = "shFBXO11.a", no = "shFBXO11.b"))
```



```r
kable(libs)
```



library_name   sequencing_effort     platform_name   platform_version   specimen_subset_external_id   cohort_no   library_qc_info   index_sequence   treatment  
-------------  --------------------  --------------  -----------------  ----------------------------  ----------  ----------------  ---------------  -----------
A54939         Karsan Lab Research   ssRNA-Seq       v1                 OCI-AML3_shSCR-1              Cohort 6    {}                CCAACA           shSCR      
A54940         Karsan Lab Research   ssRNA-Seq       v1                 OCI-AML3_shFBXO11#9-1         Cohort 6    {}                CTAGCT           shFBXO11.a 
A54941         Karsan Lab Research   ssRNA-Seq       v1                 OCI-AML3_shFBXO11#10-1        Cohort 6    {}                GATGCT           shFBXO11.b 
A54942         Karsan Lab Research   ssRNA-Seq       v1                 OCI-AML3_shSCR-2              Cohort 6    {}                TAATCG           shSCR      
A54943         Karsan Lab Research   ssRNA-Seq       v1                 OCI-AML3_shFBXO11#9-2         Cohort 6    {}                TGAATG           shFBXO11.a 
A54944         Karsan Lab Research   ssRNA-Seq       v1                 OCI-AML3_shFBXO11#10-2        Cohort 6    {}                AGTTCC           shFBXO11.b 
A54945         Karsan Lab Research   ssRNA-Seq       v1                 OCI-AML3_shSCR-3              Cohort 6    {}                CGATGT           shSCR      
A54946         Karsan Lab Research   ssRNA-Seq       v1                 OCI-AML3_shFBXO11#9-3         Cohort 6    {}                TAGCTT           shFBXO11.a 
A54947         Karsan Lab Research   ssRNA-Seq       v1                 OCI-AML3_shFBXO11#10-3        Cohort 6    {}                AACTTG           shFBXO11.b 

# Quality control checks  

### First convert row names from column "Names":  

```r
expr2 <- column_to_rownames(expr, var = "Name")
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```r
head(expr2)
```

```
## # A tibble: 6 x 9
##   A54945 A54947 A54946 A54944 A54940 A54942 A54941 A54939 A54943
##    <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
## 1  0.176  0.561  0.102  0.459  0.115  0.307  0.274 0.381   0.180
## 2  0.609  0.235  0.170  0      0.141  0.267  0.514 0.0454  0    
## 3  3.94   4.44   4.94   4.64   3.26   2.98   3.09  4.44    2.96 
## 4  0      0      0      0      0      0      0     0       0    
## 5  0      0      0      0      0      0      0     0       0    
## 6  0      0      0      0      0      0      0     0       0
```


### Dimensions:   

```r
# expression data
dim(expr2)
```

```
## [1] 195480      9
```

```r
# metadata
dim(libs)
```

```
## [1] 9 9
```


### NA values:   

```r
any(is.na(expr2))
```

```
## [1] FALSE
```


### Range:  

```r
range(expr2)
```

```
## [1]     0.00 39567.71
```


```r
hist(as.matrix(expr2))
```

![](expression_analysis_FBXO11_files/figure-html/histogram_original_data-1.png)<!-- -->
Appears not to be log transformed.    


## Data transformation / row normalization  
### Log2 transformation  

```r
# log2 transformation with constant +0.01 added to allow transformation of zeros
expr.log <- log2(expr2 + 0.01)

# Visualize log transformed data
hist(as.matrix(expr.log))
```

![](expression_analysis_FBXO11_files/figure-html/log_transf-1.png)<!-- -->
* ASK: Why does this have so many 0 values? Is this normal for RNA-seq?
    * Possible answer: Blood cells may not express that many genes highly...
* LOOK UP: Griffith lab RNAseq data analysis tutorial online -- assumptions



# DE-seq  

> NOTE: Limma analysis appears to have more power for detecting relevant gene expression changes, probably because it better estimates the variance and because it is possible to model 3 conditions to avoid losing power

> Skip to Limma section!!



```r
design.table <- select(libs, library_name, treatment) %>%
  column_to_rownames(var = "library_name")
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```r
expr.matrix <- column_to_rownames(expr, var = "Name") %>%
  select(rownames(design.table)) %>% # put columns in same order as the design matrix
  as.matrix()
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```r
# filter to "highly expressed" genes
thresh <- expr.matrix >= 1
keep <- rowSums(thresh) >= 2 # keep rows where there are > 2 samples with tpm > 1

expr.matrix.keep <- expr.matrix[keep,]

# convert to integers for DEseq
int.matrix.keep <- round(expr.matrix.keep, digits = 0)

dds <- DESeqDataSetFromMatrix(countData = int.matrix.keep, colData = design.table, design = ~treatment)
```

```
## converting counts to integer mode
```

```
## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
## design formula are characters, converting to factors
```

```r
res <- DESeq(dds)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```


## Summarise results by each treatment group  

### shFBXO11.a vs. shSCR  

```r
resultsNames(res)
```

```
## [1] "Intercept"                         
## [2] "treatment_shFBXO11.b_vs_shFBXO11.a"
## [3] "treatment_shSCR_vs_shFBXO11.a"
```

```r
result.pair <- results(res, contrast = c("treatment", "shFBXO11.a", "shSCR")) # contrast: the factor, then numerator, then denominator

kable(head(result.pair[which(result.pair$padj < 0.05),]))
```

                     baseMean   log2FoldChange       lfcSE        stat      pvalue        padj
----------------  -----------  ---------------  ----------  ----------  ----------  ----------
ENST00000373449      98.18959       -0.6872024   0.1931158   -3.558499   0.0003730   0.0164194
ENST00000396651    1376.61361       -0.3319752   0.0656212   -5.058962   0.0000004   0.0000518
ENST00000370321    1253.92354       -0.3414709   0.0754813   -4.523914   0.0000061   0.0005594
ENST00000369642      79.15805        1.2160224   0.2887543    4.211271   0.0000254   0.0017828
ENST00000369637      84.30236        0.9893354   0.2674085    3.699716   0.0002158   0.0109731
ENST00000607355      52.60446        1.9217321   0.3069893    6.259932   0.0000000   0.0000001

```r
write.csv(x = result.pair, file = "./FBXO11/DEseq/shFBXO11a_vs_shSCR.csv")
```


```r
# cutoffs: padj < 0.1, |FC| > 1.25 (permissive)
dat.a <- read.csv("./FBXO11/DEseq/shFBXO11a_vs_shSCR.csv") %>%
  mutate(Signif_p = -log10(padj) > -log10(0.1), 
         Signif_delta = abs(log2FoldChange) > log2(1.25),
         Signif = Signif_p == TRUE & Signif_delta == TRUE)

ggplot(dat.a, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(colour = Signif))
```

```
## Warning: Removed 45391 rows containing missing values (geom_point).
```

![](expression_analysis_FBXO11_files/figure-html/sh_a_volcano-1.png)<!-- -->



### shFBXO11.b vs shSCR  

```r
result.pair <- results(res, contrast = c("treatment", "shFBXO11.b", "shSCR")) # contrast: the factor, then numerator, then denominator

kable(head(result.pair[which(result.pair$padj < 0.05),]))
```

                    baseMean   log2FoldChange       lfcSE        stat     pvalue        padj
----------------  ----------  ---------------  ----------  ----------  ---------  ----------
ENST00000344843    185.20045       -1.4340034   0.1559497   -9.195293   0.00e+00   0.0000000
ENST00000376957    536.58206       -0.4561490   0.1008886   -4.521314   6.10e-06   0.0113374
ENST00000329421     92.97947       -0.9143182   0.1966855   -4.648630   3.30e-06   0.0072363
ENST00000582401     89.83154        0.8797400   0.1980195    4.442693   8.90e-06   0.0152587
ENST00000480760     48.09922        1.0812659   0.2606061    4.149043   3.34e-05   0.0417657
ENST00000367600    167.39928        0.8321301   0.1996213    4.168544   3.07e-05   0.0412676

```r
write.csv(x = result.pair, file = "./FBXO11/DEseq/shFBXO11b_vs_shSCR.csv")
```


```r
# cutoffs: padj < 0.1, |FC| > 1.25 (permissive)
dat.b <- read.csv("./FBXO11/DEseq/shFBXO11b_vs_shSCR.csv") %>%
  mutate(Signif_p = -log10(padj) > -log10(0.1), 
         Signif_delta = abs(log2FoldChange) > log2(1.25),
         Signif = Signif_p == TRUE & Signif_delta == TRUE)

ggplot(dat.b, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(colour = Signif))
```

```
## Warning: Removed 5 rows containing missing values (geom_point).
```

![](expression_analysis_FBXO11_files/figure-html/sh_b_volcano-1.png)<!-- -->



## Heatmap for differentially expressed genes in either shFBXO11 treatment  

```r
# filter to genes that are "signif" in the volcano plots above in either treatment
genes_a <- dat.a[which(dat.a$Signif == TRUE),]$X %>% as.vector()
genes_b <- dat.b[which(dat.b$Signif == TRUE),]$X %>% as.vector()

signif.genes <- c(genes_a, genes_b)

expr.signif <- expr[which(expr$Name %in% signif.genes),]

expr.signif2 <- column_to_rownames(expr.signif, var = "Name")
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```r
head(expr.signif2)
```

```
## # A tibble: 6 x 9
##   A54945 A54947 A54946 A54944 A54940 A54942 A54941 A54939 A54943
##    <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
## 1   229.   89.2  248.    75.8  259.    213.   83.9  241.   237. 
## 2   580.  451.   597.   425.   621.    572.  407.   634.   555. 
## 3   122.   71.2   94.0   62.1  101.    116.   61.0  133.    78.8
## 4   114.  119.    68.1  112.    66.0   120.  114.    95.0   74.3
## 5   152.  130.   105.   137.    72.2   129.   92.6  120.    73.4
## 6  1493. 1372.  1272.  1345.  1240.   1546. 1290.  1599.  1240.
```

Log transform, and scale rows:  

```r
# don't need to add a constant for the log scaling because all are >0
log.norm <- log2(expr.signif2) %>%
  t() %>% scale() %>% (t) %>% ## method from https://github.com/STAT540-UBC/STAT540-UBC.github.io/blob/master/seminars/seminars_winter_2019/seminar6/sm06_clustering-pca.md
  as.matrix()

log.norm <- na.omit(log.norm)
any(is.na(log.norm))
```

```
## [1] FALSE
```


```r
library(RColorBrewer)

heatmap_pallete <- colorRampPalette(brewer.pal(8, name = "RdBu"))(21) %>% rev

anno.frame <- data.frame(sample = colnames(log.norm), shRNA = libs[match(colnames(log.norm), libs$library_name),]$treatment) %>%
  column_to_rownames(var = "sample")

anno.frame$shRNA <- factor(anno.frame$shRNA, levels = c("shSCR", "shFBXO11.a", "shFBXO11.b"))

as.matrix(log.norm) %>%
pheatmap(cluster_cols = T, 
         cluster_rows = T, 
         scale = "none", 
         clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean", 
         show_rownames = F, 
         annotation_col = anno.frame, color = heatmap_pallete)
```

![](expression_analysis_FBXO11_files/figure-html/heatmap_union-1.png)<!-- -->

> Conclusion: When considering all differentially expressed genes (union), shFBXO11.b (#10) is more similar to shSCR than shFBXO11.a (#9).  

> Overall: Try Limma analysis to better estimate variance and to avoid losing statistical power when doing separate pairwise comparisons.



## Heatmap for differentially expressed genes in BOTH shFBXO11 constructs  

```r
signif.genes <- intersect(genes_a, genes_b)

expr.signif <- expr[which(expr$Name %in% signif.genes),]

expr.signif2 <- column_to_rownames(expr.signif, var = "Name")
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```r
head(expr.signif2)
```

```
## # A tibble: 6 x 9
##   A54945 A54947 A54946 A54944 A54940 A54942 A54941 A54939 A54943
##    <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
## 1  219.   324.   161.   305.   160.   229.   288.   217.   149. 
## 2   93.5   22.6   32.9   19.1   28.7   84.5   21.8   95.8   27.7
## 3   50.7   93.7   26.9   94.9   23.9   49.3   89.9   50.6   25.0
## 4  112.   191.   207.   178.   201.   107.   164.   118.   181. 
## 5 1217.   946.   795.   822.   701.  1162.   813.  1239.   716. 
## 6   53.4   22.0   27.7   22.1   31.0   54.2   22.5   49.0   29.1
```

Log transform, and scale rows:  

```r
# don't need to add a constant for the log scaling because all are >0
log.norm <- log2(expr.signif2) %>%
  t() %>% scale() %>% (t) %>% ## method from https://github.com/STAT540-UBC/STAT540-UBC.github.io/blob/master/seminars/seminars_winter_2019/seminar6/sm06_clustering-pca.md
  as.matrix()

log.norm <- na.omit(log.norm)
any(is.na(log.norm))
```

```
## [1] FALSE
```



```r
anno.frame <- data.frame(sample = colnames(log.norm), shRNA = libs[match(colnames(log.norm), libs$library_name),]$treatment) %>%
  column_to_rownames(var = "sample")

anno.frame$shRNA <- factor(anno.frame$shRNA, levels = c("shSCR", "shFBXO11.a", "shFBXO11.b"))

as.matrix(log.norm) %>%
pheatmap(cluster_cols = T, 
         cluster_rows = T, 
         scale = "none", 
         clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean", 
         show_rownames = F, 
         annotation_col = anno.frame, color = heatmap_pallete)
```

![](expression_analysis_FBXO11_files/figure-html/heatmap_intersection-1.png)<!-- -->

> Conclusion: When considering the common DE genes (intersection), the two shFBXO11 treatments cluster distinctly from shSCR.  

> Observation: All but 3 genes apprear to have the same directionality of change (up, up ; down, down) in the shFBXO11 treatments vs. shSCR.  


# Limma analysis  

> Note: Limma analysis works well for low "n" experiments because it uses a moderated t-test to better estimate variance across genes. 

> It IS possible to perform limma for 3 groups, then report results for each comparison (avoids losing power!)



Method from Limma User Guide section on RNA-seq: https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf  


```r
## Reload saved data
load("./FBXO11/expr.RData")
load("./FBXO11/libs.RData")



# Design factors
libs$treatment <- ifelse(test = grepl(pattern = "shSCR", x = libs$specimen_subset_external_id), yes = "ctrl", no = 
                           ifelse(test = grepl(pattern = "shFBXO11#9", x = libs$specimen_subset_external_id), yes = "shFBXO11.a", no = "shFBXO11.b"))

design.table <- select(libs, library_name, treatment) %>%
  column_to_rownames(var = "library_name")
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```r
design.table$treatment <- factor(design.table$treatment, levels = c("ctrl", "shFBXO11.a", "shFBXO11.b"))


# make model matrix
model.matrix <- model.matrix(~treatment, data = design.table)

  ## commented out the stuff for the makeContrasts() method
  #x <- c("shSCR-shFBXO11.b")
  #contrast.matrix <- makeContrasts(contrasts = x, levels = design.table)



# convert expr to matrix
expr.matrix <- column_to_rownames(expr, var = "Name") %>%
  select(rownames(design.table)) %>% # put columns in same order as the design matrix
  as.matrix()
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```r
# filter to "highly expressed" genes
thresh <- expr.matrix >= 5
keep <- rowSums(thresh) >= 2 # keep rows where there are > 2 samples with tpm > 5

expr.matrix.keep <- expr.matrix[keep,]
```


## Normalizing & log transformation

```r
dge <- DGEList(counts = expr.matrix.keep)

dge <- calcNormFactors(dge) # From the Limma Userguide: It is usual to apply scale normalization to RNA-seq read counts, and the TMM normalization method [34] in particular has been found to perform well in comparative studies. This can be applied to the DGEList object

kable(head(dge$counts))
```

                       A54939       A54940       A54941       A54942       A54943       A54944       A54945       A54946       A54947
----------------  -----------  -----------  -----------  -----------  -----------  -----------  -----------  -----------  -----------
ENST00000494149      7.491241     3.406867     5.204649     6.812622     4.527275     5.917108     6.312374     6.817580     5.484473
ENST00000623083      6.723140     7.166121     7.513370     6.213673     6.145557     7.718823     6.107802     6.883171     7.712641
ENST00000599771      9.128803     9.273547     9.977130    14.545309    10.078416    14.117502    13.374769    13.281399    11.277027
ENST00000416931    119.985967   165.840842   162.077902   116.268604   120.318414    94.180728    93.371523    61.124406    88.833200
ENST00000457540     58.084029    64.861148    73.284931    70.872431    73.885293    67.842562    68.317703    60.415320    67.094112
ENST00000617238    350.485344   286.532116   429.256664   433.241176   285.599971   320.779323   214.998546   175.050479   260.938612

```r
log.dge <- cpm(y = dge, normalized.lib.sizes = FALSE, log = TRUE, lib.size = NA) # using the cpm function to take the log2(count+0.01), BUT since the values are already in TPM, set all other paramters to FALSE or NA to avoid taking "per million" again.

kable(head(log.dge))
```

                     A54939     A54940     A54941     A54942     A54943     A54944     A54945     A54946     A54947
----------------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------
ENST00000494149    3.346942   2.535991   2.949503   3.242304   2.807192   3.088322   3.156834   3.249305   3.006281
ENST00000623083    3.225240   3.297230   3.350365   3.140752   3.126579   3.384192   3.120890   3.260012   3.382246
ENST00000599771    3.576493   3.595713   3.682512   4.151134   3.694735   4.114075   4.043951   4.043440   3.833244
ENST00000416931    7.030420   7.491488   7.458156   6.988757   7.034518   6.691331   6.676832   6.090700   6.607535
ENST00000457540    6.008802   6.163668   6.334243   6.290130   6.345793   6.229679   6.237169   6.074400   6.212878
ENST00000617238    8.561241   8.273117   8.852303   8.868511   8.267918   8.438079   7.862875   7.578764   8.140970

```r
kable(head(log2(dge$counts))) ## these results are fairly similar, so I think that my method using cpm() worked
```

                     A54939     A54940     A54941     A54942     A54943     A54944     A54945     A54946     A54947
----------------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------
ENST00000494149    2.905205   1.768446   2.379801   2.768210   2.178643   2.564892   2.658183   2.769260   2.455353
ENST00000623083    2.749135   2.841192   2.909460   2.635446   2.619544   2.948381   2.610653   2.783073   2.947225
ENST00000599771    3.190426   3.213121   3.318625   3.862482   3.333197   3.819413   3.741442   3.731335   3.495315
ENST00000416931    6.906722   7.373655   7.340544   6.861318   6.910714   6.557360   6.544911   5.933677   6.473027
ENST00000457540    5.860070   6.019283   6.195445   6.147153   6.207215   6.084119   6.094188   5.916842   6.068114
ENST00000617238    8.453210   8.162553   8.745697   8.759027   8.157852   8.325437   7.748183   7.451627   8.027567

## Perform limma-trend analysis

```r
fit <- lmFit(object = log.dge, design = model.matrix)  # fit with all groups to avoid losing power, but will need to extract results for each experimental group separately
fitEb <- eBayes(fit, trend=TRUE)

# results for shFBXO11.a vs. ctrl
top_a <- topTable(fitEb, coef = "treatmentshFBXO11.a", sort.by = "p", n = Inf) %>% 
  rownames_to_column(var = "Transcript") %>%
  mutate(signif = adj.P.Val < 0.05 & abs(logFC) > 1)
head(top_a)
```

```
##        Transcript     logFC  AveExpr         t      P.Value    adj.P.Val
## 1 ENST00000371732 -2.518802 3.799988 -31.19699 7.430023e-10 1.371656e-05
## 2 ENST00000397406 -1.929505 4.281757 -24.62347 5.135360e-09 4.740194e-05
## 3 ENST00000366903  1.676875 3.347192  23.04834 8.797441e-09 4.839102e-05
## 4 ENST00000319340 -1.479387 4.580358 -22.55638 1.048503e-08 4.839102e-05
## 5 ENST00000290705 -2.279134 3.713746 -20.62284 2.170724e-08 8.014746e-05
## 6 ENST00000221418 -2.043354 7.190231 -19.21486 3.848571e-08 1.184141e-04
##           B signif
## 1 11.699699   TRUE
## 2 10.562692   TRUE
## 3 10.201842   TRUE
## 4 10.080175   TRUE
## 5  9.555117   TRUE
## 6  9.119514   TRUE
```

```r
ggplot(top_a, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(color = signif))
```

![](expression_analysis_FBXO11_files/figure-html/limma_analysis-1.png)<!-- -->

```r
write.csv(x = top_a, file = "./FBXO11/Limma/Limma_result_shFBXO11-9_vs_shSCR.csv", row.names = F)

# results for shFBXO11.b vs. ctrl
top_b <- topTable(fitEb, coef = "treatmentshFBXO11.b", sort.by = "p", n = Inf) %>% 
  rownames_to_column(var = "Transcript") %>%
  mutate(signif = adj.P.Val < 0.05 & abs(logFC) > 1)
head(top_b)
```

```
##        Transcript     logFC  AveExpr         t      P.Value    adj.P.Val
## 1 ENST00000265643 -2.007431 5.456055 -21.48200 1.558758e-08 0.0001940629
## 2 ENST00000286713 -1.415996 3.725568 -20.33441 2.433265e-08 0.0001940629
## 3 ENST00000620913 -1.399649 4.792965 -19.69369 3.153614e-08 0.0001940629
## 4 ENST00000344843 -1.436316 7.506916 -18.29384 5.723279e-08 0.0002156084
## 5 ENST00000310624 -1.237790 4.349835 -18.24830 5.839564e-08 0.0002156084
## 6 ENST00000347703 -1.161754 5.208797 -17.62316 7.734843e-08 0.0002296508
##          B signif
## 1 9.646958   TRUE
## 2 9.333499   TRUE
## 3 9.145296   TRUE
## 4 8.697408   TRUE
## 5 8.681928   TRUE
## 6 8.463200   TRUE
```

```r
ggplot(top_b, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(color = signif))
```

![](expression_analysis_FBXO11_files/figure-html/limma_analysis-2.png)<!-- -->

```r
write.csv(x = top_b, file = "./FBXO11/Limma/Limma_result_shFBXO11-10_vs_shSCR.csv", row.names = F)

# List of genes with significant adj.P.Val in both shFBXO11 treatments
common_genes <- intersect(top_a[which(top_a$adj.P.Val < 0.05),]$Transcript, top_b[which(top_b$adj.P.Val < 0.05),]$Transcript)

# Summaries
nrow(top_a[which(top_a$adj.P.Val < 0.05),])
```

```
## [1] 2049
```

```r
nrow(top_b[which(top_b$adj.P.Val < 0.05),])
```

```
## [1] 1856
```

```r
length(common_genes)
```

```
## [1] 543
```

### QC Check: p-values histograms

```r
hist(top_a$P.Value)
```

![](expression_analysis_FBXO11_files/figure-html/unnamed-chunk-3-1.png)<!-- -->


```r
hist(top_b$P.Value)
```

![](expression_analysis_FBXO11_files/figure-html/unnamed-chunk-4-1.png)<!-- -->



### Identify the ~500 genes that are significantly altered in both shFBXO11 treatments

```r
signif.genes <- common_genes
```



## Heatmap for differentially expressed genes in BOTH shFBXO11 constructs  

```r
expr.signif <- expr[which(expr$Name %in% signif.genes),]

expr.signif2 <- column_to_rownames(expr.signif, var = "Name")
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```r
head(expr.signif2)
```

```
## # A tibble: 6 x 9
##   A54945 A54947 A54946 A54944 A54940 A54942 A54941 A54939 A54943
##    <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
## 1  14.6   20.7   22.8   23.8   25.7   13.3   19.4   17.7   23.2 
## 2   7.47   3.20   4.82   2.93   5.02   6.99   3.70  11.1    4.33
## 3   7.14   5.23   3.89   3.68   3.60   8.86   3.48   6.17   2.99
## 4  36.0   54.6   27.1   62.6   18.1   34.6   57.2   38.1   26.9 
## 5   8.97   4.90   5.71   4.89   3.44   7.78   4.50   9.12   3.13
## 6  64.0   54.9   53.6   52.5   54.2   66.5   49.8   66.5   49.9
```

Log transform, and scale rows:  

```r
# don't need to add a constant for the log scaling because all are >0
log.norm <- log2(expr.signif2) %>%
  t() %>% scale() %>% (t) %>% ## method from https://github.com/STAT540-UBC/STAT540-UBC.github.io/blob/master/seminars/seminars_winter_2019/seminar6/sm06_clustering-pca.md
  as.matrix()

log.norm <- na.omit(log.norm)
any(is.na(log.norm))
```

```
## [1] FALSE
```



```r
anno.frame <- data.frame(sample = colnames(log.norm), shRNA = libs[match(colnames(log.norm), libs$library_name),]$treatment) %>%
  column_to_rownames(var = "sample")

anno.frame$shRNA <- factor(anno.frame$shRNA, levels = c("ctrl", "shFBXO11.a", "shFBXO11.b"))

heatmap_pallete <- colorRampPalette(brewer.pal(8, name = "RdBu"))(21) %>% rev

out <- as.matrix(log.norm) %>%
pheatmap(cluster_cols = T, 
         cluster_rows = T, 
         scale = "none", 
         clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean", 
         show_rownames = F, 
         annotation_col = anno.frame, color = heatmap_pallete)
```

![](expression_analysis_FBXO11_files/figure-html/heatmap_int_limma-1.png)<!-- -->

> Conclusion: Compared to DEseq, Limma analysis identified more differentially expressed genes that correspond between the 2 shFBXO11 construct treatments. 

> Next step: Extract the row clusters from this heatmap to get genes that are down in both or up in both shFBXO11 construct treatments.

### Extract row clusters

```r
anno.rowclusters <- cutree(out$tree_row, k = 4) %>% as.data.frame()
colnames(anno.rowclusters) <- "cluster"

anno.rowclusters$cluster <- factor(anno.rowclusters$cluster)
levels(anno.rowclusters$cluster) <- c("D", "B", "A", "C")
anno.rowclusters$cluster <- factor(anno.rowclusters$cluster, levels = c("A", "B", "C", "D")) # to give them a logical order in the heatmap
```


## Heatmap with row clusters

```r
as.matrix(log.norm) %>%
pheatmap(cluster_cols = T, 
         cluster_rows = T, 
         scale = "none", 
         clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean", 
         show_rownames = F, 
         annotation_col = anno.frame, 
         color = heatmap_pallete, 
         annotation_row = anno.rowclusters)
```

![](expression_analysis_FBXO11_files/figure-html/heatmap_row_clusters-1.png)<!-- -->


### Clusters of interest: extract expression data in limma result for each shFBXO11 construct  

Interested in transcript clusters B (downreg in shFBXO11's) and D (upreg in shFBXO11's)

```r
anno.cl2 <- rownames_to_column(anno.rowclusters, var = "Transcript")

top_a_clust <- left_join(top_a, anno.cl2) %>%
  na.omit() %>%
  filter(cluster %in% c("B", "D"))
```

```
## Joining, by = "Transcript"
```

```r
write.csv(x = top_a_clust, file = "./FBXO11/Limma/Limma_commonEffectClusters_shFBXO11-9_vs_shSCR.csv", row.names = F)

top_b_clust <- left_join(top_b, anno.cl2) %>%
  na.omit() %>%
  filter(cluster %in% c("B", "D"))
```

```
## Joining, by = "Transcript"
```

```r
write.csv(x = top_b_clust, file = "./FBXO11/Limma/Limma_commonEffectClusters_shFBXO11-10_vs_shSCR.csv", row.names = F)
```

