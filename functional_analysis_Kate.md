# Functional analysis of Kate's RNA-seq data
Jennifer Grants  
1/24/2019  





## Load DEseq results for miR-X and miR-143  

```r
files <- list.files("./Kate/DEseq", full.names = T)

mirx <- read.csv(files[2])
mir143 <- read.csv(files[1])
```


## Load predicted targets of miR-X  
Based on exact seed sequence match.  

```r
x.target <- read.table("./Kate/exact_sponge_match/CUGUAGC.txt", sep = "\t")

colnames(x.target) <- "Target_X"
```


## Load predicted targets of miR-143  
Based on Targetscan predictions.  

```r
files2 <- list.files("./Kate/miR143_targetscan", full.names = T)

m143.target <- read.table(files2[1], sep = "\t", header = T) %>%
  select(1)

colnames(m143.target) <- "Target_143"
```



## Compare lists of miR-X and miR-143 targets to see if any shared  

```r
library(VennDiagram)
```

```
## Loading required package: grid
```

```
## Loading required package: futile.logger
```

```r
library(gridExtra)
```

```
## 
## Attaching package: 'gridExtra'
```

```
## The following object is masked from 'package:dplyr':
## 
##     combine
```

```r
venn.diagram(x = list(unlist(x.target), unlist(m143.target)), filename = "./Kate/functional_analysis_plots/Venn_Xtargets_143targets_pred.png", category.names = c("miR-X", "miR-143"), main = "Predicted miR targets", imagetype = "png", fill = c("dodgerblue", "firebrick"), cat.pos = c(0,0))
```

```
## [1] 1
```

```r
p <- venn.diagram(x = list(unlist(x.target), unlist(m143.target)), filename = NULL, category.names = c("miR-X", "miR-143"), main = "Predicted miR targets", imagetype = "png", fill = c("dodgerblue", "firebrick"), cat.pos = c(0,0))
grid.arrange(gTree(children = p))
```

![](functional_analysis_Kate_files/figure-html/compare_pred_targets-1.png)<!-- -->

> Conclusion: There is no significant overlap between targets (n = 1 common target).  


## Compare lists of significantly upregulated transcripts between miR-X SPG and miR-143 SPG-treated samples

```r
x.upreg <- mirx[which(mirx$padj < 0.05),]$X
m143.upreg <- mir143[which(mir143$padj < 0.05),]$X

venn.diagram(x = list(unlist(x.upreg), unlist(m143.upreg)), filename = "./Kate/functional_analysis_plots/Venn_Xtargets_143targets_DEseq.png", category.names = c("miR-X", "miR-143"), main = "Upregulated in DEseq (padj < 0.05)", imagetype = "png", fill = c("dodgerblue", "firebrick"), cat.pos = c(200,180))
```

```
## [1] 1
```

```r
p <- venn.diagram(x = list(unlist(x.upreg), unlist(m143.upreg)), filename = NULL, category.names = c("miR-X", "miR-143"), main = "Upregulated in DEseq (padj < 0.05)", imagetype = "png", fill = c("dodgerblue", "firebrick"), cat.pos = c(200,180))
grid.arrange(gTree(children = p))
```

![](functional_analysis_Kate_files/figure-html/compare_upreg_tx-1.png)<!-- -->

> Conclusion: There is a substantial overlap between genes upregulated between the 2 miR SPG treatments. This suggests that the presence of SPG vs. ctrl construct could be a confounding factor, or that the downstream effects of miR-143-SPG and miR-X-SPG are similar.  


## Convert ensembl transcript names to gene names

```r
library(biomaRt)

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

attr <- listAttributes(mart)

# Annotation frames for all genes: Oops, the object names shouldn't say "upreg"
anno.x.upreg <- getBM(attributes = c("ensembl_transcript_id", "external_gene_name"),
                        filters = "ensembl_transcript_id", values = mirx$X,
                        mart = mart)

anno.143.upreg <- getBM(attributes = c("ensembl_transcript_id", "external_gene_name"),
                        filters = "ensembl_transcript_id", values = mir143$X,
                        mart = mart)

detach(name = "package:biomaRt", unload = TRUE)

mirx$Gene_name <- anno.x.upreg[match(mirx$X, anno.x.upreg$ensembl_transcript_id),]$external_gene_name

mir143$Gene_name <- anno.143.upreg[match(mir143$X, anno.143.upreg$ensembl_transcript_id),]$external_gene_name
```



## Overlap upregulated genes with predicted targets  
### miR-X  

```r
venn.diagram(x = list(unlist(x.target), unlist(na.omit(mirx[which(mirx$padj < 0.05),]$Gene_name))), filename = "./Kate/functional_analysis_plots/Venn_pred_vs_upreg_mirx.png", category.names = c("Predicted targets", "Upregulated in DEseq"), main = "miR-X", imagetype = "png", fill = c("dodgerblue", "firebrick"), cat.pos = c(180, 180))
```

```
## [1] 1
```

```r
p <- venn.diagram(x = list(unlist(x.target), unlist(na.omit(mirx[which(mirx$padj < 0.05),]$Gene_name))), filename = NULL, category.names = c("Predicted targets", "Upregulated in DEseq"), main = "miR-X", imagetype = "png", fill = c("dodgerblue", "firebrick"), cat.pos = c(180, 180))
grid.arrange(gTree(children = p))
```

![](functional_analysis_Kate_files/figure-html/x_pred_vs_upreg-1.png)<!-- -->

> Conclusion: Very little overlap between upregulated genes and predicted targets, therefore upregulated genes are likely artifacts.



### miR-143  

```r
venn.diagram(x = list(unlist(m143.target), unlist(na.omit(mir143[which(mir143$padj < 0.05),]$Gene_name))), filename = "./Kate/functional_analysis_plots/Venn_pred_vs_upreg_mir143.png", category.names = c("Predicted targets", "Upregulated in DEseq"), main = "miR-143", imagetype = "png", fill = c("dodgerblue", "firebrick"), cat.pos = c(0, 0))
```

```
## [1] 1
```

```r
p <- venn.diagram(x = list(unlist(m143.target), unlist(na.omit(mir143[which(mir143$padj < 0.05),]$Gene_name))), filename = NULL, category.names = c("Predicted targets", "Upregulated in DEseq"), main = "miR-143", imagetype = "png", fill = c("dodgerblue", "firebrick"), cat.pos = c(0, 0))
grid.arrange(gTree(children = p))
```

![](functional_analysis_Kate_files/figure-html/m143_pred_vs_upreg-1.png)<!-- -->

> Conclusion: No overlap between upregulated genes and predicted targets, therefore upregulated genes are likely artifacts.
