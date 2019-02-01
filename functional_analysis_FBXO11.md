---
title: 'Functional Analysis: shFBXO11'
author: "Jennifer Grants"
date: "1/31/2019"
output: 
  html_document:
    keep_md: true
    toc: true
---



# Load common effect cluster gene Limma analysis results  

```r
files <- list.files("./FBXO11/Limma", full.names = T)
files
```

```
## [1] "./FBXO11/Limma/Limma_commonEffectClusters_shFBXO11-10_vs_shSCR.csv"
## [2] "./FBXO11/Limma/Limma_commonEffectClusters_shFBXO11-9_vs_shSCR.csv" 
## [3] "./FBXO11/Limma/Limma_result_shFBXO11-10_vs_shSCR.csv"              
## [4] "./FBXO11/Limma/Limma_result_shFBXO11-9_vs_shSCR.csv"
```

```r
a <- read.csv(files[2])
b <- read.csv(files[1])

head(a)
```

```
##        Transcript     logFC  AveExpr         t      P.Value    adj.P.Val
## 1 ENST00000371732 -3.058143 3.501391 -27.98200 1.946558e-09 3.593540e-05
## 2 ENST00000397406 -2.189460 4.101586 -24.14325 6.451672e-09 5.955216e-05
## 3 ENST00000319340 -1.626931 4.447705 -22.94618 9.741003e-09 5.994289e-05
## 4 ENST00000366903  2.038695 3.005491  21.68305 1.540036e-08 7.107653e-05
## 5 ENST00000290705 -2.817652 3.403743 -17.61194 8.219180e-08 2.528905e-04
## 6 ENST00000265643 -1.606772 5.381668 -16.43852 1.427421e-07 2.725433e-04
##           B signif cluster
## 1 11.105117   TRUE       B
## 2 10.366337   TRUE       B
## 3 10.090271   TRUE       B
## 4  9.770508   TRUE       D
## 5  8.492993   TRUE       B
## 6  8.037556   TRUE       B
```

```r
head(b)
```

```
##        Transcript      logFC  AveExpr         t      P.Value    adj.P.Val
## 1 ENST00000265643 -2.0936936 5.381668 -21.42010 1.699626e-08 0.0001553447
## 2 ENST00000620913 -1.5056111 4.682618 -19.68078 3.365900e-08 0.0001553447
## 3 ENST00000310624 -1.3619287 4.200342 -18.05436 6.735584e-08 0.0001667669
## 4 ENST00000347703 -1.2226428 5.129111 -17.89673 7.226778e-08 0.0001667669
## 5 ENST00000465384 -1.2397054 3.362522 -15.65610 2.106592e-07 0.0004195729
## 6 ENST00000409483 -0.9553173 5.116561 -13.77842 5.807585e-07 0.0005816800
##          B signif cluster
## 1 9.590293   TRUE       B
## 2 9.101342   TRUE       B
## 3 8.575900   TRUE       B
## 4 8.521017   TRUE       B
## 5 7.653411   TRUE       B
## 6 6.777927  FALSE       B
```


# Biomart conversion of transcript names to gene names

```r
human <- useMart("ensembl", "hsapiens_gene_ensembl")
attr <- listAttributes(human)
genes.of.int <- c(as.character(a$Transcript), as.character(b$Transcript)) %>% unique()
head(genes.of.int)
```

```
## [1] "ENST00000371732" "ENST00000397406" "ENST00000319340" "ENST00000366903"
## [5] "ENST00000290705" "ENST00000265643"
```


```r
anno <- getBM(attributes = c("ensembl_transcript_id", "hgnc_symbol"), filters = "ensembl_transcript_id", values = genes.of.int, mart = human)

a$Gene_name <- anno[match(a$Transcript, anno$ensembl_transcript_id),]$hgnc_symbol
b$Gene_name <- anno[match(b$Transcript, anno$ensembl_transcript_id),]$hgnc_symbol

head(a)
```

```
##        Transcript     logFC  AveExpr         t      P.Value    adj.P.Val
## 1 ENST00000371732 -3.058143 3.501391 -27.98200 1.946558e-09 3.593540e-05
## 2 ENST00000397406 -2.189460 4.101586 -24.14325 6.451672e-09 5.955216e-05
## 3 ENST00000319340 -1.626931 4.447705 -22.94618 9.741003e-09 5.994289e-05
## 4 ENST00000366903  2.038695 3.005491  21.68305 1.540036e-08 7.107653e-05
## 5 ENST00000290705 -2.817652 3.403743 -17.61194 8.219180e-08 2.528905e-04
## 6 ENST00000265643 -1.606772 5.381668 -16.43852 1.427421e-07 2.725433e-04
##           B signif cluster Gene_name
## 1 11.105117   TRUE       B     CARD9
## 2 10.366337   TRUE       B    TSPAN4
## 3 10.090271   TRUE       B    CHST13
## 4  9.770508   TRUE       D       HLX
## 5  8.492993   TRUE       B      MT1A
## 6  8.037556   TRUE       B       GAL
```

```r
head(b)
```

```
##        Transcript      logFC  AveExpr         t      P.Value    adj.P.Val
## 1 ENST00000265643 -2.0936936 5.381668 -21.42010 1.699626e-08 0.0001553447
## 2 ENST00000620913 -1.5056111 4.682618 -19.68078 3.365900e-08 0.0001553447
## 3 ENST00000310624 -1.3619287 4.200342 -18.05436 6.735584e-08 0.0001667669
## 4 ENST00000347703 -1.2226428 5.129111 -17.89673 7.226778e-08 0.0001667669
## 5 ENST00000465384 -1.2397054 3.362522 -15.65610 2.106592e-07 0.0004195729
## 6 ENST00000409483 -0.9553173 5.116561 -13.77842 5.807585e-07 0.0005816800
##          B signif cluster Gene_name
## 1 9.590293   TRUE       B       GAL
## 2 9.101342   TRUE       B   SLC38A5
## 3 8.575900   TRUE       B      NEFH
## 4 8.521017   TRUE       B     CDCA7
## 5 7.653411   TRUE       B      FUOM
## 6 6.777927  FALSE       B    PTPRN2
```


# Make GSEA expression datasets (for common effect genes only)  

> It's typical to make expression datasets for everything, and then let GSEA sort out which are significant and the same in all files. BUT, I want to focus on the things that are the same between the 2 shFBXO11 constructs, so I'm going to limit the expression data to the common effect genes.  


```r
# load the expr data and metadata
load("./FBXO11/expr.RData")
load("./FBXO11/libs.RData")

# repeat add design factors
libs$treatment <- ifelse(test = grepl(pattern = "shSCR", x = libs$specimen_subset_external_id), yes = "shSCR", no = 
                           ifelse(test = grepl(pattern = "shFBXO11#9", x = libs$specimen_subset_external_id), yes = "shFBXO11.9", no = "shFBXO11.10"))

# select expression data for common genes only
## first test if expr_a and expr_b genes are truly all same
genes_a <- dplyr::select(a, Transcript) %>% arrange(Transcript) %>% unlist()
genes_b <- dplyr::select(b, Transcript) %>% arrange(Transcript) %>% unlist()

identical(genes_a, genes_b)
```

```
## [1] TRUE
```

```r
## now limit the expr data frame to only transcripts covered in a & b (identical)
expr_common <- expr[which(expr$Name %in% a$Transcript),]
expr_common$Name <- anno[match(expr_common$Name, anno$ensembl_transcript_id),]$hgnc_symbol
expr_common$Description <- "na"

# make GSEA expression set format
expr_set <- dplyr::select(expr_common, Name, Description, everything()) %>%
  filter(Name != "")

head(expr_set)
```

```
## # A tibble: 6 x 11
##   Name  Description A54945 A54947 A54946 A54944 A54940 A54942 A54941 A54939
##   <chr> <chr>        <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
## 1 ISG15 na           14.6   20.7   22.8   23.8   25.7   13.3   19.4   17.7 
## 2 PRXL… na            7.47   3.20   4.82   2.93   5.02   6.99   3.70  11.1 
## 3 PRXL… na            7.14   5.23   3.89   3.68   3.60   8.86   3.48   6.17
## 4 RPL22 na            8.97   4.90   5.71   4.89   3.44   7.78   4.50   9.12
## 5 ACOT7 na           64.0   54.9   53.6   52.5   54.2   66.5   49.8   66.5 
## 6 VAMP3 na           19.2   26.1   28.2   26.0   29.8   20.5   27.1   20.8 
## # ... with 1 more variable: A54943 <dbl>
```


```r
write.table(x = expr_set, file = "./FBXO11/GSEA_onLimma/GSEA_expr_set_commonEffectOnly.txt", sep = "\t", row.names = F, quote = F)
```

## And .cls file for expression analysis

```r
meta <- dplyr::select(libs, library_name, treatment)

numbers <- c(nrow(meta), length(unique(meta$treatment)), 1) # 1st line of .cls file: samples, groups, (always 1)

names <- unique(meta$treatment) %>% unlist() # 2nd line of .cls file: names of groups

order <- data.frame(sample_order = colnames(expr_set[,3:ncol(expr_set)])) # make sure the order of expr_set columns and .cls file order are same!
order$treatment <- meta[match(order$sample_order, meta$library_name),]$treatment

groups <- unlist(order$treatment) # 3rd line of .cls file: group assignment of each sample (in the same order as colnames in exprset)

# combine into 1 .cls file
write(x = numbers, file = "./FBXO11/GSEA_onLimma/Phenotypes_commonEffectOnly.cls", ncolumns = length(numbers))
write(x = names, file = "./FBXO11/GSEA_onLimma/Phenotypes_commonEffectOnly.cls", ncolumns = length(names), append = TRUE)
write(x = groups, file = "./FBXO11/GSEA_onLimma/Phenotypes_commonEffectOnly.cls", ncolumns = length(groups), append = TRUE)
```

