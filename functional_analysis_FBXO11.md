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


# Make GSEA expression datasets

```r
# Create score based on logFC and how small the p-value is
a <- mutate(a, score = -log10(adj.P.Val)*logFC)
b <- mutate(b, score = -log10(adj.P.Val)*logFC)

# Make GSEA expression datasets
expr_a <- dplyr::select(a, Gene_name, score) %>%
  arrange(desc(score)) %>%
  filter(Gene_name != "")
expr_b <- dplyr::select(b, Gene_name, score) %>%
  arrange(desc(score)) %>%
  filter(Gene_name != "")

head(expr_a)
```

```
##    Gene_name    score
## 1        HLX 8.457065
## 2       FGL2 6.537112
## 3 HIST2H2AA4 5.510349
## 4      HBEGF 5.242530
## 5      BCAR1 5.232297
## 6       CD1D 5.042439
```

```r
head(expr_b)
```

```
##   Gene_name    score
## 1     THBS1 4.782465
## 2   DNAJC5B 4.622435
## 3     BCAR1 4.393856
## 4    SAMSN1 4.279557
## 5      CD1C 3.882832
## 6      TGM2 3.703369
```

```r
write.table(x = expr_a, file = "./FBXO11/GSEA_onLimma/GSEA_expr_set_shFBXO11-9_commonEffect.txt", sep = "\t", row.names = F, quote = F)
write.table(x = expr_b, file = "./FBXO11/GSEA_onLimma/GSEA_expr_set_shFBXO11-10_commonEffect.txt", sep = "\t", row.names = F, quote = F)
```

