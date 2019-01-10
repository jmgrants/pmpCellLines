# Expression analysis in MDS-L / Runx1 KO treated with Len
Jennifer Grants  
1/9/2019  





## Load saved data
`expr` is mRNA-seq TPM values, `libs` is metadata

```r
load("expr.RData")
load("libs.RData")
```


```r
libs$Genotype <- ifelse(grepl(pattern = "MDSL", x = libs$specimen_subset_external_id), yes = "WT", no = "Runx1KO")
libs$Treatment <- ifelse(grepl(pattern = "LEN", x = libs$specimen_subset_external_id), yes = "Len", no = "Ctrl")

libs$Genotype <- factor(libs$Genotype, levels = c("WT", "Runx1KO"))
libs$Treatment <- factor(libs$Treatment, levels = c("Ctrl", "Len"))
```



```r
kable(libs)
```



library_name   sequencing_effort     platform_name   platform_version   specimen_subset_external_id   cohort_no   library_qc_info   index_sequence   Genotype   Treatment 
-------------  --------------------  --------------  -----------------  ----------------------------  ----------  ----------------  ---------------  ---------  ----------
A54948         Karsan Lab Research   ssRNA-Seq       v1                 MDSL_DMSO1                    Cohort 6    {}                ACCGGC           WT         Ctrl      
A54949         Karsan Lab Research   ssRNA-Seq       v1                 MDSL_LEN1                     Cohort 6    {}                AGGCCG           WT         Len       
A54950         Karsan Lab Research   ssRNA-Seq       v1                 KO_DMSO1                      Cohort 6    {}                CAACTA           Runx1KO    Ctrl      
A54951         Karsan Lab Research   ssRNA-Seq       v1                 KO_LEN1                       Cohort 6    {}                CCACGC           Runx1KO    Len       
A54952         Karsan Lab Research   ssRNA-Seq       v1                 MDSL_DMSO2                    Cohort 6    {}                CTATAC           WT         Ctrl      
A54953         Karsan Lab Research   ssRNA-Seq       v1                 MDSL_LEN2                     Cohort 6    {}                GCAAGG           WT         Len       
A54954         Karsan Lab Research   ssRNA-Seq       v1                 KO_DMSO2                      Cohort 6    {}                TACAGC           Runx1KO    Ctrl      
A54955         Karsan Lab Research   ssRNA-Seq       v1                 KO_LEN2                       Cohort 6    {}                TGCCAT           Runx1KO    Len       
A54956         Karsan Lab Research   ssRNA-Seq       v1                 MDSL_DMSO3                    Cohort 6    {}                ATGTCA           WT         Ctrl      
A54957         Karsan Lab Research   ssRNA-Seq       v1                 MDSL_LEN3                     Cohort 6    {}                TTAGGC           WT         Len       
A54958         Karsan Lab Research   ssRNA-Seq       v1                 KO_DMSO3                      Cohort 6    {}                GGCTAC           Runx1KO    Ctrl      
A54959         Karsan Lab Research   ssRNA-Seq       v1                 KO_LEN3                       Cohort 6    {}                AAGACT           Runx1KO    Len       


## Perform DEseq2 analysis to identify differentially expressed genes in WT+Ctrl vs WT+Len only

```r
expr.wt <- expr[, which(colnames(expr) %in% c(libs[which(libs$Genotype == "WT"),]$library_name, "Name"))]

counts <- column_to_rownames(expr.wt, var = "Name") %>% apply(1:2, round) # count values must be non-negative integers for DESeq Dataset
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```r
cols <- select(libs, library_name, specimen_subset_external_id, Genotype, Treatment) %>%
  filter(Genotype == "WT")
```


```r
identical(colnames(counts), cols$library_name)
```

```
## [1] TRUE
```



```r
dds <- DESeqDataSetFromMatrix(countData = counts, colData = cols, design = ~Treatment)
```

```
## converting counts to integer mode
```

Run DEseq:

```r
# this chunk takes a long time to run; don't run it unless you update something
result <- DESeq(dds)

save(result, file = "deseqResult.RData")
```


```r
load("deseqResult.RData")
```

Run DEseq final result table:

```r
deseq.table  <- results(result)

deseq.table <- as.tibble(deseq.table) %>%
  rownames_to_column("gene_id") %>%
  na.omit()
```

```
## Warning in as.data.frame(x, row.names = NULL, optional = optional, ...):
## Arguments in '...' ignored
```

```r
range(deseq.table$log2FoldChange)
```

```
## [1] -0.6640021  1.4027618
```







