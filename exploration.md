# PMP Cell Line Data Exploration
Jennifer Grants  
1/8/2019  







## Load specimen list and identify cell lines samples of interest

```r
specimens <- amlpmpdata::db_specimens
```

Add patient metadata, and then filter out patient samples:

```r
pt_meta <- amlpmpdata::db_patients

specimens$is_patient <- pt_meta[match(specimens$patient_external_id, pt_meta$external_id),]$is_patient

kable(head(specimens))
```



patient_external_id   specimen_external_id   meta   is_patient 
--------------------  ---------------------  -----  -----------
05-53                 B31920                 NA     TRUE       
100-20                B15881                 NA     TRUE       
101-64                101-64_sample          NA     TRUE       
102-23                B16360                 NA     TRUE       
102-88                U12888                 NA     TRUE       
103-56                U12943                 NA     TRUE       




```r
cells <- specimens[which(specimens$is_patient == F),]

kable(head(cells))
```



patient_external_id   specimen_external_id   meta   is_patient 
--------------------  ---------------------  -----  -----------
AML-193               AML-193                NA     FALSE      
AML-193               AML193                 NA     FALSE      
GM12878               GM12878                NA     FALSE      
GM19240               GM19240                NA     FALSE      
HMC-1                 HMC-1                  NA     FALSE      
K562                  K562                   NA     FALSE      


```r
write.csv(cells, "./output/summary_non_patient_samples.csv", row.names = F)
```


## Investigate datasets of interest  

## Kate's miR SPG experiment  


```r
int <- filter(cells, patient_external_id == "UT-7")

kable(int)
```



patient_external_id   specimen_external_id   meta                                                                                                                         is_patient 
--------------------  ---------------------  ---------------------------------------------------------------------------------------------------------------------------  -----------
UT-7                  UT-7                   NA                                                                                                                           FALSE      
UT-7                  UT-7-GFP1              {"description": "Kate Slowski human acute megakaryocytic cell line UT-7, condition GFP, replicate 1 of 3"}                   FALSE      
UT-7                  UT-7-GFP2              {"description": "Kate Slowski human acute megakaryocytic cell line UT-7, condition GFP, replicate 2 of 3"}                   FALSE      
UT-7                  UT-7-GFP3              {"description": "Kate Slowski human acute megakaryocytic cell line UT-7, condition GFP, replicate 3 of 3"}                   FALSE      
UT-7                  UT-7-OR-SPG1           {"description": "Kate Slowski human acute megakaryocytic cell line UT-7, condition Original Sponge, replicate 1 of 3"}       FALSE      
UT-7                  UT-7-OR-SPG2           {"description": "Kate Slowski human acute megakaryocytic cell line UT-7, condition Original Sponge, replicate 2 of 3"}       FALSE      
UT-7                  UT-7-OR-SPG3           {"description": "Kate Slowski human acute megakaryocytic cell line UT-7, condition Original Sponge, replicate 3 of 3"}       FALSE      
UT-7                  UT-7-143fixSPG1        {"description": "Kate Slowski human acute megakaryocytic cell line UT-7, condition miR-143fixed sponge, replicate 1 of 3"}   FALSE      
UT-7                  UT-7-143fixSPG2        {"description": "Kate Slowski human acute megakaryocytic cell line UT-7, condition miR-143fixed sponge, replicate 2 of 3"}   FALSE      
UT-7                  UT-7-143fixSPG3        {"description": "Kate Slowski human acute megakaryocytic cell line UT-7, condition miR-143fixed sponge, replicate 3 of 3"}   FALSE      
UT-7                  UT-7-XfixSPG1          {"description": "Kate Slowski human acute megakaryocytic cell line UT-7, condition miR-Xfixed sponge, replicate 1 of 3"}     FALSE      
UT-7                  UT-7-XfixSPG2          {"description": "Kate Slowski human acute megakaryocytic cell line UT-7, condition miR-Xfixed sponge, replicate 2 of 3"}     FALSE      
UT-7                  UT-7-XfixSPG3          {"description": "Kate Slowski human acute megakaryocytic cell line UT-7, condition miR-Xfixed sponge, replicate 3 of 3"}     FALSE      


```r
write.csv(int, "./Kate/summary.csv", row.names = F)
```


> Experimental design: Kate's data  

* Original sponge (OR) that targets miR-143 and miR-X  
* miR-143-specific sponge (143fix)  
* miR-X-specific sponge (Xfix)  
* vs. GFP control or untransduced UT-7 cells

> Rationale: (Kate's thesis)  

* The original miR-143 sponge contains several potential repetitive seed sequences, besides that of miR-143 (seed sequence GAGAUGA), and so may target additional miR's  
* Among these putative miR seeds, Kate identified a complementary novel miR (miR-X, seed sequence CUGUAGC) that is: 
    * expressed in human AML  
    * predicted to target several genes whose protein products were upregulated in a SILAC experiment that used the OR sponge (greater Obs/Pred ratio of target upregulation than miR-143; still very low though)  

> Aims:  

* To identify miR-143-specific and miR-X-specific target genes  
    * based on mRNA upregulation upon specific fixed sponge treatments vs. OR sponge treatment
* To compare miR-143 and miR-X targets to Kate's SILAC results  


### Access expression data for samples of interest  

```r
libs <- amlpmpdata::db_libraries[match(int$specimen_external_id, amlpmpdata::db_libraries$specimen_subset_external_id),]
```

How to retrieve expr data tables:  
We can find all mRNA expression analysis files by querying the db_paths table. This gives a separate file for each sample. MShadbolt wrote a python script to join all these files and create one big tsv containing **tpm counts** that can be loaded as follows:  

file://isaac/klabanalysis/mshadbolt/KARSANBIO-1352_Karsan_Lab_Informatics_Support/KARSANBIO-1521_mRNA_and_miRNA_data_retrieval/report/KARSANBIO-1521_mRNA_and_miRNA_data_retrieval.html  


```r
mrnas <- read_tsv("/projects/karsanlab/PMP-AML/expression/bcbio_fastrnaseq_salmon/salmon_AML_PMP_expression_output.tsv")
```

```
## Parsed with column specification:
## cols(
##   .default = col_double(),
##   Name = col_character()
## )
```

```
## See spec(...) for full column specifications.
```


### Select only libraries of interest  

```r
Name <- mrnas$Name
expr <- mrnas[,which(colnames(mrnas) %in% libs$library_name)]
expr$Name <- Name
expr <- select(expr, Name, everything())

kable(head(expr))
```



Name                 A54974     A54963     A54964     A54966     A35856     A54971     A54972     A54973     A54967     A54969     A54968     A54965     A54970
----------------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------
ENST00000456328    0.172507   0.073746   0.218676   0.110418   0.038299   0.171210   0.156859   0.231480   0.000000   0.281069   0.138982   0.000000   0.160053
ENST00000450305    0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.003676   0.000000
ENST00000488147    2.870341   4.099234   3.408244   5.851605   1.493135   4.397266   2.983613   3.232376   4.603645   3.418479   4.786116   3.122282   2.779962
ENST00000619216    0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
ENST00000473358    0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
ENST00000469289    0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000


### Save Rdata files for downstream analysis  

```r
save(expr, file = "./Kate/expr.RData")
save(libs, file = "./Kate/libs.RData")
```



## FBXO11 knockdown  

```r
int <- filter(cells, grepl(pattern = "Linda", x = cells$meta))

kable(int)
```



patient_external_id   specimen_external_id     meta                                                                                                                                    is_patient 
--------------------  -----------------------  --------------------------------------------------------------------------------------------------------------------------------------  -----------
OCI-AML3              OCI-AML3_shSCR-1         {"description": "Linda Chang / Gerben Duns human AML cell-line OCI-AML3, transduced with control vector (shSCR). Replicate 1 of 3."}    FALSE      
OCI-AML3              OCI-AML3_shFBXO11#9-1    {"description": "Linda Chang / Gerben Duns human AML cell-line OCI-AML3, transduced with shRNA#9 against FBXO11. Replicate 1 of 3."}    FALSE      
OCI-AML3              OCI-AML3_shFBXO11#10-1   {"description": "Linda Chang / Gerben Duns human AML cell-line OCI-AML3, transduced with shRNA#10 against FBXO11. Replicate 1 of 3."}   FALSE      
OCI-AML3              OCI-AML3_shSCR-2         {"description": "Linda Chang / Gerben Duns human AML cell-line OCI-AML3, transduced with control vector (shSCR). Replicate 2 of 3."}    FALSE      
OCI-AML3              OCI-AML3_shFBXO11#9-2    {"description": "Linda Chang / Gerben Duns human AML cell-line OCI-AML3, transduced with shRNA#9 against FBXO11. Replicate 2 of 3."}    FALSE      
OCI-AML3              OCI-AML3_shFBXO11#10-2   {"description": "Linda Chang / Gerben Duns human AML cell-line OCI-AML3, transduced with shRNA#10 against FBXO11. Replicate 2 of 3."}   FALSE      
OCI-AML3              OCI-AML3_shSCR-3         {"description": "Linda Chang / Gerben Duns human AML cell-line OCI-AML3, transduced with control vector (shSCR). Replicate 3 of 3."}    FALSE      
OCI-AML3              OCI-AML3_shFBXO11#9-3    {"description": "Linda Chang / Gerben Duns human AML cell-line OCI-AML3, transduced with shRNA#9 against FBXO11. Replicate 3 of 3."}    FALSE      
OCI-AML3              OCI-AML3_shFBXO11#10-3   {"description": "Linda Chang / Gerben Duns human AML cell-line OCI-AML3, transduced with shRNA#10 against FBXO11. Replicate 3 of 3."}   FALSE      



```r
write.csv(int, "./FBXO11/summary.csv", row.names = F)
```


> Experimental design: FBXO11 knockdown data  

* FBXO11 shRNA knockdown (2 different shRNA constructs)  
* vs. sh-Scrambled controls  
* (in triplicate)  


> Rationale and Hypothesis:  

* FBXO11 is part of a E3 ubiquitin ligase complex, therefore loss of FBXO11 will affect protein expression level. These changes may in turn affect gene expression, causing changes in cellular function.  
* FBXO11 is lost in some AMLs, and FBXO11 loss in mice causes AML-like disease when combined with a predisposing mutation (AML-ETO). Therefore the functional changes evoked by FBXO11 knockdown may be causal in AML development.  
* _Hypothesis_: Loss of FBXO11 causes gene expression changes that alter cellular function to promote AML development.  


> Aims:  

* To identify processes and pathways upregulated upon knockdown of FBXO11 (changes in cellular function)  
* To compare these results to proteomics data from Linda and Gerben to determine potential FBXO11-target-protein regulators of these gene expression changes  


### Access expression data for samples of interest  

```r
libs <- amlpmpdata::db_libraries[match(int$specimen_external_id, amlpmpdata::db_libraries$specimen_subset_external_id),]
```


### Select only libraries of interest  

```r
Name <- mrnas$Name
expr <- mrnas[,which(colnames(mrnas) %in% libs$library_name)]
expr$Name <- Name
expr <- select(expr, Name, everything())

kable(head(expr))
```



Name                 A54945     A54947     A54946     A54944     A54940     A54942     A54941     A54939     A54943
----------------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------
ENST00000456328    0.176085   0.560645   0.101541   0.458976   0.115138   0.306703   0.274130   0.381055   0.179624
ENST00000450305    0.609390   0.234522   0.169559   0.000000   0.141108   0.266956   0.513587   0.045370   0.000000
ENST00000488147    3.944083   4.438525   4.941163   4.643798   3.260512   2.980054   3.091821   4.443674   2.964437
ENST00000619216    0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
ENST00000473358    0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
ENST00000469289    0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000


### Save Rdata files for downstream analysis  

```r
save(expr, file = "./FBXO11/expr.RData")
save(libs, file = "./FBXO11/libs.RData")
```


## TIRAP overexpression (mouse BM)

```r
tirap <- filter(cells, grepl(pattern = "Rawa", x = cells$meta))

kable(tirap)
```



patient_external_id   specimen_external_id   meta                                                                                                                         is_patient 
--------------------  ---------------------  ---------------------------------------------------------------------------------------------------------------------------  -----------
MIG                   MIG-1                  {"description": "Rawa Ibrahim mouse bone marrow cells. Biological replicate 1 of vector control"}                            FALSE      
MIG                   MIG-2                  {"description": "Rawa Ibrahim mouse bone marrow cells. Biological replicate 2 of vector control"}                            FALSE      
MIG                   MIG-3                  {"description": "Rawa Ibrahim mouse bone marrow cells. Biological replicate 3 of vector control"}                            FALSE      
TIRAP                 TIRAP-1                {"description": "Rawa Ibrahim mouse bone marrow cells. Biological replicate 1 of TIRAP-overexpressing experimental group"}   FALSE      
TIRAP                 TIRAP-2                {"description": "Rawa Ibrahim mouse bone marrow cells. Biological replicate 2 of TIRAP-overexpressing experimental group"}   FALSE      
TIRAP                 TIRAP-3                {"description": "Rawa Ibrahim mouse bone marrow cells. Biological replicate 3 of TIRAP-overexpressing experimental group"}   FALSE      


### Use amlpmpdata package to identify path to files

```r
# merge library names to tirap table
libraries <- amlpmpdata::db_libraries

tirap$library <- libraries[match(tirap$specimen_external_id, libraries$specimen_subset_external_id),]$library_name

# merge path to tirap table  
paths <- amlpmpdata::db_paths

tirap$path <- paths[match(tirap$library, paths$library),]$path

kable(tirap)
```



patient_external_id   specimen_external_id   meta                                                                                                                         is_patient   library   path 
--------------------  ---------------------  ---------------------------------------------------------------------------------------------------------------------------  -----------  --------  -----
MIG                   MIG-1                  {"description": "Rawa Ibrahim mouse bone marrow cells. Biological replicate 1 of vector control"}                            FALSE        A54933    NA   
MIG                   MIG-2                  {"description": "Rawa Ibrahim mouse bone marrow cells. Biological replicate 2 of vector control"}                            FALSE        A54934    NA   
MIG                   MIG-3                  {"description": "Rawa Ibrahim mouse bone marrow cells. Biological replicate 3 of vector control"}                            FALSE        A54935    NA   
TIRAP                 TIRAP-1                {"description": "Rawa Ibrahim mouse bone marrow cells. Biological replicate 1 of TIRAP-overexpressing experimental group"}   FALSE        A54936    NA   
TIRAP                 TIRAP-2                {"description": "Rawa Ibrahim mouse bone marrow cells. Biological replicate 2 of TIRAP-overexpressing experimental group"}   FALSE        A54937    NA   
TIRAP                 TIRAP-3                {"description": "Rawa Ibrahim mouse bone marrow cells. Biological replicate 3 of TIRAP-overexpressing experimental group"}   FALSE        A54938    NA   

> Paths are not available in this data package. Therefore I searched JIRA tickets and found Jenny's RPKM data for these samples, and saved to `TIRAP` / `data`.

* JIRA ticket: https://www.bcgsc.ca/jira/browse/KARSANBIO-412

### Load Jenny's pre-calculated RPKM values

```r
expr <- read_tsv(file = "./TIRAP/data/RPKM.tsv", col_names = T)
```

```
## Parsed with column specification:
## cols(
##   gene_id = col_character(),
##   chromosome = col_integer(),
##   start = col_integer(),
##   end = col_integer(),
##   strand = col_character(),
##   gene_symbol = col_character(),
##   biotype = col_character(),
##   gene_description = col_character(),
##   `A54933_normalized_coverage_(RPKM)` = col_double(),
##   `A54934_normalized_coverage_(RPKM)` = col_double(),
##   `A54935_normalized_coverage_(RPKM)` = col_double(),
##   `A54936_normalized_coverage_(RPKM)` = col_double(),
##   `A54937_normalized_coverage_(RPKM)` = col_double(),
##   `A54938_normalized_coverage_(RPKM)` = col_double()
## )
```

```
## Warning in rbind(names(probs), probs_f): number of columns of result is not
## a multiple of vector length (arg 1)
```

```
## Warning: 3141 parsing failures.
## row # A tibble: 5 x 5 col     row col        expected   actual     file                    expected   <int> <chr>      <chr>      <chr>      <chr>                   actual 1 35075 chromosome an integer GL456210.1 './TIRAP/data/RPKM.tsv' file 2 35076 chromosome an integer GL456210.1 './TIRAP/data/RPKM.tsv' row 3 35077 chromosome an integer GL456210.1 './TIRAP/data/RPKM.tsv' col 4 35078 chromosome an integer GL456210.1 './TIRAP/data/RPKM.tsv' expected 5 35079 chromosome an integer GL456210.1 './TIRAP/data/RPKM.tsv'
## ... ................. ... ................................................................ ........ ................................................................ ...... ................................................................ .... ................................................................ ... ................................................................ ... ................................................................ ........ ................................................................
## See problems(...) for more details.
```

```r
dim(expr) # should be 38215 rows according to RPKM.tsv
```

```
## [1] 38215    14
```

```r
# fix column names for samples
cols <- data.frame(cols_names = colnames(expr[,9:ncol(expr)])) %>%
  separate(col = cols_names, into = "name", remove = T)
```

```
## Warning: Expected 1 pieces. Additional pieces discarded in 6 rows [1, 2, 3,
## 4, 5, 6].
```

```r
cols2 <- data.frame(name = colnames(expr[,1:8]))
cols3 <- rbind(cols2, cols)

colnames(expr) <- cols3$name

# get expression matrix only
expr <- dplyr::select(expr, gene_id, gene_symbol, contains("A54"))
save(expr, file = "./TIRAP/expr.rda")
```

### Save metadata

```r
libs <- amlpmpdata::db_libraries[match(tirap$specimen_external_id, amlpmpdata::db_libraries$specimen_subset_external_id),] 
save(libs, file = "./TIRAP/libs.rda")
```


> Experimental design: TIRAP overexpression data

* Bone marrow transfected with TIRAP overexpression vector or MIG empty vector, transplanted to n = 3 recipient mice
* Compare gene expression in TIRAP overexpression vs. empty vector

> Hypothesis:

* TIRAP overexpression causes gene expression changes that account for the bone marrow failure phenotype observed in these mice

> Aims:

* To identify differentially expressed genes in TIRAP-oe vs. Control (Limma)
* To identify up- and downregulated pathways and processes in TIRAP-oe vs. Control (GSEA, EnrichmentMap)
