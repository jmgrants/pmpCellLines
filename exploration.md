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

> Rationale and Hypothesis: (Kate's thesis)  

* The original miR-143 sponge contains several potential repetitive seed sequences, besides that of miR-143 (seed sequence GAGAUGA), and so may target additional miR's  
* Among these putative miR seeds, Kate identified a complementary novel miR (miR-X, seed sequence CUGUAGC) that is: 
    * expressed in human AML  
    * predicted to target several genes whose protein products were upregulated in a SILAC experiment that used the OR sponge (greater Obs/Pred ratio of target upregulation than miR-143; still very low though)  
* _Hypothesis_: mi

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



