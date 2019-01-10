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


## Investigate MDS-L experiments 
> Jeremy's suggestion of most useful to our lab  
> Focus on triplicate experiment: MDS-L + DMSO vs. Len & MDS-L-Runx1-KO + DMSO vs. Len


```r
int <- filter(cells, patient_external_id == "MDS-L" & grepl(x = cells$meta, pattern = "human"))

kable(int)
```



patient_external_id   specimen_external_id   meta                                                                                                                                           is_patient 
--------------------  ---------------------  ---------------------------------------------------------------------------------------------------------------------------------------------  -----------
MDS-L                 MDSL_DMSO1             {"description": "Sergio Martinez-Hoyer human MDS cell-line MDS-L, parental cell line stimulated with DMSO control. Replicate 1 of 3."}         FALSE      
MDS-L                 MDSL_LEN1              {"description": "Sergio Martinez-Hoyer human MDS cell-line MDS-L, parental cell line stimulated with lenalidomide. Replicate 1 of 3."}         FALSE      
MDS-L                 KO_DMSO1               {"description": "Sergio Martinez-Hoyer human MDS cell-line MDS-L, RUNX1 knockout cell line stimulated with DMSO control. Replicate 1 of 3."}   FALSE      
MDS-L                 KO_LEN1                {"description": "Sergio Martinez-Hoyer human MDS cell-line MDS-L, RUNX1 knockout cell line stimulated with lenalidomide. Replicate 1 of 3."}   FALSE      
MDS-L                 MDSL_DMSO2             {"description": "Sergio Martinez-Hoyer human MDS cell-line MDS-L, parental cell line stimulated with DMSO control. Replicate 2 of 3."}         FALSE      
MDS-L                 MDSL_LEN2              {"description": "Sergio Martinez-Hoyer human MDS cell-line MDS-L, parental cell line stimulated with lenalidomide. Replicate 2 of 3."}         FALSE      
MDS-L                 KO_DMSO2               {"description": "Sergio Martinez-Hoyer human MDS cell-line MDS-L, RUNX1 knockout cell line stimulated with DMSO control. Replicate 2 of 3."}   FALSE      
MDS-L                 KO_LEN2                {"description": "Sergio Martinez-Hoyer human MDS cell-line MDS-L, RUNX1 knockout cell line stimulated with lenalidomide. Replicate 2 of 3."}   FALSE      
MDS-L                 MDSL_DMSO3             {"description": "Sergio Martinez-Hoyer human MDS cell-line MDS-L, parental cell line stimulated with DMSO control. Replicate 3 of 3."}         FALSE      
MDS-L                 MDSL_LEN3              {"description": "Sergio Martinez-Hoyer human MDS cell-line MDS-L, parental cell line stimulated with lenalidomide. Replicate 3 of 3."}         FALSE      
MDS-L                 KO_DMSO3               {"description": "Sergio Martinez-Hoyer human MDS cell-line MDS-L, RUNX1 knockout cell line stimulated with DMSO control. Replicate 3 of 3."}   FALSE      
MDS-L                 KO_LEN3                {"description": "Sergio Martinez-Hoyer human MDS cell-line MDS-L, RUNX1 knockout cell line stimulated with lenalidomide. Replicate 3 of 3."}   FALSE      


```r
write.csv(int, "./output/summary_MDSL_triplicate.csv", row.names = F)
```


> Hypothesis:  

* Loss of Runx1 reduces Lenalidomide (Len) sensitivity of MDS-L cell line by altering the gene expression response to Len.

> Aims:  

* To identify Len-responsive genes: i.e. differentially expressed DMSO vs. Len in parental MDS-L cell line.  
* To identify Runx1-dependent genes: i.e. differentially expressed in Len-treated MDS-L vs. KO (perhaps normalize to DMSO condition in both cell lines?)


## Access expression data for samples of interest

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


```r
kable(mrnas[1:5, 1:5])
```



Name                 A15372     A30673      A35631     A30490
----------------  ---------  ---------  ----------  ---------
ENST00000456328    0.266225   1.238218    0.603082   0.541529
ENST00000450305    0.000000   0.000000    0.000000   0.000000
ENST00000488147    5.534983   5.711727   10.658560   8.945182
ENST00000619216    0.000000   0.000000    0.000000   0.000000
ENST00000473358    0.000000   0.000000    0.000000   0.000000


#### Select only libraries of interest

```r
Name <- mrnas$Name
expr <- mrnas[, match(libs$library_name, colnames(mrnas))]
expr$Name <- Name
expr <- select(expr, Name, everything())

kable(head(expr))
```



Name                 A54948     A54949     A54950     A54951     A54952     A54953     A54954     A54955     A54956     A54957     A54958     A54959
----------------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------
ENST00000456328    0.074745   0.295906   0.185494   0.000000   0.231783   0.000000   0.061365   0.256001   0.200632   0.000000   0.207074   0.166321
ENST00000450305    0.000000   0.000000   0.000000   0.000000   0.000000   0.010100   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
ENST00000488147    0.684298   1.168550   1.149269   0.629464   0.981391   0.832764   0.881181   1.119525   0.746178   0.863707   0.984333   1.961214
ENST00000619216    0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
ENST00000473358    0.000000   0.000000   0.000000   0.000000   0.029053   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
ENST00000469289    0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000   0.000000







