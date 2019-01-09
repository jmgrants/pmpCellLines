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

head(specimens)
```

```
## # A tibble: 6 x 4
##   patient_external_id specimen_external_id meta  is_patient
##   <chr>               <chr>                <chr> <lgl>     
## 1 05-53               B31920               <NA>  TRUE      
## 2 100-20              B15881               <NA>  TRUE      
## 3 101-64              101-64_sample        <NA>  TRUE      
## 4 102-23              B16360               <NA>  TRUE      
## 5 102-88              U12888               <NA>  TRUE      
## 6 103-56              U12943               <NA>  TRUE
```


```r
cells <- specimens[which(specimens$is_patient == F),]

head(cells)
```

```
## # A tibble: 6 x 4
##   patient_external_id specimen_external_id meta  is_patient
##   <chr>               <chr>                <chr> <lgl>     
## 1 AML-193             AML-193              <NA>  FALSE     
## 2 AML-193             AML193               <NA>  FALSE     
## 3 GM12878             GM12878              <NA>  FALSE     
## 4 GM19240             GM19240              <NA>  FALSE     
## 5 HMC-1               HMC-1                <NA>  FALSE     
## 6 K562                K562                 <NA>  FALSE
```


```r
write.csv(cells, "./output/summary_non_patient_samples.csv", row.names = F)
```

