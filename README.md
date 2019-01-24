# pmpCellLines  
Exploring the PMP cell line data  

## Purpose  
The purpose of this repository is to explore the PMP cell line data using the skills gained in the STAT540 course.  

## Workflow
* _exploration.Rmd_: Identify cell line samples of interest, develop hypotheses 

    * (1) Kate's miR SPG data: (folder: "./Kate/")  
    - save: `expr` (mRNA-seq tpms) and `libs` (metadata)  
    
* _expression_analysis.Rmd_: Starting from saved `expr` and `libs`, perform DEseq2 analysis  
    - save: `result` (running DEseq(dds) -- takes a long time to run!; file: deseqResult.RData)  
