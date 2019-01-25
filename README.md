# pmpCellLines  
Exploring the PMP cell line data  

## Purpose  
The purpose of this repository is to explore the PMP cell line data using the skills gained in the STAT540 course.  

## Workflow
* _exploration.Rmd_: Identify cell line samples of interest, develop hypotheses 

    * (1) Kate's miR SPG data: (folder: "./Kate/")  
    - save: `expr` (mRNA-seq tpms) and `libs` (metadata)  
    
* _expression_analysis.Rmd_: Starting from saved `expr` and `libs`, perform DEseq2 analysis  
    - save: Results of DEseq comparisons for each construct (construct/GFP ctrl) as .csv files in `DEseq` folder  
    
* _functional_analysis.Rms_: Starting from files in the `DEseq` folder, perform functional analysis  


## Specifics of functional analysis
1. Kate's data:

* Overlap genes significantly upregulated upon miR-X SPG treatment with  miR-X exact seed match genes  
* Overlap genes significantly upregulated upon miR-143 SPG treatment with miR-143 targets predicted in TargetScan  
