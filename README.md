# pmpCellLines  
Exploring the PMP cell line data  

## Purpose  
The purpose of this repository is to explore the PMP cell line data using the skills gained in the STAT540 course.  

## Workflow
* _exploration_: Identify cell line samples of interest, develop hypotheses 
    - save: `expr` (mRNA-seq tpms) and `libs` (metadata)  

* _expression_analysis_: Starting from saved `expr` and `libs`, perform DEseq2 analysis  
    - save: Results of DEseq comparisons for each condition (condition/ctrl) as .csv files in `DEseq` folder  
    
* _functional_analysis_: Starting from files in the `DEseq` folder, perform functional analysis  

    * Specific analyses:
    1. Kate's miR SPG data: folder = "./Kate/" 
    2. FBXO11 knockdown: folder = "./FBXO11/"


## Specifics of functional analysis
1. Kate's data:

* Overlap genes significantly upregulated upon miR-X SPG treatment with  miR-X exact seed match genes  
* Overlap genes significantly upregulated upon miR-143 SPG treatment with miR-143 targets predicted in TargetScan 


2. FBXO11 knockdown data:
