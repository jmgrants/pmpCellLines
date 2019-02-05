# pmpCellLines  
Exploring the PMP cell line data  

## Purpose  
The purpose of this repository is to explore the PMP cell line data using the skills gained in the STAT540 course.  

## Specific analyses:

1. Kate's miR SPG data: folder = Kate
2. FBXO11 knockdown: folder = FBXO11


## Workflow
* _exploration_: Identify cell line samples of interest, develop hypotheses 
    - save: `expr` (mRNA-seq tpms) and `libs` (metadata)  

* _expression_analysis_: Starting from saved `expr` and `libs`, perform DEseq2 analysis  
    - save: Results of DEseq comparisons for each condition (condition/ctrl), as `.csv` files in `DEseq` folder  
        - NOTE: Do not use DEseq results, as pairwise comparisons reduce power!
    - save: Results of Limma analysis (full model with 3 shRNA treatment groups), exported as separate files showing fold change for each shFBXO11 construct, as `.csv` files in `Limma` folder
    - save: Filtered results of Limma analysis, with results for gene clusters that had common effect in both shFBXO11 constructs, as `commonEffectClusters ... .csv` files in `Limma` folder
    - save: List of gene names in common effect
  
    
* _functional_analysis_: Starting from files in the `Limma` folder, perform functional analysis  
    - save: Results are in `GSEA_onLimma` folder

* _comparison_dataset_: Starting from proteomics data from Linda/Gerben/Angela in the `Proteomics` folder, compile ubiquitination targets list and binding targets list
    - save: `ub.target.list` from diGly results, and `binding.targets` from coIP results
    - save: .csv versions of the above (for doing STRING analysis) are in `Results_summaries` folders
    
* _comparison_results_: Make comparisons between `commonEffectClusters` genes (`Limma`) and `ub.target.list`, `binding.targets` (`Proteomics`)
    - save: `STRING_onProteomics` and `STRING_onImmuneLE` folders
    - save: `Hypothesis` .md file for hypothesis generated from comparing diGly Ub target network to Immune Leading Edge gene network  

    
### Specifics of functional/comparison analyses
1. Kate's data:

* Overlap genes significantly upregulated upon miR-X SPG treatment with  miR-X exact seed match genes  
* Overlap genes significantly upregulated upon miR-143 SPG treatment with miR-143 targets predicted in TargetScan 


2. FBXO11 knockdown data:

* Using `commonEffectClusters` genes, prep data and perform GSEA analysis, EnrichmentMap network generation
* From Linda/Gerben/Angela's proteomics data, identify ubiquitination targets (diGly experiment) and binding targets list (Flag-FBXO11 coIP/MS experiment)
* Compare diGly STRING network (contains MAVS signalosome cluster) to GO-term "Immune" leading edge genes STRING network, to determine

