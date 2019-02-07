# pmpCellLines  
Exploring the PMP cell line data  

## Purpose  
The purpose of this repository is to explore the PMP cell line data using the skills gained in the STAT540 course.  

## Specific analyses:

1. Kate's miR SPG data: folder = Kate
2. FBXO11 knockdown: folder = FBXO11
3. TIRAP overexpression: folder = TIRAP


## Workflow pipeline
* _exploration_: Identify cell line samples of interest, develop hypotheses 
    - save: `expr` (mRNA-seq as tpm values unless otherwise stated) and `libs` (metadata)  
    - note: in TIRAP, `expr` values are RPKM  

* _expression_analysis_: Starting from saved `expr` and `libs`, perform differential expression analysis  
    - save: Results of DEseq comparisons for each condition (condition/ctrl), as `.csv` files in `DEseq` folder  
        - __NOTE: Do not use DEseq results, as splitting to pairwise comparisons reduce power!__
    - save: Results of Limma analysis (full model), exported as separate files showing fold change for each group comparison to ctrl, as `.csv` files in `Limma` folder  
    
* _functional_analysis_: Starting from files in the `Limma` folder, perform functional analysis  
    - save: Results are in `GSEA_onLimma` folder  

* _comparison_dataset_: An alternative method of looking at same effect  
    
* _comparison_results_: Compare comparison dataset with RNA-seq data to see if results agree  
    



### Specific notes on functional analyses  
1. Kate's data:

* Overlap genes significantly upregulated upon miR-X SPG treatment with  miR-X exact seed match genes  
* Overlap genes significantly upregulated upon miR-143 SPG treatment with miR-143 targets predicted in TargetScan 


2. FBXO11 knockdown data:

  * save: Filtered results of Limma analysis, with results for gene clusters that had common effect in both shFBXO11 constructs, as `commonEffectClusters ... .csv` files in `Limma` folder
  * save: List of gene names in common effect

* Using `commonEffectClusters` genes, prep data and perform GSEA analysis, EnrichmentMap network generation

3.  TIRAP overexpression data:

* Using `expr.matrix.keep` genes, which are filtered for 'higher' expressed genes (RPKM >= 0.1 in at least 3 samples)  
  * save: GSEA expression set and .cls file in `GSEA_fromLimma` folder  

### Specific notes on comparison datasets 
2.  FBXO11 knockdown data:

* Starting from proteomics data from Linda/Gerben/Angela in the `Proteomics` folder, compile ubiquitination targets list and binding targets list
  * save: `ub.target.list` from diGly results, and `binding.targets` from coIP results
  * save: .csv versions of the above (for doing STRING analysis) are in `Results_summaries` folders

* i.e. identify ubiquitination targets (diGly experiment) and binding targets list (Flag-FBXO11 coIP/MS experiment)

### Specific notes on comparison results  
2.  FBXO11 knockdown data:

* Make comparisons between `commonEffectClusters` genes (`Limma`) and `ub.target.list`, `binding.targets` (`Proteomics`)
  * save: `STRING_onProteomics` and `STRING_onImmuneLE` folders
  * save: `Hypothesis` .md file for hypothesis generated from comparing diGly Ub target network to Immune Leading Edge gene network  

* i.e. ompare diGly STRING network (contains MAVS signalosome cluster) to GO-term "Immune" leading edge genes STRING network, to determine

