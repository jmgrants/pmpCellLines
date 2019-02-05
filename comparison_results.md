# Compare FBXO11 transcriptome and proteomics results
Jennifer Grants  
2/4/2019  



# Predictions  

* Proteins that are ubiquitination targets may be detected as FBXO11 binding targets.  

* The network of FBXO11 protein targets may share common functions with the network of genes upregulated in shFBXO11 cells. i.e. Loss of FBXO11 produces gene expression changes due to increased abundance of FBXO11 protein targets.  



# Ubiquitination/binding overlap

```r
load("./FBXO11/Proteomics/rda/ub.target.list.rda")
load("./FBXO11/Proteomics/rda/binding.target.rda")
```


```r
overlap <- intersect(ub.target.list, binding.target)
overlap
```

```
## character(0)
```

> Conclusion: There is no overlap between the Ub target list and the binding target list. They are very small lists, and Ub interactions could be too short-lived/weak to be detected in the coIP.  


# Prepare proteomics lists for STRING analysis  

```r
ub.targets <- data.frame(Target = ub.target.list)
write.csv(x = ub.targets, file = "./FBXO11/Proteomics/Results_summaries/diGly_result_ub_targets.csv", row.names = F)

bd.targets <- data.frame(Binding = binding.target)
write.csv(x = bd.targets, file = "./FBXO11/Proteomics/Results_summaries/coIP_result_binding_partners.csv", row.names = F)
```

> Performed STRING analysis on protein-protein interactions, then Autoannotate based on MCL clusters analysis (4 Biggest Words) of STRING "description".

*  See: `STRING_onProteomics` folder.
