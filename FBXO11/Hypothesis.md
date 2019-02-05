# Hypothesis generated from comparing diGly Ub target network to Immune Leading Edge gene network  

## Observations and evidence

1.  MAVS signalosome is a (putative) Ub target of FBXO11 E3 Ub ligase complex. See: `STRING_onProteomics` folder.  __BUT it's probably not as straightforward as that, because UBXN1 is a negative regulator of MAVS activity... RETHINK!!!!__
2.  MAVS signalosome activates NF-kB (literature). Therefore, upregulation of MAVS signalosome activity upon loss of FBXO11 may cause NF-kB hyperactivation.  
3.  Custom GSEA with NF-KB target genes (https://www.bu.edu/nf-kb/gene-resources/target-genes/) shows upregulation in shFBXO11 RNA-seq (p = 0.05). See: `GSEA_onLimma` folder.  
4.  GO-term EnrichmentMap of RNA-seq shows upregulation of immune response networks in in shFBXO11 RNA-seq. See: `GSEA_onLimma` > `EnrichmentMap` folder.  
5.  To get a better picture of what is regulating "immune" genes, I took the leading edge genes from all GO terms with "Immune" that were significantly upregulated in GSEA analysis of shFBXO11 vs. shSCR RNA-seq (5 GO terms). STRING analysis showed a large network of interacting proteins, many of which were NFKB target genes (http://amp.pharm.mssm.edu/Harmonizome/gene_set/NFKB1/JASPAR+Predicted+Transcription+Factor+Targets; http://amp.pharm.mssm.edu/Harmonizome/gene_set/RELA/JASPAR+Predicted+Transcription+Factor+Targets). Among these were 2 NFKB-target TF's (PPARG and EGR1), suggesting that additional TFs could contribute to changes in gene expression signatures.  

## Hypothesis  
Loss of FBXO11 may lead to upregulation of MAVS signalosome (__MAY BE WRONG...rethink!!__), and NF-kB activation. NF-kB activation may lead to activation of immune response genes, as well as additional transcription factors including PPARG and EGR1 that may also contribute to activation of immune response genes and other gene expression responses (e.g. lipid metabolism, cell growth).