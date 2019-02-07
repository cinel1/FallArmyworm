# FallArmyworm
Analysis pipeline for RNAseq read QC, assembly, annotation, differential expression, and functional interpretation analyses on brain tissue from chronically bat call-exposed male fall armyworm moths

Abridged pipeline:
  1) Raw RNA read QC and trimming via FastQC and Trimmomatic
  2) Read in silico normalization in Trinity followed by assembly.
  3) Transcriptome annotation via Annocript
  4) Sample-specific reads mapped and quantified against reference with Kallisto
  5) Differential expression analysis in R using EdgeR, limma, voom, and Ebayes
  6) Gene Ontology enrichment analysis using BiNGO in the Cytoscape software suite
  7) KEGG pathway analysis using GhostKOALA
