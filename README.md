# FallArmyworm
Analysis pipeline for RNAseq read QC, assembly, annotation, differential expression, and functional interpretation analyses on brain tissue from chronically bat call-exposed male fall armyworm moths

Abridged pipeline:
  1) Raw RNA read QC and trimming via FastQC and Trimmomatic
  2) Read in silico normalization in Trinity followed by multiple assembly using Combined De Novo Transcriptome Assembly pipeline provided by the National Center for Genome Analysis Support (see https://ncgas.org/files/PAG%20Transcriptome%20Pipeline%20Demo-final.pdf for more info).
  3) Composite, non-redundant reference transcriptome generated using EviGenes (see http://arthropods.eugenes.org/EvidentialGene/trassembly.html for more info)
  4) Composite transcriptome annotation via Annocript
  5) Sample-specific reads mapped and quantified against reference with Kallisto
  6) Differential expression analysis in R using EdgeR, limma, voom, and Ebayes
  7) Gene Ontology enrichment analysis using BiNGO in the Cytoscape software suite
  8) KEGG pathway analysis using GhostKOALA
