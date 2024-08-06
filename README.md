# ATAC-sequencing

This repository contains the a ATAC-seq analysis as described in [Schmidbaur et al., nature communications 2022](https://www.nature.com/articles/s41467-022-29694-7). 
It is largely based on the [Harvard Informatics pipeline](https://github.com/harvardinformatics/ATAC-seq) and [Thomas Caroll's pipeline](https://rockefelleruniversity.github.io/RU_ATAC_Workshop.html). 

# Background and aim

The large scale genome rearrangement in cephalopods leads to the question how the reorganization and the emergence of novel syntenies contributed to the evolution of innovations and gene regulatory network. A method to address this question is the assay for transposase accessible chromatin sequencing (ATAC-Seq) ([Buenrostro et al., 2013](https://www.nature.com/articles/nmeth.2688)). This methods allows the epigenomic profiling by sequencing open chromatin regions in the genome. A hyperactive Tn5 transposase cleaves and tags accessible regions of the genome with sequencing adapters. ATAC-Seq data of the hawaiian bobtail squid, E. scolopes was collected by A. Kawaguchi and H. Schmidbaur. The aim of this project was to first analyse the generated data by testing different protocols and pipelines.

Aims:
1. Analysis of the ATAC-Seq data: Finding a suitable pipeline for downstream analysis.
2. Functional assessment of open chromatin regions

**Stages:** 

<p align="center">
    <img src="https://github.com/user-attachments/assets/a5ad4b84-fcc2-47e1-9e5b-0d4a54bf5b27" width="30%">
</p>

# The Analysis workflow 
<p align="center">
  <img src="https://github.com/user-attachments/assets/6e1c0533-1331-4eb4-b6df-c0b722d3ac94" width="70%">
</p>

1. Pre-processing: 
   - Adapter trimming and quality control
   - [ngmerge](https://github.com/jsh58/NGmerge), [bbduk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/), [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
2. Genome alignment: 
   - against Euprymna sc. chromosomal scale genome
   - [Bowtie2](https://github.com/BenLangmead/bowtie2)
3. Peak calling: 
   - comparison of peak calling tools
   - [Genrich](https://github.com/jsh58/Genrich), [MACS2](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html), [HOMER](http://homer.ucsd.edu/homer/ngs/peaks.html)
4. Genome reannotation:
   - Genome reannotation adding: Transcription start sites (TSS), intron-exon boundaries, 3' and 5' untranslated regions (3'/5' UTR)
   - [PASA](https://github.com/PASApipeline/PASApipeline/blob/master/docs/index.asciidoc)
5. Downstream analysis:
   - Gene ontology analysis: [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html)
   - Motif enrichment: [HOMER](http://homer.ucsd.edu/homer/ngs/peaks.html)
  
   
  
  

