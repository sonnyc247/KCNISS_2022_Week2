# KCNI Summer School 2022
## Project 2 - Single Cell Transcriptomics

**Integrative analysis of mouse and human single-cell gene expression data.** What are conserved cell types in the brain and what are their characteristics?

### What’s this project about? 

**Main idea:** perform integrative analysis of mouse and human (and possibly other species, including macaques and marmosets) neocortical cell types according to transcriptomics and possibly intrinsic electrophysiology. These analyses will directly feed into cell and circuit models of mouse and human circuits, led by Etay’s group.

**Key questions:**

1. Can we identify orthologous cell types between species?
2. What genes / features distinguish cell types? 
    1. Are there aspects of these features that seem relevant to computational processing of cells and circuits?
    2. What characteristics of these features can be modelled in the context of cell and circuit models?

**What (dataset) resources are available to help answer this question?**
[Allen Institute for Brain Sciences Cell Types database](https://celltypes.brain-map.org/)

**Transcriptomics TAs:** Sonny Chen, Mel Davie

## Resources

Dataset information:
* [Region mapping between species](https://github.com/sonnyc247/KCNISS_2022_Week2/blob/main/Data/Region_Mapping.csv)
    * Includes number of cells available within each region for full human and mouse datasets

Seurat tutorials:
* [Dataset integration workflow](https://satijalab.org/seurat/articles/integration_introduction.html)
* [Differential expression analysis](https://satijalab.org/seurat/articles/de_vignette.html)
* [Data visualization methods](https://satijalab.org/seurat/articles/visualization_vignette.html)

Getting started:
1. Preprocessing completed using code [here](https://github.com/sonnyc247/KCNISS_2022_Week2/blob/main/Code/Preprocessing.R)
    * Finding homologous genes between species, downsampling large datasets, filtering out undetected genes
2. [Initial Seurat analysis](https://github.com/sonnyc247/KCNISS_2022_Week2/blob/main/Code/Processing.R)
    * Dataset integration, *de novo* clustering, visualization
