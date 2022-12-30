# PATH

*Please note that this software is still in development and being updated*  

### Introduction
Here we introduce **PATH** -- the **P**hylogenetic **A**nalysis of **T**ranscriptional 
**H**eritability -- an analytical framework for quantifying cell state heritability
versus plasticity and inferring cell state transition dynamics, in the context of
*somatic evolution*, using single-cell lineage trees or phylogenies. For more information,
see our [bioRxiv preprint]().

Single-cell phylogenies represent the ancestral
relatedness of individual cells, and if annotated with additional (single-cell) 
measurements, such as transcriptional state (for instance, as provided by 
scRNAseq), present a unique opportunity to study somatic evolutionary dynamics. 
We can apply PATH's analytical tools to these types of data to, 
* (i) measure cell state heritability versus plasticity, 
* (ii) infer phenotypic transition dynamics, and further, 
* (iii) identify heritable gene modules or pathways.  

Specifically, by using PATH, we can quantify the phylogenetic dependency of a 
single-cell measurement, such as cellular phenotype or state, 
broadly defined  (*e.g.*, cell type, location, or transcriptional profile), by
calculating *phylogenetic correlations*. 

For categorical cellular phenotypes or states, such as cell type, 
we can further use PATH to infer cell state transition dynamics -- 
the probabilities of transitioning between cell states -- 
from phylogenetic correlations. This type of analysis,
in the context of measuring cell types,
can be used to map ontogenetic relationships or differentiation hierarchies, or
in the context of cellular location, to infer metastatic seeding routes. 

### Vignettes
For examples of how to apply PATH, see the vignettes listed below.

In 
[PATH analysis of idealized phylogenies](https://htmlpreview.github.io/?https://github.com/landau-lab/PATH/blob/main/vignettes/Idealized_phylogenies.html)
and 
[PATH analysis of somatic evolution](https://htmlpreview.github.io/?https://github.com/landau-lab/PATH/blob/main/vignettes/Somatic_evolution.html) 
are
demonstrations on how to use PATH to measure phylogenetic correlations
and to infer cell state transition dynamics on simulated phylogenies,
in two contexts, using *idealized* phylogenies, or phylogenies simulated
by a *sampled somatic evolutionary process*. 

Finally, going beyond the analysis of predefined cellular states, 
we can measure the phylogenetic correlations of
gene transcription,  
to search for highly heritable or phylogenetically dependent
gene modules or pathways, in an unbiased fashion.
In [PATH analysis of glioblastoma]() is a demonstration
on how to apply PATH to a published single-cell human patient-derived
glioblastoma phylogeny 
([Chaligne R, Gaiti F, Silverbush D, and Schiffman JS *et. al. Nature Genetics*. 2021.](https://doi.org/10.1038/s41588-021-00927-7)), 
to measure phylogenetic correlations, infer
transition dynamics, and to discover heritable gene modules.

