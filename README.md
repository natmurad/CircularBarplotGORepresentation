# Circular barplot representing Gene Ontology Enrichment Analysis

This code creates a circular barplot representation of Gene
Ontology data from topGO enrichment analysis. This code was
adquired from https://www.r-graph-gallery.com/circular-barplot.html
and adapted for our data. This analysis is available on
Silva-Brandao et al (2020) Transcriptome differential co-expression 
reveals distinct molecular response of fall-armyworm strains to DIMBOA.
DOI: 10.1002/ps.6051 Journal: Pest Management Science.

### File Description

This folder contains the following files:


```circularbarplotGO.R``` - script for plotting the graphs

```exp{num}.csv``` -  a table with 4 columns: valueUP (significance of GO for UPregulated), representativeUP (class of GO), valueDOWN (significance of GO for DOWNregulated), representativeDOWN (class of GO) for experiment comparison {num}.

```exp{num}.tiff``` - circular barplot for experiment comparison {num}
