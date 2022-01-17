# svalflux
Microbial community responses to increasing methane flux in Arctic marine sediments at Storfjordrenna, Svalbard    

Authors: Scott A. Klasek, WeiLi Hong, Marta E. Torres, Stella Ross, Katelyn Hostetler, Alexey Portnov, Friederike Gründger, and Frederick S. Colwell    
[Published in Nature Communications, November 2021](https://www.nature.com/articles/s41467-021-26549-5)   

**Reproducible workflow:**    
Raw 16S rRNA sequence reads are publicly available as [NCBI BioProject PRJNA533183](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA533183). Additional metadata data are available [here](https://github.com/sklasek/svalflux/tree/master/data).   

[Our amplicon sequence processing workflow is available here](https://github.com/sklasek/svalflux/blob/master/markdowns/1_sequence_processing.md).   

[Phyloseq](https://joey711.github.io/phyloseq/) is a convienient way to work with amplicon sequence data in R. ASVs, taxonomic annotations, phylogenetic trees, and sample information can be easily stored together in one object. If you're cool with how we've processed sequences, removed unwanted taxa, decontaminated, and pruned samples with too few reads, you can begin working with the phyloseq object [ps.frdp](https://github.com/sklasek/svalflux/blob/master/data/ps.frdp). If you're interested, the phyloseq object [ps.f](https://github.com/sklasek/svalflux/blob/master/data/ps.f) contains those potential contaminants as well as taxa and samples we've removed.  

For those interested in how we've analyzed and presented data shown in the main text of the manuscript, see [community composition bubble plots](https://github.com/sklasek/svalflux/blob/master/markdowns/2_bubble_plots.md), [panel figures showing porewater, 16S community, and ddPCR data](https://github.com/sklasek/svalflux/blob/master/markdowns/3_core_info_and_panel_figures.md), and [alpha diversity and ordinations](https://github.com/sklasek/svalflux/blob/master/markdowns/4_community_analysis.md). We also include [Supplemental figures and tables](https://github.com/sklasek/svalflux/blob/master/markdowns/6_supplemental_figures.md) and a [co-occurrence network of ASVs](https://github.com/sklasek/svalflux/blob/master/markdowns/7_network.md). 

Codes have been archived [as a DOI on zenodo here.](https://doi.org/10.5281/zenodo.5347746)  

[Licensing info, if you want to use the code here.](https://github.com/sklasek/svalflux/blob/master/markdowns/mit_license.md)