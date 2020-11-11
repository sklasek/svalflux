# svalflux
**Microbial community responses to increasing methane flux in Arctic marine sediments at Storfjordrenna, Svalbard**
*Authors:* Scott A. Klasek, WeiLi Hong, Marta E. Torres, Stella Ross, Katelyn Hostetler, Alexey Portnov, Friederike Gründger, and Frederick S. Colwell
*Publication:* not published yet.

*Reproducible workflow:*
If you would like to reproduce this workflow, starting from raw 16S reads, they are publicly available as NCBI BioProject PRJNA533183. 
If you're cool with the way we've processed our sequences, you can skip the NCBI download and begin working with the phyloseq object ps.f in /data.
If you're cool with how we've removed unwanted sequences, decontaminated, and pruned samples with too few reads, you can begin working with the phyloseq object ps.frdp in /data.
Porewater geochemistry, modeled AOM rates, and ddPCR gene count data are available in the /data folder. 
R code is available in /markdowns that can be followed in numerical order (1 through 4, from sequence processing to community analysis). 
If you're just here to see the figures without the code, they are in documents numbered 5 and 6 (main figures and supplemental figures). 

Please reach out if you have questions!
-Scott