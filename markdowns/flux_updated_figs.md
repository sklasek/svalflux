Storfjordrenna Manuscript – Updated Figures
================
Scott Klasek
9-2-20

## Data repository

I’ve uploaded this data to my Github page. This repository is private,
but (fingers crossed) upon publication, I can set it to public so that
anyone who wants to reproduce this work or make similar figures can look
at the code, import the data, and easily do so.

## The map is unchanged:

![Figure 1. Bathymetric map of Storfjordrenna gas hydrate mounds and
modeling approach. (A) Storfjordrenna is located south of the Svalbard
Archipelago. Locations of cores and gas hydrate mounds (GHMs) are shown.
White polygons indicate areas of seafloor gas release observed at the
time of the cruise. All cores were collected from GHMs, with the
exception of the reference core GC1048, which was sampled 400 m west of
GHM5. The schematic in (B) depicts sulfate and methane concentrations
throughout a sediment column at a steady-state condition. As methane
flux increases, (C), SR-AOM is stimulated at shallower depths and
sulfate profiles show a concave-up curvature. After decades of steadily
increasing methane flux at a particular area, (D), reactive-transport
modeling can be used to estimate how quickly the methane front diffused
the distance between the current and prior sulfate depletion depths
(indicated by the orange bracket).](figures/F1.png)

## Panel figures

Porewater, microbial community, and gene count data from the seep site
(core PC1029). AOM rate data has been
removed.

<img src="flux_updated_figs_files/figure-gfm/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

Figure 2. Geochemical, microbial community, and gene abundance data from
an active seep site. Push core PC1029 is located at the seep in the
center of GHM3. (A) shows porewater sulfate, sulfide, and alkalinity
concentrations in addition to frenulate siboglinid tubeworms and gas
hydrate nodules several cm in diameter recovered in a replicate core at
40-50 cm below seafloor. (B) depicts percent abundances of dominant
Bacterial and Archaeal classes within the microbial community (left
panel), dominant anaerobic methanotrophic archaeal (ANME) families
(center panel) and sulfate-reducing bacterial (SRB) genera (right
panel). (C) shows log10 copy numbers of mcrA and dsrAB genes per gram
bulk sediment.

And from cores where methane flux is increasing (GC1045 and
GC1081):

<img src="flux_updated_figs_files/figure-gfm/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

Figure 3. Geochemical, numerical, microbial community, and gene
abundance data from two sites showing sulfate-methane dynamics
suggestive of recent increases in methane flux. Gravity cores GC1045 and
GC1081 are located at GHMs 3 and 4, respectively. Sulfate-methane
transition depths are indicated by dashed lines. (A) shows porewater
sulfate, sulfide, and alkalinity, and (B) the temporal progression of
modeled AOM rates from 10 years ago to up to 2 years after sampling. (C)
indicates percent abundances of dominant bacterial and archaeal classes,
dominant anaerobic methanotrophic archaeal (ANME) familes, and
sulfate-reducing bacterial (SRB) genera. (D) shows copy numbers of mcrA
and dsrAB genes per gram bulk sediment.

And from four steady-state
cores:

<img src="flux_updated_figs_files/figure-gfm/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

Figure 4. Geochemical, numerical, microbial community, and gene
abundance data from four sites showing steady-state sulfate-methane
dynamics. GC1068–1070 are from GHM5, and reference core GC1048 is
located to the west of GHM5. Sulfate-methane transition depths are
indicated by dashed lines. (A) shows porewater sulfate, sulfide, and
alkalinity, and (B) present-day modeled AOM rates. (C) indicates percent
abundances of dominant bacterial and archaeal classes, dominant
anaerobic methanotrophic archaeal (ANME) familes, and sulfate-reducing
bacterial (SRB) genera.(D) shows copy numbers of mcrA and dsrAB genes
per gram bulk sediment, with values below the detectable limit (10^3
g-1) along the margin of the panel. Macroscopic translucent-to-yellow
biofilms, shown as yellow symbols in panel (A), were observed at
approximate SMT depths in cores GC1048 and GC1070 (symbol size not to
scale with depth axis).

## Differentially abundant ASVs across zones and methane flux types:

This plot shows differentially abundant ASVs, replacing the old figure
S2. Each point represents an ASV whose relative abundance is different
between communities in cores of different methane stages (more abundant
in one stage as compared to the other two combined). These stages are
further subset into above-SMT and below-SMT communities and facetted
into taxonomic Families (for ANMEs), Genera (for Deltaproteobacteria),
or classes (for all others). NAs at each taxonomic category of plotting
were omitted from this graph. Note there were no samples taken below
SMTZ at the seep.

As before, we still observe ASV-scale niche differences among ANME and
SRB subpopulations, but this plot allows us to see a few other things:
1) ASVs associated with seepage are often dominant within communities,
and include JS1 Atribacteria, several Gammaproteobacteria and
Campylobacteria, and many genera of Deltaproteobacteria  
2\) Calditrichia may be associated with recent upwards methane
migration  
3\) Many ASVs of Aminicenantia and Dehalococcoidia are indicative of
steady-state sulfate reduction zones  
4\) Far fewer biomarkers distinguish stages below the SMTZ (and none are
unique to below-SMTZ) but Lokiarchaea may be one key taxon associated
with steady-state methanogenic
zones.

<img src="flux_updated_figs_files/figure-gfm/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

Figure 5. Differentially abundant ASVs between three discrete stages of
methane dynamics, either above or below the sulfate-methane transition
zone. Higher differential abundance values represent stronger
associations between an ASV and a particular group. In addition to
position relative to the SMTZ, vertical panels convey class-level
taxonomic annotations for all taxa except ANME (at family level) and SRB
(at genus level). Bubble sizes represent mean percent abundances within
each group. Differential abundance was inferred using an alpha of 0.05
and a Benjamini-Hochberg correction for multiple comparisons.

### Alpha-diversity (diversity within communities)

This figure is essentially the same, but I’m wondering if it’s fair to
include PC1029 samples in A if we can’t really estimate a peak AOM rate…
Think about this some
more.  
<img src="flux_updated_figs_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

Figure 6. Microbial community diversity patterns across peak modeled AOM
depths. Shannon diversity indices of microbial communities for
individual samples plotted by their distance above (positive) or below
(negative) the depths corresponding to present-day maximum AOM rates
across all cores. Cores are divided by panel based on whether methane
flux is (A) high, or (B) low, and colored according to different stages
of methane dynamics. Dotted horizontal lines show the distance interval
corresponding to high methane flux samples. Multiple R^2 and slope
p-values are shown for linear regressions of points within these
intervals.

### Beta-diversity (diversity between communities)

The two most helpful categories that explain the most of the variance
within these microbial communities are methane stage and redox zone.
There is a significant interaction among them
too.

<img src="flux_updated_figs_files/figure-gfm/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

Figure 7. Principal Coordinates Analysis (PCoA) ordination of weighted
Unifrac distances between all communities, colored according to stages
of methane dynamics (A), and redox zone within the sediment column (B).
A Hellinger transformation was applied to ASV count tables before
calculating the distance matrix. PERMANOVA tests verify distinct
community structures present across stages and redox zones, with R^2 and
p-values shown in corresponding plots.

### An updated conceptual figure

![](figures/flux.conceptual.2.png)

Figure 8. Conceptual depiction of microbial community changes concurrent
with distinct stages of methane dynamics at Storfjordrenna GHMs. Methane
and sulfate profiles are shown in gray and blue lines, respectively,
with microbial community changes indicated by blowup circles. ANME and
SRB represent red and green ovals, with other bacteria and archaea in
gray (shapes representing the diversity of other taxa). (A) Methane
seepage, potentially driven by dissociation of gas hydrates at the upper
limit of stability, stimulates high rates of AOM and densities of
ANME/SRB. Sulfide fluxes from AOM are sufficient to support frenulate
siboglinid tubeworms, which may irrigate sulfate over the upper several
cm of sediment, increase redox gradients in underlying sediments, and
further ANME/SRB growth. (B) In the absence of bioirrigation, a
concave-up bend in porewater sulfate suggests recent methane migration
into the sulfate-reduction zone that may have resulted from a pulse of
methane beginning around 290 years ago. Methane travels upward
throughout the sediment column, and ANME/SRB growth follows with less
than a year of lag time, driving down alpha diversity. (C) Steady-state
sulfate profiles suggest no recent methane influx, and low rates of AOM
are observed at the SMTZ. Stable conditions may allow for higher
microbial diversity and the growth of macroscopic biofilms.

### Supplemental Figures

I’m not sure about my values for these flux calculations… but basically
increasing methane flux \>\> steady-state

![](flux_updated_figs_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->
Figure S1. Methane fluxes across each core, integrated from modeled AOM
rates derived from data at the time of sampling. Cores are colored by
methane stage. PC1029 was omitted due to high uncertainty in modeling
AOM rates (explain more?)

![](flux_updated_figs_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Figure S2. Regression of mcrA gene copy numbers (per gram bulk sediment)
to modeled AOM rates shows a linear relationship across samples from all
cores (log-log transformation). Samples are colored according to stages
of methane dynamics, and those that did not contain detectable mcrA were
omitted.

    ## Warning: Removed 5 rows containing missing values (geom_point).

![](flux_updated_figs_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Figure S3. Empirically-measured sulfate porewater profiles from cores
taken from steady-state areas, shown with modeled sulfate profiles
assuming a non-steady state scenario where methane flux is increasing.
The discrepancies in the shapes of these two profiles for these four
cores is evidence that they are not experiencing a recent increase in
methane flux, and are thus characterized as steady-state.

### Supplemental Tables

    ##     Core   Latitude  Longitude Water depth (m) Core recovery (cm)
    ## 1 PC1029  76 06.398  15 58.151             381                 27
    ## 2 GC1045  76 06.347  15 57.959             387                130
    ## 3 GC1081 76 07.022  16 02.593              369                102
    ## 4 GC1048 76 06.737  15 59.845              387                335
    ## 5 GC1068 76 06.739  16 00.311              384                295
    ## 6 GC1069 76 06.719  16 00.334              383                227
    ## 7 GC1070 76 06.703  16 00.162              385                326
    ##   SMT depth (cm)
    ## 1             NA
    ## 2             82
    ## 3             56
    ## 4            320
    ## 5            108
    ## 6            138
    ## 7             69

Table S1. Latitude, longitude, water depth, core recovery, and
sulfate-methane transition depth of all cores analyzed in this study.

In this table below, I realize that the total-core flux values do not
match those that I showed in the first supplemental figure… so I must
have made a
    mistake?

    ##      Core Year before sampling CH4 flux mols m^-2 yr ^-1 Peak AOM depth (cm)
    ## 1  GC1045                   21                     3.674               280.0
    ## 2  GC1045                   20                     3.744               270.0
    ## 3  GC1045                   19                     3.808               257.5
    ## 4  GC1045                   18                     3.874               247.5
    ## 5  GC1045                   17                     3.944               237.5
    ## 6  GC1045                   16                     4.021               227.5
    ## 7  GC1045                   15                     4.101               217.5
    ## 8  GC1045                   14                     4.184               207.5
    ## 9  GC1045                   13                     4.268               197.5
    ## 10 GC1045                   12                     4.351               187.5
    ## 11 GC1045                   11                     4.433               177.5
    ## 12 GC1045                   10                     4.514               167.5
    ## 13 GC1045                    9                     4.593               157.5
    ## 14 GC1045                    8                     4.671               147.5
    ## 15 GC1045                    7                     4.749               137.5
    ## 16 GC1045                    6                     4.826               127.5
    ## 17 GC1045                    5                     4.902               117.5
    ## 18 GC1045                    4                     4.978               107.5
    ## 19 GC1045                    3                     5.054                97.5
    ## 20 GC1045                    0                     5.279                67.5
    ## 21 GC1081                   22                     2.930               275.0
    ## 22 GC1081                   17                     3.253               225.0
    ## 23 GC1081                   12                     3.589               175.0
    ## 24 GC1081                    7                     3.899               125.0
    ## 25 GC1081                    2                     4.207                75.0
    ## 26 GC1081                    1                     4.267                65.0
    ## 27 GC1081                    0                     4.325                55.0
    ## 28 GC1081                   -1                     4.379                45.0
    ## 29 GC1081                   -2                     4.430                35.0

Table S2. Increases in methane flux over the past two decades for cores
GC1045 and GC1081 and corresponding depths of modeled peak AOM rates.
Fluxes are integrated from AOM rate data, using cell widths of 2.5 cm.
(Peak AOM depths also at 2.5 cm resolution).