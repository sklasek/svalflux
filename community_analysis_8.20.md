community\_analysis\_8.20
================
Scott Klasek
8/26/2020

## load necessary libraries and data

``` r
library(phyloseq)
library(tidyverse)
```

    ## ── Attaching packages ────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.2     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.4
    ## ✓ tidyr   1.0.2     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## Warning: package 'ggplot2' was built under R version 3.6.2

    ## ── Conflicts ───────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(egg)
```

    ## Loading required package: gridExtra

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
library(lemon)
```

    ## Warning: package 'lemon' was built under R version 3.6.2

    ## 
    ## Attaching package: 'lemon'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     %||%

    ## The following objects are masked from 'package:ggplot2':
    ## 
    ##     CoordCartesian, element_render

``` r
library(gridExtra)
library(grid)
library(patchwork)
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-6

``` r
library(DESeq2)
```

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     sampleNames

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## Loading required package: BiocParallel

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following object is masked from 'package:purrr':
    ## 
    ##     simplify

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply, rowsum

``` r
library(here)
```

    ## here() starts at /Users/scottklasek/Desktop/svalflux

``` r
sessioninfo <- sessionInfo()
ps.frdp <- readRDS(file="ps.frdp") # imports the final phyloseq object
```

## alpha diversity

``` r
# subset by high flux vs low flux 
ps.frdp.ss <- subset_samples(ps.frdp, sample_data(ps.frdp)$stage=="steadystate") # subsetsall steady-state samples
ps.frdp.nss.seep <- subset_samples(ps.frdp, sample_data(ps.frdp)$stage!="steadystate") # subsets the non-steady-state samples and the seep

alphadiv.ss <- estimate_richness(ps.frdp.ss, split = TRUE, measures=c("Shannon","InvSimpson","Observed")) # calculate tables of alpha diversity
```

    ## Warning in estimate_richness(ps.frdp.ss, split = TRUE, measures = c("Shannon", : The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

``` r
alphadiv.nss <- estimate_richness(ps.frdp.nss.seep, split = TRUE, measures=c("Shannon","InvSimpson","Observed"))
```

    ## Warning in estimate_richness(ps.frdp.nss.seep, split = TRUE, measures = c("Shannon", : The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

``` r
alphadiv.ss$depth <- substring(rownames(alphadiv.ss), 8,12) # pull out depth values from the sample names
alphadiv.ss$depth <- as.numeric(alphadiv.ss$depth) # convert them to numeric
alphadiv.nss$depth <- substring(rownames(alphadiv.nss), 8,12)
alphadiv.nss$depth <- as.numeric(alphadiv.nss$depth)

alphadiv.ss$core <- substring(rownames(alphadiv.ss), 1,6) # pull out core names from sample names
alphadiv.nss$core <- substring(rownames(alphadiv.nss), 1,6)

alphadiv.nss$peakaom <- NA # creates a NA value for peak AOM depth
alphadiv.ss$peakaom <- NA
alphadiv.nss[which(alphadiv.nss$core=="GC1081"),6] <- 55 # add in peak AOM depths for all cores
alphadiv.nss[which(alphadiv.nss$core=="GC1045"),6] <- 67.5
alphadiv.nss[which(alphadiv.nss$core=="PC1029"),6] <- 13
alphadiv.ss[which(alphadiv.ss$core=="GC1068"),6] <- 74
alphadiv.ss[which(alphadiv.ss$core=="GC1069"),6] <- 106
alphadiv.ss[which(alphadiv.ss$core=="GC1070"),6] <- 68
alphadiv.ss[which(alphadiv.ss$core=="GC1048"),6] <- 302.5
alphadiv.ss$dist_above_peakaom <- alphadiv.ss$peakaom-alphadiv.ss$depth # calculate distance above/below peak AOM
alphadiv.nss$dist_above_peakaom <- alphadiv.nss$peakaom-alphadiv.nss$depth 

alphadiv.ss.range <- subset(alphadiv.ss, dist_above_peakaom < 63.5 & dist_above_peakaom > -42.5) # subset the dataframe to include steady-state samples from depths across the SMT that correspond to depths of non-steady-state samples

nss.lm <- lm(alphadiv.nss$Shannon~alphadiv.nss$dist_above_peakaom) # linear model of shannon alpha diversity over distance to SMT for non-steady-state cores
ss.lm <- lm(alphadiv.ss$Shannon~alphadiv.ss$dist_above_peakaom) # linear model of shannon alpha diversity over distance to SMT for non-steady-state cores
ss.lm.range <- lm(alphadiv.ss.range$Shannon~alphadiv.ss.range$dist_above_peakaom)
summary(nss.lm) # summarizes the linear regression: slope is 0.029618, intercept is 4.167367, both are significant. Multiple R^2 is 0.5341, p value of slope is 1.5e-05
```

    ## 
    ## Call:
    ## lm(formula = alphadiv.nss$Shannon ~ alphadiv.nss$dist_above_peakaom)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.1531 -0.3537  0.2309  0.4064  0.9437 
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                     4.167367   0.115310  36.141  < 2e-16 ***
    ## alphadiv.nss$dist_above_peakaom 0.029618   0.005532   5.353  1.5e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5948 on 25 degrees of freedom
    ## Multiple R-squared:  0.5341, Adjusted R-squared:  0.5155 
    ## F-statistic: 28.66 on 1 and 25 DF,  p-value: 1.497e-05

``` r
summary(ss.lm) # slope is significant considering all the points across depth 
```

    ## 
    ## Call:
    ## lm(formula = alphadiv.ss$Shannon ~ alphadiv.ss$dist_above_peakaom)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1816 -0.5484  0.1525  0.4978  1.5101 
    ## 
    ## Coefficients:
    ##                                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                    4.036265   0.128180  31.489  < 2e-16 ***
    ## alphadiv.ss$dist_above_peakaom 0.003526   0.001070   3.294  0.00188 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.885 on 47 degrees of freedom
    ## Multiple R-squared:  0.1876, Adjusted R-squared:  0.1703 
    ## F-statistic: 10.85 on 1 and 47 DF,  p-value: 0.001881

``` r
summary(ss.lm.range) # slope is 0.002780, intercept is 3.831772, slope not significant across the SMT range. Multiple R^2 is 0.003901, p value of slope is 0.762
```

    ## 
    ## Call:
    ## lm(formula = alphadiv.ss.range$Shannon ~ alphadiv.ss.range$dist_above_peakaom)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -2.97898 -0.51723 -0.07824  0.57861  1.71682 
    ## 
    ## Coefficients:
    ##                                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                          3.831772   0.219516  17.456 3.83e-15 ***
    ## alphadiv.ss.range$dist_above_peakaom 0.002780   0.009067   0.307    0.762    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.05 on 24 degrees of freedom
    ## Multiple R-squared:  0.003901,   Adjusted R-squared:  -0.0376 
    ## F-statistic: 0.09399 on 1 and 24 DF,  p-value: 0.7618

``` r
# assign facet grid labels and methane flux stages
alphadiv.nss$core_flowtype_label <- "High Methane Flux"
alphadiv.ss$core_flowtype_label <- "Low Methane Flux"
alphadiv.nss$stage <- "increasing methane flux"
alphadiv.nss[which(alphadiv.nss$core=="PC1029"),9] <- "active methane seepage"
alphadiv.ss$stage <- "steady-state"

plotrichness.ss <- ggplot(alphadiv.ss, aes(dist_above_peakaom, Shannon, color=stage, shape=core)) 
ssplot <- plotrichness.ss+
  geom_point(size=2)+
  facet_grid(~core_flowtype_label)+
  scale_y_continuous("Shannon diversity index", breaks=c(1,2,3,4,5,6), limits = c(0.8,6.5))+
  scale_x_continuous("Distance above present-day peak AOM rate (cm)", limits = c(-250,306))+
  scale_color_manual("Methane stage", values=c("#8da0cb"))+
  scale_shape_manual("Core",values = c(4:7))+
  geom_vline(aes(xintercept=-42.5), linetype="dashed", size=0.8)+
  geom_vline(aes(xintercept=63.5), linetype="dashed", size=0.8)+
  geom_segment(x=-42.5,y=3.713622,xend=63.5,yend=4.008302, color="black", size=1)+
  theme_bw()+ # now for some reason you need to put the theme_bw() ahead of the theme() or nothing in theme() will work
  theme(strip.text = element_text(size = 11))+
  annotate("text", x = -150, y = 4.9, label = "paste(italic(R) ^ 2, \" = 0.0039\")", parse = TRUE)+
  annotate("text", x = -150, y = 5.6, label = "paste(italic(p), \" = 0.762\")", parse = TRUE)+
  annotate("rect", xmin = -225, xmax = -75, ymin = 4.5, ymax = 6, alpha = .2)

plotrichness.nss <- ggplot(alphadiv.nss, aes(dist_above_peakaom, Shannon, color=stage, shape=core))
nssplot <- plotrichness.nss+
  geom_point(size=2)+
  facet_grid(~core_flowtype_label)+
  scale_y_continuous("Shannon diversity index", breaks=c(1,2,3,4,5,6), limits = c(0.8,6.5))+
  scale_x_continuous("Distance above present-day peak AOM rate (cm)", limits = c(-250,306))+
  scale_color_manual("Methane stage", values=c("#fc8d62", "#66c2a5"))+
  scale_shape_manual("Core", values = c(1:3))+
  geom_vline(aes(xintercept=-42.5), linetype="dashed", size=0.8)+
  geom_vline(aes(xintercept=63.5), linetype="dashed", size=0.8)+
  geom_segment(x=-42.5,y=2.908602,xend=63.5,yend=6.04811, color="black", size=1)+
  theme_bw()+ # now for some reason you need to put the theme_bw() ahead of the theme() or nothing in theme() will work
  theme(strip.text = element_text(size = 11))+
  annotate("text", x = -150, y = 4.9, label = "paste(italic(R) ^ 2, \" = 0.5341\")", parse = TRUE)+
  annotate("text", x = -150, y = 5.6, label = "paste(italic(p), \" = 1.5e-05\")", parse = TRUE)+
  annotate("rect", xmin = -225, xmax = -75, ymin = 4.5, ymax = 6, alpha = .2)

div.fig <- ggarrange(nssplot,ssplot, heights = c(1,1), ncol=1, nrow=2, labels = c("A", "B")) # Figure 6
```

![](community_analysis_8.20_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
div.fig
```

![](community_analysis_8.20_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
div.fig <- saveRDS(div.fig, "figures/figure6") # export alpha diversity figure
```

## beta diversity

#### transforming ASV table

``` r
# Hellinger transformation
otu.hel <- otu_table(decostand(otu_table(ps.frdp), method = "hellinger"), taxa_are_rows=FALSE)
ps.hel <- phyloseq(tax_table(ps.frdp),
                    sample_data(ps.frdp),
                    otu_table(otu.hel),
                    phy_tree(ps.frdp),
                    refseq(ps.frdp)) 
```

#### permanova tests

Testing the influence of stage (seep, non-steady-state, and
steady-state) and geochemical zone (linear sulfate-reduction zone,
non-steady-state reduction zone, and
below-sulfate-methane-transition-zone) on community
composition

``` r
sample_data(ps.hel)$smtpos <- ifelse(sample_data(ps.hel)$geochem_zone == "linear SR zone","above SMTZ",
                              ifelse(sample_data(ps.hel)$geochem_zone == "nonlinear SR zone","above SMTZ",
                              ifelse(sample_data(ps.hel)$geochem_zone == "below SMTZ", "below SMTZ", NA))) # writing above-below smtz as another categorical variable
metadata.ps.hel <- as(sample_data(ps.hel), "data.frame") #  create sample data dataframe

# binary jaccard
dm.jac <- phyloseq::distance(ps.frdp, method = "jaccard", binary = TRUE) # makes binary Jaccard distance matrix (no need to transform)
adonis(dm.jac ~ geochem_zone*stage, data=metadata.ps.hel)
```

    ## 
    ## Call:
    ## adonis(formula = dm.jac ~ geochem_zone * stage, data = metadata.ps.hel) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                    Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## geochem_zone        2     1.899 0.94929  2.4343 0.05845  0.001 ***
    ## stage               2     2.485 1.24263  3.1865 0.07651  0.001 ***
    ## geochem_zone:stage  2     1.189 0.59469  1.5250 0.03662  0.001 ***
    ## Residuals          69    26.908 0.38997         0.82842           
    ## Total              75    32.481                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Bray-Curtis
dm.bc <- phyloseq::distance(ps.hel, method = "bray") # makes bray-curtis distance matrix
adonis(dm.bc ~ geochem_zone*stage, data=metadata.ps.hel)
```

    ## 
    ## Call:
    ## adonis(formula = dm.bc ~ geochem_zone * stage, data = metadata.ps.hel) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                    Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## geochem_zone        2    2.4565 1.22825  4.3089 0.09240  0.001 ***
    ## stage               2    3.3416 1.67078  5.8614 0.12570  0.001 ***
    ## geochem_zone:stage  2    1.1178 0.55891  1.9607 0.04205  0.001 ***
    ## Residuals          69   19.6684 0.28505         0.73985           
    ## Total              75   26.5843                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# unweighted Unifrac
set.seed(1)
dm.hel.unifrac <- UniFrac(ps.hel, weighted = FALSE, normalized = TRUE,  parallel = FALSE, fast = TRUE) # create weighted Unifrac distance matrix
```

    ## Warning in UniFrac(ps.hel, weighted = FALSE, normalized = TRUE, parallel =
    ## FALSE, : Randomly assigning root as -- ASV4966 -- in the phylogenetic tree in
    ## the data you provided.

``` r
adonis(dm.hel.unifrac ~ geochem_zone*stage, data=metadata.ps.hel) 
```

    ## 
    ## Call:
    ## adonis(formula = dm.hel.unifrac ~ geochem_zone * stage, data = metadata.ps.hel) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                    Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## geochem_zone        2    1.5197 0.75983  2.7836 0.06519  0.001 ***
    ## stage               2    2.0504 1.02522  3.7559 0.08796  0.001 ***
    ## geochem_zone:stage  2    0.9067 0.45337  1.6609 0.03890  0.002 ** 
    ## Residuals          69   18.8346 0.27297         0.80796           
    ## Total              75   23.3115                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# weighted Unifrac 
set.seed(1)
dm.hel.wunifrac <- UniFrac(ps.hel, weighted = TRUE, normalized = TRUE,  parallel = FALSE, fast = TRUE) # create weighted Unifrac distance matrix
```

    ## Warning in UniFrac(ps.hel, weighted = TRUE, normalized = TRUE, parallel =
    ## FALSE, : Randomly assigning root as -- ASV4966 -- in the phylogenetic tree in
    ## the data you provided.

``` r
adonis(dm.hel.wunifrac ~ geochem_zone*stage, data=metadata.ps.hel) # geochem zone is 9.6%
```

    ## 
    ## Call:
    ## adonis(formula = dm.hel.wunifrac ~ geochem_zone * stage, data = metadata.ps.hel) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                    Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
    ## geochem_zone        2    0.3328 0.166388  4.8776 0.09650  0.001 ***
    ## stage               2    0.5932 0.296596  8.6945 0.17202  0.001 ***
    ## geochem_zone:stage  2    0.1686 0.084296  2.4711 0.04889  0.002 ** 
    ## Residuals          69    2.3538 0.034113         0.68258           
    ## Total              75    3.4483                  1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# other variables are all significant, but possibly redundant:
adonis(dm.hel.wunifrac ~ stage, data=metadata.ps.hel) #  stage alone is 19.9%
```

    ## 
    ## Call:
    ## adonis(formula = dm.hel.wunifrac ~ stage, data = metadata.ps.hel) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## stage      2    0.6872 0.34359  9.0837 0.19928  0.001 ***
    ## Residuals 73    2.7612 0.03782         0.80072           
    ## Total     75    3.4483                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
adonis(dm.hel.wunifrac ~ core_flowtype, data=metadata.ps.hel) # significant, but at 13% is less than stage (which it is redundant with)
```

    ## 
    ## Call:
    ## adonis(formula = dm.hel.wunifrac ~ core_flowtype, data = metadata.ps.hel) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## core_flowtype  1    0.4786 0.47859  11.925 0.13879  0.001 ***
    ## Residuals     74    2.9698 0.04013         0.86121           
    ## Total         75    3.4483                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
adonis(dm.hel.wunifrac ~ SO4_mM, data=metadata.ps.hel) # SO4 more than geochem zone (maybe more helpful, too? but redundant)
```

    ## 
    ## Call:
    ## adonis(formula = dm.hel.wunifrac ~ SO4_mM, data = metadata.ps.hel) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## SO4_mM     1    0.3746 0.37464  9.0196 0.10864  0.001 ***
    ## Residuals 74    3.0737 0.04154         0.89136           
    ## Total     75    3.4483                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
adonis(dm.hel.wunifrac ~ HS_mM, data=metadata.ps.hel) # significant but redundant
```

    ## 
    ## Call:
    ## adonis(formula = dm.hel.wunifrac ~ HS_mM, data = metadata.ps.hel) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
    ## HS_mM      1    0.2088 0.208804  4.7696 0.06055  0.001 ***
    ## Residuals 74    3.2395 0.043778         0.93945           
    ## Total     75    3.4483                  1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
adonis(dm.hel.wunifrac ~ Alk_mM, data=metadata.ps.hel) # significant but redundant
```

    ## 
    ## Call:
    ## adonis(formula = dm.hel.wunifrac ~ Alk_mM, data = metadata.ps.hel) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
    ## Alk_mM     1    0.2735 0.273485  6.3744 0.07931  0.001 ***
    ## Residuals 74    3.1749 0.042904         0.92069           
    ## Total     75    3.4483                  1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
adonis(dm.hel.wunifrac ~ core, data=metadata.ps.hel) # core matters more than any other, but "core" itself is a meaningless variable
```

    ## 
    ## Call:
    ## adonis(formula = dm.hel.wunifrac ~ core, data = metadata.ps.hel) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs  MeanSqs F.Model     R2 Pr(>F)    
    ## core       6    1.1273 0.187880  5.5852 0.3269  0.001 ***
    ## Residuals 69    2.3211 0.033639         0.6731           
    ## Total     75    3.4483                  1.0000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
adonis(dm.hel.wunifrac ~ pingo, data=metadata.ps.hel) # gas hydrate mound the cores are sampled from matters more as well, but also a meaningless variable
```

    ## 
    ## Call:
    ## adonis(formula = dm.hel.wunifrac ~ pingo, data = metadata.ps.hel) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs  MeanSqs F.Model     R2 Pr(>F)    
    ## pingo      3    0.7893 0.263114  7.1245 0.2289  0.001 ***
    ## Residuals 72    2.6590 0.036931         0.7711           
    ## Total     75    3.4483                  1.0000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
adonis(dm.hel.wunifrac ~ smtzposition, data=metadata.ps.hel) # above/below SMT explains only 5%
```

    ## 
    ## Call:
    ## adonis(formula = dm.hel.wunifrac ~ smtzposition, data = metadata.ps.hel) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##              Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
    ## smtzposition  1    0.1840 0.183970  4.1704 0.05335  0.001 ***
    ## Residuals    74    3.2644 0.044113         0.94665           
    ## Total        75    3.4483                  1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
adonis(dm.hel.wunifrac ~ geochem_zone*core, data=metadata.ps.hel) # no significant interaction between core and geochem zone
```

    ## 
    ## Call:
    ## adonis(formula = dm.hel.wunifrac ~ geochem_zone * core, data = metadata.ps.hel) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                   Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
    ## geochem_zone       2    0.3328 0.166388  5.3888 0.09650  0.001 ***
    ## core               6    0.9888 0.164803  5.3375 0.28675  0.001 ***
    ## geochem_zone:core  5    0.2124 0.042484  1.3759 0.06160  0.071 .  
    ## Residuals         62    1.9143 0.030876         0.55515           
    ## Total             75    3.4483                  1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Using the weighted Unifrac distance, both geochem\_zone & stage were
significant (stage explaining 16-17% of variance, while zone 8-9%).
Interactions were also significant at around 4.9%. This did not vary by
the transformation type (I tried variance-stabilizing transformation as
well and got nearly the same results). Bray-Curtis, unweighted Unifrac,
and Jaccard distance metrics showed the same patterns, though they
respectively explained decreasing amounts of variance overall. So it
seems relative abundances are a bit more influential than taxonomy alone
in explaining differences between these communities, though they both
help.

Pingo (which Gas Hydrate Mound a sample was collected from) explains a
bit more variance than does stage, and Core does even more. But it is a
confounding variable (each pingo is at a separate methane flux stage),
and which particular pingo a sample comes from offers no real helpful
information. Core also confounds with pingo, core\_flowtype, stage, etc.
Even though core was the single variable that explained more of the
variance than any other, geochem\_zone\*stage explained slightly more.
Furthermore, core and geochem\_zone do not have a significant
interaction.

#### ordinations

weighted Unifrac with Hellinger transformation

``` r
# make better labels for plotting
sample_data(ps.hel)$stage <- ifelse(sample_data(ps.hel)$stage == "seep", "active methane seepage",
                        ifelse(sample_data(ps.hel)$stage == "fluxincreasing", "increasing methane flux",
                        ifelse(sample_data(ps.hel)$stage == "steadystate", "steady-state", NA)))

sample_data(ps.hel)$geochem_zone <- ifelse(sample_data(ps.hel)$geochem_zone == "lin", "linear SR zone",
                        ifelse(sample_data(ps.hel)$geochem_zone == "nss", "nonlinear SR zone",
                        ifelse(sample_data(ps.hel)$geochem_zone == "below", "below SMTZ", NA)))

sample_data(ps.hel)$geochem_zone <- factor(sample_data(ps.hel)$geochem_zone, levels = c("linear SR zone", "nonlinear SR zone", "below SMTZ")) # reorder 

# trying NMDS and PCoAs
set.seed(1)
ord.ps.hel.wuni.nmds <- ordinate(ps.hel, "NMDS", "unifrac", weighted=TRUE) # Stress = 0.135
```

    ## Warning in UniFrac(physeq, ...): Randomly assigning root as -- ASV4966 -- in the
    ## phylogenetic tree in the data you provided.

    ## Run 0 stress 0.1306678 
    ## Run 1 stress 0.1246761 
    ## ... New best solution
    ## ... Procrustes: rmse 0.02366616  max resid 0.1580323 
    ## Run 2 stress 0.134905 
    ## Run 3 stress 0.1246751 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001182949  max resid 0.0007144867 
    ## ... Similar to previous best
    ## Run 4 stress 0.1306665 
    ## Run 5 stress 0.1413206 
    ## Run 6 stress 0.1246746 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0008163292  max resid 0.005657734 
    ## ... Similar to previous best
    ## Run 7 stress 0.1349034 
    ## Run 8 stress 0.1365192 
    ## Run 9 stress 0.1306541 
    ## Run 10 stress 0.1246775 
    ## ... Procrustes: rmse 0.0003118379  max resid 0.002190585 
    ## ... Similar to previous best
    ## Run 11 stress 0.1585576 
    ## Run 12 stress 0.1246767 
    ## ... Procrustes: rmse 0.0009869864  max resid 0.006751837 
    ## ... Similar to previous best
    ## Run 13 stress 0.134899 
    ## Run 14 stress 0.124674 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001310897  max resid 0.0008783968 
    ## ... Similar to previous best
    ## Run 15 stress 0.134901 
    ## Run 16 stress 0.1702664 
    ## Run 17 stress 0.1306596 
    ## Run 18 stress 0.1349008 
    ## Run 19 stress 0.1389394 
    ## Run 20 stress 0.1349014 
    ## *** Solution reached

``` r
set.seed(1)
ord.ps.hel.wuni.pcoa <- ordinate(ps.hel, "PCoA", "unifrac", weighted=TRUE)
```

    ## Warning in UniFrac(physeq, ...): Randomly assigning root as -- ASV4966 -- in the
    ## phylogenetic tree in the data you provided.

``` r
# plotting ordinations by stage 
lb1 <- paste("R^2 == 0.197")
ord1.p.stage <- plot_ordination(ps.hel, ord.ps.hel.wuni.pcoa, color = "stage")+
  stat_ellipse()+
  scale_color_manual("",values = c("#fc8d62", "#66c2a5", "#8da0cb"))+
  ggtitle("A", subtitle = "Cores classified by methane stage")+
  annotate(geom="text", x=-0.1, y=0.28, label="p < 0.001") +
  annotate(geom="text", x=-0.1, y=0.34, label=lb1, parse=TRUE)+
  theme_bw()+
  theme(legend.position = "bottom")+
  guides(col = guide_legend(ncol = 1))

# plotting ordinations by redox zone 
lb2 <- paste("R^2 == 0.096")
ord1.p.zone<- plot_ordination(ps.hel, ord.ps.hel.wuni.pcoa, color = "geochem_zone")+
  stat_ellipse()+
  scale_color_manual("",values = c("#e78ac3","#a6d854","#ffd92f"))+
  ggtitle("B", subtitle = "Samples classified by redox zone")+
  annotate(geom="text", x=-0.15, y=0.17, label="p < 0.001", parse=TRUE)+
  annotate(geom="text", x=-0.15, y=0.23, label=lb2, parse=TRUE)+
  theme_bw()+
  theme(legend.position = "bottom")+
  guides(col = guide_legend(ncol = 1))

# showing them both together
ord.fig <- ord1.p.stage + ord1.p.zone
ord.fig
```

![](community_analysis_8.20_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ord1.p.stage <- saveRDS(ord1.p.stage, "figures/ord1.p.stage")
ord1.p.zone <- saveRDS(ord1.p.zone, "figures/ord1.p.zone")
```

The PCoA has fewer outliers, and explains 25% x 17.4% variance on the
two major axes, so I’ll stick with that over NMDS.

Maybe don’t even show the plot by geochem zone… it’s just a mess that
doesn’t suggest any clear grouping pattern. Just showing it because it
tested as significant doesn’t seem a good enough reason…

Can also add the sections “more permanovas” from lines 844-920, and
possibly beta-diversity from lines 938-974. Just to be thorough…

## follow-up tests on beta-diversity

``` r
library(phyloseq)
# omitting below-SMT samples, is there still a difference by stage?
nsamples(ps.hel)
ps.hel.above <- subset_samples(ps.hel, geochem_zone!="below SMTZ") # subset
nsamples(ps.hel.above)
set.seed(1)
dm.hel.above.wunifrac <- UniFrac(ps.hel.above, weighted = TRUE, normalized = TRUE,  parallel = FALSE, fast = TRUE) # run unifrac
metadata.ps.hel.above <- as(sample_data(ps.hel.above), "data.frame") # write metadata data frame
adonis(dm.hel.above.wunifrac ~ stage, data=metadata.ps.hel.above) # yes there is a difference: R2, p are 0.21521,  0.001 ***

# considering only below-SMT samples, is there still a difference by stage?
ps.hel.below <- subset_samples(ps.hel, geochem_zone=="below SMTZ") # subset
nsamples(ps.hel.below)
set.seed(1)
dm.hel.below.wunifrac <- UniFrac(ps.hel.below, weighted = TRUE, normalized = TRUE,  parallel = FALSE, fast = TRUE) # run unifrac
metadata.ps.hel.below <- as(sample_data(ps.hel.below), "data.frame") # write metadata data frame
adonis(dm.hel.below.wunifrac ~ stage, data=metadata.ps.hel.below) # yes there is a difference: R2, p are 0.18597  0.003 **

# differences between the two stages steady-state and fluxincreasing? 
ps.hel.gc <- subset_samples(ps.hel, stage!="active methane seepage") # subset
nsamples(ps.hel.gc)
set.seed(1)
dm.hel.gc.wunifrac <- UniFrac(ps.hel.gc, weighted = TRUE, normalized = TRUE,  parallel = FALSE, fast = TRUE) # run unifrac
metadata.ps.hel.gc <- as(sample_data(ps.hel.gc), "data.frame") # write metadata data frame
adonis(dm.hel.gc.wunifrac ~ stage, data=metadata.ps.hel.gc) # there is a difference: R2, p are 0.06055  0.001 ***

# differences between the two stages (steady-state and fluxincreasing) when only considering above-SMT samples? 
ps.hel.above.gc <- subset_samples(ps.hel.above, stage!="active methane seepage") # subset
nsamples(ps.hel.above.gc)
set.seed(11)
dm.hel.above.gc.wunifrac <- UniFrac(ps.hel.above.gc, weighted = TRUE, normalized = TRUE,  parallel = FALSE, fast = TRUE) # run unifrac
metadata.ps.hel.above.gc <- as(sample_data(ps.hel.above.gc), "data.frame") # write metadata data frame
adonis(dm.hel.above.gc.wunifrac ~ stage, data=metadata.ps.hel.above.gc) #  there is barely a difference: R2, p are 0.04448  0.056 (Honestly, if you really like p < 0.05, significance depends on which set.seed you use)

# differences between linear SR and non-steady-state SR across both seep and fluxincreasing?
ps.hel.above.flux <- subset_samples(ps.hel, stage!="steady-state" & geochem_zone!="below SMTZ") # subset
nsamples(ps.hel.above.flux)
set.seed(1)
dm.hel.above.flux.wunifrac <- UniFrac(ps.hel.above.flux, weighted = TRUE, normalized = TRUE,  parallel = FALSE, fast = TRUE) # run unifrac
metadata.ps.hel.above.flux <- as(sample_data(ps.hel.above.flux), "data.frame") # write metadata data frame
adonis(dm.hel.above.flux.wunifrac ~ geochem_zone, data=metadata.ps.hel.above.flux) #  there is a difference: R2, p are 0.16792  0.004 **

# differences between linear SR and non-steady-state SR across fluxincreasing only?
ps.hel.above.fluxinc <- subset_samples(ps.hel.above.flux, stage=="increasing methane flux") # subset
nsamples(ps.hel.above.fluxinc) # there are only 8 samples
set.seed(1)
dm.hel.above.fluxinc.wunifrac <- UniFrac(ps.hel.above.fluxinc, weighted = TRUE, normalized = TRUE,  parallel = FALSE, fast = TRUE) # run unifrac
metadata.ps.hel.above.fluxinc <- as(sample_data(ps.hel.above.fluxinc), "data.frame") # write metadata data frame
adonis(dm.hel.above.fluxinc.wunifrac ~ geochem_zone, data=metadata.ps.hel.above.fluxinc) #  there is barely a difference: R2, p are 0.3008  0.052 (Honestly, if you really like p < 0.05, significance might depend on which set.seed you use)

# differences between linear SR and non-steady-state SR across seep only?
ps.hel.above.fluxseep <- subset_samples(ps.hel.above.flux, stage=="active methane seepage") # subset
nsamples(ps.hel.above.fluxseep)
set.seed(1)
dm.hel.above.fluxseep.wunifrac <- UniFrac(ps.hel.above.fluxseep, weighted = TRUE, normalized = TRUE,  parallel = FALSE, fast = TRUE) # run unifrac
metadata.ps.hel.above.fluxseep <- as(sample_data(ps.hel.above.fluxseep), "data.frame") # write metadata data frame
adonis(dm.hel.above.fluxseep.wunifrac ~ geochem_zone, data=metadata.ps.hel.above.fluxseep) #  there is a difference: R2, p are 0.34069  0.002 **

# differences between non-steady-state SR and below-SMT across fluxincreasing only?
ps.hel.fluxinc.lower <- subset_samples(ps.hel, stage=="increasing methane flux" & geochem_zone!="lin") # subset
nsamples(ps.hel.fluxinc.lower)# only 12 samples
set.seed(1)
dm.hel.fluxinc.lower.wunifrac <- UniFrac(ps.hel.fluxinc.lower, weighted = TRUE, normalized = TRUE,  parallel = FALSE, fast = TRUE) # run unifrac
metadata.ps.hel.fluxinc.lower <- as(sample_data(ps.hel.fluxinc.lower), "data.frame") # write metadata data frame
adonis(dm.hel.fluxinc.lower.wunifrac ~ geochem_zone, data=metadata.ps.hel.fluxinc.lower) #  there is a difference: R2, p are 0.36624  0.026 * 

# following this up with an ordination of fluxincreasing GC samples and color by geochem_zone
ps.hel.fluxinc <- subset_samples(ps.hel, stage=="increasing methane flux") # subset
nsamples(ps.hel.fluxinc)
ord.ps.helfluxinc.wuni.pcoa <- ordinate(ps.hel.fluxinc, "PCoA", "unifrac", weighted=TRUE) # ordinate
plot_ordination(ps.hel.fluxinc, ord.ps.helfluxinc.wuni.pcoa, color = "geochem_zone")+stat_ellipse() # from the very few samples we have here, the idea that the linear SR zone is the outlier is supported

# t-tests of beta-diversity between stages
ps.hel.seep <- subset_samples(ps.hel, stage=="active methane seepage") # subset
ps.hel.fluxinc
ps.hel.ss <- subset_samples(ps.hel, stage=="steady-state")

set.seed(1)
dm.seep <- UniFrac(ps.hel.seep, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE) # calculate distance matrices for beta-diversity calculations
set.seed(1)
dm.fluxinc <- UniFrac(ps.hel.fluxinc, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
set.seed(1)
dm.ss <- UniFrac(ps.hel.ss, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)

bdiv.seep <- as.vector(dm.seep) # rewrite the triangle distance matrices as vectors
bdiv.fluxinc <- as.vector(dm.fluxinc) 
bdiv.ss <- as.vector(dm.ss) 

t.test(bdiv.seep, bdiv.ss) 
t.test(bdiv.fluxinc, bdiv.ss) 
t.test(bdiv.fluxinc, bdiv.seep) 
```

Conclusions: all stages are distinctly different from one another. This
is true even when considering only above-SMTZ or below-SMTZ samples,
though remember that seep stages did not have below-SMTZ samples. The
weakest difference (just not quite statistically significant) was
between above-SMTZ samples when comparing fluxincreasing and
steady-state directly.

There appear to be differences between steady-state and non-steady-state
regions of sulfate-reduction zones, though these difference only tested
significant in PC1029 samples. In fluxincreasing cores, samples from the
nonlinear SR zone clustered closer with below SMTZ samples than those
from linear SR zones, suggesting they have begun to change towards
below-SMTZ communities. However, there are only 12 samples here and we’d
need more to make this conclusion for sure.

There is also evidence for lower beta-diversity among seep samples than
among the other two stages, suggesting there is some aspect of community
convergence (not that surprising). But remember, this stage is
represented by only one core, and we saw a lot of inter-core
variance.

## differentially abundant ASVs across methane flux stages and above/below the SMTZ

We have 76 samples: 15 from seep, 49 from steady-state, and 12 from
flux-increasing stages. There are no seep samples from below-SMTZ, 16
from steady-state, and 4 from flux-increasing. In an ideal world, we
would have a more balanced sample design for assessing ASVs that are
differentially abundant between these three stages (above the SMTZ) and
between the latter two stages (below SMTZ). However, DESeq2 is able to
handle imbalances in group numbers fairly robustly:
<https://support.bioconductor.org/p/115433/>

``` r
# enter logical metadata categories corresponding to whether sample stage is seep or increasing flux or steady-state
sample_data(ps.frdp)$is.seep <- FALSE
sample_data(ps.frdp)[which(sample_data(ps.frdp)$stage=="seep"),22] <- TRUE
sample_data(ps.frdp)$is.flux <- FALSE
sample_data(ps.frdp)[which(sample_data(ps.frdp)$stage=="fluxincreasing"),23] <- TRUE
sample_data(ps.frdp)$is.ss <- FALSE
sample_data(ps.frdp)[which(sample_data(ps.frdp)$stage=="steadystate"),24] <- TRUE

# subset the phyloseq object by above/below SMTZ
ps.above <- subset_samples(ps.frdp, smtzposition=="above") 
ps.below <- subset_samples(ps.frdp, smtzposition=="below") 

# make the same phyloseq object, but otu table with relative abundance instead of raw counts
ps.above.ra <- transform_sample_counts(ps.above, function(OTU) OTU/sum(OTU))
ps.below.ra <- transform_sample_counts(ps.below, function(OTU) OTU/sum(OTU))

# define a function which will make a table out of deseq biomarker results, and return taxonomy and relative abundance info 
table.stage.from.deseq <- function(deseq, relabundps, group){
  de.results <- DESeq(deseq, test = "Wald", fitType = "parametric", sfType = "poscounts")
  de.results.table <- results(de.results, cooksCutoff = FALSE)
  de.results.table <- de.results.table[which(de.results.table$padj < 0.05), ]
  meanabund.vector <- vector("numeric", length(rownames(de.results.table)))
  for (i in rownames(de.results.table)) {meanabund.vector[[i]] <- mean(otu_table(relabundps)[which(sample_data(relabundps)$stage==group),i])}
  meanabund.vector <- meanabund.vector[(0.5*length(meanabund.vector)+1):length(meanabund.vector)] 
  de.results.table <- cbind(as(meanabund.vector, "matrix"), as(de.results.table, "data.frame"), as(tax_table(relabundps)[rownames(de.results.table), ], "matrix"))
  colnames(de.results.table)[1] <- "relabund"
  return(de.results.table)
}

# then convert phyloseq object to deseq, and run the function for each stage of methane dynamics
de.above.is.seep <- phyloseq_to_deseq2(ps.above, ~ is.seep)
```

    ## converting counts to integer mode

``` r
de.results.table.above.seep <- table.stage.from.deseq(deseq = de.above.is.seep, relabundps = ps.above.ra, group = "seep") 
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 5102 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
de.above.is.flux <- phyloseq_to_deseq2(ps.above, ~ is.flux)
```

    ## converting counts to integer mode

``` r
de.results.table.above.flux <- table.stage.from.deseq(deseq = de.above.is.flux, relabundps = ps.above.ra, group = "fluxincreasing") 
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 5301 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
de.above.is.ss <- phyloseq_to_deseq2(ps.above, ~ is.ss)
```

    ## converting counts to integer mode

``` r
de.results.table.above.ss <- table.stage.from.deseq(deseq = de.above.is.ss, relabundps = ps.above.ra, group = "steadystate") 
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 5536 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
de.below.is.flux <- phyloseq_to_deseq2(ps.below, ~ is.flux)
```

    ## converting counts to integer mode

``` r
de.results.table.below.flux <- table.stage.from.deseq(deseq = de.below.is.flux, relabundps = ps.below.ra, group = "fluxincreasing") 
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 1312 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
de.below.is.ss <- phyloseq_to_deseq2(ps.below, ~ is.ss)
```

    ## converting counts to integer mode

``` r
de.results.table.below.ss <- table.stage.from.deseq(deseq = de.below.is.ss, relabundps = ps.below.ra, group = "steadystate") 
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 1312 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
# assign category for each de.results.table
de.results.table.above.seep$de.group <- "above.seep"
de.results.table.above.flux$de.group <- "above.flux"
de.results.table.above.ss$de.group <- "above.ss"
de.results.table.below.flux$de.group <- "below.flux"
de.results.table.below.ss$de.group <- "below.ss"

# rewrite ASV numbers so rbind will not rename duplicates
de.results.table.above.seep$asv <- rownames(de.results.table.above.seep)
de.results.table.above.flux$asv <- rownames(de.results.table.above.flux)
de.results.table.above.ss$asv <- rownames(de.results.table.above.ss)
de.results.table.below.flux$asv <- rownames(de.results.table.below.flux)
de.results.table.below.ss$asv <- rownames(de.results.table.below.ss)

# merge the biomarkers from each block into a single data frame
all.biomarkers <- rbind(as(de.results.table.above.seep, "data.frame"), 
                        as(de.results.table.above.flux, "data.frame"), 
                        as(de.results.table.above.ss, "data.frame"),
                        as(de.results.table.below.flux, "data.frame"),
                        as(de.results.table.below.ss, "data.frame")) 

# select only the ones more abundant in the blocks (as opposed to less)
all.biomarkers <- subset(all.biomarkers, log2FoldChange>0) 

# add in metadata about SMT position only, and stage only
all.biomarkers$smtpos <- "above"
all.biomarkers[which(all.biomarkers$de.group=="below.flux"),16] <- "below"
all.biomarkers[which(all.biomarkers$de.group=="below.ss"),16] <- "below"
all.biomarkers$stage <- "increasing methane flux"
all.biomarkers[which(all.biomarkers$de.group=="above.seep"),17] <- "active methane seepage"
all.biomarkers[which(all.biomarkers$de.group=="above.ss"),17] <- "steady-state"
all.biomarkers[which(all.biomarkers$de.group=="below.ss"),17] <- "steady-state"

# subset all.biomarkers dataframe by categories we want to graph at
most.biomarkers <- all.biomarkers %>% filter(Class!="Deltaproteobacteria" & Class!="Methanomicrobia") # remove Deltaproteobacteria & ANMEs because we'll graph them at higher taxonomic resolution, also removes NAs
anme.biomarkers <- all.biomarkers %>% filter(Class == "Methanomicrobia") # subset ANMEs only (all Methanomicrobial biomarkers are ANME in this dataset)
delta.biomarkers <- all.biomarkers %>% filter(Class == "Deltaproteobacteria" & is.na(Genus)==FALSE) # subset Deltaproteobacteria only

# when plotting, facet by more than two stages: Plotting most biomarkers at the class level, but ANMEs at Family and Deltas at Genus
most.biomarkers$plotlevel <- "Community"
anme.biomarkers$plotlevel <- "ANME"
delta.biomarkers$plotlevel <- "SRB"
most.biomarkers$taxlevel <- most.biomarkers$Class
anme.biomarkers$taxlevel <- anme.biomarkers$Family
delta.biomarkers$taxlevel <- delta.biomarkers$Genus
biomarkers.to.plot <- rbind(most.biomarkers, anme.biomarkers, delta.biomarkers) # combine all subsetted dataframes for plotting

gg.all.biom <- ggplot(biomarkers.to.plot, aes(x=factor(taxlevel, levels = rev(levels(factor(taxlevel)))), y=log2FoldChange, color=Kingdom)) + 
  geom_point(aes(size = relabund*100, fill=Kingdom), color="black", pch=21) +
  scale_size_area("Percent abundance", max_size = max(all.biomarkers$relabund)*400)+
  scale_fill_discrete("Domain")+
  scale_x_discrete("")+
  scale_y_continuous("log2 Differential Abundance Change", position = "right", limits = c(0,27)) +
  facet_grid(smtpos+plotlevel~stage, scales = "free", space = "free") +
  coord_flip() +
  theme_bw()+
  theme(strip.text.x = element_text(size=11), strip.text.y = element_text(size=8),
        axis.text.x = element_text(size=10))
gg.all.biom
```

![](community_analysis_8.20_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
gg.all.biom <- saveRDS(gg.all.biom, "figures/figure5") # export figure
```
