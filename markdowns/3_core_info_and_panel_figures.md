3\_core\_info\_and\_panel\_figures
================
Scott Klasek
6/9/2021

Figures 3, 4, 5, and S3 that show porewater data, microbial community
composition, and dsrAB & mcrA gene concentrations for cores sorted by
methane states.

## load necessary libraries

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.3     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.0     ✓ dplyr   1.0.5
    ## ✓ tidyr   1.1.3     ✓ stringr 1.4.0
    ## ✓ readr   1.4.0     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(phyloseq)
library(png)
library(grid)
library(egg)
```

    ## Loading required package: gridExtra

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
library(here)
```

    ## here() starts at /Users/scottklasek/Desktop/svalflux

``` r
sessioninfo <- sessionInfo()
```

## import and graph porewater data

``` r
# import the porewater data
porewater <- read.csv(file="/Users/scottklasek/Desktop/svalflux/data/porewater.csv")

# add in +/- 4% error bars on sulfate measurements
porewater$maxerror <- 0
porewater$minerror <- 0
porewater <- porewater %>% mutate(maxerror = ifelse(species=="sulfate (mM)", value*1.04, ifelse(species!="sulfate (mM)", value*1, maxerror)))
porewater <- porewater %>% mutate(minerror = ifelse(species=="sulfate (mM)", value*0.96, ifelse(species!="sulfate (mM)", value*1, minerror)))

# subset by core flowtype 
pore.1029 <- subset(porewater, porewater$stage=="seep") 
pore.fi <- subset(porewater, porewater$stage=="inc") 
pore.ss <- subset(porewater, porewater$stage=="ss")

# seep site
hydrate <- readPNG("/Users/scottklasek/Desktop/svalflux/figures/hydrates_symbol.png") # import hydrate symbol
hydrateg <- rasterGrob(hydrate, interpolate=TRUE) # raster-grob the hydrate symbol
worms <- readPNG("/Users/scottklasek/Desktop/svalflux/figures/pingo_tubeworms.png") # import tubeworm symbol
wormsg <- rasterGrob(worms, interpolate=TRUE) # raster-grob the tubeworm symbol
ggp.1029 <- ggplot(pore.1029,aes(depth,value,color=species))
gg.porewater.1029 <- ggp.1029+geom_line(size=1.5)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin = minerror, ymax = maxerror))+
  coord_flip()+
  ggtitle("Porewater Geochemistry")+
  scale_y_continuous("",limits = c(0,40),position = "right")+
  scale_x_reverse("Depth (cmbsf)",breaks=c(0,5,10,15,20,25,30,35),limits=c(35,0))+
  scale_color_manual("",values=c("orangered1","black","dodgerblue2","gold3"),guide=guide_legend())+
  annotation_custom(hydrateg, xmin = -37, xmax = -27, ymin = 0, ymax = 15)+
  annotation_custom(wormsg, xmin = -11, xmax = Inf, ymin = 13, ymax = 25)+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "vertical",axis.text=element_text(size=12),legend.text=element_text(size=11),axis.title=element_text(size=14))

# non-steady-state
ggp.fi <- ggplot(pore.fi,aes(depth,value,color=species))
gg.porewater.nss <- ggp.fi+geom_line(size=1.5)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin = minerror, ymax = maxerror))+
  geom_vline(aes(xintercept=smt),linetype="dashed",size=0.8)+
  coord_flip()+
  facet_grid(core~.)+
  scale_y_continuous("",limits = c(0,40),position = "right")+
  scale_x_reverse("Depth (cmbsf)",breaks=c(0,25,50,75,100,125,150),limits=c(150,0))+
  scale_color_manual("",values=c("orangered1","black","dodgerblue2","gold3"))+
  labs(title="Porewater\nGeochemistry")+
  theme_bw()+
  theme(axis.text=element_text(size=13),
        strip.text.y = element_text(size=13),
        title=element_text(size=13),
        axis.title=element_text(size=15),
        legend.text = element_text(size=13),
        legend.direction = "vertical",legend.position = "bottom",plot.margin=margin(0,0,0,0.25,"cm"))
# steady-state
biofilm <- readPNG("/Users/scottklasek/Desktop/svalflux/figures/biofilm_symbol.png") # import biofilm symbol
# to make sure this image is plotted only within facets corresponding to cores biofilms were found in, use this function found here: https://stackoverflow.com/questions/44688623/adding-custom-images-to-ggplot-facets
annotation_custom2 <- 
function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){ layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))} # then add this to the ggplot call, with data corresponding to the facets the image should appear in
ggp.ss <- ggplot(pore.ss,aes(depth,value,color=species))
gg.porewater.ss <- ggp.ss+geom_line(size=1.5)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin = minerror, ymax = maxerror))+
  geom_vline(aes(xintercept=smt),linetype="dashed",size=0.8)+
  coord_flip()+
  facet_grid(core~.)+
  scale_y_continuous("",limits = c(0,40),position = "right")+
  scale_x_reverse("Depth (cmbsf)",breaks=c(0,50,100,150,200,250,300,350),limits=c(350,0))+
  scale_color_manual("",values=c("orangered1","black","dodgerblue2","gold3"))+
  labs(title="Porewater\nGeochemistry")+
  annotation_custom2(rasterGrob(biofilm, interpolate=TRUE), xmin=-105, xmax=-20, ymin=0, ymax=40, data=pore.ss[1,])+ # for GC1070
  annotation_custom2(rasterGrob(biofilm, interpolate=TRUE), xmin=-350, xmax=-260, ymin=0, ymax=40, data=pore.ss[75,])+ # for GC1048
  theme_bw()+
  theme(axis.text=element_text(size=13),
        strip.text.y = element_text(size=13),
        title=element_text(size=13),
        axis.title=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        legend.direction = "vertical",legend.position = "bottom",plot.margin=margin(0,0,0,0.25,"cm"))
```

## import and graph aom rate data

``` r
aomrates <- read.csv(file="/Users/scottklasek/Desktop/svalflux/data/aomrates2020update.csv") # import rate data

# subset by methane flux condition
aomrates.inc <- subset(aomrates, aomrates$core=="GC1045"|aomrates$core=="GC1081")
aomrates.ss <- subset(aomrates, aomrates$core=="GC1068"|aomrates$core=="GC1069"|aomrates$core=="GC1070"|aomrates$core=="GC1048")
aomrates.seep <- subset(aomrates, aomrates$core=="PC1029")

# graph seep rates
aom1029 <- ggplot(aomrates.seep, aes(depth,aom))
raom1029 <- aom1029+geom_line(size=1.5)+
  coord_flip()+
  scale_y_continuous("",position="right", breaks=c(0,2000,4000,6000,8000))+
  scale_x_reverse("",breaks=c(0,5,10,15,20,25,30,35),limits=c(35,0))+
  labs(title=expression("AOM rate (nmol cm"^{-3}*" d"^{-1}*")"))+
  theme_bw()+
  theme(axis.text=element_text(size=12),legend.text=element_text(size=11),axis.title=element_text(size=14))

# graph non-steady-state rates 
aomrates.inc$yrsbp <- factor(aomrates.inc$yrsbp,levels=c(-2,-1,0,1,2,3,4,5,6,7,8,9,10)) 
aom.nss <- ggplot(aomrates.inc,aes(depth,aom,color=yrsbp))
raom.nss <- aom.nss+geom_line(size=1)+
  scale_color_manual("years \nbefore \nsampling",values=c("dodgerblue3","dodgerblue1","black","#662506","#993404","#cc4c02","#ec7014","#fe9929","#fec44f","grey75","grey65","grey55","grey45"))+
  geom_vline(aes(xintercept=smt),linetype="dashed",size=0.8)+
  coord_flip()+
  facet_grid(core~.)+
  labs(title=expression("Modeled \nAOM rate (nmol cm"^{-3}*" d"^{-1}*")"))+
  scale_y_continuous("",breaks=c(0,100,200),limits = c(-5,220),position = "right")+
  scale_x_reverse("",breaks=c(0,25,50,75,100,125,150),limits=c(150,0))+
  theme_bw()+
  theme(axis.text=element_text(size=13),
        strip.text.y = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        title=element_text(size=13),
        legend.title = element_text(size=13),
        legend.text = element_text(size=13),
        plot.margin=margin(1,0,0,-0.25,"cm"))

# graph steady-state rates
aom.ss <- ggplot(aomrates.ss, aes(depth,aom))
raom.ss <- aom.ss+geom_line(size=1)+
  geom_vline(aes(xintercept=smt),linetype="dashed",size=0.8)+
  coord_flip()+
  facet_grid(core~.)+
  scale_y_continuous("",limits = c(0,20),position = "right")+
  scale_x_reverse("",breaks=c(0,50,100,150,200,250,300,350),limits=c(350,0))+
  theme_bw()+
  labs(title=expression("AOM rate (nmol cm"^{-3}*" d"^{-1}*")"))+
  theme(axis.text=element_text(size=13),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        axis.text.y=element_blank(),
        title=element_text(size=13),
        legend.title = element_text(size=13),
        legend.text = element_text(size=13),
        plot.margin=margin(1,0,0,-0.25,"cm"))
```

## import bubble plots

``` r
bp.seep <- readRDS(file = "/Users/scottklasek/Desktop/svalflux/figures/bp.seep") 
bp.fluxincreasing <- readRDS(file = "/Users/scottklasek/Desktop/svalflux/figures/bp.fluxincreasing") 
bp.steadystate <- readRDS(file = "/Users/scottklasek/Desktop/svalflux/figures/bp.steadystate") 
```

## graph ddpcr data (already in phyloseq object)

``` r
# import phyloseq object
ps.frdp <- readRDS(file="/Users/scottklasek/Desktop/svalflux/data/ps.frdp") # imports the final phyloseq object

# make a dataframe with gene counts in long format
core <- rep(sample_data(ps.frdp)$core, times=2) # repeat core data for each gene
depth <- rep(sample_data(ps.frdp)$depth, times=2) # repeat depth 
smt <- rep(sample_data(ps.frdp)$smt, times=2) # repeat smt
genecount <- c(sample_data(ps.frdp)$dsrab, sample_data(ps.frdp)$mcra) # combine gene counts into one vector, mcra first
gene <- c(rep("dsrAB", each=76), rep("mcrA", each=76)) # create vector with type of gene being counted
gene <- factor(gene, levels = c("mcrA", "dsrAB"))
stage <- rep(sample_data(ps.frdp)$stage, times=2) # repeat stage
ddpcr.all <- data.frame(core, depth, smt, stage, genecount, gene) # make data frame
ddpcr.all$logcopies <- log10(ddpcr.all$genecount) # transform gene counts to log

# subset by stage of methane flux
dd.pcr.fi <- subset(ddpcr.all, ddpcr.all$stage=="fluxincreasing")
dd.pcr.ss <- subset(ddpcr.all, ddpcr.all$stage=="steadystate")
dd.pcr.seep <- subset(ddpcr.all, ddpcr.all$stage=="seep")

# seep ddPCR plot
dd1029 <- ggplot(dd.pcr.seep, aes(depth, logcopies, color=gene))
dd.seep <- dd1029+geom_line(size=1.5)+
  geom_point(size=3)+
  coord_flip()+
  scale_y_continuous("",breaks=c(4,5,6,7,8),limits = c(3.7,8.5),position="right")+
  scale_x_reverse("",breaks=c(0,5,10,15,20,25,30,35),limits=c(35,0))+
  scale_color_discrete("")+
  theme_bw()+
  labs(title=expression(paste("log"[10]," gene copies")))+
  theme(legend.position = "top",axis.text=element_text(size=12),legend.text=element_text(size=11),
        axis.title=element_text(size=8),axis.text.y = element_blank())

# non-steady-state ddPCR plot
ddfi<- ggplot(dd.pcr.fi, aes(depth, logcopies, color=gene))
dd.fi <- ddfi+geom_line(size=1.5)+
  geom_point(size=3)+
  coord_flip()+
  scale_y_continuous("", breaks=c(4,5,6,7), limits = c(3.7,7.4), position="right")+
  scale_x_reverse("", breaks=c(0,25,50,75,100,125,150), limits=c(150,0))+
  scale_color_discrete("")+
  geom_vline(aes(xintercept=smt),linetype="dashed",size=0.8)+
  facet_grid(core~.)+
  theme_bw()+
  labs(title=expression(paste("log"[10]," gene copies")))+
  theme(legend.position = "top",
        axis.title=element_text(size=13),
        axis.text = element_text(size=13),
        legend.text = element_text(size=13),
        strip.text = element_text(size=13))

# steady-state ddPCR plot
ddss<- ggplot(dd.pcr.ss, aes(depth, logcopies, color=gene))
dd.ss <- ddss+geom_line(size=1.5)+
  geom_point(size=3)+
  coord_flip()+
  scale_y_continuous("", breaks=c(4,5,6,7), limits = c(3.7,7.6), position="right")+
  scale_x_reverse("", breaks=c(0,50,100,150,200,250,300,350), limits=c(350,0))+
  geom_vline(aes(xintercept=smt), linetype="dashed",size=0.8)+
  facet_grid(core~.)+
  theme_bw()+
  labs(title=expression(paste("log"[10]," gene copies")))+
  theme(legend.position = "top",
        legend.text=element_text(size=13),
        strip.text = element_text(size=13),
        axis.title=element_text(size=13), 
        axis.text.x = element_text(size=13),
        axis.text.y = element_blank())
```

## create panel figures

``` r
# Figure 4 (Non-steady-state sites) export dimensions at 6 x 12
fig.nss <- ggarrange(gg.porewater.nss, bp.fluxincreasing, dd.fi, raom.nss, ncol = 4, nrow = 1,
  widths = c(2.4, 5.4, 1.8, 2.1),labels=c("A","B","C","D"))
```

    ## Warning: Removed 62 row(s) containing missing values (geom_path).

![](3_core_info_and_panel_figures_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# Figure 5 (Seep) export dimensions at 9.8 x 4.9
fig1029 <- ggarrange(gg.porewater.1029, bp.seep, dd.seep, widths=c(1.1,2.5,0.9),ncol = 3,nrow = 1,labels=c("A","B","C"))
```

    ## Warning: Removed 1 row(s) containing missing values (geom_path).

    ## Warning: Removed 3 rows containing missing values (geom_point).

![](3_core_info_and_panel_figures_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

Panel figure of steady-state without GC1048

``` r
# porewater plot
gg.porewater.ss.no48 <- ggplot(pore.ss %>% filter(core!="GC1048"), aes(depth,value,color=species))+
  geom_line(size=1.5)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin = minerror, ymax = maxerror))+
  geom_vline(aes(xintercept=smt),linetype="dashed",size=0.8)+
  coord_flip()+
  facet_grid(core~.)+
  scale_y_continuous("",limits = c(0,40),position = "right")+
  scale_x_reverse("Depth (cmbsf)",breaks=c(0,50,100,150,200,250,300,350),limits=c(350,0))+
  scale_color_manual("",values=c("orangered1","black","dodgerblue2","gold3"))+
  labs(title="Porewater\nGeochemistry")+
  annotation_custom2(rasterGrob(biofilm, interpolate=TRUE), xmin=-105, xmax=-20, ymin=0, ymax=40, data=pore.ss[1,])+ # for GC1070
  theme_bw()+
  theme(axis.text=element_text(size=13),
        strip.text.y = element_text(size=13),
        title=element_text(size=13),
        axis.title=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        legend.direction = "vertical",legend.position = "bottom",plot.margin=margin(0,0,0,0.25,"cm"))

# import steady-state bubble plot without GC1048
bp.steadystate.no48 <- readRDS(file = "/Users/scottklasek/Desktop/svalflux/figures/bp.steadystate.no48")

# ddPCR plot
dd.ss.no48 <- ggplot(dd.pcr.ss %>% filter(core!="GC1048"), aes(depth, logcopies, color=gene))+
  geom_line(size=1.5)+
  geom_point(size=3)+
  coord_flip()+
  scale_y_continuous("", breaks=c(4,5,6,7), limits = c(3.7,7.6), position="right")+
  scale_x_reverse("", breaks=c(0,50,100,150,200,250,300,350), limits=c(350,0))+
  geom_vline(aes(xintercept=smt), linetype="dashed",size=0.8)+
  facet_grid(core~.)+
  theme_bw()+
  labs(title=expression(paste("log"[10]," gene copies")))+
  theme(legend.position = "top",
        legend.text=element_text(size=13),
        strip.text = element_text(size=13),
        axis.title=element_text(size=13), 
        axis.text.x = element_text(size=13),
        axis.text.y = element_blank())

# Fig 3 panel figure
fig.ss.no48 <- ggarrange(gg.porewater.ss.no48, bp.steadystate.no48, dd.ss.no48,
                        ncol = 3, nrow = 1,widths = c(2.4,4.8,2),labels=c("A","B","C"))
```

    ## Warning: Removed 1 row(s) containing missing values (geom_path).

    ## Warning: Removed 4 rows containing missing values (geom_point).

![](3_core_info_and_panel_figures_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Panel figure of only GC1048 for supplementals

``` r
# porewater plot
gg.porewater.ss.only48 <- ggplot(pore.ss %>% filter(core=="GC1048"), aes(depth,value,color=species))+
  geom_line(size=1.5)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin = minerror, ymax = maxerror))+
  geom_vline(aes(xintercept=smt),linetype="dashed",size=0.8)+
  coord_flip()+
  facet_grid(core~.)+
  scale_y_continuous("",limits = c(0,40),position = "right")+
  scale_x_reverse("Depth (cmbsf)",breaks=c(0,50,100,150,200,250,300,350),limits=c(350,0))+
  scale_color_manual("",values=c("orangered1","black","dodgerblue2","gold3"))+
  labs(title="Porewater\nGeochemistry")+
  annotation_custom2(rasterGrob(biofilm, interpolate=TRUE), xmin=-350, xmax=-260, ymin=0, ymax=40, data=pore.ss[75,])+ # for GC1048
  theme_bw()+
  theme(axis.text=element_text(size=13),
        strip.text.y = element_text(size=13),
        title=element_text(size=13),
        axis.title=element_text(size=15),
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        legend.direction = "vertical",legend.position = "bottom",plot.margin=margin(0,0,0,0.25,"cm"))

# import steady-state bubble plot without GC1048
bp.steadystate.only48 <- readRDS(file = "/Users/scottklasek/Desktop/svalflux/figures/bp.steadystate.only48")

# ddPCR plot
dd.ss.only48 <- ggplot(dd.pcr.ss %>% filter(core=="GC1048"), aes(depth, logcopies, color=gene))+
  geom_line(size=1.5)+
  geom_point(size=3)+
  coord_flip()+
  scale_y_continuous("", breaks=c(4,5,6,7), limits = c(3.7,7.6), position="right")+
  scale_x_reverse("", breaks=c(0,50,100,150,200,250,300,350), limits=c(350,0))+
  geom_vline(aes(xintercept=smt), linetype="dashed",size=0.8)+
  facet_grid(core~.)+
  theme_bw()+
  labs(title=expression(paste("log"[10]," gene copies")))+
  theme(legend.position = "top",
        legend.text=element_text(size=13),
        strip.text = element_text(size=13),
        axis.title=element_text(size=13), 
        axis.text.x = element_text(size=13),
        axis.text.y = element_blank())

# make panel figure
fig.ss.only48 <- ggarrange(gg.porewater.ss.only48, bp.steadystate.only48, dd.ss.only48,
                        ncol = 3, nrow = 1,widths = c(2.4,4.8,2),labels=c("A","B","C"))
```

    ## Warning: Removed 1 rows containing missing values (geom_point).

![](3_core_info_and_panel_figures_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->
