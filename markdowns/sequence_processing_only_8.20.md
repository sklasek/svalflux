sequence\_processing\_only\_8.20
================
Scott Klasek
9/3/2020

## 16S sequence processing and decontamination

Analyses are in the documents:  
bubble\_plots (making community composition bubble plots)  
core\_info\_and\_panel\_figures (adding in geochemistry, modeled AOM
rates, and gene count data along with community composition)  
community\_analysis (alpha and beta diversity, ordinations, biomarkers
between different methane flux
    stages)

## load necessary libraries

``` r
library(dada2)
```

    ## Loading required package: Rcpp

``` r
library(tidyverse)
```

    ## ── Attaching packages ───────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.2     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.4
    ## ✓ tidyr   1.0.2     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## Warning: package 'ggplot2' was built under R version 3.6.2

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(phyloseq)
library(decontam)
library(DECIPHER)
```

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

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

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

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

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'XVector'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     compact

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: RSQLite

``` r
library(phangorn)
```

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     complement

``` r
library(here)
```

    ## here() starts at /Users/scottklasek/Desktop/svalflux

``` r
sessioninfo <- sessionInfo()
```

## run dada2 workflow to get ASVs

Estimated time to run this chunk on my laptop computer (2019 macbook
pro): an hour or two. Taxonomy was the longest single step for me.

``` r
fqpath <- "/Users/scottklasek/Desktop/flux_resubmission/fastq"
original_file_names <- list.files(fqpath)
new_file_names <- list.files(fqpath)
fnFs.2 <- sort(list.files(fqpath, pattern="_R1.fastq", full.names = TRUE))
fnRs.2 <- sort(list.files(fqpath, pattern="_R2.fastq", full.names = TRUE))
sample.names.2 <- sapply(strsplit(basename(fnFs.2), "_"), `[`, 1)

# you can plot profiles of the read quality like this, but I'm omitting it in the markdown
# plotQualityProfile(fnFs.2[30:35])
# plotQualityProfile(fnRs.2[30:35])

# problem: some of our samples are sequenced twice. which fastq files are better quality?
dup_samples <- c("GC1048-305.0A","GC1048-305.0B","GC1070-063.0A","GC1070-063.0B","GC1070-068.0A","GC1070-068.0B")
dup_vectors <- match(dup_samples,sample.names.2)
# plotQualityProfile(fnFs.2[dup_vectors])
# plotQualityProfile(fnRs.2[dup_vectors])
filtFs.2 <- file.path(fqpath, "filtered", paste0(sample.names.2, "_F_filt.fastq"))
filtRs.2 <- file.path(fqpath, "filtered", paste0(sample.names.2, "_R_filt.fastq"))
names(filtFs.2) <- sample.names.2
names(filtRs.2) <- sample.names.2
out.2 <- filterAndTrim(fnFs.2, filtFs.2, fnRs.2, filtRs.2, truncLen=c(230,230), trimLeft=c(19,20),
                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=TRUE)
errF.2 <- learnErrors(filtFs.2, multithread=TRUE)
```

    ## 106606273 total bases in 505243 reads from 13 samples will be used for learning the error rates.

``` r
errR.2 <- learnErrors(filtRs.2, multithread=TRUE)
```

    ## 106101030 total bases in 505243 reads from 13 samples will be used for learning the error rates.

``` r
plotErrors(errF.2, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](sequence_processing_only_8.20_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
plotErrors(errR.2, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](sequence_processing_only_8.20_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
dadaFs.2 <- dada(filtFs.2, err=errF.2, multithread=TRUE)
```

    ## Sample 1 - 10185 reads in 2417 unique sequences.
    ## Sample 2 - 552 reads in 337 unique sequences.
    ## Sample 3 - 62504 reads in 31923 unique sequences.
    ## Sample 4 - 32003 reads in 15973 unique sequences.
    ## Sample 5 - 32920 reads in 14814 unique sequences.
    ## Sample 6 - 41408 reads in 13437 unique sequences.
    ## Sample 7 - 39234 reads in 9062 unique sequences.
    ## Sample 8 - 37067 reads in 12295 unique sequences.
    ## Sample 9 - 29608 reads in 10283 unique sequences.
    ## Sample 10 - 77475 reads in 36288 unique sequences.
    ## Sample 11 - 47942 reads in 23240 unique sequences.
    ## Sample 12 - 43258 reads in 16809 unique sequences.
    ## Sample 13 - 51087 reads in 21121 unique sequences.
    ## Sample 14 - 65629 reads in 16672 unique sequences.
    ## Sample 15 - 65909 reads in 15472 unique sequences.
    ## Sample 16 - 112 reads in 108 unique sequences.
    ## Sample 17 - 46570 reads in 13283 unique sequences.
    ## Sample 18 - 51853 reads in 12452 unique sequences.
    ## Sample 19 - 35860 reads in 9356 unique sequences.
    ## Sample 20 - 53928 reads in 17665 unique sequences.
    ## Sample 21 - 26824 reads in 7721 unique sequences.
    ## Sample 22 - 32912 reads in 7106 unique sequences.
    ## Sample 23 - 49762 reads in 5250 unique sequences.
    ## Sample 24 - 3490 reads in 368 unique sequences.
    ## Sample 25 - 1150 reads in 335 unique sequences.
    ## Sample 26 - 177 reads in 165 unique sequences.
    ## Sample 27 - 148 reads in 134 unique sequences.
    ## Sample 28 - 24995 reads in 7860 unique sequences.
    ## Sample 29 - 37086 reads in 12530 unique sequences.
    ## Sample 30 - 33520 reads in 10261 unique sequences.
    ## Sample 31 - 24261 reads in 7133 unique sequences.
    ## Sample 32 - 33107 reads in 10418 unique sequences.
    ## Sample 33 - 42225 reads in 13360 unique sequences.
    ## Sample 34 - 40202 reads in 11432 unique sequences.
    ## Sample 35 - 215 reads in 199 unique sequences.
    ## Sample 36 - 1578 reads in 554 unique sequences.
    ## Sample 37 - 58540 reads in 16009 unique sequences.
    ## Sample 38 - 43474 reads in 9641 unique sequences.
    ## Sample 39 - 11053 reads in 3158 unique sequences.
    ## Sample 40 - 4310 reads in 1006 unique sequences.
    ## Sample 41 - 55935 reads in 23144 unique sequences.
    ## Sample 42 - 38606 reads in 14095 unique sequences.
    ## Sample 43 - 31065 reads in 8692 unique sequences.
    ## Sample 44 - 99396 reads in 32454 unique sequences.
    ## Sample 45 - 118309 reads in 36428 unique sequences.
    ## Sample 46 - 111599 reads in 32739 unique sequences.
    ## Sample 47 - 89632 reads in 26877 unique sequences.
    ## Sample 48 - 35287 reads in 10555 unique sequences.
    ## Sample 49 - 27880 reads in 8379 unique sequences.
    ## Sample 50 - 39234 reads in 11148 unique sequences.
    ## Sample 51 - 77005 reads in 15004 unique sequences.
    ## Sample 52 - 60498 reads in 28543 unique sequences.
    ## Sample 53 - 30306 reads in 10582 unique sequences.
    ## Sample 54 - 52678 reads in 16062 unique sequences.
    ## Sample 55 - 29046 reads in 7174 unique sequences.
    ## Sample 56 - 35118 reads in 8491 unique sequences.
    ## Sample 57 - 24768 reads in 5682 unique sequences.
    ## Sample 58 - 32089 reads in 8279 unique sequences.
    ## Sample 59 - 3044 reads in 956 unique sequences.
    ## Sample 60 - 22006 reads in 5663 unique sequences.
    ## Sample 61 - 973 reads in 276 unique sequences.
    ## Sample 62 - 49974 reads in 10170 unique sequences.
    ## Sample 63 - 51856 reads in 11470 unique sequences.
    ## Sample 64 - 91493 reads in 27560 unique sequences.
    ## Sample 65 - 28051 reads in 8427 unique sequences.
    ## Sample 66 - 42234 reads in 9389 unique sequences.
    ## Sample 67 - 47478 reads in 9855 unique sequences.
    ## Sample 68 - 72706 reads in 20010 unique sequences.
    ## Sample 69 - 37210 reads in 13125 unique sequences.
    ## Sample 70 - 31861 reads in 8132 unique sequences.
    ## Sample 71 - 32974 reads in 7533 unique sequences.
    ## Sample 72 - 33059 reads in 8122 unique sequences.
    ## Sample 73 - 33842 reads in 7945 unique sequences.
    ## Sample 74 - 50705 reads in 19997 unique sequences.
    ## Sample 75 - 105483 reads in 34929 unique sequences.
    ## Sample 76 - 71950 reads in 28129 unique sequences.
    ## Sample 77 - 34489 reads in 12426 unique sequences.
    ## Sample 78 - 25886 reads in 9845 unique sequences.
    ## Sample 79 - 79701 reads in 28205 unique sequences.
    ## Sample 80 - 32952 reads in 11322 unique sequences.
    ## Sample 81 - 48745 reads in 17084 unique sequences.
    ## Sample 82 - 33738 reads in 9974 unique sequences.
    ## Sample 83 - 36914 reads in 14115 unique sequences.
    ## Sample 84 - 40970 reads in 13891 unique sequences.
    ## Sample 85 - 42330 reads in 14676 unique sequences.
    ## Sample 86 - 34791 reads in 12486 unique sequences.
    ## Sample 87 - 35853 reads in 12932 unique sequences.
    ## Sample 88 - 73919 reads in 19032 unique sequences.

``` r
dadaRs.2 <- dada(filtRs.2, err=errR.2, multithread=TRUE)
```

    ## Sample 1 - 10185 reads in 2684 unique sequences.
    ## Sample 2 - 552 reads in 371 unique sequences.
    ## Sample 3 - 62504 reads in 32404 unique sequences.
    ## Sample 4 - 32003 reads in 16522 unique sequences.
    ## Sample 5 - 32920 reads in 15177 unique sequences.
    ## Sample 6 - 41408 reads in 13573 unique sequences.
    ## Sample 7 - 39234 reads in 10415 unique sequences.
    ## Sample 8 - 37067 reads in 13285 unique sequences.
    ## Sample 9 - 29608 reads in 11447 unique sequences.
    ## Sample 10 - 77475 reads in 37196 unique sequences.
    ## Sample 11 - 47942 reads in 23241 unique sequences.
    ## Sample 12 - 43258 reads in 17788 unique sequences.
    ## Sample 13 - 51087 reads in 23527 unique sequences.
    ## Sample 14 - 65629 reads in 17512 unique sequences.
    ## Sample 15 - 65909 reads in 15894 unique sequences.
    ## Sample 16 - 112 reads in 112 unique sequences.
    ## Sample 17 - 46570 reads in 14251 unique sequences.
    ## Sample 18 - 51853 reads in 14927 unique sequences.
    ## Sample 19 - 35860 reads in 10961 unique sequences.
    ## Sample 20 - 53928 reads in 17733 unique sequences.
    ## Sample 21 - 26824 reads in 9054 unique sequences.
    ## Sample 22 - 32912 reads in 8875 unique sequences.
    ## Sample 23 - 49762 reads in 7393 unique sequences.
    ## Sample 24 - 3490 reads in 1908 unique sequences.
    ## Sample 25 - 1150 reads in 950 unique sequences.
    ## Sample 26 - 177 reads in 171 unique sequences.
    ## Sample 27 - 148 reads in 147 unique sequences.
    ## Sample 28 - 24995 reads in 9069 unique sequences.
    ## Sample 29 - 37086 reads in 12572 unique sequences.
    ## Sample 30 - 33520 reads in 11225 unique sequences.
    ## Sample 31 - 24261 reads in 7571 unique sequences.
    ## Sample 32 - 33107 reads in 10826 unique sequences.
    ## Sample 33 - 42225 reads in 14846 unique sequences.
    ## Sample 34 - 40202 reads in 12890 unique sequences.
    ## Sample 35 - 215 reads in 208 unique sequences.
    ## Sample 36 - 1578 reads in 574 unique sequences.
    ## Sample 37 - 58540 reads in 18069 unique sequences.
    ## Sample 38 - 43474 reads in 11386 unique sequences.
    ## Sample 39 - 11053 reads in 2986 unique sequences.
    ## Sample 40 - 4310 reads in 1156 unique sequences.
    ## Sample 41 - 55935 reads in 28454 unique sequences.
    ## Sample 42 - 38606 reads in 19259 unique sequences.
    ## Sample 43 - 31065 reads in 13175 unique sequences.
    ## Sample 44 - 99396 reads in 32316 unique sequences.
    ## Sample 45 - 118309 reads in 38381 unique sequences.
    ## Sample 46 - 111599 reads in 33507 unique sequences.
    ## Sample 47 - 89632 reads in 28710 unique sequences.
    ## Sample 48 - 35287 reads in 11825 unique sequences.
    ## Sample 49 - 27880 reads in 9436 unique sequences.
    ## Sample 50 - 39234 reads in 10835 unique sequences.
    ## Sample 51 - 77005 reads in 21556 unique sequences.
    ## Sample 52 - 60498 reads in 34919 unique sequences.
    ## Sample 53 - 30306 reads in 11237 unique sequences.
    ## Sample 54 - 52678 reads in 19397 unique sequences.
    ## Sample 55 - 29046 reads in 7606 unique sequences.
    ## Sample 56 - 35118 reads in 9031 unique sequences.
    ## Sample 57 - 24768 reads in 6251 unique sequences.
    ## Sample 58 - 32089 reads in 10252 unique sequences.
    ## Sample 59 - 3044 reads in 2342 unique sequences.
    ## Sample 60 - 22006 reads in 7453 unique sequences.
    ## Sample 61 - 973 reads in 780 unique sequences.
    ## Sample 62 - 49974 reads in 13814 unique sequences.
    ## Sample 63 - 51856 reads in 17613 unique sequences.
    ## Sample 64 - 91493 reads in 38917 unique sequences.
    ## Sample 65 - 28051 reads in 8740 unique sequences.
    ## Sample 66 - 42234 reads in 9851 unique sequences.
    ## Sample 67 - 47478 reads in 10131 unique sequences.
    ## Sample 68 - 72706 reads in 20334 unique sequences.
    ## Sample 69 - 37210 reads in 13748 unique sequences.
    ## Sample 70 - 31861 reads in 8421 unique sequences.
    ## Sample 71 - 32974 reads in 8581 unique sequences.
    ## Sample 72 - 33059 reads in 9892 unique sequences.
    ## Sample 73 - 33842 reads in 9220 unique sequences.
    ## Sample 74 - 50705 reads in 19055 unique sequences.
    ## Sample 75 - 105483 reads in 34445 unique sequences.
    ## Sample 76 - 71950 reads in 26998 unique sequences.
    ## Sample 77 - 34489 reads in 11668 unique sequences.
    ## Sample 78 - 25886 reads in 9557 unique sequences.
    ## Sample 79 - 79701 reads in 27096 unique sequences.
    ## Sample 80 - 32952 reads in 11885 unique sequences.
    ## Sample 81 - 48745 reads in 18915 unique sequences.
    ## Sample 82 - 33738 reads in 10655 unique sequences.
    ## Sample 83 - 36914 reads in 14157 unique sequences.
    ## Sample 84 - 40970 reads in 13801 unique sequences.
    ## Sample 85 - 42330 reads in 14100 unique sequences.
    ## Sample 86 - 34791 reads in 12731 unique sequences.
    ## Sample 87 - 35853 reads in 13907 unique sequences.
    ## Sample 88 - 73919 reads in 21216 unique sequences.

``` r
dadaFs.2[[75]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 1136 sequence variants were inferred from 34929 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

``` r
dadaRs.2[[75]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 1067 sequence variants were inferred from 34445 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

``` r
mergers.2 <- mergePairs(dadaFs.2, filtFs.2, dadaRs.2, filtRs.2, verbose=TRUE)
```

    ## 9972 paired-reads (in 15 unique pairings) successfully merged out of 10018 (in 26 pairings) input.

    ## 350 paired-reads (in 26 unique pairings) successfully merged out of 379 (in 40 pairings) input.

    ## 52186 paired-reads (in 1419 unique pairings) successfully merged out of 56307 (in 2512 pairings) input.

    ## 26002 paired-reads (in 852 unique pairings) successfully merged out of 28255 (in 1454 pairings) input.

    ## 28379 paired-reads (in 749 unique pairings) successfully merged out of 30143 (in 1235 pairings) input.

    ## 38148 paired-reads (in 591 unique pairings) successfully merged out of 39552 (in 1028 pairings) input.

    ## 37379 paired-reads (in 222 unique pairings) successfully merged out of 38272 (in 419 pairings) input.

    ## 33414 paired-reads (in 597 unique pairings) successfully merged out of 35185 (in 987 pairings) input.

    ## 26349 paired-reads (in 522 unique pairings) successfully merged out of 27776 (in 844 pairings) input.

    ## 67069 paired-reads (in 1543 unique pairings) successfully merged out of 71443 (in 2647 pairings) input.

    ## 40272 paired-reads (in 1039 unique pairings) successfully merged out of 44093 (in 2102 pairings) input.

    ## 38459 paired-reads (in 845 unique pairings) successfully merged out of 40611 (in 1392 pairings) input.

    ## 43304 paired-reads (in 935 unique pairings) successfully merged out of 46627 (in 1714 pairings) input.

    ## 62289 paired-reads (in 464 unique pairings) successfully merged out of 63937 (in 834 pairings) input.

    ## 64655 paired-reads (in 374 unique pairings) successfully merged out of 65332 (in 473 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 3 (in 1 pairings) input.

    ## 43588 paired-reads (in 782 unique pairings) successfully merged out of 45065 (in 1174 pairings) input.

    ## 49268 paired-reads (in 651 unique pairings) successfully merged out of 50437 (in 941 pairings) input.

    ## 35278 paired-reads (in 475 unique pairings) successfully merged out of 35545 (in 566 pairings) input.

    ## 52776 paired-reads (in 682 unique pairings) successfully merged out of 53414 (in 871 pairings) input.

    ## 26437 paired-reads (in 136 unique pairings) successfully merged out of 26693 (in 182 pairings) input.

    ## 32304 paired-reads (in 87 unique pairings) successfully merged out of 32733 (in 147 pairings) input.

    ## 49572 paired-reads (in 29 unique pairings) successfully merged out of 49603 (in 49 pairings) input.

    ## 3448 paired-reads (in 3 unique pairings) successfully merged out of 3474 (in 6 pairings) input.

    ## 1036 paired-reads (in 8 unique pairings) successfully merged out of 1057 (in 12 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 48 (in 5 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 27 (in 1 pairings) input.

    ## 23045 paired-reads (in 327 unique pairings) successfully merged out of 23840 (in 551 pairings) input.

    ## 34372 paired-reads (in 487 unique pairings) successfully merged out of 35432 (in 762 pairings) input.

    ## 30975 paired-reads (in 329 unique pairings) successfully merged out of 32110 (in 623 pairings) input.

    ## 22510 paired-reads (in 219 unique pairings) successfully merged out of 23097 (in 401 pairings) input.

    ## 31187 paired-reads (in 405 unique pairings) successfully merged out of 31946 (in 637 pairings) input.

    ## 39258 paired-reads (in 535 unique pairings) successfully merged out of 40601 (in 827 pairings) input.

    ## 38114 paired-reads (in 513 unique pairings) successfully merged out of 39194 (in 803 pairings) input.

    ## 35 paired-reads (in 1 unique pairings) successfully merged out of 39 (in 2 pairings) input.

    ## 1444 paired-reads (in 19 unique pairings) successfully merged out of 1477 (in 29 pairings) input.

    ## 56366 paired-reads (in 554 unique pairings) successfully merged out of 57688 (in 886 pairings) input.

    ## 41758 paired-reads (in 448 unique pairings) successfully merged out of 42743 (in 628 pairings) input.

    ## 10900 paired-reads (in 120 unique pairings) successfully merged out of 10941 (in 137 pairings) input.

    ## 4082 paired-reads (in 19 unique pairings) successfully merged out of 4158 (in 29 pairings) input.

    ## 47257 paired-reads (in 843 unique pairings) successfully merged out of 52137 (in 1925 pairings) input.

    ## 32221 paired-reads (in 589 unique pairings) successfully merged out of 35442 (in 1296 pairings) input.

    ## 27982 paired-reads (in 282 unique pairings) successfully merged out of 29552 (in 594 pairings) input.

    ## 91536 paired-reads (in 1032 unique pairings) successfully merged out of 95309 (in 1918 pairings) input.

    ## 110426 paired-reads (in 1811 unique pairings) successfully merged out of 114444 (in 2857 pairings) input.

    ## 105164 paired-reads (in 1451 unique pairings) successfully merged out of 108258 (in 2330 pairings) input.

    ## 84462 paired-reads (in 1055 unique pairings) successfully merged out of 86987 (in 1674 pairings) input.

    ## 33588 paired-reads (in 504 unique pairings) successfully merged out of 34277 (in 688 pairings) input.

    ## 26139 paired-reads (in 353 unique pairings) successfully merged out of 26923 (in 539 pairings) input.

    ## 37546 paired-reads (in 586 unique pairings) successfully merged out of 38310 (in 798 pairings) input.

    ## 74680 paired-reads (in 509 unique pairings) successfully merged out of 76101 (in 825 pairings) input.

    ## 49298 paired-reads (in 932 unique pairings) successfully merged out of 55291 (in 2294 pairings) input.

    ## 27211 paired-reads (in 642 unique pairings) successfully merged out of 28375 (in 943 pairings) input.

    ## 48885 paired-reads (in 646 unique pairings) successfully merged out of 50851 (in 1180 pairings) input.

    ## 27629 paired-reads (in 321 unique pairings) successfully merged out of 28223 (in 468 pairings) input.

    ## 33784 paired-reads (in 340 unique pairings) successfully merged out of 34405 (in 484 pairings) input.

    ## 23823 paired-reads (in 242 unique pairings) successfully merged out of 24160 (in 327 pairings) input.

    ## 30649 paired-reads (in 359 unique pairings) successfully merged out of 31355 (in 484 pairings) input.

    ## 2269 paired-reads (in 12 unique pairings) successfully merged out of 2728 (in 34 pairings) input.

    ## 21093 paired-reads (in 175 unique pairings) successfully merged out of 21542 (in 243 pairings) input.

    ## 868 paired-reads (in 7 unique pairings) successfully merged out of 875 (in 8 pairings) input.

    ## 47518 paired-reads (in 397 unique pairings) successfully merged out of 48784 (in 673 pairings) input.

    ## 48974 paired-reads (in 487 unique pairings) successfully merged out of 50474 (in 812 pairings) input.

    ## 84564 paired-reads (in 750 unique pairings) successfully merged out of 88759 (in 1571 pairings) input.

    ## 26730 paired-reads (in 427 unique pairings) successfully merged out of 27334 (in 556 pairings) input.

    ## 40076 paired-reads (in 651 unique pairings) successfully merged out of 40936 (in 892 pairings) input.

    ## 46163 paired-reads (in 490 unique pairings) successfully merged out of 46964 (in 630 pairings) input.

    ## 68786 paired-reads (in 753 unique pairings) successfully merged out of 70814 (in 1305 pairings) input.

    ## 34431 paired-reads (in 497 unique pairings) successfully merged out of 35393 (in 752 pairings) input.

    ## 30768 paired-reads (in 190 unique pairings) successfully merged out of 31251 (in 306 pairings) input.

    ## 32149 paired-reads (in 169 unique pairings) successfully merged out of 32548 (in 234 pairings) input.

    ## 31991 paired-reads (in 222 unique pairings) successfully merged out of 32405 (in 294 pairings) input.

    ## 32947 paired-reads (in 197 unique pairings) successfully merged out of 33373 (in 272 pairings) input.

    ## 46063 paired-reads (in 795 unique pairings) successfully merged out of 47760 (in 1272 pairings) input.

    ## 97799 paired-reads (in 955 unique pairings) successfully merged out of 101613 (in 1957 pairings) input.

    ## 64442 paired-reads (in 888 unique pairings) successfully merged out of 67971 (in 1908 pairings) input.

    ## 31248 paired-reads (in 576 unique pairings) successfully merged out of 32667 (in 1008 pairings) input.

    ## 23478 paired-reads (in 422 unique pairings) successfully merged out of 24470 (in 730 pairings) input.

    ## 71742 paired-reads (in 945 unique pairings) successfully merged out of 75749 (in 1930 pairings) input.

    ## 30054 paired-reads (in 498 unique pairings) successfully merged out of 31278 (in 789 pairings) input.

    ## 46082 paired-reads (in 576 unique pairings) successfully merged out of 47370 (in 922 pairings) input.

    ## 31063 paired-reads (in 391 unique pairings) successfully merged out of 32036 (in 662 pairings) input.

    ## 32743 paired-reads (in 470 unique pairings) successfully merged out of 34596 (in 1040 pairings) input.

    ## 37476 paired-reads (in 471 unique pairings) successfully merged out of 39006 (in 913 pairings) input.

    ## 38571 paired-reads (in 536 unique pairings) successfully merged out of 40303 (in 1119 pairings) input.

    ## 31388 paired-reads (in 400 unique pairings) successfully merged out of 32788 (in 868 pairings) input.

    ## 32516 paired-reads (in 433 unique pairings) successfully merged out of 33879 (in 792 pairings) input.

    ## 68774 paired-reads (in 464 unique pairings) successfully merged out of 71446 (in 1158 pairings) input.

``` r
# quick trim
seqtab.2 <- makeSequenceTable(mergers.2)
# can inspect the distribution of sequences like this:
# table(nchar(getSequences(seqtab.2)))
seqtab.trim.2 <- seqtab.2[,nchar(colnames(seqtab.2)) %in% 212:217] # trim the sequences to 212-217 bp

# remove chimeras and graph reads through the workflow
seqtab.nochim.2 <- removeBimeraDenovo(seqtab.trim.2, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 705 bimeras out of 17795 input sequences.

``` r
sum(seqtab.nochim.2)/sum(seqtab.trim.2) # 97.9% of sequences are not chimeras
```

    ## [1] 0.9794688

``` r
getN <- function(x) sum(getUniques(x))
track.2 <- cbind(out.2, sapply(dadaFs.2, getN), sapply(dadaRs.2, getN), sapply(mergers.2, getN), rowSums(seqtab.nochim.2))
colnames(track.2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track.2) <- sample.names.2
tdf.2 <- as.data.frame(track.2)
tdf.2$filtered_out <- tdf.2$input-tdf.2$filtered
tdf.2$noise_removed <- tdf.2$filtered-with(tdf.2, pmin(denoisedF, denoisedR))
tdf.2$unmerged <- (tdf.2$filtered-tdf.2$noise_removed)-tdf.2$merged
tdf.2$chimeras <- tdf.2$merged-tdf.2$nonchim
tdf.2 <- data.frame(sample = row.names(tdf.2), tdf.2)
# select only the columns we want to plot:
tdfs.2 <- tdf.2[,c(1,7,8,9,10,11)]
tdfl.2 <- gather(tdfs.2, step, reads, nonchim:chimeras, factor_key=FALSE)

# reorder the steps to plot them in an order that makes sense:
tdfl.2$step <- factor(tdfl.2$step, levels = c("filtered_out","noise_removed","unmerged","chimeras", "nonchim"))
track.reads.2 <- ggplot(tdfl.2,aes(sample,reads,fill=step))
track.reads.plot.2 <- track.reads.2+
  geom_bar(stat="identity")+
  scale_y_continuous(breaks = c(50000,100000,150000,200000,250000))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
track.reads.plot.2
```

![](sequence_processing_only_8.20_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

``` r
# assign taxonomy (change this path to wherever your database is!)
taxa.2 <- assignTaxonomy(seqtab.nochim.2, "~/Desktop/flux_resubmission/silva_nr_v132_train_set.fa", multithread=TRUE)
taxa.print.2 <- taxa.2
rownames(taxa.print.2) <- NULL

# print percentages of reads identified to each taxonomic level
message("% phyla identified is ", (1-(sum(is.na(taxa.print.2[,2]))/length(taxa.print.2[,2])))*100)
```

    ## % phyla identified is 88.7419543592744

``` r
message("% classes identified is ", (1-(sum(is.na(taxa.print.2[,3]))/length(taxa.print.2[,3])))*100)
```

    ## % classes identified is 76.7349327091867

``` r
message("% orders identified is ", (1-(sum(is.na(taxa.print.2[,4]))/length(taxa.print.2[,4])))*100)
```

    ## % orders identified is 52.0830895260386

``` r
message("% families identified is ", (1-(sum(is.na(taxa.print.2[,5]))/length(taxa.print.2[,5])))*100)
```

    ## % families identified is 31.6676418958455

``` r
message("% genera identified is ", (1-(sum(is.na(taxa.print.2[,6]))/length(taxa.print.2[,6])))*100)
```

    ## % genera identified is 13.6102984201287

## make a phylogenetic tree

This took a few hours on my laptop.

``` r
seqs <- getSequences(seqtab.nochim.2)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA) # this took ~50 minutes
```

    ## Determining distance matrix based on shared 8-mers:
    ## ================================================================================
    ## 
    ## Time difference of 2191.93 secs
    ## 
    ## Clustering into groups by similarity:
    ## ================================================================================
    ## 
    ## Time difference of 90.32 secs
    ## 
    ## Aligning Sequences:
    ## ================================================================================
    ## 
    ## Time difference of 157.86 secs
    ## 
    ## Iteration 1 of 2:
    ## 
    ## Determining distance matrix based on alignment:
    ## ================================================================================
    ## 
    ## Time difference of 183.66 secs
    ## 
    ## Reclustering into groups by similarity:
    ## ================================================================================
    ## 
    ## Time difference of 51.9 secs
    ## 
    ## Realigning Sequences:
    ## ================================================================================
    ## 
    ## Time difference of 72.89 secs
    ## 
    ## Iteration 2 of 2:
    ## 
    ## Determining distance matrix based on alignment:
    ## ================================================================================
    ## 
    ## Time difference of 187.12 secs
    ## 
    ## Reclustering into groups by similarity:
    ## ================================================================================
    ## 
    ## Time difference of 45.15 secs
    ## 
    ## Realigning Sequences:
    ## ================================================================================
    ## 
    ## Time difference of 16.32 secs
    ## 
    ## Refining the alignment:
    ## ================================================================================
    ## 
    ## Time difference of 1.28 secs

``` r
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align) # started 3:05 pm, finished around an hour later? 
treeNJ <- NJ(dm) # Note, tip order != sequence order ## this takes at least an hour... started 6:01 pm
fit <- pml(treeNJ, data=phang.align)
```

    ## negative edges length changed to 0!

``` r
fitGTR <- update(fit, k=4, inv=0.2)
```

## Add in metadata and create phyloseq object

Very quick to run. The metadata file contains sample name, core, depth
below seafloor, pingo, sulfate, sulfide, alkalinity concentrations,
core\_flowtype (increasing “inc” or steady state “ss”), geochem\_zone
(linear “lin”, nss “non-steady-state”, or below sulfate-methane
transition zone “below”), whether the read quality was ok, stage
(fluxincreasing, steadystate, or seep), sample or control, library
concentration (ng/ul), who extracted the DNA, and mcrA and dsrAB ddPCR
gene counts in copies per gram.

``` r
metadata.2 <- read.csv(file="~/Desktop/svalflux/metadata.csv")
rownames(metadata.2) <- sample.names.2
ps.f <- phyloseq(tax_table(taxa.2),
                 sample_data(metadata.2),
                 phy_tree(fitGTR$tree),
                 otu_table(seqtab.nochim.2, taxa_are_rows = FALSE))
dna <- Biostrings::DNAStringSet(taxa_names(ps.f))
names(dna) <- taxa_names(ps.f)
ps.f <- merge_phyloseq(ps.f, dna)
taxa_names(ps.f) <- paste0("ASV", seq(ntaxa(ps.f)))
ps.f # the first phyloseq object
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 17090 taxa and 88 samples ]
    ## sample_data() Sample Data:       [ 88 samples by 16 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 17090 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 17090 tips and 17088 internal nodes ]
    ## refseq()      DNAStringSet:      [ 17090 reference sequences ]

## Removal of taxa, pruning, and decontamination steps

These chunks only took several minutes, and most steps could be run
interactively. The decontam steps were the longest ones for me.

``` r
# first make a table of taxa that I may or may not want to keep:
ASV_classifications <- c("Bacteria","Archaea","Eukaryota", "Unclassified", "Chloroplast","Mitochondria")
num_ASVs <- c(sum(tax_table(ps.f)[,1]=="Bacteria", na.rm = TRUE),
               sum(tax_table(ps.f)[,1]=="Archaea", na.rm=TRUE),
               sum(tax_table(ps.f)[,1]=="Eukaryota", na.rm=TRUE),
               sum(is.na(tax_table(ps.f)[,1])),
               sum(tax_table(ps.f)[,4]=="Chloroplast", na.rm=TRUE), 
               sum(tax_table(ps.f)[,5]=="Mitochondria", na.rm=TRUE))
num_ASVcounts <- c(sum(otu_table(ps.f)[,which(tax_table(ps.f)[,1]=="Bacteria")], na.rm = TRUE),
              sum(otu_table(ps.f)[,which(tax_table(ps.f)[,1]=="Archaea")], na.rm = TRUE),
              sum(otu_table(ps.f)[,which(tax_table(ps.f)[,1]=="Eukaryota")], na.rm = TRUE), 
              sum(is.na(tax_table(ps.f)[,1])),
              sum(otu_table(ps.f)[,which(tax_table(ps.f)[,4]=="Chloroplast")], na.rm = TRUE), 
              sum(otu_table(ps.f)[,which(tax_table(ps.f)[,5]=="Mitochondria")], na.rm = TRUE))
asv.table <- cbind.data.frame(ASV_classifications, num_ASVs, num_ASVcounts)
asv.table
```

    ##   ASV_classifications num_ASVs num_ASVcounts
    ## 1            Bacteria    14032       2565941
    ## 2             Archaea     2607        604967
    ## 3           Eukaryota       73          3851
    ## 4        Unclassified      378           378
    ## 5         Chloroplast       54          3447
    ## 6        Mitochondria       11           304

remove Eukarya, Chloroplasts,
Mitochondria

``` r
ps.fr <- subset_taxa(ps.f, (Kingdom!="Eukaryota")) # unknown and eukaryote ASVs removed
ps.fr <- subset_taxa(ps.fr, (Order!="Chloroplast") | is.na(Order)) # chloroplasts removed
ps.fr <- subset_taxa(ps.fr, (Family!="Mitochondria") | is.na(Family)) # mitochondria removed
ntaxa(ps.f)-ntaxa(ps.fr) # number of taxa removed 
```

    ## [1] 516

``` r
ps.fr # phyloseq object with these sequences removed
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 16574 taxa and 88 samples ]
    ## sample_data() Sample Data:       [ 88 samples by 16 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 16574 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 16574 tips and 16572 internal nodes ]
    ## refseq()      DNAStringSet:      [ 16574 reference sequences ]

run decontam using the combined approach, which combines the scores from
prevalence and frequency using Fisher’s
method

``` r
sample_data(ps.fr)$is.neg <- sample_data(ps.fr)$Sample_or_Control == "Control Sample"
contamdf.flux.comb <- isContaminant(ps.fr, method="combined", neg="is.neg", conc="quant_reading")
```

    ## Warning in isContaminant(ps.fr, method = "combined", neg = "is.neg", conc =
    ## "quant_reading"): Removed 3 samples with zero total counts (or frequency).

``` r
table(contamdf.flux.comb$contaminant) # prints the number of contaminants
```

    ## 
    ## FALSE  TRUE 
    ## 16493    81

``` r
flux.comb.contam <- (which(contamdf.flux.comb$contaminant))
flux.comb.contam.mx <- subset(tax_table(ps.fr)[flux.comb.contam,])
```

plot library sizes of true samples and blanks

``` r
df.flux <- as.data.frame(sample_data(ps.fr)) 
df.flux$LibrarySize <- sample_sums(ps.fr)
df.flux <- df.flux[order(df.flux$LibrarySize),]
df.flux$Index <- seq(nrow(df.flux))
flux.sample.depth.samples.and.controls <- ggplot(data=df.flux, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point() + ggtitle("Before decontamination")
flux.sample.depth.samples.and.controls 
```

![](sequence_processing_only_8.20_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

Remove the contaminants

``` r
contam <- rownames(flux.comb.contam.mx) # the 81 contaminant names
allnames <- rownames(tax_table(ps.fr)) # the taxa names
flux.uncontam <- allnames[!allnames %in% contam] # gives us the noncontaminant names
ps.frd <- prune_taxa(flux.uncontam, ps.fr) # creates new phyloseq object without contaminants
ps.frd # decontaminated phyloseq object
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 16493 taxa and 88 samples ]
    ## sample_data() Sample Data:       [ 88 samples by 17 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 16493 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 16493 tips and 16491 internal nodes ]
    ## refseq()      DNAStringSet:      [ 16493 reference sequences ]

I remembered there were Micrococcus in some groups that I had removed
from my previous OTU dataset. Some of them required manual removal. Just
to double check, I plotted library sizes
again:

``` r
which(flux.comb.contam.mx[,6]=="Micrococcus") # no, decontam did not
```

    ## integer(0)

``` r
sum(otu_table(ps.fr)[,which(tax_table(ps.fr)[,6]=="Micrococcus")])/sum(otu_table(ps.fr))*100 # 0.67% of sequences
```

    ## [1] 0.6690227

``` r
ps.blanks <- subset_samples(ps.fr, Sample_or_Control=="Control Sample") # phyloseq object with only the blanks 
blanks.top100 <- names(sort(taxa_sums(ps.blanks), decreasing=TRUE))[1:100]
blanks.ps.ra <- transform_sample_counts(ps.blanks, function(OTU) OTU/sum(OTU))
blanks.ps.top100 <- prune_taxa(blanks.top100, blanks.ps.ra)
blanks_barplot <- plot_bar(blanks.ps.top100, fill="Class")
blanks_barplot # yes indeed these blanks are different
```

![](sequence_processing_only_8.20_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
df.flux["EBkati-000.0A","LibrarySize"] # over 9k reads
```

    ##               LibrarySize
    ## EBkati-000.0A        9972

``` r
EBkati.ra <- get_taxa(blanks.ps.ra, "EBkati-000.0A") # extract relative abundances from one of our blanks
EBkati.ra <- subset(EBkati.ra, EBkati.ra>0) # omit ASVs not detected in the blank
EBkati.ra <- sort(EBkati.ra, decreasing = TRUE) # sort the list of ASVs by descending order
tax_table(ps.blanks)[names(EBkati.ra)[EBkati.ra>0],] # lists taxonomies of ASVs in the sample by descending relative abundance. Here we see the two most abundant ASVs belong to micrococcus. What is the abundance of Micrococcus across all samples?
```

    ## Taxonomy Table:     [15 taxa by 6 taxonomic ranks]:
    ##         Kingdom    Phylum           Class                
    ## ASV27   "Bacteria" "Actinobacteria" "Actinobacteria"     
    ## ASV95   "Bacteria" "Actinobacteria" "Actinobacteria"     
    ## ASV175  "Bacteria" "Actinobacteria" "Actinobacteria"     
    ## ASV111  "Bacteria" "Proteobacteria" "Gammaproteobacteria"
    ## ASV717  "Bacteria" "Actinobacteria" "Actinobacteria"     
    ## ASV1065 "Bacteria" "Proteobacteria" "Alphaproteobacteria"
    ## ASV1269 "Bacteria" "Proteobacteria" "Gammaproteobacteria"
    ## ASV1383 "Bacteria" "Bacteroidetes"  "Bacteroidia"        
    ## ASV588  "Bacteria" "Firmicutes"     "Bacilli"            
    ## ASV2099 "Bacteria" "Actinobacteria" "Actinobacteria"     
    ## ASV2505 "Bacteria" "Proteobacteria" "Gammaproteobacteria"
    ## ASV139  "Bacteria" "Proteobacteria" "Gammaproteobacteria"
    ## ASV1084 "Archaea"  "Euryarchaeota"  "Methanomicrobia"    
    ## ASV6554 "Bacteria" "Proteobacteria" "Deltaproteobacteria"
    ## ASV457  "Bacteria" "Proteobacteria" "Gammaproteobacteria"
    ##         Order                   Family                Genus                 
    ## ASV27   "Micrococcales"         "Micrococcaceae"      "Micrococcus"         
    ## ASV95   "Micrococcales"         "Micrococcaceae"      "Micrococcus"         
    ## ASV175  "Micrococcales"         "Dermabacteraceae"    "Brachybacterium"     
    ## ASV111  "Pseudomonadales"       "Moraxellaceae"       "Enhydrobacter"       
    ## ASV717  "Micrococcales"         "Micrococcaceae"      "Kocuria"             
    ## ASV1065 "Sphingomonadales"      "Sphingomonadaceae"   "Qipengyuania"        
    ## ASV1269 "Pseudomonadales"       "Pseudomonadaceae"    "Pseudomonas"         
    ## ASV1383 "Flavobacteriales"      "Flavobacteriaceae"   "Flavobacterium"      
    ## ASV588  "Lactobacillales"       "Streptococcaceae"    "Streptococcus"       
    ## ASV2099 "Micrococcales"         "Micrococcaceae"      "Micrococcus"         
    ## ASV2505 "Enterobacteriales"     "Enterobacteriaceae"  "Klebsiella"          
    ## ASV139  "Enterobacteriales"     "Enterobacteriaceae"  "Escherichia/Shigella"
    ## ASV1084 "Methanosarcinales"     "Methanosarcinaceae"  "Methanococcoides"    
    ## ASV6554 "Desulfuromonadales"    "Desulfuromonadaceae" "Desulfuromonas"      
    ## ASV457  "Betaproteobacteriales" "Burkholderiaceae"    "Variovorax"

``` r
ps.fr.ra <- transform_sample_counts(ps.fr, function(OTU) OTU/sum(OTU)) # first transform ps.fr to relative abundances
micrococcus <- subset_taxa(ps.fr.ra, Genus=="Micrococcus")
plot_bar(micrococcus) # only the 1048 201, 260, 261 depths also have Micrococcus in significant amounts, ~20%. They all have >25k reads:
```

    ## Warning: Removed 12 rows containing missing values (position_stack).

![](sequence_processing_only_8.20_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

``` r
# plot library sizes of true samples and blanks after decontamination
df.dc.flux <- as.data.frame(sample_data(ps.frd)) 
df.dc.flux$LibrarySize <- sample_sums(ps.frd)
df.dc.flux <- df.dc.flux[order(df.dc.flux$LibrarySize),]
df.dc.flux$Index <- seq(nrow(df.dc.flux))
flux.dc.sample.depth.samples.and.controls <- ggplot(data=df.dc.flux, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point() + ggtitle("After decontamination")
flux.dc.sample.depth.samples.and.controls # now the other negative control is lower
```

![](sequence_processing_only_8.20_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->

``` r
ps.frd <- subset_taxa(ps.frd, (Genus!="Micrococcus") | is.na(Genus)) # removed 4 Micrococcus ASVs from ps.frd
ps.frd
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 16489 taxa and 88 samples ]
    ## sample_data() Sample Data:       [ 88 samples by 17 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 16489 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 16489 tips and 16487 internal nodes ]
    ## refseq()      DNAStringSet:      [ 16489 reference sequences ]

``` r
# of the 85 contaminants, what % of the reads in the dataset?
sum(otu_table(ps.fr)[,contam])/sum(otu_table(ps.fr))*100 + sum(otu_table(ps.fr)[,which(tax_table(ps.fr)[,6]=="Micrococcus")])/sum(otu_table(ps.fr))*100 # 1.05%
```

    ## [1] 1.045701

prune samples with low-abundance reads to create the final phyloseq
object

``` r
df.dc.flux[1:13,18] # here are our lowest 13 libraries
```

    ##               LibrarySize
    ## GC1048-121.0A           0
    ## GC1068-021.0A           0
    ## GC1068-041.0A           0
    ## GC1068-201.0A          35
    ## EBstel-000.0A         306
    ## GC1070-068.0B         848
    ## GC1048-320.0A        1036
    ## GC1068-226.0A        1050
    ## GC1070-063.0B        2267
    ## GC1048-305.0B        3426
    ## GC1068-282.0A        3816
    ## GC1068-255.0A        8931
    ## EBkati-000.0A        9697

``` r
ps.frdp <- prune_samples(sample_sums(ps.frd)>=8931, ps.frd) # 12 samples with less than 8931 reads removed (max library size removed was 3816)
ps.frdp # PhyloSeq object Flux, Removed unwanted taxa, Decontaminated, Pruned
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 16489 taxa and 76 samples ]
    ## sample_data() Sample Data:       [ 76 samples by 17 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 16489 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 16489 tips and 16487 internal nodes ]
    ## refseq()      DNAStringSet:      [ 16489 reference sequences ]

add more metadata (SMT depth, peakAOM depth, distance above/below
peakAOM depth, and classification of samples as above or below SMT) and
save
it

``` r
sample_data(ps.frdp)[which(sample_data(ps.frdp)$core=="GC1048"),8] <- "ss" # overwrites "inc" to "ss" for core GC1048
sample_data(ps.frdp)$smt <- NA # write a new metadata column for SMT depth
sample_data(ps.frdp)[which(sample_data(ps.frdp)$core=="GC1045"),18] <- 82 # assign SMT depths by core
sample_data(ps.frdp)[which(sample_data(ps.frdp)$core=="GC1048"),18] <- 320
sample_data(ps.frdp)[which(sample_data(ps.frdp)$core=="GC1068"),18] <- 108
sample_data(ps.frdp)[which(sample_data(ps.frdp)$core=="GC1069"),18] <- 138
sample_data(ps.frdp)[which(sample_data(ps.frdp)$core=="GC1070"),18] <- 69
sample_data(ps.frdp)[which(sample_data(ps.frdp)$core=="GC1081"),18] <- 56

sample_data(ps.frdp)$peakaom <- NA # creates a new metada column for peak AOM depth (CHECK THIS)
sample_data(ps.frdp)[which(sample_data(ps.frdp)$core=="GC1045"),19] <- 65 # assign peak AOM depths by core
sample_data(ps.frdp)[which(sample_data(ps.frdp)$core=="GC1048"),19] <- 313.5
sample_data(ps.frdp)[which(sample_data(ps.frdp)$core=="GC1068"),19] <- 108
sample_data(ps.frdp)[which(sample_data(ps.frdp)$core=="GC1069"),19] <- 137
sample_data(ps.frdp)[which(sample_data(ps.frdp)$core=="GC1070"),19] <- 67
sample_data(ps.frdp)[which(sample_data(ps.frdp)$core=="GC1081"),19] <- 56

sample_data(ps.frdp)$smtzposition <- "above" # creates default value "above" for samples 
sample_data(ps.frdp)[which(sample_data(ps.frdp)$geochem_zone=="below"),20] <- "below" # overwrites as "below" for samples where geochem_zone = below
sample_data(ps.frdp)$dist_above_peakaom <- sample_data(ps.frdp)$peakaom-sample_data(ps.frdp)$depth # creates a new variable, distance above peak AOM
sample_data(ps.frdp)[which(sample_data(ps.frdp)$core=="GC1048"&sample_data(ps.frdp)$depth<310),9] <- "lin" # GC1048 is a steady-state core, so geochem zone should be 
saveRDS(ps.frdp, "ps.frdp") # saved phyloseq object, pushed to Github (need to update with metadata)
```
