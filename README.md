
# Code and data OTU-PITI

This repository provides R code and data to reproduce results and figures from the paper:

> Keck F., Vasselon V., Rimet F., Bouchez A. & Kahlert M. Boosting DNA metabarcoding for biomonitoring with phylogenetic estimation of OTUs' ecological profiles.

## Analyses

All the analyses are included in the file `main_script.R`. Be aware that running the main script can take some time (few hours on a standard laptop).
If you just can't wait, an image of the R session after all the calculations were completed is available in the results folder.
A report of the main results and figures can be generated from `report.Rmd`.

## Dependencies

Install R packages from CRAN

    install.packages(c("ape", "adephylo", "stringr",
                       "phylobase", "phylosignal",
                       "Rphylopars", "phytools", "phangorn"))


Install R packages from GitHub

    devtools::install_github("fkeck/diatobc")

Other programs

`MUSCLE`, `RAxML 8.2.10`, `PATHd8 1.0`

You need to install these programs and edit their paths in `main_script.R`.

## R Session Info

    R version 3.3.1 (2016-06-21)
    Platform: x86_64-w64-mingw32/x64 (64-bit)
    Running under: Windows >= 8 x64 (build 9200)
    
    locale:
    [1] LC_COLLATE=English_Ireland.1252  LC_CTYPE=English_Ireland.1252   
    [3] LC_MONETARY=English_Ireland.1252 LC_NUMERIC=C                    
    [5] LC_TIME=English_Ireland.1252    
    
    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods  
    [7] base     
    
    other attached packages:
    [1] diatobc_1.0      Rphylopars_0.2.9 phylosignal_1.2  phylobase_0.8.2 
    [5] stringr_1.2.0    adephylo_1.1-6   ade4_1.7-4       ape_4.1         
    
    loaded via a namespace (and not attached):
      [1] colorspace_1.2-6        seqinr_3.3-1           
      [3] deldir_0.1-12           rsconnect_0.4.3        
      [5] rprojroot_1.2           roxygen2_6.0.1         
      [7] mvtnorm_1.0-5           xml2_1.0.0             
      [9] codetools_0.2-15        splines_3.3.1          
     [11] mnormt_1.5-4            bold_0.3.5             
     [13] doBy_4.5-15             knitr_1.14             
     [15] ips_0.0-7               jsonlite_1.1           
     [17] phylolm_2.5             cluster_2.0.4          
     [19] mvnmle_0.1-11           shiny_1.0.5            
     [21] rentrez_1.0.2           httr_1.2.1             
     [23] backports_1.1.1         assertthat_0.1         
     [25] Matrix_1.2-6            lazyeval_0.2.0         
     [27] rotl_3.0.1              htmltools_0.3.5        
     [29] tools_3.3.1             bindrcpp_0.2           
     [31] igraph_1.1.2            coda_0.18-1            
     [33] gtable_0.2.0            glue_1.1.1             
     [35] taxize_0.7.9            reshape2_1.4.2         
     [37] clusterGeneration_1.3.4 dplyr_0.7.4            
     [39] maps_3.1.1              gmodels_2.16.2         
     [41] fastmatch_1.1-0         Rcpp_0.12.14           
     [43] msm_1.6.4               spdep_0.6-8            
     [45] gdata_2.17.0            nlme_3.1-128           
     [47] iterators_1.0.8         mime_0.5               
     [49] phangorn_2.2.0          gtools_3.5.0           
     [51] devtools_1.13.4         XML_3.98-1.5           
     [53] geiger_2.0.6            LearnBayes_2.15        
     [55] MASS_7.3-45             scales_0.5.0           
     [57] subplex_1.4-1           parallel_3.3.1         
     [59] expm_0.999-2            animation_2.5          
     [61] yaml_2.1.14             memoise_1.0.0          
     [63] ggplot2_2.2.1           reshape_0.8.6          
     [65] stringi_1.1.2           foreach_1.4.3          
     [67] plotrix_3.6-6           permute_0.9-4          
     [69] phytools_0.6-20         boot_1.3-18            
     [71] chron_2.3-47            rlang_0.1.6            
     [73] pkgconfig_2.0.1         commonmark_1.4         
     [75] rncl_0.6.0              evaluate_0.10          
     [77] lattice_0.20-33         purrr_0.2.4            
     [79] bindr_0.1               rredlist_0.1.0         
     [81] deSolve_1.20            plyr_1.8.4             
     [83] magrittr_1.5            R6_2.2.2               
     [85] combinat_0.0-8          DBI_0.5-1              
     [87] pillar_1.0.1            withr_1.0.2            
     [89] mgcv_1.8-12             survival_2.40-1        
     [91] scatterplot3d_0.3-40    sp_1.2-3               
     [93] tibble_1.4.1            uuid_0.1-2             
     [95] rmarkdown_1.6           RNeXML_2.0.7           
     [97] adegenet_2.0.1          grid_3.3.1             
     [99] data.table_1.9.6        vegan_2.4-1            
    [101] digest_0.6.13           xtable_1.8-2           
    [103] tidyr_0.7.2             httpuv_1.3.5           
    [105] numDeriv_2016.8-1       munsell_0.4.3          
    [107] quadprog_1.5-5