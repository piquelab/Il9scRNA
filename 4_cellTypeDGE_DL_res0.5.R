
#############################################################
### cluster marker identification
### 
#############################################################

library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
library(harmony)

sc <- readRDS("1_soupx_doubletfinder_chrM/sc.NormByLibrary.cellclassify_newfilter-res0.5.2023-12-27.rds")

outFolder="~/rs_grp_il9/Il9scRNA/4_celltypeDGE_IL9_DL_res0.5/"
system(paste0("mkdir -p ", outFolder))
##setwd(outFolder)

dim(sc)

table(sc$Library)

table(sc$Location) 


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 60 * 1024 ^ 3)



markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

opfn <- paste0(outFolder,"markersredo0.5.rds") 
write_rds(markers, opfn)
opfn <- paste0(outFolder,"markersredo0.5.tsv") 
write_tsv(markers,opfn)

m2 <- markers %>% arrange(cluster,-avg_log2FC) %>% group_by(cluster)


top20 <- m2 %>% top_n(n = 20, wt = avg_log2FC)

fname=paste0(outFolder,"ClusterDEG0.5.tsv");
write_tsv(m2,fname)
m2<-read_tsv(fname)
write.csv(m2,file=paste0(outFolder,"ClusterDEG0.5.csv"))

fname=paste0(outFolder,"ClusterDEGtop200.5.tsv");
write_tsv(top20,fname)
top20<-read_tsv(fname)
write.csv(top20,file=paste0(outFolder,"ClusterDEGtop200.5.csv"))


