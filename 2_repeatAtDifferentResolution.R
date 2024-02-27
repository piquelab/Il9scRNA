
######################################################################
# 1. preprocessing and clustering steps
######################################################################
library(Seurat)
library(Matrix)
library(future)
library(tidyverse)


outFolder="./1_soupx_doubletfinder_chrM/"
system(paste0("mkdir -p ", outFolder))

future::plan(strategy = 'multicore', workers = 12)
options(future.globals.maxSize = 30 * 1024 ^ 3)


opfn <- paste0(outFolder,"seuratObj-before-clustering.","2023-10-09",".rds") 

sc <- read_rds(opfn)

sc <- FindClusters(sc, verbose = TRUE,resolution=0.2)

opfn <- paste0(outFolder,"sc.NormByLibrary.cellclassify_newfilter-res0.2.",Sys.Date(),".rds") 
write_rds(sc, opfn)

fname=paste0(outFolder,"/UMAP_Harmony-res0.2.png");
png(fname,width=1000,height=1000)
DimPlot(sc, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6) #+ NoLegend()
dev.off()




####sc<-read_rds("4_harmony_cellClass_soupx_doubletfinder_chrM/sc.NormByLibrary.cellclassify_newfilter-res0.4.2021-06-28.rds")
##fname=paste0("5_harmony_cellClass_soupx-doubletfinder_chrM_plots/UMAP_Harmony-res0.4.png");
##png(fname,width=1000,height=1000)##

