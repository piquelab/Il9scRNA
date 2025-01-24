######################################################################
# 1. Cleaning up 
######################################################################
library(Seurat)
library(Matrix)
library(future)
library(tidyverse)


outFolder="./2_soupx_doubletfinder_Cleaning/"
system(paste0("mkdir -p ", outFolder))

ifn <- paste0("./1_soupx_doubletfinder_chrM/seuratObj-SoupX-Doublet-merge.2024-03-18.rds") 
sc <- read_rds(ifn)


#################################################################
# Quality control
#################################################################


####################################################
# percentage of reads mapping to mitochondrial
####################################################

cmd <- paste0("cat ","/wsu/home/groups/piquelab/data/refGenome10x/refdata-gex-mm10-2020-A/genes/genes.gtf",
              " | awk '$3~/gene/'",
              " | sed 's/gene_id //;s/;.*gene_type/\t/;s/; gene_name /\t/;s/; level.*//'")##
##cmd <- paste0("cat ","/rs/rs_grp_il9/novogene/ref/genes.gtf", 
##               " | awk '$3~/gene/'",
##               " | sed 's/gene_id //;s/;.*gene_type/\t/;s/; gene_name /\t/;s/; level.*//'")
cat(cmd,"\n")

aux <- data.table::fread(cmd=cmd) %>% mutate(TSS=ifelse(V7=="+",V4,V5))%>% select(Chr=V1,Min=V4,Max=V5,kbid=V9,TSS,Strand=V7,Type=V10,gene_name=V11) 

anno <- tibble(gene_name=rownames(sc),rs=rowSums(sc)) %>% filter(rs>0) %>% left_join(aux) %>% filter(!is.na(Chr))

table(is.na(anno$Chr))

table(anno$Chr)

table(anno$Type)

head(anno)

sc <- sc[anno$gene_name,]


sc[["percent.mt"]] <- PercentageFeatureSet(sc,features=anno[anno$Chr=="chrM",]$gene_name)

sc[["percent.mt"]] %>% summary()

mean(sc[["percent.mt"]]>20)


dens <- density(sc[["percent.mt"]]$percent.mt)
# # plot density
plot(dens, frame = FALSE, col = "steelblue", 
     main = "mitochondrial genes ") 


# ## Filter sc for things matching genotype. chrM or number of RNAs.  
# 
sc[["percent.Y"]] <- PercentageFeatureSet(sc,features=anno[anno$Chr=="chrY",]$gene_name)
sc[["percent.Y"]] %>% summary()
# 
anno[anno$gene_name=="Xist",]$kbid
sc[["percent.Xist"]] <- PercentageFeatureSet(sc,features=anno[anno$gene_name=="Xist",]$gene_name )
sc[["percent.Xist"]] %>% summary()


## VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## plot1 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "percent.mt")
## plot2 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
## CombinePlots(plots = list(plot1, plot2))

sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mt < 20)

opfn <- paste0(outFolder,"seuratObj-after-mt-filtering.",Sys.Date(),".rds") 
write_rds(sc, opfn)



##sc<-read_rds("4_harmony_cellClass_soupx_doubletfinder_chrM/seuratObj-after-mt-filtering.2021-06-28.rds")##
##sc <- readRDS("/rs/rs_grp_il9/Il9scRNA/1_soupx_doubletfinder_chrM/seuratObj-after-mt-filtering.2023-10-09.rds")


############################################################
## Clustering
############################################################

future::plan(strategy = 'multicore', workers = 12)
options(future.globals.maxSize = 8 * 1024 ^ 3)


sc <- NormalizeData(sc, verbose=TRUE)

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)

sc <- ScaleData(sc, verbose = TRUE)

sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 100, verbose = TRUE)

ElbowPlot(sc, ndims = 100, reduction = "pca")

library(harmony)

sc <- RunHarmony(sc,c("Library"))#,reduction="pca")

sc <- RunUMAP(sc,reduction = "harmony", dims = 1:30)

sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30, verbose = TRUE)


#before clustering
opfn <- paste0(outFolder,"seuratObj-before-clustering.",Sys.Date(),".rds") 
write_rds(sc, opfn)

opfn <- paste0(outFolder,"seuratObj-before-clustering.","2024-03-09",".rds") 
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

###########
###########


sc <- FindClusters(sc, verbose = TRUE,resolution=0.1)

opfn <- paste0(outFolder,"sc.NormByLibrary.cellclassify_newfilter-res0.1.",Sys.Date(),".rds") 
write_rds(sc, opfn)

fname=paste0(outFolder,"/UMAP_Harmony-res0.1.png");
png(fname,width=1000,height=1000)
DimPlot(sc, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6) #+ NoLegend()
dev.off()




####sc<-read_rds("4_harmony_cellClass_soupx_doubletfinder_chrM/sc.NormByLibrary.cellclassify_newfilter-res0.4.2021-06-28.rds")
##fname=paste0("5_harmony_cellClass_soupx-doubletfinder_chrM_plots/UMAP_Harmony-res0.4.png");
##png(fname,width=1000,height=1000)##

