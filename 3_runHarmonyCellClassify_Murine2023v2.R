library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
library(harmony)
library(SingleR)

########################################################
# cell type classification using SingleR 
# reference: GSE200289_Murine2023
########################################################


future::plan(strategy = 'multicore', workers = 12)
options(future.globals.maxSize = 30 * 1024 ^ 3)

outFolder="./3_harmony_cellClass_Murine2023/"
system(paste0("mkdir -p ", outFolder))

#######################################################
# seurat R objects

# uploaded in grid: 
sc1 <- read_rds("../murine2023ref/GSE200289_SeuratObject-GEO.rds")

sc2 <- read_rds("1_soupx_doubletfinder_chrM/sc.NormByLibrary.cellclassify_newfilter-res0.2.2023-12-27.rds")

length(rownames(sc1))
length(rownames(sc2))
cgenes <- intersect(rownames(sc1),rownames(sc2))
length(cgenes)


md1 <- sc1@meta.data %>% 
  mutate(BARCODE=paste0("R_",BARCODE)) %>%
    select(BARCODE,Library,Location,seurat_clusters,celltype)
    

md2 <- sc2@meta.data %>% rownames_to_column("BARCODE") %>% 
  mutate(BARCODE=paste0("Q_",BARCODE)) %>%
  mutate(celltype=NA) %>% 
  select(BARCODE,Library,Location,seurat_clusters,celltype)


md <- rbind(md1,md2)  %>% mutate(BARCODE2=BARCODE) %>% column_to_rownames("BARCODE")

sc1count <- GetAssayData(sc1, slot="counts", assay="RNA")   
sc2count <- GetAssayData(sc2, slot="counts", assay="RNA")   

sc1count <- sc1count[cgenes,]
colnames(sc1count) <- paste0("R_",colnames(sc1count))
sc2count <- sc2count[cgenes,]
colnames(sc2count) <- paste0("Q_",colnames(sc2count))


identical(md1$BARCODE,colnames(sc1count))
identical(md2$BARCODE,colnames(sc2count))


scc <- cbind(sc1count,sc2count)


length(intersect(rownames(md),colnames(scc)))
identical(rownames(md),colnames(scc))

sc <- CreateSeuratObject(counts=scc,meta.data=md)

sc <- NormalizeData(sc, verbose=TRUE)

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)

sc <- ScaleData(sc, verbose = TRUE)

sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 100, verbose = TRUE)


#loc<-names(sc$Location)[which(is.na(sc$Location))]

tissues<-sc1@meta.data$celltype  ## Change tissues by celltype maybe in other places. 
names(tissues)<-rownames(sc1@meta.data)
names(tissues)<-paste0(names(tissues),"_1")
sc$Location[which( colnames(sc) %in% loc)]<- tissues[loc ] #sc1$celltype [ colnames(sc)[which(colnames(sc) %in% loc)]]

sc <- RunHarmony(sc,c("Library"),reduction.use = "pca")



sc <- RunUMAP(sc,reduction = "harmony", dims = 1:30)

###### Cluster

sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30, verbose = TRUE)

sc <- FindClusters(sc, verbose = TRUE, resolution=0.2)





he <- t(sc@reductions$harmony@cell.embeddings[,1:30])

query.he <- he[,is.na(sc@meta.data$celltype)]
ref.he <- he[,!is.na(sc@meta.data$celltype)]
ref.labels <- sc@meta.data$celltype[!is.na(sc@meta.data$celltype)]


pred.labels <- SingleR(test = query.he, ref = ref.he, labels = ref.labels)

##table(pred.labels)

table(pred.labels$pruned.labels)

sum(is.na(pred.labels$pruned.labels))

fname=paste0(outFolder,"sc_harmony_cellClass_Murine2023.res0.2.rds")
write_rds(pred.labels,fname)

md <- pred.labels %>% as.data.frame() %>% 
  rownames_to_column("BARCODES") %>%
  left_join(sc2@meta.data %>% rownames_to_column("BARCODES"))


fname=paste0(outFolder,"3_harmony_cellClass_Murine2023.res0.2.csv")
write_csv(md,fname)

## save object.
fname=paste0("4_FullIntegrated.refmurine2023.Harmony.res0.2.rds")
write_rds(sc,fname)


pred.labels<-read_rds(paste0("4_FullIntegrated.refmurine2023.Harmony.res0.2.rds"))



