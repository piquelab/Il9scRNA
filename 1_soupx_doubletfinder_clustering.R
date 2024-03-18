######################################################################
# 1. preprocessing and clustering steps
######################################################################
library(Seurat)
library(dplyr)
library(Matrix)
library(harmony)
library(future)
library(readr)
library(DoubletFinder)
library(SoupX)
library(diem)
library(tidyverse)


outFolder="./1_soupx_doubletfinder_chrM/"
system(paste0("mkdir -p ", outFolder))

future::plan(strategy = 'multicore', workers = 12)
options(future.globals.maxSize = 30 * 1024 ^ 3)


##samples<-read_tsv("Bacteria_inducedPTLproject-Details.txt")


#Load in the required files (barcodes.tsv, genes.tsv, and matrix.mtx), from the raw data files (raw_feature_bc_matrix).

cellrangerFolder="../demux_cellranger_IL9mice/"

folders <- list.files(cellrangerFolder,pattern="^M_*")

samples <- folders %>% as_tibble() %>% select(Sample=value) %>% 
  separate(into=c(NA,"SampleindexID","Tissue","Condition"),
            col=Sample,remove=FALSE)

sc_list<-sapply(folders, function(x){
  #x<-forders[i]
  print(x)
  #mouse.data <- Read10X(data.dir = paste0("../counts_BIPTL_2021-05-01/",x,"/outs/filtered_feature_bc_matrix"))
  mouse.data<-load10X(paste0(cellrangerFolder,x,"/"))

  
  ##############################################################################
  # SoupX automatic way
  ##############################################################################
  
  sc = autoEstCont(mouse.data)
  out = adjustCounts(sc)
##  mouse.data<-out
  sc4 <- out
  #################################################################################
  # creating seurat object
  #################################################################################
  
  sc4 <- CreateSeuratObject(counts = mouse.data$toc, project = "Mouse",min.cells = 0, min.features=200)
  sc4@meta.data$Library<-rep(x,nrow(sc4@meta.data))
  sc4@meta.data$Location<-rep(samples %>% filter(Sample==x) %>%select(Tissue)%>% unlist, nrow(sc4@meta.data))
  sc4@meta.data$Condition<-rep(samples %>% filter(Sample==x) %>%select(Condition)%>% unlist, nrow(sc4@meta.data))
  sc4@meta.data$SampleindexID<-rep(samples %>% filter(Sample==x) %>%select(SampleindexID)%>% unlist, nrow(sc4@meta.data))
  
  sc4 <- NormalizeData(sc4, verbose=TRUE)
  # 
  sc4 <- FindVariableFeatures(sc4, selection.method = "vst", nfeatures = 3000)
  # 
  sc4 <- ScaleData(sc4, verbose = TRUE)
  sc4 <- RunPCA(sc4)
  # sc4 <- RunUMAP(sc4, dims = 1:10)

  #################################################################################
  # doublet removal
  #################################################################################
  
  nExp_poi <- round(0.075*nrow(sc4@meta.data))  ## Assuming 7.5% doublet formation rate
  
  #as maxima in BCmvn distributions
  sc_dbf <- doubletFinder(sc4, PCs = 1:75, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE)
  
  DF.name = colnames(sc_dbf@meta.data)[grepl("DF.classification", colnames(sc_dbf@meta.data))]
  sc_dbf = sc_dbf[, sc_dbf@meta.data[, DF.name] == "Singlet"]
  sc_dbf
##sc4 
  } )


sc <- merge(sc_list[[1]],sc_list[-1], 
            add.cell.ids= samples$Sample,
            project="mouse-IL9")


sc[["RNA"]] <-JoinLayers(sc[["RNA"]])


opfn <- paste0(outFolder,"seuratObj-SoupX-Doublet-merge.",Sys.Date(),".rds") 
write_rds(sc, opfn)


######################################
### FINISH and continue.  #####

