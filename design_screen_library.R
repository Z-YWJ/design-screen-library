
library(scales)
library(UpSetR)
library(ggplot2)
library(reshape2)  
library(patchwork)
library(RColorBrewer)
suppressMessages(library(dplyr)) 
suppressMessages(library(tidyr))
suppressMessages(library(enrichplot))
library(org.Hs.eg.db)
suppressMessages(library(clusterProfiler))
options(warn = 1)  # 显示警告信息设置为1   不显示则设置为-1
setwd("/home/zhaoyp/design-screen-library/")
# options(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS=TRUE)   让Bioconductor通过联网去验证版本
# options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) 
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# options("download.file.method"="libcurl")             配置R联网环境
# options("url.method"="libcurl")

########################################## 参数设置区域 ######################################### 

GeneList <- NULL
StepAnnotationOfGeneList <- NULL

# 选择细胞类型 可供选择的有  1. CD4 T-cell     2. macrophage     3. CD8 T-cell    4. naive B-cell    5. NK-cell  6. myeloid DC
CellType <- "myeloid DC"

# 根据表达量高低来选取基因数量
GeneNumberByExpreesion <- 1500

# 根据 nTPM 的大小来筛去一些基因  nTPM > nTPMcutoff 的才会保留
nTPMcutoff <- 1

# 选取肿瘤组织中特定细胞类型的基因
IncludingGeneSpecificInCancer_Celltype <- TRUE; CancerCellType <- c( "Dendritic cells"); Least_sample_number <- 10
# 可供选择的细胞类型 "Dendritic cells",  "CD4+ effector T cells", "CD4+ naive T cells", "NK cells", "Macrophages", "Monocytes"
#                "CD8+ effector T cells", "CD8+ naive T cells", "Activated B cells", "Naive B cells", "Plasma cells", "Neutrophils"(如果用这个 下面要改), 

# 是否包含转录因子  1639个
IncludingTFs <- TRUE

# 是否包含具有特定免疫功能的基因家族 共17种
functional_family <- c("Antigen_Processing_and_Presentation", "BCRSignalingPathway", "Chemokines", "Cytokines", "Interferons", "Interleukins_Receptor", "TCRsignalingPathway", 
                       "TGFb_Family_Member_Receptor", "TNF_Family_Members_Receptors","Antimicrobials", "Chemokine_Receptors", "Cytokine_Receptors", "Interferon_Receptor", 
                       "Interleukins", "NaturalKiller_Cell_Cytotoxicity", "TGFb_Family_Member", "TNF_Family_Members")
Including_functional_family <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)

### 对于CD8 或 CD4细胞   是否包含 Activated naive cd4/8 cells 相比于 naive cd4/8 cells 上升的基因，如若包含的话，包含19380基因中的前百分之几
ActivateVsNaiveTcellGene <- FALSE; AVNTcellGenePercent <- 0.02   # 共19380genes

### 对于macrophage   是否包含 normal tissue(主要是Lung)中的 High-staining 的蛋白基因
Including_highexpr_genes_in_LungMacrophage <- FALSE

# 对于macrophage   是否包含肝中 macrophage 的高表达基因
Including_highexpr_genes_in_liver_macrophage <- FALSE; livermacroPercent <- 0.15

# 对于macrophage   是否包含 Kupffercells 中的 High-expr 基因  ( 单细胞数据 )
Including_highexpr_genes_in_Kupffercells  <- FALSE;  KupfferExprGenePercent <- 0.15

### 对于NK-cells  是否包含单细胞数据中的high-expr基因
Including_highexpr_genes_in_NKcells <- FALSE; NKcellsExprGenePercent <- 0.01

### 对于B-cells 是否包含单细胞数据中的high-expr基因
Including_highexpr_genes_in_Bcells <- FALSE; BcellsExprGenePercent <- 0.01

### 对于DC 是否包含单细胞数据中的high-expr基因
Including_highexpr_genes_in_DC <- TRUE; DCExprGenePercent <- 0.15

##############################################################################################
# 只传入ensembl编号  注释基因加入的步骤
step_ann <- function(OldAnndf, NewAnndf, step_description){
    for (row in NewAnndf){
        if (!row %in% OldAnndf[,"ENSEMBL"]){
            row <- c(row,step_description)
            OldAnndf <- rbind(OldAnndf, row)
        }else{
            idx <- which(OldAnndf[,"ENSEMBL"] == row )
            OldAnndf[idx,"Step"] <- paste0(OldAnndf[idx,"Step"],",",step_description)
        }
    }
    return(OldAnndf)
}



#####################################    读取文件    ###########################################

hpa_files_path <- "/home/zhaoyp/design-screen-library/HPA"

### DB - Human Protein Atlas
normal_tissue <- read.table(paste0(hpa_files_path,"/1.normal_tissue.tsv"), header=T, sep="\t", check.name=F,quote="", comment.char="")               # 15318 genes, 64 tissue, 145 cell type,  1194479行
# macrophages      
# Tissue        macrophages 
# endometrium         120
# lung              13441    
# pathology <- read.table(paste0(hpa_files_path,"/2.pathology.tsv"), header=T,sep="\t",check.name=F)                                # 20090 genes, 20 cancers  考虑基因的致癌性
rna_sg_celltype <- read.table(paste0(hpa_files_path,"/8.rna_single_cell_type.tsv"), header=T,sep="\t",check.name=F)                 # 20090 genes, 79 cell types from 30 datasets  (nTPM) 
# B-cells   dendritic cells    Kupffer cells    monocytes     NK-cells    T-cells  
rna_sg_celltype_tissue <- read.table(paste0(hpa_files_path,"/9.rna_single_cell_type_tissue.tsv"), header=T,sep="\t",check.name=F)

HPA_immune_cell_gene <- read.table(paste0(hpa_files_path,"/18.rna_immune_cell.tsv"), header=T,sep="\t",check.name=F,comment.char="")                # 20090 genes, 19 immune cell types   重要
rna_immune_cell_monaco <- read.table(paste0(hpa_files_path,"/20.rna_immune_cell_monaco.tsv"), header=T,sep="\t",check.name=F,comment.char="")       # 20090 genes, 30 immune cell types 
rna_immune_cell_schmiedel <- read.table(paste0(hpa_files_path,"/21.rna_immune_cell_schmiedel.tsv"), header=T,sep="\t",
                                               check.name=F, comment.char="", stringsAsFactors=FALSE) # 19380 genes, 15 immune cell types
protein_atlas <- read.table(paste0(hpa_files_path,"/26.proteinatlas.tsv"), header=T,sep="\t",check.name=F,fill=T,
                                               quote="", comment.char="", stringsAsFactors=FALSE)          # 20090 genes annotation

### DB - cancerSCEM
celltype_vec <- c("Dendritic cells", "CD4+ effector T cells", "CD4+ naive T cells", "NK cells", "Macrophages", "Monocytes", "CD8+ effector T cells", 
                 "CD8+ naive T cells", "Activated B cells", "Naive B cells", "Neutrophils", "Plasma cells")
celltype_list <- list( `Dendritic cells` = c(), `CD4+ effector T cells` = c(), `CD4+ naive T cells` = c(), `NK cells`= c(), `Macrophages` = c(), `Monocytes` =c(),
       `CD8+ effector T cells` = c(), `CD8+ naive T cells` = c(), `Activated B cells`=c(), `Naive B cells` = c(), `Neutrophils` = c(), `Plasma cells`=c())
cancer_celltype_list <- list(`AA`=celltype_list, `AML`=celltype_list, `ATC`=celltype_list, `BCC`=celltype_list,`CRC`=celltype_list, `DCIS` =celltype_list, `GBM`=celltype_list,
                             `HCC`=celltype_list, `LUAD`=celltype_list, `LUSC`=celltype_list, `MCC`=celltype_list,`MIUBC`=celltype_list, `MPAL`=celltype_list, `NBL` =celltype_list,
                             `NSCLC`=celltype_list, `OV`=celltype_list, `PDAC`=celltype_list,`STAD`=celltype_list,`TNBC`=celltype_list,`UCEC`=celltype_list)
SCEM_files <- list.files("/home/zhaoyp/design-screen-library/cancerSCEM/", pattern="all.celltypes.tsv", full.names=T)
for (file in SCEM_files) {
    sample <- strsplit(file, "/")[[1]][length(strsplit(file, "/")[[1]])];  
    cancer <- strsplit(sample, "-")[[1]][1]
    tmp_cancer_file <- read.table(file, header=T, sep="\t", check.names=F, comment.char="", stringsAsFactors=FALSE)[,c(6,7)]
    tmp_cancer_file <- subset(tmp_cancer_file, cluster %in% celltype_vec)
    if (length(tmp_cancer_file[,1]) == 0) {next}
    for (row in 1:nrow(tmp_cancer_file)) {
        celltype <- tmp_cancer_file[row,1];  gene <- tmp_cancer_file[row,2]
        cancer_celltype_list[[cancer]][[celltype]] <- c(cancer_celltype_list[[cancer]][[celltype]], gene)
    }
}

### DB - Immport
ImmPort_genelist <- read.table("/home/zhaoyp/design-screen-library/Immport/GeneList.txt", header=T, sep="\t", check.names=F,comment.char="", stringsAsFactors=FALSE)

### DB - The Human Transcription Factors 2834
TF_ann <- read.csv("/home/zhaoyp/design-screen-library/HTF_DatabaseExtract.csv")

## essential genes list from addgene  共 2449 个基因
essential_genelist <- read.csv("/home/zhaoyp/design-screen-library/essential_genes_list_from_addgene.csv")
tmp_essential_genelist <- AnnotationDbi::select(org.Hs.eg.db, keys=essential_genelist$gene, columns="ENSEMBL", keytype="SYMBOL", multiVals="first")
tmp_essential_genelist <- tmp_essential_genelist[!duplicated(tmp_essential_genelist$SYMBOL), ]
na_tmp_essential_genelist <- tmp_essential_genelist[is.na(tmp_essential_genelist$ENSEMBL),"SYMBOL"]
na_tmp_essential_genelist <- AnnotationDbi::select(org.Hs.eg.db, keys=na_tmp_essential_genelist, columns="ENSEMBL", keytype="ALIAS", multiVals="first")
na_tmp_essential_genelist <- na_tmp_essential_genelist[!duplicated(na_tmp_essential_genelist$ALIAS), ]; 
colnames(tmp_essential_genelist) <- c("gene","ENSEMBL"); colnames(na_tmp_essential_genelist) <- c("gene","ENSEMBL");
essential_genelist <- rbind(na.omit(tmp_essential_genelist), na_tmp_essential_genelist)

### DB - HouseKeeping Transcript Atlas
HKGs <- read.csv("/home/zhaoyp/design-screen-library/Housekeeping_GenesHuman.csv",sep=";")            # 转录组水平 2833     Gene水平只有 2176 
HKGs <- as.data.frame(lapply(HKGs, as.character), stringsAsFactors = FALSE)
# HKGsMostStableTwenty <- read.csv("/home/zhaoyp/design-screen-library/MostStable.csv",sep=";")
# 增加 ENSG 的信息 ;     此处不运行 之后用symbol去除
# symbol2ensembl <- AnnotationDbi::select(org.Hs.eg.db, keys=HKGs[,"Gene.name"], columns="ENSEMBL", keytype="SYMBOL")  # many to many
# symbol2ensembl <- symbol2ensembl[!duplicated(symbol2ensembl[,1]) ,]  #
# HKGs <- merge(HKGs, symbol2ensembl, by.x="Gene.name", by.y="SYMBOL")

### addgenes - https://www.addgene.org/pooled-library/bassik-human-crispr-knockout/   # 20533     3
Gene_func_addgene <- read.table("/home/zhaoyp/design-screen-library/gene_function_annotation_from_addgene.txt", sep="_", check.names=F, stringsAsFactors=FALSE)[,c(1,2,3)]
Gene_func_addgene <- Gene_func_addgene[!duplicated(Gene_func_addgene[,1]),]; Gene_func_addgene <- Gene_func_addgene[-c(1,2),]; colnames(Gene_func_addgene) <- c("ENSEMBL","gene","function")



### 按照细胞类型特异的基因表达量排序 ###
if ( CellType=="CD8 T-cell" ){
    HPA_immune <- c("naive CD8 T-cell"); monaco_immune <- c("naive CD8 T-cell"); schmiedel_immune <- c("naive CD8 T-cell","Naive CD8 T-cell activated")
}else if(CellType=="CD4 T-cell") {
    HPA_immune <- c("naive CD4 T-cell"); monaco_immune <- c("naive CD4 T-cell"); schmiedel_immune <- c("naive CD4 T-cell","Naive CD4 T-cell activated")
}else if(CellType=="macrophage"){  HPA_immune <- c("classical monocyte"); 
    # monaco_immune <- c("classical monocate"); schmiedel_immune <- c("classical monocate") 应该用不到
}else if(CellType=="naive B-cell") { HPA_immune <- c("naive B-cell"); 
}else if(CellType=="NK-cell")      { HPA_immune <- c("NK-cell"); 
}else if(CellType=="myeloid DC")   { HPA_immune <- c("myeloid DC")}

HPA_genes_subset_celltype <- HPA_immune_cell_gene[which(HPA_immune_cell_gene$`Immune cell` %in% HPA_immune  ),]     # 合理的写法 每种 cell type 都对应20090种基因
HPA_genes_subset_expre <- HPA_genes_subset_celltype[order(HPA_genes_subset_celltype$nTPM, decreasing=T),]   # 最初始的 gene list
HPA_genes_subset_expre <- HPA_genes_subset_expre[1:GeneNumberByExpreesion, c("Gene name", "Gene")]

# !!
writeLines('==================Selecting genes with top expression=================')
if(CellType=="macrophage"){writeLines("*** Annotation: for macrophage, we select the top-expr genes from classical monocyte ***") }
writeLines(paste("Before the selection, the length of the GeneList: ", 0))
GeneList <- HPA_genes_subset_expre
colnames(GeneList) <- c("Gene","ENSEMBL")
writeLines(paste("After  the selection, the length of the GeneList: ", nrow(GeneList),"\n"))

StepAnnotationOfGeneList <- data.frame(ENSEMBL=GeneList[,2], Step = "TopGeneInFirstStep" )  
# writeLines('==========================================================\n')



### 添加细胞特异的肿瘤高表达基因 ###
if (IncludingGeneSpecificInCancer_Celltype ){
    genes_dup_set <- list()
    for (cell in CancerCellType) {
        genes_dup_set[[cell]] <- c()
        for (index in 1:length(cancer_celltype_list)){
            genes_dup_set[[cell]] <- c(genes_dup_set[[cell]], cancer_celltype_list[[index]][[cell]])
        }
        genes_dup_set[[cell]] <- sort( table(genes_dup_set[[cell]] ) )
    }
    tmp_df <- data.frame(unlist(genes_dup_set))

    if (CellType != "macrophage"){
        tmp_df$cell_type <- factor(sapply(rownames(tmp_df), function(x) strsplit(x, "cells\\.")[[1]][1] ));
        tmp_df$genes <- sapply(rownames(tmp_df), function(x) strsplit(x, "cells\\.")[[1]][2] ) 
    }else{
        tmp_df$cell_type <- factor(sapply(rownames(tmp_df), function(x) strsplit(x, "s\\.")[[1]][1] ));
        tmp_df$genes <- sapply(rownames(tmp_df), function(x) strsplit(x, "s\\.")[[1]][2]) }
    colnames(tmp_df) <- c("Number", "CellType", "GeneName")
    
    genesname_filterd_list <- list()
    label <- "\n"
    for (cell in CancerCellType){
        genesname_filterd_list[[cell]] <- names(which(genes_dup_set[[cell]]>Least_sample_number))
        label <- paste0(label, cell, " - ", "all DEGs: ", length(genes_dup_set[[cell]]), "; filtered DEGs: ", sum(genes_dup_set[[cell]]>Least_sample_number ),"\n" )
    } 
    #genesname_intersect <- Reduce(intersect, genesname_filterd_list)
    genesname_union <- Reduce(union, genesname_filterd_list)
    label <- paste0(label, "\n", "union DEGs: ", length(genesname_union))

    ggplot(tmp_df, aes(x=CellType, y=Number, color=CellType))+
      geom_boxplot() + geom_jitter(position=position_jitter(0.2)) +
      theme_classic() + geom_hline(aes( yintercept = Least_sample_number), color="blue",size=1) +
      annotate("text", x=1, y=50,size=3, label=label , color="blue") + ggtitle(CellType," DEGsforCelltypesInCancer")
    ggsave(paste0("/home/zhaoyp/design-screen-library/library_design_output/",CellType,"DEGsforCelltypesInCancer.png"), width = 8, height = 8, dpi = 300 )

    # !!
    tmpGeneList <- AnnotationDbi::select(org.Hs.eg.db, keys=genesname_union, columns="ENSEMBL", keytype="SYMBOL", multiVals="first")
    tmpGeneList <- tmpGeneList[!duplicated(tmpGeneList$SYMBOL), ]  
    writeLines('===========Selecting high-expr genes in specific cell type of pan-cancer ===========')
    writeLines("*** Annotation: here one photo is output which can help with the selection of cutoff ***")
    writeLines(paste("The number of Union DEGs: ", nrow(tmpGeneList)))
    writeLines(paste("Before the selection, the length of the GeneList: ", nrow(GeneList)))
    colnames(tmpGeneList) <- c("Gene","ENSEMBL")
    GeneList <- rbind(GeneList, tmpGeneList)
    writeLines(paste("After  the selection, the length of the GeneList: ", nrow(GeneList),"\n")) 
    #writeLines('=================================================================================')

    StepAnnotationOfGeneList <- step_ann(StepAnnotationOfGeneList, tmpGeneList[,"ENSEMBL"],"CancerDExpr" )
}


### 添加转录因子 ###
if (IncludingTFs){
    TF_ann_subset <- subset(TF_ann, `Is.TF.`  == "Yes" , select=c("Ensembl.ID", "HGNC.symbol", "EntrezGene.ID"))   # 1639  3
    
    # !!
    writeLines('===========Selecting genes encoding transcription factors===========')
    writeLines(paste("The number of TFs: ", nrow(TF_ann_subset)))
    writeLines(paste("Before the selection, the length of the GeneList: ", nrow(GeneList)))
    tmpGeneList <- TF_ann_subset[, c("HGNC.symbol", "Ensembl.ID")]
    colnames(tmpGeneList) <- c("Gene","ENSEMBL")
    GeneList <- rbind(GeneList, tmpGeneList)
    writeLines(paste("After  the selection, the length of the GeneList: ", nrow(GeneList),"\n"))
    #writeLines('===================================================================')

    StepAnnotationOfGeneList <- step_ann(StepAnnotationOfGeneList, tmpGeneList[,"ENSEMBL"],"TF" )

}


### optional 增加 Immport 注释的功能基因 ###
tmp_Ensembl <- AnnotationDbi::select(org.Hs.eg.db, keys=as.character(ImmPort_genelist$ID), columns="ENSEMBL", keytype="ENTREZID", multiVals="first") 
tmp_Ensembl <- tmp_Ensembl[!duplicated(tmp_Ensembl$ENTREZID), ]         # 有44个基因找不到ensembl id，此处就先直接省去吧
ImmPort_genelist_ensembl  <- merge(ImmPort_genelist, tmp_Ensembl, by.x="ID", by.y="ENTREZID",all.x=TRUE)   
ImmPort_genelist_ensembl <- na.omit(ImmPort_genelist_ensembl); ImmPort_genelist_ensembl <- ImmPort_genelist_ensembl[,c("ID","Symbol","Category","ENSEMBL")]


# !!
writeLines('===================Selecting genes by its function====================')
writeLines(paste("Before the selection, the length of the GeneList: ", nrow(GeneList)))
tmpGeneList <- NULL
for (i in 1:length(Including_functional_family)){
    if (Including_functional_family[i]){ 
        tmpGeneList1 <- subset(ImmPort_genelist_ensembl, Category=Including_functional_family[i], select=c("Symbol","ENSEMBL"))
        colnames(tmpGeneList1) <- c("Gene","ENSEMBL")
        tmpGeneList <- rbind(tmpGeneList, tmpGeneList1)
}}
tmpGeneList <- tmpGeneList[!duplicated(tmpGeneList$Gene),]
GeneList <- rbind(GeneList, tmpGeneList) 
writeLines(paste("After  the selection, the length of the GeneList: ", nrow(GeneList),"\n"))

StepAnnotationOfGeneList <- step_ann(StepAnnotationOfGeneList, tmpGeneList[,"ENSEMBL"],"FunctionGene" )
#writeLines('==========================================================')



### optional 增加 naive T-cell 中的上调基因 ###
if (ActivateVsNaiveTcellGene && HPA_immune %in% c("naive CD8 T-cell", "naive CD4 T-cell")){
    ActNaiveData <- subset(rna_immune_cell_schmiedel, `Immune cell` %in% schmiedel_immune)

    ggplot(ActNaiveData, aes(x=TPM, fill=`Immune cell`)) + theme_classic() +
      geom_histogram(data=ActNaiveData[which(ActNaiveData$`Immune cell` %in% c("naive CD8 T-cell","naive CD4 T-cell")),], aes(y=after_stat(count * (-1)))) + 
      geom_histogram(data=ActNaiveData[which(ActNaiveData$`Immune cell` %in% c("Naive CD8 T-cell activated","Naive CD4 T-cell activated")),], aes(y=after_stat(count))) +     
      scale_x_continuous(expand=c(0,0), trans = log10_trans() )  +  xlab("log10(nTPM)") +
      scale_y_continuous(expand=c(0,0), limits = c(-1500,1500), breaks=seq(-1500,1500,100), labels=abs(seq(-1500,1500,100)))   
    ggsave(paste0("/home/zhaoyp/design-screen-library/library_design_output/",CellType,"ActVsNaiveTcellTPM.png"), width = 6, height = 8, dpi = 300 )
 
    ActNaiveDataWide <- spread(ActNaiveData, key=`Immune cell`, value=TPM)
    if(HPA_immune == "naive CD8 T-cell"){
    ActNaiveDataWide$Dexp <- (ActNaiveDataWide$`Naive CD8 T-cell activated` - ActNaiveDataWide$`naive CD8 T-cell`)/(ActNaiveDataWide$`naive CD8 T-cell`+1); ann_xy <-c(100,600)}
    if(HPA_immune == "naive CD4 T-cell"){
    ActNaiveDataWide$Dexp <- (ActNaiveDataWide$`Naive CD4 T-cell activated` - ActNaiveDataWide$`naive CD4 T-cell`)/(ActNaiveDataWide$`naive CD4 T-cell`+1); ann_xy <-c(600,600)}
    ActNaiveDataWide <- ActNaiveDataWide[order(ActNaiveDataWide$Dexp, decreasing=T),]

    cutoff <- min(ActNaiveDataWide[1:(nrow(ActNaiveDataWide)*AVNTcellGenePercent),"Dexp"])
    ggplot(ActNaiveDataWide, aes(x=Dexp)) + theme_classic() + 
      geom_histogram(fill="grey") + scale_x_continuous(expand=c(0,0), trans = log10_trans() )  + xlab("log10(Dexp_TPM)")  + ylab("count") +
      geom_vline(xintercept = cutoff, linewidth = 1, colour = "#FF3721")  + ggtitle("Activation Vs Naive Tcell Dexpression TPM ") +
      annotate("text", x=ann_xy[1], y=ann_xy[2], label=paste0("Percent: ",AVNTcellGenePercent*100," %\n Number: ",round(AVNTcellGenePercent * nrow(ActNaiveDataWide),0), 
        "\nLowest DexpTPM: ",round(cutoff,2)))
    ggsave(paste0("/home/zhaoyp/design-screen-library/library_design_output/",CellType,"ActVsNaiveTcellTPMExp.png"), width = 6, height = 8, dpi = 300 )

    # !!
    writeLines('==========Selecting genes with differential expression in T cell=========')
    writeLines("*** Annotation: here two photos are output which can help with the selection of cutoff ***")
    writeLines(paste("Before the selection, the length of the GeneList: ", nrow(GeneList)))
    tmpGeneList <- ActNaiveDataWide[1:round(nrow(ActNaiveDataWide)*AVNTcellGenePercent, 0), c("Gene name", "Gene")]
    colnames(tmpGeneList) <- c("Gene","ENSEMBL")
    GeneList <- rbind(GeneList, tmpGeneList)
    writeLines(paste("After  the selection, the length of the GeneList: ", nrow(GeneList),"\n"))
    #writeLines('===================================================================')

    StepAnnotationOfGeneList <- step_ann(StepAnnotationOfGeneList, tmpGeneList[,"ENSEMBL"],"TcellActVsNaive" )
}

########################################## 巨噬细胞 ######################################### 

### optional  增加normal tissue(主要是Lung)中的 High-staining 的蛋白基因  ###
if(Including_highexpr_genes_in_LungMacrophage && CellType=="macrophage"){
    macrophage_in_lung_tissue <- subset(normal_tissue, `Cell type` == "macrophages" & `Tissue` == "lung")  # 13441
    macrophage_in_lung_tissue <- subset(macrophage_in_lung_tissue, `Level` == "High" & `Reliability` %in% c("Enhanced","Supported")) # 前两个置信程度是693 前三个置信程度是1370

    writeLines('=========Selecting genes with high-expr in lung-macrophages========')
    writeLines(paste0("According to HPA database, the number of High-credibility-high-expr-genes in lung: ",nrow(macrophage_in_lung_tissue)))
    writeLines(paste("Before the selection, the length of the GeneList: ", nrow(GeneList)))
    tmpGeneList <- macrophage_in_lung_tissue[,c("Gene name", "Gene")]
    colnames(tmpGeneList) <- c("Gene","ENSEMBL")
    GeneList <- rbind(GeneList, tmpGeneList)
    writeLines(paste("After  the selection, the length of the GeneList: ", nrow(GeneList),"\n"))
    #writeLines('===================================================================')
    StepAnnotationOfGeneList <- step_ann(StepAnnotationOfGeneList, tmpGeneList[,"ENSEMBL"],"LungMacrophageProtein" )
}


### optional  增加肝中巨噬细胞的高表达基因   ###
if(Including_highexpr_genes_in_liver_macrophage && CellType=="macrophage"){
    Expr_in_liverMacro <- subset(rna_sg_celltype_tissue, `Cell type` == "macrophages" & `Tissue` == "liver" ) 
    # 有两类巨噬细胞 可以取两者的均值
    Expr_in_liverMacro_tmp <- data.frame(Expr_in_liverMacro %>% group_by(`Gene`)%>% summarise(Mean_nTPM=mean(nTPM)))
    Expr_in_liverMacro <- merge(Expr_in_liverMacro_tmp, Expr_in_liverMacro, by="Gene", all.x=TRUE)[,c("Gene","Mean_nTPM","Gene name")]
    Expr_in_liverMacro <- Expr_in_liverMacro[order(Expr_in_liverMacro$Mean_nTPM, decreasing=T), ]
    Expr_in_liverMacro <- Expr_in_liverMacro[!duplicated(Expr_in_liverMacro$Gene),]

    cutoff <- min(Expr_in_liverMacro[1:(livermacroPercent* 20090),"Mean_nTPM"])

    ggplot(Expr_in_liverMacro, aes(x=Mean_nTPM)) + theme_classic() +
      geom_histogram(fill="cadetblue2") + scale_x_continuous(expand=c(0,0), trans = log10_trans() )  + ylab("count") +
      geom_vline(xintercept = cutoff, linewidth = 1, colour = "#FF3721")  + xlab("log10(mean nTPM)") + ggtitle("liver_macrophage_mean_nTPM - two clusters") +
      annotate("text", x=13000, y=1500, label=paste0("Percent: ",livermacroPercent*100," %\n Number: ",round(livermacroPercent * 20090,0), "\nLowest nTPM: ",cutoff))
    ggsave("/home/zhaoyp/design-screen-library/library_design_output/macrophage_liver_mean_nTPM.png", width = 6, height = 8, dpi = 300 )
    
    writeLines('===========Selecting genes with high-expr in liver macrophage==========')
    writeLines("*** Annotation: here one photos are output which can help the selection of cutoff ***")
    writeLines(paste("Before the selection, the length of the GeneList: ", nrow(GeneList)))
    tmpGeneList <- Expr_in_liverMacro[1:(livermacroPercent * 20090), c("Gene name","Gene")]
    colnames(tmpGeneList) <- c("Gene","ENSEMBL")
    GeneList <- rbind(GeneList, tmpGeneList)
    writeLines(paste("After  the selection, the length of the GeneList: ", nrow(GeneList),"\n"))
    #writeLines('===================================================================')

    StepAnnotationOfGeneList <- step_ann(StepAnnotationOfGeneList, tmpGeneList[,"ENSEMBL"],"LiverMacrophage" )
}



### optional  增加 Kupffercells 中的 High-expr 基因 ###
if(Including_highexpr_genes_in_Kupffercells && CellType=="macrophage"){
    Expr_in_KupfferCells <- subset(rna_sg_celltype, `Cell type` == "Kupffer cells" )  # 20090
    Expr_in_KupfferCells <- Expr_in_KupfferCells[order(Expr_in_KupfferCells$nTPM, decreasing=T),]
    
    cutoff <- min(Expr_in_KupfferCells[1:(KupfferExprGenePercent * 20090),"nTPM"])

    ggplot(Expr_in_KupfferCells, aes(x=nTPM)) + theme_classic() +
      geom_histogram(fill="cadetblue2") + scale_x_continuous(expand=c(0,0), trans = log10_trans() )  +
      geom_vline(xintercept = cutoff, linewidth = 1, colour = "#FF3721")  + xlab("log10(nTPM)") + ggtitle("KupfferCellsExpr_nTPM") +
      annotate("text", x=8000, y=1000, label=paste0("Percent: ",KupfferExprGenePercent*100," %\n Number: ",round(KupfferExprGenePercent * 20090,0), "\nLowest nTPM: ",cutoff))
    ggsave("/home/zhaoyp/design-screen-library/library_design_output/KupfferCellsExprnTPM.png", width = 6, height = 8, dpi = 300 )
    
    writeLines('===========Selecting genes with high-expr in Kupffer Cells==========')
    writeLines("*** Annotation: here one photos are output which can help the selection of cutoff ***")
    writeLines(paste("Before the selection, the length of the GeneList: ", nrow(GeneList)))
    tmpGeneList <- Expr_in_KupfferCells[1:(KupfferExprGenePercent * 20090), c("Gene name","Gene")]
    colnames(tmpGeneList) <- c("Gene","ENSEMBL")
    GeneList <- rbind(GeneList, tmpGeneList)
    writeLines(paste("After  the selection, the length of the GeneList: ", nrow(GeneList),"\n"))
    #writeLines('===================================================================')

    StepAnnotationOfGeneList <- step_ann(StepAnnotationOfGeneList, tmpGeneList[,"ENSEMBL"],"KupfferCells" )
}



### optional  增加 NK-cells 单细胞数据中 High-expr 基因 ###
if(Including_highexpr_genes_in_NKcells && CellType=="NK-cell"){
    Expr <- subset(rna_sg_celltype, `Cell type` == "NK-cells" )  # 20090
    Expr <- Expr[order(Expr$nTPM, decreasing=T),]
    
    cutoff <- min(Expr[1:(NKcellsExprGenePercent * 20090),"nTPM"])

    ggplot(Expr, aes(x=nTPM)) + theme_classic() +
      geom_histogram(fill="cadetblue2") + scale_x_continuous(expand=c(0,0), trans = log10_trans() )  +
      geom_vline(xintercept = cutoff, linewidth = 1, colour = "#FF3721")  + xlab("log10(nTPM)") + ggtitle("NK-CellsExpr_nTPM") +
      annotate("text", x=5000, y=1000, label=paste0("Percent: ",NKcellsExprGenePercent*100," %\n Number: ",round(NKcellsExprGenePercent * 20090,0), "\nLowest nTPM: ",cutoff))
    ggsave("/home/zhaoyp/design-screen-library/library_design_output/NK-CellsExprnTPM.png", width = 6, height = 8, dpi = 300 )
    
    writeLines('===========Selecting genes with high-expr in NK Cells==========')
    writeLines("*** Annotation: here one photo are output which can help the selection of cutoff ***")
    writeLines(paste("Before the selection, the length of the GeneList: ", nrow(GeneList)))
    tmpGeneList <- Expr[1:(NKcellsExprGenePercent * 20090), c("Gene name","Gene")]
    colnames(tmpGeneList) <- c("Gene","ENSEMBL")
    GeneList <- rbind(GeneList, tmpGeneList)
    writeLines(paste("After  the selection, the length of the GeneList: ", nrow(GeneList),"\n"))
    #writeLines('===================================================================')

    StepAnnotationOfGeneList <- step_ann(StepAnnotationOfGeneList, tmpGeneList[,"ENSEMBL"],"NKHighExpr" )
}



### 对于B细胞 是否包含单细胞数据中的high-expr基因
if(Including_highexpr_genes_in_Bcells && CellType=="naive B-cell"){
    Expr <- subset(rna_sg_celltype, `Cell type` == "B-cells" )  # 20090
    Expr <- Expr[order(Expr$nTPM, decreasing=T),]
    cutoff <- min(Expr[1:(BcellsExprGenePercent * 20090),"nTPM"])

    ggplot(Expr, aes(x=nTPM)) + theme_classic() +
      geom_histogram(fill="cadetblue2") + scale_x_continuous(expand=c(0,0), trans = log10_trans() )  +
      geom_vline(xintercept = cutoff, linewidth = 1, colour = "#FF3721")  + xlab("log10(nTPM)") + ggtitle("B-CellsExpr_nTPM") +
      annotate("text", x=5000, y=1000, label=paste0("Percent: ",BcellsExprGenePercent*100," %\n Number: ",round(BcellsExprGenePercent * 20090,0), "\nLowest nTPM: ",cutoff))
    ggsave("/home/zhaoyp/design-screen-library/library_design_output/naive-B-CellsExprnTPM.png", width = 6, height = 8, dpi = 300 )
    
    writeLines('===========Selecting genes with high-expr in B Cells==========')
    writeLines("*** Annotation: here one photo is output which can help the selection of cutoff ***")
    writeLines(paste("Before the selection, the length of the GeneList: ", nrow(GeneList)))
    tmpGeneList <- Expr[1:(BcellsExprGenePercent * 20090), c("Gene name","Gene")]
    colnames(tmpGeneList) <- c("Gene","ENSEMBL")
    GeneList <- rbind(GeneList, tmpGeneList)
    writeLines(paste("After  the selection, the length of the GeneList: ", nrow(GeneList),"\n"))
    #writeLines('===================================================================')
    StepAnnotationOfGeneList <- step_ann(StepAnnotationOfGeneList, tmpGeneList[,"ENSEMBL"],"BCellHighExpr" )
}



### 对于DC细胞 是否包含单细胞数据中的high-expr基因
if(Including_highexpr_genes_in_DC && CellType=="myeloid DC"){
    Expr <- subset(rna_sg_celltype, `Cell type` == "dendritic cells" )  # 20090
    Expr <- Expr[order(Expr$nTPM, decreasing=T),]
    cutoff <- min(Expr[1:(DCExprGenePercent * 20090),"nTPM"])

    ggplot(Expr, aes(x=nTPM)) + theme_classic() +
      geom_histogram(fill="cadetblue2") + scale_x_continuous(expand=c(0,0), trans = log10_trans() )  +
      geom_vline(xintercept = cutoff, linewidth = 1, colour = "#FF3721")  + xlab("log10(nTPM)") + ggtitle("DCExpr_nTPM") +
      annotate("text", x=1000, y=800, label=paste0("Percent: ",DCExprGenePercent*100," %\n Number: ",round(DCExprGenePercent * 20090,0), "\nLowest nTPM: ",cutoff))
    ggsave("/home/zhaoyp/design-screen-library/library_design_output/myeloid DCExprnTPM.png", width = 6, height = 8, dpi = 300 )
    
    writeLines('===========Selecting genes with high-expr in B Cells==========')
    writeLines("*** Annotation: here one photo is output which can help the selection of cutoff ***")
    writeLines(paste("Before the selection, the length of the GeneList: ", nrow(GeneList)))
    tmpGeneList <- Expr[1:(DCExprGenePercent * 20090), c("Gene name","Gene")]
    colnames(tmpGeneList) <- c("Gene","ENSEMBL")
    GeneList <- rbind(GeneList, tmpGeneList)
    writeLines(paste("After  the selection, the length of the GeneList: ", nrow(GeneList),"\n"))
    #writeLines('===================================================================')
    StepAnnotationOfGeneList <- step_ann(StepAnnotationOfGeneList, tmpGeneList[,"ENSEMBL"],"DCHighExpr" )
}


############################################################################################## 

writeLines('========================= De duplication ===========================')
GeneList <- GeneList[!duplicated(GeneList$ENSEMBL),]
writeLines(paste("After  the de-duplication, the length of the GeneList: ", nrow(GeneList),"\n"))

############################################################################################## 
#### 去除表达量低的基因 ####

if ( CellType=="CD8 T-cell" )     { rna_alias <- "naive CD8 T-cell"; sgrna_alias <- "T-cells"
}else if(CellType=="CD4 T-cell")  { rna_alias <- "naive CD4 T-cell"; sgrna_alias <- "T-cells"
}else if(CellType=="naive B-cell"){ rna_alias <- "naive B-cell"; sgrna_alias <- "B-cells"
}else if(CellType=="NK-cell")     { rna_alias <- "NK-cell"; sgrna_alias <- "NK-cells"
}else if(CellType=="macrophage")  { rna_alias <- "classical monocyte"; sgrna_alias <- "Kupffer cells"; sgrna_liver_alias <- "liver macrophages"
}else if(CellType=="myeloid DC")  { rna_alias <- "myeloid DC"; sgrna_alias <- "dendritic cells"}


subset_HPA_immune_cell_gene <- subset(HPA_immune_cell_gene, `Immune cell` == rna_alias)
subset_rna_sg_celltype <- subset(rna_sg_celltype, `Cell type` == sgrna_alias)
gene_set_1 <- subset_HPA_immune_cell_gene[subset_HPA_immune_cell_gene$nTPM > nTPMcutoff, c("Gene name","Gene")]
gene_set_2 <- subset_rna_sg_celltype[subset_rna_sg_celltype$nTPM > nTPMcutoff, c("Gene name","Gene")]

p1 <- ggplot(subset_HPA_immune_cell_gene, aes(x=nTPM)) + theme_classic() +
        geom_histogram() + scale_x_continuous(expand=c(0,0), trans = log10_trans() )  +
        geom_vline(xintercept = nTPMcutoff, linewidth = 1, colour = "red")  + xlab("log10(nTPM)") + ggtitle(paste0(rna_alias," bulkRNA-nTPM")) +      
        annotate("text", x=1000, y=800,size=2,color="red",
                label=paste0("All genes: ",nrow(subset_HPA_immune_cell_gene),"\nLowest nTPM: ",nTPMcutoff, "\nFilterd genes: ",sum(subset_HPA_immune_cell_gene$nTPM > nTPMcutoff)))
p2 <- ggplot(subset_rna_sg_celltype, aes(x=nTPM)) + theme_classic() +
        geom_histogram() + scale_x_continuous(expand=c(0,0), trans = log10_trans() )  +
        geom_vline(xintercept = nTPMcutoff, linewidth = 1, colour = "red")  + xlab("log10(nTPM)") + ggtitle(paste0(sgrna_alias," sgRNA-nTPM")) +      
        annotate("text", x=1000, y=800,size=2, color="red",
                 label=paste0("All genes: ",nrow(subset_rna_sg_celltype),"\nLowest nTPM: ",nTPMcutoff, "\nFilterd genes: ",sum(subset_rna_sg_celltype$nTPM > nTPMcutoff)))
if (CellType != "macrophage"){
    gene_set <- intersect(gene_set_1$Gene, gene_set_2$Gene)
    p1|p2
    ggsave(paste0("/home/zhaoyp/design-screen-library/library_design_output/",CellType,"ExprnTPMFilter.png"), width = 10, height = 6, dpi = 300 )   
}else{
    gene_set_3 <- Expr_in_liverMacro[Expr_in_liverMacro$Mean_nTPM > nTPMcutoff, c("Gene name","Gene")]
    gene_set <- Reduce(intersect, list(gene_set_1$Gene, gene_set_2$Gene, gene_set_3$Gene))
    p3 <- ggplot(Expr_in_liverMacro, aes(x=Mean_nTPM)) + theme_classic() +
        geom_histogram() + scale_x_continuous(expand=c(0,0), trans = log10_trans() )  +
        geom_vline(xintercept = nTPMcutoff, linewidth = 1, colour = "red")  + xlab("log10(nTPM)") + ggtitle(paste0(sgrna_liver_alias," Mean-nTPM")) +      
        annotate("text", x=1000, y=800,size=2, color="red",
        label=paste0("All genes: ",nrow(Expr_in_liverMacro),"\nLowest nTPM: ",nTPMcutoff, "\nNumber of genes: ",sum(Expr_in_liverMacro$Mean_nTPM > nTPMcutoff)))
    p1|p2|p3
    ggsave(paste0("/home/zhaoyp/design-screen-library/library_design_output/",CellType,"ExprnTPMFilter.png"), width = 15, height = 6, dpi = 300 )
}

writeLines('========================= De low expr genes ===========================')
writeLines("*** Annotation: here one photo is output which can help the selection of cutoff ***")
writeLines(paste("The length of the Intersect-non-low-expr-genes: ", length(gene_set)))
writeLines(paste("Before the selection, the length of the GeneList: ", nrow(GeneList)))
GeneList <- GeneList[which( GeneList$ENSEMBL %in% gene_set) , ]
writeLines(paste("After  the de-essential genes, the length of the GeneList: ", nrow(GeneList),"\n"))


############################################################################################## 

### 去除 GeneList 中的 essential genes ###
writeLines('========================= De essential genes ===========================')
writeLines(paste("The length of the essential genes list: ", nrow(essential_genelist )))
GeneList <- GeneList[which(! GeneList$ENSEMBL %in% essential_genelist$ENSEMBL) , ]
writeLines(paste("After  the de-essential genes, the length of the GeneList: ", nrow(GeneList),"\n"))

### 去除 GeneList 中的 housekeeping genes ###
### 基因列表编号注释  通过SYMBOL 和 ALIAS 去除genes ###
whole_gene_list_ann <- protein_atlas[,c(1,2,3,4,8,12,13)]
EnsemblID <- as.character(whole_gene_list_ann[,"Ensembl"])
ens2entrez <- AnnotationDbi::select(org.Hs.eg.db, keys=EnsemblID, columns="ENTREZID", keytype="ENSEMBL")  # 1 to many
ens2entrez <- ens2entrez[!duplicated(ens2entrez[,1]) ,]  #20092
whole_gene_list_ann <- merge(whole_gene_list_ann, ens2entrez, by.x="Ensembl", by.y="ENSEMBL")
whole_gene_list_ann <- as.data.frame(lapply(whole_gene_list_ann, as.character), stringsAsFactors = FALSE)
whole_gene_list_ann[,"X.Gene.synonym."] <- sapply(whole_gene_list_ann[,"X.Gene.synonym."], function(x) {x<-paste0(" ",gsub("\"","",x),",")} )

### 去除 housekeeping genes ###
#DB https://housekeeping.unicamp.br/?homePageGlobal 
# 下面这行 重新加载时 要再次运行
remove_index <- na.omit(sapply(HKGs[,"Gene.name"], function(x) grep(paste0(" ",x,","),whole_gene_list_ann[,"X.Gene.synonym."])[1]))  # 通过alias筛选 
GenesListRemoveHKGs <- whole_gene_list_ann[-remove_index ,]                                            # 剩余20002
GenesListRemoveHKGs <- GenesListRemoveHKGs[!(GenesListRemoveHKGs[,"Gene"] %in% HKGs[,"Gene.name"]), ]  #通过symbol筛选   剩余 17871

writeLines('========================= De housekeeping genes ===========================')
writeLines(paste("The length of the housekeeping genes list: ", nrow(HKGs)))
GeneList <- merge(GenesListRemoveHKGs, GeneList, by.x="Ensembl", by.y="ENSEMBL")
writeLines(paste("After  the de-housekeeping genes, the length of the GeneList: ", nrow(GeneList),"\n"))


### 最后输出并保存 ###
colnames(GeneList) <- c("Ensembl", "Gene", "Synonym", "Description", "Protein.class", "Evidence", "HPA.evidence", "ENTREZID", "Gene.name.HPA")
GeneList <- GeneList[order(GeneList$Ensembl),]
GeneList <- merge(GeneList, StepAnnotationOfGeneList, by.x="Ensembl",by.y="ENSEMBL")
colnames(GeneList)[ncol(GeneList)] <- paste0(CellType,"Step")
write.table(GeneList, file=paste0("/home/zhaoyp/design-screen-library/library_design_output/",CellType,"_GeneList.txt"), row.names=FALSE, quote=FALSE,sep="\t")


########  对输出的基因的功能进行作图  addgene  #########
function_of_GeneList <- GeneList[,c("Ensembl", "Gene")]
function_of_GeneList <- merge(function_of_GeneList, Gene_func_addgene, by.x="Ensembl", by.y="ENSEMBL")

df <- table(function_of_GeneList[,"function"]);  labs <- paste0(names(df),"\n", df," (",round(df/sum(df)*100,1),"%)")
legend.text <- data.frame(n1=names(df), n2=c("Apoptosis and cancer","Drug targets, kinases, phosphatases","Gene Expression","Membrane proteins",
                                      "Proteostasis","Trafficking, mitochodrial, motility","Unassigned1","Unassigned2","Unassigned3"))
pdf(paste0("/home/zhaoyp/design-screen-library/library_design_output/",CellType,"_Genes_function_annotation_from_addgene.pdf"),height=5,width=10)
pie(df, labels=labs, main=paste(nrow(GeneList), "Genes function annotation from addgene"), col=brewer.pal(9,"Set3"),cex=0.5)
legend("bottomleft", legend=paste0( legend.text$n1," - ", legend.text$n2 ), fill=brewer.pal(9,"Set3"), cex=0.5)
dev.off()


########  对输出的基因进行GO分析并作图  GO annotation  #########
# 进行Go分析
# go <- enrichGO(function_of_GeneList[,"Ensembl"], OrgDb = org.Hs.eg.db, ont = "all", pvalueCutoff = 0.5, qvalueCutoff = 0.1, readable = TRUE)
# png(paste0("/home/zhaoyp/design-screen-library/libraary_design_output",CellType,"_Genes_GO_analysis"), height=700, width=500)
# barplot(go, showCategory = 15, title = paste(nrow(GeneList), "Genes GO annotation"), split="ONTOLOGY" ) + facet_grid(ONTOLOGY~., scales="free")
# dev.off()


########  对输出的基因的功能进行作图  HPA  #########
function_of_GeneList <- GeneList[,c("Gene", "Protein.class")]

# 向量转为存在的矩阵 HPA
function_mat <- matrix(0, nrow=nrow(function_of_GeneList), ncol=25)
rownames(function_mat) <- function_of_GeneList$Gene
colnames(function_mat) <- c("Blood group antigen proteins","Cancer-related genes","Candidate cardiovascular disease genes", "CD markers", "Citric acid cycle related proteins", 
                            "Disease related genes","Enzymes", "FDA approved drug targets", "G-protein coupled receptors", "Human disease related genes", "Immunoglobulin genes",
                            "Metabolic proteins", "Nuclear receptors", "Plasma proteins", "Potential drug targets", "Predicted intracellular proteins","Predicted membrane proteins", 
                            "Predicted secreted proteins",  "RAS pathway related proteins", "Ribosomal proteins",
                            "RNA polymerase related proteins", "T-cell receptor genes", "Transcription factors", "Transporters", "Voltage-gated ion channels")
function_of_GeneList <- sapply(function_of_GeneList[,"Protein.class"], function(x) { x<-gsub("\"", "",x); strsplit(x, ", ")})

# 生成 0 - 1 矩阵
for ( i in 1:length(function_of_GeneList) ) {          # copilot 太强了吧
  for ( j in 1:length(function_of_GeneList[[i]]) ) {
    function_mat[i, function_of_GeneList[[i]][[j]]] <- 1
}}
function_mat <- data.frame(function_mat); function_mat$Gene <- rownames(function_mat)

pdf(paste0("/home/zhaoyp/design-screen-library/library_design_output/",CellType,"_Genes_function_annotation_from_HPA_0.pdf"),height=6,width=10)
upset(function_mat, nsets=50, sets.x.label="Protein class", mb.ratio=c(0.55,0.45), mainbar.y.label = paste0(CellType,"\n All Gene Number: ",nrow(GeneList)),
      main.bar.color = "steelblue", sets.bar.color = "steelblue")
dev.off()


pdf(paste0("/home/zhaoyp/design-screen-library/library_design_output/",CellType,"_Genes_function_annotation_from_HPA_1.pdf"),height=6,width=12)
upset(function_mat, nsets=50, sets.x.label="Protein class", mb.ratio=c(0.55,0.45), mainbar.y.label = paste0(CellType,"\n All Gene Number: ",nrow(GeneList)),
      main.bar.color = "steelblue", sets.bar.color = "steelblue", order.by=c("freq"), decreasing=c(TRUE),nintersects = 80 )
dev.off()

##############################################
writeLines('*** Annotation: Finally, three photos are output here to understand gene function ***')





