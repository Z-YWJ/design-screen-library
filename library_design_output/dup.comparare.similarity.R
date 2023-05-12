

library(UpSetR)
library(ggplot2)


setwd("/home/zhaoyp/design-screen-library/library_design_output/")
table_cd8 <- read.table("./CD8 T-cell_GeneList.txt",header=T,sep="\t",quote="")
table_cd4 <- read.table("./CD4 T-cell_GeneList.txt",header=T,sep="\t",quote="")
table_macro <- read.table("./macrophage_GeneList.txt",header=T,sep="\t",quote="")

df <- data.frame(table_cd8[,c("Ensembl","Gene")]) #1655
df <- rbind(df,table_macro[,c("Ensembl","Gene")]) #1980
df <- rbind(df,table_cd4[,c("Ensembl","Gene")]) #1980
dup_df <- df[duplicated(df[,"Ensembl"]),] 
dup_df <- dup_df[duplicated(dup_df[,"Ensembl"]),] #重复项1007

step_df <- table_cd8[,c(1,ncol(table_cd8))]; step_df <- merge(step_df, table_cd4[,c(1,ncol(table_cd4))]); step_df <- merge(step_df, table_macro[,c(1,ncol(table_macro))] )

dup_df <- merge(dup_df,table_macro[,-ncol(table_macro)])
dup_df <- merge(dup_df,step_df,by.x="Ensembl",by.y="Ensembl")
write.table(dup_df, file="./dup_cd8_cd4_macro.txt", row.names=FALSE, quote=FALSE,sep="\t")




function_of_GeneList <- dup_df[,c("Gene", "Protein.class")]

function_mat <- matrix(0, nrow=nrow(function_of_GeneList), ncol=25)
rownames(function_mat) <- function_of_GeneList$Gene
colnames(function_mat) <- c("Blood group antigen proteins","Cancer-related genes","Candidate cardiovascular disease genes", "CD markers", "Citric acid cycle related proteins", 
                            "Disease related genes","Enzymes", "FDA approved drug targets", "G-protein coupled receptors", "Human disease related genes", "Immunoglobulin genes",
                            "Metabolic proteins", "Nuclear receptors", "Plasma proteins", "Potential drug targets", "Predicted intracellular proteins","Predicted membrane proteins", 
                            "Predicted secreted proteins",  "RAS pathway related proteins", "Ribosomal proteins",
                            "RNA polymerase related proteins", "T-cell receptor genes", "Transcription factors", "Transporters", "Voltage-gated ion channels")
function_of_GeneList <- sapply(function_of_GeneList[,"Protein.class"], function(x) { x<-gsub("\"", "",x); strsplit(x, ", ")})

for ( i in 1:length(function_of_GeneList) ) {          # copilot 太强了吧
  for ( j in 1:length(function_of_GeneList[[i]]) ) {
    function_mat[i, function_of_GeneList[[i]][[j]]] <- 1
}}
function_mat <- data.frame(function_mat); function_mat$Gene <- rownames(function_mat)

pdf("./dup_cd8_cd4_macrophage.pdf",height=6,width=12)
upset(function_mat, nsets=50, sets.x.label="Protein class", mb.ratio=c(0.55,0.45), 
      main.bar.color = "steelblue", sets.bar.color = "steelblue", order.by=c("freq"), decreasing=c(TRUE),nintersects = 80 )
dev.off()


dup_ensembel <- dup_df[,1]
 
table_cd8[! table_cd8[,1] %in% dup_ensembel, 1]


table_macro[! table_macro[,1] %in% dup_ensembel, 1 ]



