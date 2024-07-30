rm(list=ls())
setwd("/Users/chenzifa/Desktop/analyst/zhongqian/20210707-heatmap/figure7/")
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)

## HaCaT 和 Fibroblast s
gene <- read.csv("./article simplified.csv")
gene$Term <- as.factor(gene$Term)
gene$SASP_class <- as.factor(gene$SASP_class)

mapping0 <- subset(gene,select = c(Term,SASP_class,genes))
rownames(mapping0) <- mapping0$genes
data <- subset(gene,select=-c(Term,SASP_class))
rownames(data) <- data$genes
data <- data[,-1]
data0 <- data




## HaCaT
df <- data0[,1:3]
data <- df
data[is.na(data)] <- 0
data <- data[!(rowSums(data)==0),]
names(data) <- c("M.osloensis","C.acnes","S.epidermidis")
all <- subset(data,select = c(c(C.acnes,M.osloensis,S.epidermidis)))

direction <- readxl::read_xlsx("./基因-衰老对应关系(1).xlsx",sheet = 4)
names(direction) <- c("genes","direction","ratio")
direction <- direction[-1,]
direction <- as.data.frame(direction)
rownames(direction) <- direction$genes
direction <- direction[,-3]
mapping <- merge(mapping0,direction,by="genes")
rownames(mapping) <- mapping$genes
all$genes <- rownames(all)
data <- merge(mapping,all,by="genes")
data <- data[order(data$direction) ,]

  
## FB
df <- data0[,4:6]
data <- df
data[is.na(data)] <- 0
data <- data[!(rowSums(data)==0),]
names(data) <- c("M.osloensis","C.acnes","S.epidermidis")
all <- subset(data,select = c(c(C.acnes,M.osloensis,S.epidermidis)))

direction <- readxl::read_xlsx("./基因-衰老对应关系(1).xlsx",sheet = 2)
names(direction) <- c("genes","direction","ratio")
direction <- direction[-1,]
direction <- as.data.frame(direction)
rownames(direction) <- direction$genes
direction <- direction[,-3]
mapping <- merge(mapping0,direction,by="genes")
rownames(mapping) <- mapping$genes
all$genes <- rownames(all)
data <- merge(mapping,all,by="genes")
data <- data[order(data$direction) ,]



rownames(data) <- data$genes
df <- as.matrix(data[,5:7])
mapping <- data[,1:4]
range(df)
brew9 <- brewer.pal(9, "Set1")
col = colorRamp2(c(-1,0,17), c(brew9[2],'white',brew9[1]))  ## Fb
col = colorRamp2(c(-1,0,2), c(brew9[2],'white',brew9[1]))  ## HACaT

col2 = col3 = brewer.pal(9,"Oranges")
annotation = rowAnnotation(Term = mapping$Term,
                            SASP_class = mapping$SASP_class,
                           Direction = mapping$direction,
                            col = list( Term=c("Cell cycle" = col2[9], 
                                               "Deregulated Metabolic Profile" = col2[8], 
                                               "Macromolecular Damage"=col2[7],
                                               "SASP"=col2[6]),
                                        SASP_class=c("Interleukins" = col3[2], 
                                                     "Growth factors; regulators" = col3[1], 
                                                     "Proteases and regulators"=col3[3],
                                                     "Receptors;ligands"= col3[4],
                                                     "none"=col3[5]),
                                        Direction=c("U" = brew9[1], 
                                               "D" = brew9[2])
                                        ))






annotation_row = data.frame(
  GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6)))
)
rownames(annotation_row) = paste("Gene", 1:20, sep = "")
head(annotation_row)

pdf("h.pdf",height = 10,width = 6)
  heatmap(df, name = "Log2(FC)", 
           cluster_rows = F,   
           cluster_columns = F,
           row_order = rownames(df),
           column_order = colnames(df),
           show_row_names = T,  
           show_column_names = T
          )

dev.off()
heatmap(df)

