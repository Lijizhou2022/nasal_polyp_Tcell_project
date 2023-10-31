#for LCM-seq data analysis
## Step1 load salmon results into Deseq2
library(tximport)
library("DESeq2")
sample_info <- read.table(file = info_file_dir,header = T)
files <- file.path(work_dir,samples$file,"quant.genes.sf")
txi <- tximport(files,type = "salmon",txOut = TRUE)
# deseq2 analysis
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
dds<-DESeq(ddsTxi)
## Step2 plot PCA
vsd <- vst(ddsTxi, blind=FALSE)
library("ggthemes")
library("RColorBrewer")
pcaData <- plotPCA(vsd, intgroup="condition", returnData=T)
percentVar <- round(100 * attr(pcaData, "percentVar"))
library(ggplot2)
library(ggrepel)
ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(aes(color = condition),size=3) + 
  scale_discrete_manual(values=colors,
    aesthetics = 'colour')+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title=element_blank())+
  labs(title="tnasa 857Land 912 CM",x="PC1",y="PC2")+
  xlab(paste0("PC1: ",round(percentVar[1],digits = 2),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2],digits = 2),"% variance"))

## Step3 gene diffent expression analysis
sampleA <- "GZMK_P"
sampleB <- "GZMK_N"
res <- results(dds,contrast=c( "condition",sampleA, sampleB))