getwd()
t_type <- read.delim("../../gene_type",header = F)
colnames(t_type) <- c("GeneID","Type")
g <- read.csv("../gene_count_matrix.csv",row.names = "gene_id")
countdata <- as.matrix(g)
library("DESeq2")
library("ggplot2")
coldata <- data.frame("genotype"=c(rep(c("caf1a-1 caf1b-3 ddm1","ccr4a-1 ddm1"),each=2),
                                   rep(c("ccr4b-1 ddm1","ddm1"),each=2)))
colnames(countdata) <- c("caf1a-1 caf1b-3 ddm1-L1","caf1a-1 caf1b-3 ddm1-L2",
                         "ccr4a-1 ddm1-L1","ccr4a-1 ddm1-L2",
                         "ccr4b-1 ddm1-L1","ccr4b-1 ddm1-L2",
                         "ddm1-L1","ddm1-L2")
rownames(coldata) <- colnames(countdata)
all(rownames(coldata) %in% colnames(countdata))
dds <- DESeqDataSetFromMatrix(countData = countdata,
                               colData = coldata, 
                               design = ~genotype)
dds
dds <- DESeq(dds)

###Drawing sample matrix and PCA
rld <- rlog(dds)
sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
library( pheatmap )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap( sampleDistMatrix,
          clustering_distance_rows = sampleDists, 
          clustering_distance_cols = sampleDists,
          col = colours, border_color = NA )
rld$genotype<-factor(rld$genotype,levels = c("ddm1","ccr4a-1 ddm1","ccr4b-1 ddm1","caf1a-1 caf1b-3 ddm1"))
plotPCA(rld,intgroup=c("genotype"))

##ccr4a volcano plot
colData_4A<- data.frame("genotype"=c(rep(c("ccr4a-1 ddm1"),each=2),rep(c("ddm1"),each=2)))
colData_4A
dds_4A <- DESeqDataSetFromMatrix(countData = countData[,c(3,4,7,8)], 
                                 colData = colData_4A, 
                                 design = ~genotype)
dds_4A$genotype
dds_4A$genotype<- factor(dds_4A$genotype,levels = c("ddm1","ccr4a-1 ddm1"))
dds_4A$genotype
ds_4A <- DESeq(dds_4A)
res_4A <- results(ds_4A)
summary(res_4A)
sum(res_4A$padj < 0.05, na.rm=TRUE)
re<- as.data.frame(res_4A)
re$Significant  <- "No"
re[is.na(re$padj),]$Significant  <- "Not test"
re[re$padj <= 0.05 & re$log2FoldChange>=1 &!is.na(re$padj),]$Significant <-"Up"
re[re$padj <= 0.05 &re$log2FoldChange<=-1&!is.na(re$padj),]$Significant  <-"Down"
re$Significant  <- factor(re$Significant,levels = c("Up","Down","No"))
table(re[,c(7)])
#v <- ggplot(re,aes(log2FoldChange,-log10(padj)))
#v + geom_point(aes(colour=Significant))+geom_vline(xintercept=c(-1,1),linetype=4,colour="grey")+geom_hline(yintercept=-log10(0.05),linetype=4,colour="grey")+xlim(-17, 17)
re$GeneID <- rownames(re)
re <- re[,c(8,1:7)]
te_name <- read.delim("../../../../genome/TE_famile_type.txt")
t_type_name <- merge(t_type,te_name,by="GeneID",all = TRUE)
res4a <- merge(re,t_type_name,by="GeneID")
table(res4a[,c(8:9)])
#v <- ggplot(res4a,aes(log2FoldChange,-log10(padj)))
#v + geom_point(aes(colour=Significant,shape=Type))
#+geom_vline(xintercept=c(-1,1),linetype=4,colour="grey")+geom_hline(yintercept=-log10(0.05),linetype=4,colour="grey")+xlim(-17, 17)
library("ggrepel")
ggplot(res4a,aes(log2FoldChange,-log10(padj))) + geom_point(aes(colour=Type,shape=Type))+
  geom_vline(xintercept=c(-1,1),linetype=4,colour="black")+
  geom_hline(yintercept=-log10(0.05),linetype=4,colour="black")+
  xlim(-13, 13)+ylim(0,70)+theme_bw()+
  theme(axis.text=element_text(size=12,color="black"))+ 
  scale_color_manual(values=c("grey","grey"))+
  geom_point(data = res4a[res4a$Type=="Transposon"&res4a$Significant!="No",], 
             aes(x=log2FoldChange,-log10(padj),shape=Type), color = 'red',alpha=0.3)+
  labs(x="log2 FC (ccr4a-1 ddm1/ddm1)",y="-Log10 FDR")
  
write.table(ccr4a,"ccr4a_ddm1_DEG.txt",sep="\t",quote=F,row.names = F)


###ccr4b volcano plot
colData_4B<- data.frame("genotype"=c(rep(c("ccr4b-1"),each=2),rep(c("ddm1"),each=2)))
colData_4B
dds_4B <- DESeqDataSetFromMatrix(countData = countData[,c(5:8)], 
                                 colData = colData_4B, 
                                 design = ~genotype)
dds_4B$genotype
dds_4B$genotype<- factor(dds_4B$genotype,levels = c("ddm1","ccr4b-1"))
dds_4B$genotype
ds_4B <- DESeq(dds_4B)
res_4B <- results(ds_4B)
summary(res_4B)
sum(res_4B$padj < 0.05, na.rm=TRUE)
re<- as.data.frame(res_4B)
re$Significant  <- "No"
re[is.na(re$padj),]$Significant  <- "Not test"
re[re$padj <= 0.05 & re$log2FoldChange>=1 &!is.na(re$padj),]$Significant <-"Up"
re[re$padj <= 0.05 &re$log2FoldChange<=-1&!is.na(re$padj),]$Significant  <-"Down"
re$Significant  <- factor(re$Significant,levels = c("Up","Down","No"))
table(re[,c(7)])
v <- ggplot(re,aes(log2FoldChange,-log10(padj)))
v + geom_point(aes(colour=Significant))+geom_vline(xintercept=c(-1,1),linetype=4,colour="grey")+geom_hline(yintercept=-log10(0.05),linetype=4,colour="grey")+xlim(-17, 17)
re$GeneID <- rownames(re)
re <- re[,c(8,1:7)]
res4b <- merge(re,t_type_name,by="GeneID")
table(res4b[,c(8:9)])
v <- ggplot(res4b,aes(log2FoldChange,-log10(padj)))
v + geom_point(aes(colour=Type,shape=Type))+
  geom_vline(xintercept=c(-1,1),linetype=4,colour="black")+
  geom_hline(yintercept=-log10(0.05),linetype=4,colour="black")+
  xlim(-13, 13)+ylim(0,40)+theme_bw()+
  theme(axis.text=element_text(size=12,color="black"))+
  scale_color_manual(values=c("grey","grey"))+
  geom_point(data = res4b[res4b$Type=="Transposon"&res4b$Significant!="No",], aes(x=log2FoldChange,-log10(padj),shape=Type), color = 'red',alpha=0.3)

write.table(res4b,"ccr4B_ddm1_DEG.txt",sep="\t",quote=F,row.names = F)


###caf1a-1caf1b-3 volcano plot
colData_4AB<- data.frame("genotype"=c(rep(c("cafAB"),each=2),rep(c("ddm1"),each=2)))
colData_4AB
dds_4AB <- DESeqDataSetFromMatrix(countData = countData[,c(1:2,7:8)], 
                                  colData = colData_4AB, 
                                  design = ~genotype)
dds_4AB$genotype
dds_4AB$genotype<- factor(dds_4AB$genotype,levels = c("ddm1","cafAB"))
dds_4AB$genotype
ds_4AB <- DESeq(dds_4AB)
res_4AB <- results(ds_4AB)
summary(res_4AB)
sum(res_4AB$padj < 0.05, na.rm=TRUE)
re<- as.data.frame(res_4AB)
re$Significant  <- "No"
re[is.na(re$padj),]$Significant  <- "Not test"
re[re$padj <= 0.05 & re$log2FoldChange>=1 &!is.na(re$padj),]$Significant <-"Up"
re[re$padj <= 0.05 &re$log2FoldChange<=-1&!is.na(re$padj),]$Significant  <-"Down"
re$Significant  <- factor(re$Significant,levels = c("Up","Down","No"))
table(re[,c(7)])
#v <- ggplot(re,aes(log2FoldChange,-log10(padj)))
#v + geom_point(aes(colour=Significant))+geom_vline(xintercept=c(-1,1),linetype=4,colour="grey")+geom_hline(yintercept=-log10(0.05),linetype=4,colour="grey")+xlim(-17, 17)
re$GeneID <- rownames(re)
re <- re[,c(8,1:7)]
res4ab <- merge(re,t_type_name,by="GeneID")
table(res4ab[,c(8:9)])
v <- ggplot(res4ab,aes(log2FoldChange,-log10(padj)))
v + geom_point(aes(colour=Type,shape=Type))+
  geom_vline(xintercept=c(-1,1),linetype=4,colour="black")+
  geom_hline(yintercept=-log10(0.05),linetype=4,colour="black")+
  xlim(-13, 13)+ylim(0,60)+theme_bw()+
  theme(axis.text=element_text(size=12,color="black"))+
  scale_color_manual(values=c("grey","grey"))+
  geom_point(data = res4ab[res4ab$Type=="Transposon"&res4ab$Significant!="No",], 
             aes(x=log2FoldChange,-log10(padj),shape=Type), color = 'red',alpha=0.3)+
  labs(x="log2 FC (caf1a−1 caf1b−3 ddm1/ddm1)",y="-Log10 FDR")

write.table(res4ab,"cafAcafB_ddm1_DEG.txt",sep="\t",quote=F,row.names = F)


##venn upregulated transposon
library(VennDiagram)
library("scales")
overrideTriple=F
data4A<-read.delim("ccr4a_ddm1_DEG.txt",header = T)
data4Aup <- data4A[data4A$Significant=="Up" & data4A$Type=="Transposon" & !is.na(data4A$baseMean),]
data4Aup<-data4Aup[!is.na(data4Aup$GeneID),]
data4B<-read.delim("ccr4B_ddm1_DEG.txt",header = T)
data4Bup <- data4B[data4B$Significant=="Up" & data4B$Type=="Transposon" & !is.na(data4B$baseMean),]
data4Bup<-data4Bup[!is.na(data4Bup$GeneID),]
data4AB<-read.delim("cafAcafB_ddm1_DEG.txt",header = T)
data4ABup <- data4AB[data4AB$Significant=="Up" & data4AB$Type=="Transposon" & !is.na(data4AB$baseMean),]
data4ABup<-data4ABup[!is.na(data4ABup$GeneID),]
x <- list("ccr4a-1"= data4Aup$GeneID,"ccr4b-1"= data4Bup$GeneID,"caf1a-1 caf1b-3"= data4ABup$GeneID)
my_plot <- 
  venn.diagram(x,filename=NULL,
             height = 400, 
             width = 400, scaled = FALSE,
             resolution = 300,
             compression = "lzw",
             lwd = 1.5,
             col=c("#e41a1c","#4daf4a","#377eb8"),
             fill = c(alpha('#e41a1c',0.3), alpha('#4daf4a',0.3), alpha('#377eb8',0.3)),
             cex = 0.9,
             fontfamily = "Arial",
             cat.cex = 0.9,
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "Arial",
             cat.fontface = "italic",
             rotation = 1,
     main="Up regulated TEs",
     main.fontfamily = "Arial",
     main.cex = 0.9,
)

ggsave(my_plot, file="up_te.svg", device = "svg",width =4, height=4)
library("svglite")#install.packages("svglite")


###TE family braplot
te_family <- read.delim("../../../../genome/TE_famile_type.txt",header=T)
te <- as.data.frame(table(te_family$Transposon_Super_Family))
te$class <- "All TEs"
te<-te[c(2:5,8:12),]  ###DNA/Pogo  1       LINE? 1 DNA4 in reference were not calculated.
te4ab<-as.data.frame(table(data4ABup$Transposon_Super_Family))
te4ab$class <- "Upregulated in caf1a-1 caf1b-3"
te4a<-as.data.frame(table(data4Aup$Transposon_Super_Family))
te4a$class <- "Upregulated in ccr4a-1"
te4b<-as.data.frame(table(data4Bup$Transposon_Super_Family))
te4b$class <- "Upregulated in ccr4b-1"
rdr6tar <- read.delim("RDR6-targeted.txt")
ddm1tar<- read.delim("ddm1_target_TEs.txt")
rdr6target <- as.data.frame(table(te_family[te_family$GeneID%in%rdr6tar$GeneID,]$Transposon_Super_Family))
rdr6target$class <- "Targeted by RDR6"
rdr6target <- rdr6target[2:10,]  #DNA 1 was not calculated
ddm1target <- as.data.frame(table(te_family[te_family$GeneID%in%ddm1tar$GeneID,]$Transposon_Super_Family))
ddm1target$class <- "Targeted by ddm1"

tefamily <- rbind(te,te4a,te4b,te4ab,rdr6target,ddm1target)
tapply(tefamily[,2],tefamily[,c(3)],sum)
write.table(tefamily,"te_family.txt",sep="\t",quote=FALSE,row.names = F)
te_class<-read.delim("te_family.txt")
colnames(te_class)[1] <- c("Type")
te_class$class <- factor(te_class$class,
                  levels = c("All TEs",
                             "Targeted by ddm1",
                             "Targeted by RDR6",
                  "Upregulated in ccr4a-1",
                  "Upregulated in ccr4b-1",
                  "Upregulated in caf1a-1 caf1b-3"))
te_class$Fraction <- te_class$Freq/te_class$sum*100

library("RColorBrewer")
mcolors <- brewer.pal(9, "Set1")
ggplot(te_class, aes(x=class, y=Fraction, fill=Type))+
  geom_bar(width=0.7,alpha=0.8,stat = "identity")+theme_bw()+scale_fill_manual(values=mcolors)+
  theme(panel.grid=element_blank())+  #,axis.line=element_line(size=0.5)
  theme(axis.text=element_text(size=12,color="black"),axis.text.x = element_text(angle=30,hjust=1,vjust=1))+ 
  labs(x="",y="Fraction (%)")


#up-down-gene-number

df <- read.delim("ddm1_2replicates_ccr4b_2replicates_up_down.txt")
df$genotype <- factor(df$genotype,levels=c("ccr4a-1","ccr4b-1","caf1a-1 caf1b-3"))
df$coding <-factor(df$coding,levels=c("Gene","TE"))
df$type <- factor(df$type,levels=c("Up","Down"))
library("ggmap")#install.packages("ggmap")
# Create a barplot
colors=c("Black","red","grey","pink")
ggplot(df, aes(x=genotype, y=number, fill=interaction(coding,type)))+
  geom_bar(width = 0.5, stat = "identity")+
  facet_wrap(~coding)+
  geom_text(aes(label=number), 
            position = position_dodge2(width = 0.8),#,preserve = 'single'),
            vjust = -0.3, hjust = 0.15)+ scale_fill_manual(values=colors)+ 
  theme_bw()+#geom_hline(yintercept = 0,color="black")+
  theme(panel.grid=element_blank())+  #,axis.line=element_line(size=0.5)
  theme(axis.text=element_text(size=12,color="black"),axis.text.x = element_text(angle=30,hjust=1,vjust=1,face="italic"))+ 
  scale_y_continuous(breaks =c(-1000,-800,-600,-400,-200,0,200,400,600))+
  labs(x="",y="Number")


# Fold change plot
res4ab$genotype <- c("caf1a-1 caf1b-3")
res4a$genotype <- c("ccr4a-1")
res4b$genotype <- c("ccr4b-1")
ccr4a_cafab <- merge(res4a,res4ab,by="GeneID")
ccr4a_cafab$Type.x <- factor(ccr4a_cafab$Type.x,levels=c("Transposon","Gene"))
df_layer_1 <- ccr4a_cafab[ccr4a_cafab$Type.x=="Gene"&!is.na(ccr4a_cafab$log2FoldChange.x)&!is.na(ccr4a_cafab$log2FoldChange.y),]
df_layer_2 <- ccr4a_cafab[ccr4a_cafab$Type.x=="Transposon",]
df_layer_2_up <- df_layer_2[df_layer_2$Significant.x=="Up"&df_layer_2$Significant.y=="Up",]
df_layer_2_upp <- df_layer_2_up[!is.na(df_layer_2_up$GeneID),]

ggplot() + geom_hline(yintercept = 0,color="black" )+
  geom_vline(xintercept = 0,color="black")+
  geom_point(
    data=df_layer_2_upp, 
    aes(x=log2FoldChange.y, 
        y=log2FoldChange.x), 
   shape=17,alpha=0.3)+scale_color_manual(values=c("#377eb8"))+
  theme_bw()+xlim(0,10)+ylim(0,10)+geom_abline(slope=1,linetype ="dashed",col="red" )+
  theme(panel.grid=element_blank())+
  theme(axis.text=element_text(size=10,color="black"))+ 
  labs(x="Log 2 FC (caf1a-1 caf1b-3 ddm1/ddm1)",y="Log 2 FC (ccr4a-1 ddm1/ddm1)")

table(df_layer_2[,c(8,23)])
library("pheatmap")
#pheatmap(df_layer_2_upp[,c(3,18)],scale="row",show_rownames = F)
fpkm <- read.delim("../all_sample_gene_fpkm")
fpkm_te <- merge(fpkm,t_type_name,by="GeneID")
fpkm_ete <- fpkm_te[fpkm_te$Type=="Transposon",] 
fpkm_exte  <- fpkm_ete[fpkm_ete$cafA1cafB3.ddm1.1>=1 |
                         fpkm_ete$cafA1cafB3.ddm1.2>=1|
                         fpkm_ete$ccr4A.ddm1.1>=1|
                         fpkm_ete$ccr4A.ddm1.2>=1|
                         fpkm_ete$ccr4B.ddm1.1>=1|
                         fpkm_ete$ccr4B.ddm1.2>=1|
                         fpkm_ete$ddm1.1>=1|
                         fpkm_ete$ddm1.2>=1,] 
fpkm_exte  <- fpkm_ete[fpkm_ete$GeneID%in%df_layer_2_upp$GeneID,]
colnames(fpkm_exte)[c(8,10,12,13,14,15,17,18)] <-
  c("caf1a-1 caf1b-3 ddm1-L1","caf1a-1 caf1b-3 ddm1-L2","ccr4a-1 ddm1-L1","ccr4a-1 ddm1-L2","ccr4b-1 ddm1-L1","ccr4b-1 ddm1-L2","ddm1-L1","ddm1-L2")
pheatmap(fpkm_exte[,c(17,18,12:15,8,10)],scale="row",
         show_rownames = F,cluster_cols = F)
               

##Genome TE number

ccra_te <- read.delim("/Volumes/ChoLab/Ling/circos-0.69-9/ccr4a_te_num.txt",header = F)
rdr6_te <- read.delim("/Volumes/ChoLab/Ling/circos-0.69-9/rdr6_te_num.txt",header = F)
ddm_te <- read.delim("/Volumes/ChoLab/Ling/circos-0.69-9/ddm1_te_num.txt",header = F)

max(ccra_te$V4)
max(rdr6_te$V4)
max(ddm_te$V4)
plot(1:1193,ccra_te$V4,
     ylim=c(-32,15),type="l",lwd=1.5,xlab="Coordinates (Kb)",ylab="TEs Number",col="red")
lines(1:1193,-rdr6_te$V4,lwd=1.5)
lines(1:1193,-ddm_te$V4,lwd=1.5,col="grey")
tair<- read.delim("/Volumes/ChoLab/Ling/circos-0.69-9/tair.genome.kar.txt")
abline(v=c(0,tair$cvend/100000))
abline(v=tair$cvpericenstart/100000,col="grey")
abline(v=tair$cvpericenend/100000,col="grey")


plot(1:2385,ccra_te$V4,
     ylim=c(-40,40),type="l",lwd=1.5,xlab="Coordinates (Kb)",ylab="TEs Number",col="red")
lines(1:2385,-rdr6_te$V4,lwd=1.5)
lines(1:2385,ddm_te$V4,lwd=1.5,col="grey")
lines(1:2385,-ddm_te$V4,lwd=1.5,col="grey")
abline(v=c(0,tair$cvend/50000))
abline(v=tair$cvpericenstart/50000,col="grey")
abline(v=tair$cvpericenend/50000,col="grey")



plot(1:2385,ddm_te$V4,
     ylim=c(-40,40),type="l",lwd=1.5,xlab="Coordinates (Kb)",
     ylab="TEs Number",col="grey")
lines(1:2385,-ddm_te$V4,lwd=1.5,col="grey")
lines(1:2385,ccra_te$V4,lwd=1.5,col="red")
lines(1:2385,-rdr6_te$V4,lwd=1.5)
abline(v=c(0,tair$cvend/50000))
abline(v=tair$cvpericenstart/50000,col="grey")
abline(v=tair$cvpericenend/50000,col="grey")
