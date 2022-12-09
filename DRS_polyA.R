library("ggplot2")
library(ggpubr)
library("scales")
t_type <- read.delim("~/genome/gene_transid_te_assemble_type_for_r",header=F)
colnames(t_type) <- c("GeneID","TransID","Type")

ddm1_trans_ngs <- read.delim("../6.stringtie/ddm1-2_trans_fpkm",header = F)
ccr4a_trans_ngs <- read.delim("../6.stringtie/ccr4A-ddm1-2_trans_fpkm",header = F)

ddm1_drs <- read.delim("../4.nanocount/ddm1_merge.counts.tsv")
ccr4a_drs <- read.delim("../4.nanocount/ccr4a_merge.counts.tsv")

ddm1_expression <- merge(ddm1_drs,ddm1_trans_ngs,by.x = "transcript_name",by.y="V1")
colnames(ddm1_expression)
colnames(ddm1_expression)[6]<- "FPKM"
cor(ddm1_expression$tpm,ddm1_expression$FPKM) #0.8466633
cor.test(ddm1_expression$tpm,ddm1_expression$FPKM)
ddm1_expression <- merge(ddm1_expression,t_type,by.x ="transcript_name",by.y="V2")
ddm1_gene_expression<- aggregate(ddm1_expression$tpm, by=list(GeneID=ddm1_expression$V1),sum)
###plot
plot(log2(ddm1_expression$FPKM+1),log2(ddm1_expression$tpm+1),
     pch=20,col=alpha("black", 0.3),
     xlab="NGS: log2 (ddm1-L2 FPKM +1)",ylab="DRS: log2 (ddm1-L2 TPM +1)")
abline(0,1,col="red")
text(1,14,expression(paste(italic(R)^italic('2'),' = ','0.85')))
text(1.5,13,expression(paste(italic(P),' < ','2.2 e -16')))
points(log2(ddm1_expression[ddm1_expression$transcript_name=="AT5G17125.1",]$FPKM+1),
       log2(ddm1_expression[ddm1_expression$transcript_name=="AT5G17125.1",]$tpm+1),
       col="blue",pch=20)



ccr4a_expression <- merge(ccr4a_drs,ccr4a_trans_ngs,by.x = "transcript_name",by.y="V1")
colnames(ccr4a_expression)
colnames(ccr4a_expression)[6]<- "FPKM"
cor(ccr4a_expression$tpm,ccr4a_expression$FPKM)  #0.8152016 
cor.test(ccr4a_expression$tpm,ccr4a_expression$FPKM)
ccr4a_expression <- merge(ccr4a_expression,t_type,by.x ="transcript_name",by.y="V2" )
ccr4a_gene_expression <- aggregate(ccr4a_expression$tpm, by=list(GeneID=ccr4a_expression$V9),sum)
plot(ccr4a_expression$FPKM,ccr4a_expression$tpm,pch=20)
table(ddm1_expression[ddm1_expression$tpm>=1,]$Type)
#expressed (TPM>=1) ddm1
# #25968     # 24451       1517   5.84% TE
table(ccr4a_expression[ccr4a_expression$tpm>=1,]$Type)
#expressed (TPM>=1) ccr4a
#26408      #24475       1933    7.31% similar pattern with NGS
###plot
plot(log2(ccr4a_expression$FPKM+1),log2(ccr4a_expression$tpm+1),pch=20,col=alpha("black", 0.3),
     xlab="NGS: log2 (ccr4a-1 ddm1-L2 FPKM +1)",ylab="DRS: log2 (ccr4a-1 ddm1-L2 TPM +1)")
abline(0,1,col="red")
text(1,14,expression(paste(italic(R)^italic('2'),' = ','0.82')))
text(1.5,13,expression(paste(italic(P),' < ','2.2 e -16')))
points(log2(ccr4a_expression[ccr4a_expression$transcript_name=="AT5G17125.1",]$FPKM+1),
       log2(ccr4a_expression[ccr4a_expression$transcript_name=="AT5G17125.1",]$tpm+1),
       col="blue",pch=20)

ddm1_expression <- ddm1_expression[,c(1,3:8)]
colnames(ddm1_expression)
colnames(ddm1_expression)[c(2:3,5:6)] <- c("ddm1 DRS est_count","ddm1 DRS TPM","ddm1 NGS FPKM","GeneID")
ccr4a_expression <- ccr4a_expression[,c(1,3:4,6)]
colnames(ccr4a_expression)
colnames(ccr4a_expression)[c(2:4)] <- c("ccr4a DRS est_count","ccr4a DRS TPM","ccr4a NGS FPKM")
ex <- merge(ddm1_expression,ccr4a_expression,by= "transcript_name")
cor(ex$`ddm1 NGS FPKM`,ex$`ccr4a NGS FPKM`)#0.941816
cor(ex$`ddm1 DRS TPM`,ex$`ddm1 NGS FPKM`)#0.8466633    first time0.8187243
cor(ex$`ccr4a DRS TPM`,ex$`ccr4a NGS FPKM`)# 0.8152016  first time0.7859564
cor(ex$`ddm1 DRS TPM`,ex$`ccr4a DRS TPM`)#0.990699  first time 0.9928805
cor.test(ex[ex$Type=="Gene",]$`ddm1 DRS TPM`,ex[ex$Type=="Gene",]$`ccr4a DRS TPM`) #0.9914149 
cor.test(ex[ex$Type=="Transposon",]$`ccr4a DRS TPM`,ex[ex$Type=="Transposon",]$`ddm1 DRS TPM`)#0.8751833 


###plot
plot(log2(ex$`ddm1 DRS TPM`+1),log2(ex$`ccr4a DRS TPM`+1),pch=20,col=alpha("black", 0.3),
     xlab="DRS: log2 (ddm1-L2 TPM +1)",ylab="DRS: log2 (ccr4a-1 ddm1-L2 TPM +1)")
abline(0,1,col="red")
text(1,14,expression(paste(italic(R)^italic('2'),' = ','0.99')))
text(1.5,13,expression(paste(italic(P),' < ','2.2 e -16')))

plot(log2(ex[ex$Type=="Transposon",]$`ddm1 DRS TPM`+1),
     log2(ex[ex$Type=="Transposon",]$`ccr4a DRS TPM`+1),
     col=alpha("black", 0.3),pch=20,
     xlab="DRS: log2 (ddm1-L2 TPM +1)",ylab="DRS: log2 (ccr4a-1 ddm1-L2 TPM +1)")
abline(0,1,col="red")
text(1,10,expression(paste(italic(R)^italic('2'),' = ','0.88')))
text(1.5,9,expression(paste(italic(P),' < ','2.2 e -16')))

ccr4a_polya <- read.delim("../3.polya/ccr4a_assembly_pass.tsv")
ddm1_polya <- read.delim("../3.polya/ddm1_assembly_pass.tsv")
length(table(ccr4a_polya$contig)) #34231  first time28878
length(table(ddm1_polya$contig)) #33091   first time30471
ccr4a_reads <- as.data.frame(table(ccr4a_polya$contig))
ddm1_reads <- as.data.frame(table(ddm1_polya$contig))
ccr4a_polya_median <- aggregate(ccr4a_polya$polya_length,by=list(TransID=ccr4a_polya$contig),median)
ddm1_polya_median <- aggregate(ddm1_polya$polya_length,by=list(TransID=ddm1_polya$contig),median)
colnames(ccr4a_polya_median)[2]<-"median"
colnames(ddm1_polya_median)[2]<-"median"
ccr4a_polya_mean <- aggregate(ccr4a_polya$polya_length,by=list(TransID=ccr4a_polya$contig),mean)
ddm1_polya_mean <- aggregate(ddm1_polya$polya_length,by=list(TransID=ddm1_polya$contig),mean)
colnames(ccr4a_polya_mean)[2]<-"mean"
colnames(ddm1_polya_mean)[2]<-"mean"
ccr4a_polya_median_reads <- merge(ccr4a_polya_median,ccr4a_reads,
                                  by.x="TransID",by.y="Var1")
ccr4a_polya_median_reads <- merge(ccr4a_polya_median_reads,ccr4a_polya_mean,
                                  by="TransID")
ddm1_polya_median_reads<- merge(ddm1_polya_median,ddm1_reads,
                                by.x="TransID",by.y="Var1")
ddm1_polya_median_reads<- merge(ddm1_polya_median_reads,ddm1_polya_mean,
                                by="TransID")
ccr4a_polya_median_expression <-  merge(ccr4a_polya_median_reads,ccr4a_drs,
                                        by.x="TransID",by.y="transcript_name")
ccr4a_polya_type <- merge(ccr4a_polya_median_expression,t_type,by="TransID")

ddm1_polya_median_expression <-  merge(ddm1_polya_median_reads,ddm1_drs,
                                       by.x="TransID",by.y="transcript_name")
ddm1_polya_type <- merge(ddm1_polya_median_expression,t_type,by="TransID")
write.table(ddm1_polya_type,"ddm1-L2_median_ploya_length.txt",sep="\t",row.names = F,quote=F)
write.table(ccr4a_polya_type,"ccr4a-1_ddm1-L2_median_ploya_length.txt",sep="\t",row.names = F,quote=F)

table(ccr4a_polya_type$Type)
#Gene 30714   TE    3408  TE 9.987691%  
table(ddm1_polya_type$Type)
#Gene 30310   TE    2685  TE 8.137597%   
hist((ccr4a_polya_type[ccr4a_polya_type$Type=="Gene",]$median),xlim=c(0,600),breaks=120)
hist((ddm1_polya_type[ddm1_polya_type$Type=="Gene",]$median),xlim=c(0,600),breaks=120)

hist((ccr4a_polya_type[ccr4a_polya_type$Type=="Transposon",]$median),xlim=c(0,600),breaks=600)
hist((ddm1_polya_type[ddm1_polya_type$Type=="Transposon",]$median),xlim=c(0,600),breaks=600)


cor.test(ddm1_polya_type$mean,ddm1_polya_type$median) ##0.9741043 p-value < 2.2e-16
cor.test(ccr4a_polya_type$mean,ccr4a_polya_type$median)##0.9738196 p-value < 2.2e-16


###ccr4a-1/ddm1-L2 mutants
ch_15 <- ccr4a_polya[ccr4a_polya$contig %in% ccr4a_polya_type[ccr4a_polya_type$Freq>=15,]$TransID,]
###transposon set frequency >=5
ch_5_te <- ccr4a_polya[ccr4a_polya$contig %in% ccr4a_polya_type[ccr4a_polya_type$Freq>=5&ccr4a_polya_type$Type=="Transposon",]$TransID,]
length(unique(ch_5_te$readname))
ch_5_te$genotype <- c("ccr4a-1 ddm1-L2")
summary(ch_5_te$polya_length)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.44   64.77  106.92  119.43  159.74  548.00 
###transposon all unique reads 5nt bin
ch_5_te_unique <- unique(ch_5_te[,c(1,4:10)])
ch_5_te_unique$genotype <- c("ccr4a-1 ddm1-L2")
ch <- hist(ch_5_te_unique$polya_length,breaks = 500)
View(ch$breaks)
chr <- as.data.frame(cbind(ch$breaks[-1],ch$density))
###all gene all reads 1nt bin
length(unique(ch_15$readname))
ch_15_unique <- unique(ch_15[,c(1,4:10)])
ch_15_unique$genotype <- c("ccr4a-1 ddm1-L2")
summary(ch_15_unique$polya_length)
ch <- hist(ch_15_unique$polya_length,breaks = 702)
chr <- as.data.frame(cbind(ch$breaks[-1],ch$density))
library("pheatmap")
pheatmap(t(chr$V2),cluster_rows = F,cluster_cols = F)
summary(ccr4a_polya$polya_length)
summary(ddm1_polya$polya_length)


###ddm1-L2 mutants
dh_15 <- ddm1_polya[ddm1_polya$contig %in% ddm1_polya_type[ddm1_polya_type$Freq>=15,]$TransID,]
###transposon set frequency >=5
dh_5_te <- ddm1_polya[ddm1_polya$contig %in% ddm1_polya_type[ddm1_polya_type$Freq>=5&ddm1_polya_type$Type=="Transposon",]$TransID,]
length(unique(dh_5_te$readname))
###transposon all unique reads 5nt bin
dh_5_te <- unique(dh_5_te[,c(1,4:10)])
dh_5_te$genotype <- c("ddm1-L2")
summary(dh_5_te$polya_length)
dh<- hist(dh_5_te$polya_length,breaks =600)
View(dh)
dhr <- as.data.frame(cbind(dh$breaks[-1],dh$density))
pheatmap(t(dhr$V2),cluster_rows = F,cluster_cols = F)
###draw all unique reads 1nt bin
dh_15_unique <- unique(dh_15[,c(1,4:10)])
summary(dh_15_unique$polya_length)
dh<- hist(dh_15_unique$polya_length,breaks = 669)
dhr <- as.data.frame(cbind(dh$breaks[-1],dh$density))
pheatmap(t(chr$V2),cluster_rows = F,cluster_cols = F)
hr <- merge(dhr,chr,by="V1")
colnames(hr)[2:3] <- c("ddm1-L2","ccr4a-1 ddm1-L2")
write.table(hr,"~/hr.txt",sep="\t",row.names = F,quote=F)
hr<-read.delim("~/hr.txt")
pheatmap(t(hr[1:250,2:3]),cluster_rows = F,cluster_cols = F,
         cellheight = 60,border_color = NA,show_colnames = F)
#color = colorRampPalette(c("white","yellow", "red"))(100))
mr <- merge(dmr,cmr,by="V1",all=T)
mr[is.na(mr$`cc4a-1 ddm1-L2`),]$`cc4a-1 ddm1-L2` <- c(0)
colnames(mr)[2:3] <- c("ddm1-L2","cc4a-1 ddm1-L2")
write.table(mr,"~/mr.txt",sep="\t",row.names = F,quote=F)
mr<-read.delim("~/mr.txt")
pheatmap(t(mr[,2:3]),cluster_rows = F,cluster_cols = F,
         cellheight = 60)

a <- merge(hr,mr,by="V1",all=T)
pheatmap(t(a[,2:5]),cluster_rows = F,cluster_cols = F,
         cellheight = 60)
plot(density(dh_15_unique$polya_length),col="green",xlim=c(0,250))
lines(density(ch_15_unique$polya_length))
summary(polya$`ddm1 median`)
summary(polya$`ccr4a median`)

###draw all unique reads hist+density plot
dh_15_unique <- unique(dh_15[,c(1,4:10)])
dh_15_unique$genotype <- "ddm1-L2"
ch_15_unique <- unique(ch_15[,c(1,4:10)])
ch_15_unique$genotype <- c("ccr4a-1 ddm1-L2")
hr <- rbind(ch_15_unique,dh_15_unique)
library("ggplot2")
ggplot(hr, aes(x =polya_length,fill=genotype,color=genotype)) + 
  scale_fill_manual(values=c("red", "blue"))+
  scale_color_manual(values=c("#e41a1c", "#253494"))+
  geom_histogram(aes(y = ..density..),position="identity",
                 binwidth=8,alpha=0.2,colour = "white")+
  geom_density(alpha=0.1)+xlim(0,250)+labs(x="All reads ployA length")+
  theme_classic()+theme(legend.position = c(0.8, 0.8))

###draw all TE unique reads hist+density plot
dh_5_te <- ddm1_polya[ddm1_polya$contig %in% ddm1_polya_type[ddm1_polya_type$Freq>=5&ddm1_polya_type$Type=="Transposon",]$TransID,]
length(unique(dh_5_te$readname))
length(unique(dh_5_te$contig))
dh_5_te$genotype <- c("ddm1-L2")
ch_5_te <- ccr4a_polya[ccr4a_polya$contig %in% ccr4a_polya_type[ccr4a_polya_type$Freq>=5&ccr4a_polya_type$Type=="Transposon",]$TransID,]
length(unique(ch_5_te$readname))
length(unique(ch_5_te$contig))
wilcox.test(hr[hr$genotype=="ccr4a-1 ddm1-L2",]$polya_length,hr[hr$genotype=="ddm1-L2",]$polya_length)

dh_5_unique <- unique(dh_5_te[,c(1,4:10)])
dh_5_unique$genotype <- "ddm1-L2"
ch_5_unique <- unique(ch_5_te[,c(1,4:10)])
ch_5_unique$genotype <- c("ccr4a-1 ddm1-L2")
hr <- rbind(ch_5_unique,dh_5_unique)
library("ggplot2")
ggplot(hr, aes(x =polya_length,fill=genotype,color=genotype)) + 
  scale_fill_manual(values=c("red", "blue"))+
  scale_color_manual(values=c("#e41a1c", "#253494"))+
  geom_histogram(aes(y = ..density..),position="identity",
                 binwidth=8,alpha=0.2,colour = "white")+
  geom_density(alpha=0.1)+xlim(0,250)+labs(x="All TE reads ployA length")+
  theme_classic()+theme(legend.position = c(0.8, 0.8))
wilcox.test(hr[hr$genotype=="ccr4a-1 ddm1-L2",]$polya_length,hr[hr$genotype=="ddm1-L2",]$polya_length)


###calculate p values for transcripts 
polya$polya.pvalue <- "0"
#pvalue for reads polya length
for (i in 1:29146) {
  m <- polya$TransID[i]
  polya[i,]$polya.pvalue <-wilcox.test(ccr4a_polya[ccr4a_polya$contig==m,]$polya_length,
                                                        ddm1_polya[ddm1_polya$contig==m,]$polya_length,paired=FALSE)$p.value
}
#multiple FDR adjust p-value
polya$polya.fdr =p.adjust(polya$polya.pvalue,method = "BH")
##One gene chose one highest reads number transcripts
c <- ccr4a_polya_type[ccr4a_polya_type$TransID %in% polya$TransID,]
chose1 <- c[order(c$GeneID,c$Freq,decreasing = T),]
cc <- chose1[!duplicated(chose1$GeneID),]
sel_polya<- polya[polya$TransID%in%cc$TransID,]
write.table(sel_polya,"ccr4a_ddm1_median_polya_length_with_pvalue_highest_expression_transcrips_in_ccr4a_selected.txt",sep="\t",row.names = F,quote=F)
