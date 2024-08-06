# https://rockefelleruniversity.github.io/RU_ATAC_Workshop.html
library(ChIPQC)
library(GenomicFeatures)
library(ChIPseeker)
library(dplyr)
library(genomation)
library(soGGi)
library(ggupset)
library(limma)
library(tidyr)
library(data.table)
library(RColorBrewer)

#We want to annotate for each cluster the peaks
#where are the peaks within the clusters 
#in which genomic regions and how many peaks are in those regions within the clusters 

#read the file
#this transcriptome was generated with PASA, the longest reads were taken for the transcriptome 
escgff=read.table("/Users/pui/Documents/Lab_Vienna/ATAC/Peakfiles/PASA_transcript/GRanges_escgff")
names=c("chromosome","start","end","strand","type","transcript_id","gene_id","exon_id")
colnames(escgff)<-names
head(escgff,n=10)

#reformating into a granges object for the programmes, into a gene finding formate, annotations with gene information
GRanges_escgff_ex<-makeGRangesFromDataFrame(escgff,keep.extra.columns=TRUE)
#Granges formate is useful because it has start and stop in one column

#formating into a database 
TxDB_escsynttad.db<-makeTxDbFromGRanges(GRanges_escgff_ex)

#open peak file and will format to Granges object
#wrote a python scrip IDR_renames to rename lachesis_group to chr
Stage20=readNarrowPeak("/Users/pui/Documents/Lab_Vienna/ATAC/IDR/stage_20_idr.dms.short")
Stage25=readNarrowPeak("/Users/pui/Documents/Lab_Vienna/ATAC/IDR/stage_25_idr.dms.short")
Stage29=readNarrowPeak("/Users/pui/Documents/Lab_Vienna/ATAC/IDR/stage_29_idr.dms.short")

#promoter region should be bigger than the default, default is 3000
#10kb promoter region 
Stage20_peak <- annotatePeak(Stage20,tssRegion = c(-10000,10000), TxDb = TxDB_escsynttad.db)
Stage25_peak <- annotatePeak(Stage25,tssRegion = c(-10000,10000), TxDb = TxDB_escsynttad.db)
Stage29_peak <- annotatePeak(Stage29,tssRegion = c(-10000,10000), TxDb = TxDB_escsynttad.db)

#Pie plot with peaks annotated to different regions of the genome 
plotAnnoPie(Stage20_peak)
plotAnnoPie(Stage25_peak)
plotAnnoPie(Stage29_peak)


#Bar plot with peaks annotated to different regions of the genome 
plotAnnoBar(Stage20_peak)
plotAnnoBar(Stage25_peak)
plotAnnoBar(Stage29_peak)


#Plot in which you can see the peaks that overlap with different regions of the genome 
upsetplot(Stage20_peak)
upsetplot(Stage25_peak)
upsetplot(Stage29_peak)

#make a list in which we have the geneID and to which synteny it belongs to
#from a list to a matrix
Ceph=scan("/Users/pui/Documents/Lab_Vienna/ATAC/Cluster/ceph_list.txt", what=character(1),sep = ",")
Ceph<-as.data.frame(Ceph)
#adds Ceph as additional column 
Ceph$syn_type<-c("Ceph")
colnames(Ceph)<-c("geneId","syn_type")
Ceph

Meta=scan("/Users/pui/Documents/Lab_Vienna/ATAC/Cluster/meta_list.txt",what=character(1),sep=",")
Meta<-as.data.frame(Meta)
Meta$syn_type<-c("Meta")
colnames(Meta)<-c("geneId","syn_type")
Meta

All=scan("/Users/pui/Documents/Lab_Vienna/ATAC/Cluster/allgenes_list.txt",what = character(1), sep = ",")
All<-as.data.frame(All)
All$syn_type=c("non-syntenic")
colnames(All)<-c("geneId","syn_type")
All

#rbind to append Meta to Ceph list 
Synteny<-rbind(Ceph,Meta)
#adding gene ID from all genes, that are not already in the Meta and Ceph list 
#[!duplicated(geneId)] this will skip all genIDs that are alsready in the table 
#no duplicated geneIDs will be in the new synteny list 
Synteny<-rbindlist(list(Synteny,All))[!duplicated(geneId)]
tail(Synteny)

#use merge to add synteny to each GRanges annotations of each stage 
#in GRanges object all the informations are given 
#all.x keeps all the rows
GRanges_Stage20_annotations<-as.GRanges(Stage20_peak)
GRanges_Stage20_annotations<-merge(GRanges_Stage20_annotations,Synteny,by='geneId', all.x=TRUE)
#removes the rest of exon... information which is unneccessary
GRanges_Stage20_annotations$annotation<-sub("\\(.*\\)","",GRanges_Stage20_annotations$annotation)


GRanges_Stage25_annotations<-as.GRanges(Stage25_peak)
GRanges_Stage25_annotations<-merge(GRanges_Stage25_annotations,Synteny,by='geneId',all.x=TRUE)
GRanges_Stage25_annotations$annotation<-sub("\\(.*\\)","",GRanges_Stage25_annotations$annotation)

GRanges_Stage29_annotations<-as.GRanges(Stage29_peak)
GRanges_Stage29_annotations<-merge(GRanges_Stage29_annotations,Synteny,by='geneId',all.x=TRUE)
GRanges_Stage29_annotations$annotation<-sub("\\(.*\\)","",GRanges_Stage29_annotations$annotation)

#only load plyr when needed, problems when dplyr is also loaded 
library(plyr)

#counts the frequency of each annotation to each synteny type 
counts_Stage20<-ddply(GRanges_Stage20_annotations,.(GRanges_Stage20_annotations$syn_type,GRanges_Stage20_annotations$annotation),nrow)
colnames(counts_Stage20)<-c("Synteny","annotation","frequency")
counts_Stage20

counts_Stage25<-ddply(GRanges_Stage25_annotations,.(GRanges_Stage25_annotations$syn_type,GRanges_Stage25_annotations$annotation),nrow)
colnames(counts_Stage25)<-c("Synteny","annotation","frequency")
counts_Stage25

counts_Stage29<-ddply(GRanges_Stage29_annotations,.(GRanges_Stage29_annotations$syn_type,GRanges_Stage29_annotations$annotation),nrow)
colnames(counts_Stage29)<-c("Synteny", "annotation","frequency")
counts_Stage29

# to unload plyr package: detach(package:"plyr",unload=TRUE)
unloadNamespace("plyr")
library("dplyr")
#group.by generates a table by synteny types  
#mutates adds a new variable and preserves existing ones 
counts_Stage20<-counts_Stage20%>%
  group_by(Synteny)%>%
  mutate(percentage= frequency/sum(frequency))

counts_Stage25<-counts_Stage25%>%
  group_by(Synteny)%>%
  mutate(percentage=frequency/sum(frequency))

counts_Stage29<-counts_Stage29%>%
  group_by(Synteny)%>%
  mutate(percentage=frequency/sum(frequency))

#Grey colour pallette for the plots
#colors<-brewer.pal(7,"Greys"), "#000000"



#names(group.colors) <- levels(factor(c(levels(counts_Stage20$annotation), levels(as.factor(counts_Stage25$annotation)))) # Extract all levels of both data
#my_scale <- scale_fill_manual(name = annotation, values = colors)               # Specify scale_fill_manual
#Annotations of peaks that are present in both replicates 
#distribution of annotations within the different clusters 
#counts_Stage20$annotation= as.factor(counts_Stage20$annotation)

theme_templ=   theme(
  # Hide panel borders and remove grid lines
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Change axis line
  axis.line = element_line(colour = "black"),
  axis.text.x = element_text(colour = "black", size = 20),
  axis.text.y = element_text(colour = "black", size= 20, margin = margin(t = 3, r = 50, b = 0, l = 0)),
  legend.title = element_blank(),
  text = element_text(size = 20))

group.colors<-setNames(c('#a84848','#306078', '#78a890',"#9370DB","#7890a8","#a8a8c0",'#c0c0c0'), levels(as.factor(counts_Stage25$annotation)))  


counts_Stage20_plot <- ggplot(counts_Stage20, aes(x=Synteny,y=percentage,fill=as.factor(annotation))) +   
  geom_bar(position= "dodge", stat="identity",width = 0.9)+
  #geom_col(width = 0.8) +
  #scale_fill_n +
  xlab("\nSynteny type") + 
  scale_fill_manual(values = group.colors) +
  theme_bw()+ #black and white theme 
  theme_templ +
  scale_y_continuous(name='Read annotation \n to genomic regions (%)\n ', expand = c(0,0),
                    limits = c(0,1)) +
  ggtitle("Stage 20")

counts_Stage20_plot


counts_Stage25_plot<-ggplot(counts_Stage25, aes(x=Synteny,y=percentage,fill=as.factor(annotation))) +   
  geom_bar(position= "dodge", stat="identity",width = 0.9)+
  #geom_col(width = 0.8) +
  #scale_fill_n +
  xlab("\nSynteny type") + 
  scale_fill_manual(values = group.colors) +
  theme_bw()+ #black and white theme 
  theme_templ +
  scale_y_continuous(name='Read annotation \n to genomic regions (%)\n ', expand = c(0,0),
                     limits = c(0,1)) +
  ggtitle("Stage 25")

counts_Stage25_plot


counts_Stage29_plot<-ggplot(counts_Stage29, aes(x=Synteny,y=percentage,fill=as.factor(annotation))) +   
  geom_bar(position= "dodge", stat="identity",width = 0.9)+
  #geom_col(width = 0.8) +
  #scale_fill_n +
  xlab("\nSynteny type") + 
  scale_fill_manual(values = group.colors) +
  theme_bw()+ #black and white theme 
  theme_templ +
  scale_y_continuous(name='Read annotation \n to genomic regions (%)\n ', expand = c(0,0),
                     limits = c(0,1)) +
  ggtitle("Stage 29")

counts_Stage29_plot


ggsave(file="/Users/pui/Documents/Lab_Vienna/ATAC/IDR/figures/Genrich_Stage20_ClusterAnnotation.svg", plot=counts_Stage20_plot, width = 9, height = 6)
ggsave(file="/Users/pui/Documents/Lab_Vienna/ATAC/IDR/figures/Genrich_Stage25_ClusterAnnotation.svg", plot=counts_Stage25_plot, width = 9, height = 6)
ggsave(file="/Users/pui/Documents/Lab_Vienna/ATAC/IDR/figures/Genrich_Stage29_ClusterAnnotation.svg", plot=counts_Stage29_plot, width = 9, height = 6)

#to test significance in the differences between the synteny classes we test with Fishers test 
#Long format, wide format
#we need to convert the long format to wide formate 
counts_Stage20_wide<-spread(counts_Stage20[,1:3],annotation,frequency)
#replaces NA with a 0
counts_Stage20_wide[is.na(counts_Stage20_wide)]<-0
counts_Stage20_wide2<-as.data.frame(counts_Stage20_wide[,2:7])
row.names(counts_Stage20_wide2)<-counts_Stage20_wide$Synteny

----
#To test signifance we use the Fishers test, with this we can only compare two different parameters
#intron promoter
fisher.test(counts_Stage20_wide2[1:2,c(5,6)])

Distal_Intron<-counts_Stage20_wide2[1:2,c(2,5)]
#Transposase to check with different matrix 
fish<-fisher.test(t(Distal_Intron))
chisq.test(Distal_Intron)
#Nullhypothesis proportion between the variables are equal 
Meta_nonsyn<-counts_Stage20_wide2[2:3,c(3,6)]
fisher.test(Meta_nonsyn)
chisq.test(Meta_nonsyn)

Ceph_nonsyn<-counts_Stage20_wide2[c(1,3),c(3,6)]
fisher.test(Ceph_nonsyn)
chisq.test(Ceph_nonsyn)
-------
  
  #make tables for all stages between the three different synteny types 
  
  #Stage 20 
  #convert counts table to a wide formate, use only first three columns Synteny, annotation and frequency, the percentage is not needed for fishers test 
  #spread by annotation and frequency
  counts_Stage20_wide<-spread(counts_Stage20[,1:3],annotation,frequency)
#replaces NA with a 0
counts_Stage20_wide[is.na(counts_Stage20_wide)]<-0
counts_Stage20_wide2<-as.data.frame(counts_Stage20_wide[,2:7])
row.names(counts_Stage20_wide2)<-counts_Stage20_wide$Synteny
counts_Stage20_wide2

#To test signifance we use the Fishers test, with this we can only compare two different parameters
#we will only test intron, distal because of biggest difference in graph 
#Stage 20: Ceph to Meta
UTR5_Distal<-counts_Stage20_wide2[1:2,c(1,2)]
UTR5_Downstream<-counts_Stage20_wide2[1:2,c(1,3)]
UTR5_Exon<-counts_Stage20_wide2[1:2,c(1,4)]
UTR5_Intron<-counts_Stage20_wide2[1:2,c(1,5)]
UTR5_Promoter<-counts_Stage20_wide2[1:2,c(1,6)]
Distal_Downstream<-counts_Stage20_wide2[1:2,c(2,3)]
Distal_Exon<-counts_Stage20_wide2[1:2,c(2,4)]
Distal_Intron<-counts_Stage20_wide2[1:2,c(2,5)]
Distal_Promoter<-counts_Stage20_wide2[1:2,c(2,6)]
Downstream_Exon<-counts_Stage20_wide2[1:2,c(3,4)]
Downstream_Intron<-counts_Stage20_wide2[1:2,c(3,5)]
Downstream_Promoter<-counts_Stage20_wide2[1:2,c(3,6)]
Exon_Intron<-counts_Stage20_wide2[1:2,c(4,5)]
Exon_Promoter<-counts_Stage20_wide2[1:2,c(4,6)]
Intron_Promoter<-counts_Stage20_wide2[1:2,c(5,6)]

Stage20_CephMeta_fisher<-cbind(c(NA,fisher.test(UTR5_Distal)$p.value,fisher.test(UTR5_Downstream)$p.value,fisher.test(UTR5_Exon)$p.value,fisher.test(UTR5_Intron)$p.value,fisher.test(UTR5_Promoter)$p.value),
                               c(fisher.test(UTR5_Distal)$p.value,NA,fisher.test(Distal_Downstream)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Distal_Promoter)$p.value),
                               c(fisher.test(UTR5_Downstream)$p.value,fisher.test(Distal_Downstream)$p.value,NA,fisher.test(Downstream_Exon)$p.value,fisher.test(Downstream_Intron)$p.value,fisher.test(Downstream_Promoter)$p.value),
                               c(fisher.test(UTR5_Exon)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Downstream_Exon)$p.value, NA,fisher.test(Exon_Intron)$p.value,fisher.test(Exon_Promoter)$p.value),
                               c(fisher.test(UTR5_Intron)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Downstream_Intron)$p.value,fisher.test(Exon_Intron)$p.value,NA,fisher.test(Intron_Promoter)$p.value),
                               c(fisher.test(UTR5_Promoter)$p.value,fisher.test(Distal_Promoter)$p.value,fisher.test(Downstream_Promoter)$p.value,fisher.test(Exon_Promoter)$p.value,fisher.test(Intron_Promoter)$p.value,NA))
#format to e number so the table is smaller in size
Stage20_CephMeta_fisher<-as.data.frame(round(Stage20_CephMeta_fisher,digits = 4))
colnames(Stage20_CephMeta_fisher)<-c("5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
rownames(Stage20_CephMeta_fisher)<-c("5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
Stage20_CephMeta_fisher

#Stage20: Meta and non-syntenic
counts_Stage20_wide2
UTR5_Distal<-counts_Stage20_wide2[2:3,c(1,2)]
UTR5_Downstream<-counts_Stage20_wide2[2:3,c(1,3)]
UTR5_Exon<-counts_Stage20_wide2[2:3,c(1,4)]
UTR5_Intron<-counts_Stage20_wide2[2:3,c(1,5)]
UTR5_Promoter<-counts_Stage20_wide2[2:3,c(1,6)]
Distal_Downstream<-counts_Stage20_wide2[2:3,c(2,3)]
Distal_Exon<-counts_Stage20_wide2[2:3,c(2,4)]
Distal_Intron<-counts_Stage20_wide2[2:3,c(2,5)]
Distal_Promoter<-counts_Stage20_wide2[2:3,c(2,6)]
Downstream_Exon<-counts_Stage20_wide2[2:3,c(3,4)]
Downstream_Intron<-counts_Stage20_wide2[2:3,c(3,5)]
Downstream_Promoter<-counts_Stage20_wide2[2:3,c(3,6)]
Exon_Intron<-counts_Stage20_wide2[2:3,c(4,5)]
Exon_Promoter<-counts_Stage20_wide2[2:3,c(4,6)]
Intron_Promoter<-counts_Stage20_wide2[2:3,c(5,6)]


Stage20_MetaNonSyn_fisher<-cbind(c(NA,fisher.test(UTR5_Distal)$p.value,fisher.test(UTR5_Downstream)$p.value,fisher.test(UTR5_Exon)$p.value,fisher.test(UTR5_Intron)$p.value,fisher.test(UTR5_Promoter)$p.value),
                                 c(fisher.test(UTR5_Distal)$p.value,NA,fisher.test(Distal_Downstream)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Distal_Promoter)$p.value),
                                 c(fisher.test(UTR5_Downstream)$p.value,fisher.test(Distal_Downstream)$p.value,NA,fisher.test(Downstream_Exon)$p.value,fisher.test(Downstream_Intron)$p.value,fisher.test(Downstream_Promoter)$p.value),
                                 c(fisher.test(UTR5_Exon)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Downstream_Exon)$p.value, NA,fisher.test(Exon_Intron)$p.value,fisher.test(Exon_Promoter)$p.value),
                                 c(fisher.test(UTR5_Intron)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Downstream_Intron)$p.value,fisher.test(Exon_Intron)$p.value,NA,fisher.test(Intron_Promoter)$p.value),
                                 c(fisher.test(UTR5_Promoter)$p.value,fisher.test(Distal_Promoter)$p.value,fisher.test(Downstream_Promoter)$p.value,fisher.test(Exon_Promoter)$p.value,fisher.test(Intron_Promoter)$p.value,NA))
#format to e number so the table is smaller in size
Stage20_MetaNonSyn_fisher<-as.data.frame(round(Stage20_MetaNonSyn_fisher, digits=4))
colnames(Stage20_MetaNonSyn_fisher)<-c("5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
rownames(Stage20_MetaNonSyn_fisher)<-c("5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
Stage20_MetaNonSyn_fisher

#Stage20: Ceph to non-syntenic
counts_Stage20_wide2
UTR5_Distal<-counts_Stage20_wide2[c(1,3),c(1,2)]
UTR5_Downstream<-counts_Stage20_wide2[c(1,3),c(1,3)]
UTR5_Exon<-counts_Stage20_wide2[c(1,3),c(1,4)]
UTR5_Intron<-counts_Stage20_wide2[c(1,3),c(1,5)]
UTR5_Promoter<-counts_Stage20_wide2[c(1,3),c(1,6)]
Distal_Downstream<-counts_Stage20_wide2[c(1,3),c(2,3)]
Distal_Exon<-counts_Stage20_wide2[c(1,3),c(2,4)]
Distal_Intron<-counts_Stage20_wide2[c(1,3),c(2,5)]
Distal_Promoter<-counts_Stage20_wide2[c(1,3),c(2,6)]
Downstream_Exon<-counts_Stage20_wide2[c(1,3),c(3,4)]
Downstream_Intron<-counts_Stage20_wide2[c(1,3),c(3,5)]
Downstream_Promoter<-counts_Stage20_wide2[c(1,3),c(3,6)]
Exon_Intron<-counts_Stage20_wide2[c(1,3),c(4,5)]
Exon_Promoter<-counts_Stage20_wide2[c(1,3),c(4,6)]
Intron_Promoter<-counts_Stage20_wide2[c(1,3),c(5,6)]

Stage20_CephNonSyn_fisher<-cbind(c(NA,fisher.test(UTR5_Distal)$p.value,fisher.test(UTR5_Downstream)$p.value,fisher.test(UTR5_Exon)$p.value,fisher.test(UTR5_Intron)$p.value,fisher.test(UTR5_Promoter)$p.value),
                                 c(fisher.test(UTR5_Distal)$p.value,NA,fisher.test(Distal_Downstream)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Distal_Promoter)$p.value),
                                 c(fisher.test(UTR5_Downstream)$p.value,fisher.test(Distal_Downstream)$p.value,NA,fisher.test(Downstream_Exon)$p.value,fisher.test(Downstream_Intron)$p.value,fisher.test(Downstream_Promoter)$p.value),
                                 c(fisher.test(UTR5_Exon)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Downstream_Exon)$p.value, NA,fisher.test(Exon_Intron)$p.value,fisher.test(Exon_Promoter)$p.value),
                                 c(fisher.test(UTR5_Intron)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Downstream_Intron)$p.value,fisher.test(Exon_Intron)$p.value,NA,fisher.test(Intron_Promoter)$p.value),
                                 c(fisher.test(UTR5_Promoter)$p.value,fisher.test(Distal_Promoter)$p.value,fisher.test(Downstream_Promoter)$p.value,fisher.test(Exon_Promoter)$p.value,fisher.test(Intron_Promoter)$p.value,NA))
#format to e number so the table is smaller in size
Stage20_CephNonSyn_fisher<-as.data.frame(round(Stage20_CephNonSyn_fisher, digits=4))
colnames(Stage20_CephNonSyn_fisher)<-c("5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
rownames(Stage20_CephNonSyn_fisher)<-c("5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
Stage20_CephNonSyn_fisher

#Stage 25
#get the wide format with only frequencies and sort by annotations
counts_Stage25_wide<-spread(counts_Stage25[,1:3],annotation,frequency)
#replacing NA with a 0 
counts_Stage25_wide[is.na(counts_Stage25_wide)]<-0
#to remove first column 
counts_Stage25_wide2<-as.data.frame(counts_Stage25_wide[,2:8])
#to add rownames of column Synteny of counts_Stage20_wide table 
row.names(counts_Stage25_wide2)<-counts_Stage25_wide$Synteny
counts_Stage25_wide2

#Stage 25:Ceph to Meta 
UTR3_UTR5<-counts_Stage25_wide2[1:2,c(1,2)]
UTR3_Distal<-counts_Stage25_wide2[1:2,c(1,3)]
UTR3_Downstream<-counts_Stage25_wide2[1:2,c(1,4)]
UTR3_Exon<-counts_Stage25_wide2[1:2,c(1,5)]
UTR3_Intron<-counts_Stage25_wide2[1:2,c(1,6)]
UTR3_Promoter<-counts_Stage25_wide2[1:2,c(1,7)]
UTR5_Distal<-counts_Stage25_wide2[1:2,c(2,3)]
UTR5_Downstream<-counts_Stage25_wide2[1:2,c(2,4)]
UTR5_Exon<-counts_Stage25_wide2[1:2,c(2,5)]
UTR5_Intron<-counts_Stage25_wide2[1:2,c(2,6)]
UTR5_Promoter<-counts_Stage25_wide2[1:2,c(2,7)]
Distal_Downstream<-counts_Stage25_wide2[1:2,c(3,4)]
Distal_Exon<-counts_Stage25_wide2[1:2,c(3,5)]
Distal_Intron<-counts_Stage25_wide2[1:2,c(3,6)]
Distal_Promoter<-counts_Stage25_wide2[1:2,c(3,7)]
Downstream_Exon<-counts_Stage25_wide2[1:2,c(4,5)]
Downstream_Intron<-counts_Stage25_wide2[1:2,c(4,6)]
Downstream_Promoter<-counts_Stage25_wide2[1:2,c(4,7)]
Exon_Intron<-counts_Stage25_wide2[1:2,c(5,6)]
Exon_Promoter<-counts_Stage25_wide2[1:2,c(5,7)]
Intron_Promoter<-counts_Stage25_wide2[1:2,c(6,7)]


#fisher test 
Stage25_CephMeta_fisher<-cbind(c(NA,fisher.test(UTR3_UTR5)$p.value,fisher.test(UTR3_Distal)$p.value,fisher.test(UTR3_Downstream)$p.value,fisher.test(UTR3_Exon)$p.value,fisher.test(UTR3_Intron)$p.value,fisher.test(UTR3_Promoter)$p.value),
                               c(fisher.test(UTR3_UTR5)$p.value,NA,fisher.test(UTR5_Distal)$p.value,fisher.test(UTR5_Downstream)$p.value,fisher.test(UTR5_Exon)$p.value,fisher.test(UTR5_Intron)$p.value,fisher.test(UTR5_Promoter)$p.value),
                               c(fisher.test(UTR3_Distal)$p.value,fisher.test(UTR5_Distal)$p.value,NA,fisher.test(Distal_Downstream)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Distal_Promoter)$p.value),
                               c(fisher.test(UTR3_Downstream)$p.value,fisher.test(UTR5_Downstream)$p.value,fisher.test(Distal_Downstream)$p.value,NA,fisher.test(Downstream_Exon)$p.value,fisher.test(Downstream_Intron)$p.value,fisher.test(Downstream_Promoter)$p.value),
                               c(fisher.test(UTR3_Exon)$p.value,fisher.test(UTR5_Exon)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Downstream_Exon)$p.value,NA,fisher.test(Exon_Intron)$p.value,fisher.test(Exon_Promoter)$p.value),
                               c(fisher.test(UTR3_Intron)$p.value,fisher.test(UTR5_Intron)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Downstream_Intron)$p.value, fisher.test(Exon_Intron)$p.value,NA,fisher.test(Intron_Promoter)$p.value),
                               c(fisher.test(UTR3_Promoter)$p.value,fisher.test(UTR5_Promoter)$p.value,fisher.test(Distal_Promoter)$p.value,fisher.test(Downstream_Promoter)$p.value,fisher.test(Exon_Promoter)$p.value,fisher.test(Intron_Promoter)$p.value,NA))
#format to e number so the table is smaller in size
Stage25_CephMeta_fisher<-as.data.frame(round(Stage25_CephMeta_fisher, digits=4))
colnames(Stage25_CephMeta_fisher)<-c("3'UTR","5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
rownames(Stage25_CephMeta_fisher)<-c("3'UTR","5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
Stage25_CephMeta_fisher

#Stage 25: Meta and non-syntenic
counts_Stage25_wide2
UTR3_UTR5<-counts_Stage25_wide2[2:3,c(1,2)]
UTR3_Distal<-counts_Stage25_wide2[2:3,c(1,3)]
UTR3_Downstream<-counts_Stage25_wide2[2:3,c(1,4)]
UTR3_Exon<-counts_Stage25_wide2[2:3,c(1,5)]
UTR3_Intron<-counts_Stage25_wide2[2:3,c(1,6)]
UTR3_Promoter<-counts_Stage25_wide2[2:3,c(1,7)]
UTR5_Distal<-counts_Stage25_wide2[2:3,c(2,3)]
UTR5_Downstream<-counts_Stage25_wide2[2:3,c(2,4)]
UTR5_Exon<-counts_Stage25_wide2[2:3,c(2,5)]
UTR5_Intron<-counts_Stage25_wide2[2:3,c(2,6)]
UTR5_Promoter<-counts_Stage25_wide2[2:3,c(2,7)]
Distal_Downstream<-counts_Stage25_wide2[2:3,c(3,4)]
Distal_Exon<-counts_Stage25_wide2[2:3,c(3,5)]
Distal_Intron<-counts_Stage25_wide2[2:3,c(3,6)]
Distal_Promoter<-counts_Stage25_wide2[2:3,c(3,7)]
Downstream_Exon<-counts_Stage25_wide2[2:3,c(4,5)]
Downstream_Intron<-counts_Stage25_wide2[2:3,c(4,6)]
Downstream_Promoter<-counts_Stage25_wide2[2:3,c(4,7)]
Exon_Intron<-counts_Stage25_wide2[2:3,c(5,6)]
Exon_Promoter<-counts_Stage25_wide2[2:3,c(5,7)]
Intron_Promoter<-counts_Stage25_wide2[2:3,c(6,7)]

Stage25_MetaNonSyn_fisher<-cbind(c(NA,fisher.test(UTR3_UTR5)$p.value,fisher.test(UTR3_Distal)$p.value,fisher.test(UTR3_Downstream)$p.value,fisher.test(UTR3_Exon)$p.value,fisher.test(UTR3_Intron)$p.value,fisher.test(UTR3_Promoter)$p.value),
                                 c(fisher.test(UTR3_UTR5)$p.value,NA,fisher.test(UTR5_Distal)$p.value,fisher.test(UTR5_Downstream)$p.value,fisher.test(UTR5_Exon)$p.value,fisher.test(UTR5_Intron)$p.value,fisher.test(UTR5_Promoter)$p.value),
                                 c(fisher.test(UTR3_Distal)$p.value,fisher.test(UTR5_Distal)$p.value,NA,fisher.test(Distal_Downstream)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Distal_Promoter)$p.value),
                                 c(fisher.test(UTR3_Downstream)$p.value,fisher.test(UTR5_Downstream)$p.value,fisher.test(Distal_Downstream)$p.value,NA,fisher.test(Downstream_Exon)$p.value,fisher.test(Downstream_Intron)$p.value,fisher.test(Downstream_Promoter)$p.value),
                                 c(fisher.test(UTR3_Exon)$p.value,fisher.test(UTR5_Exon)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Downstream_Exon)$p.value,NA,fisher.test(Exon_Intron)$p.value,fisher.test(Exon_Promoter)$p.value),
                                 c(fisher.test(UTR3_Intron)$p.value,fisher.test(UTR5_Intron)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Downstream_Intron)$p.value, fisher.test(Exon_Intron)$p.value,NA,fisher.test(Intron_Promoter)$p.value),
                                 c(fisher.test(UTR3_Promoter)$p.value,fisher.test(UTR5_Promoter)$p.value,fisher.test(Distal_Promoter)$p.value,fisher.test(Downstream_Promoter)$p.value,fisher.test(Exon_Promoter)$p.value,fisher.test(Intron_Promoter)$p.value,NA))
#format to e number so the table is smaller in size
Stage25_MetaNonSyn_fisher<-as.data.frame(round(Stage25_MetaNonSyn_fisher, digits=4))
colnames(Stage25_MetaNonSyn_fisher)<-c("3'UTR","5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
rownames(Stage25_MetaNonSyn_fisher)<-c("3'UTR","5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
Stage25_MetaNonSyn_fisher

#Stage25: Ceph to non-syntenic
counts_Stage25_wide2
UTR3_UTR5<-counts_Stage25_wide2[c(1,3),c(1,2)]
UTR3_Distal<-counts_Stage25_wide2[c(1,3),c(1,3)]
UTR3_Downstream<-counts_Stage25_wide2[c(1,3),c(1,4)]
UTR3_Exon<-counts_Stage25_wide2[c(1,3),c(1,5)]
UTR3_Intron<-counts_Stage25_wide2[c(1,3),c(1,6)]
UTR3_Promoter<-counts_Stage25_wide2[c(1,3),c(1,7)]
UTR5_Distal<-counts_Stage25_wide2[c(1,3),c(2,3)]
UTR5_Downstream<-counts_Stage25_wide2[c(1,3),c(2,4)]
UTR5_Exon<-counts_Stage25_wide2[c(1,3),c(2,5)]
UTR5_Intron<-counts_Stage25_wide2[c(1,3),c(2,6)]
UTR5_Promoter<-counts_Stage25_wide2[c(1,3),c(2,7)]
Distal_Downstream<-counts_Stage25_wide2[c(1,3),c(3,4)]
Distal_Exon<-counts_Stage25_wide2[c(1,3),c(3,5)]
Distal_Intron<-counts_Stage25_wide2[c(1,3),c(3,6)]
Distal_Promoter<-counts_Stage25_wide2[c(1,3),c(3,7)]
Downstream_Exon<-counts_Stage25_wide2[c(1,3),c(4,5)]
Downstream_Intron<-counts_Stage25_wide2[c(1,3),c(4,6)]
Downstream_Promoter<-counts_Stage25_wide2[c(1,3),c(4,7)]
Exon_Intron<-counts_Stage25_wide2[c(1,3),c(5,6)]
Exon_Promoter<-counts_Stage25_wide2[c(1,3),c(5,7)]
Intron_Promoter<-counts_Stage25_wide2[c(1,3),c(6,7)]

Stage25_CephNonSyn_fisher<-cbind(c(NA,fisher.test(UTR3_UTR5)$p.value,fisher.test(UTR3_Distal)$p.value,fisher.test(UTR3_Downstream)$p.value,fisher.test(UTR3_Exon)$p.value,fisher.test(UTR3_Intron)$p.value,fisher.test(UTR3_Promoter)$p.value),
                                 c(fisher.test(UTR3_UTR5)$p.value,NA,fisher.test(UTR5_Distal)$p.value,fisher.test(UTR5_Downstream)$p.value,fisher.test(UTR5_Exon)$p.value,fisher.test(UTR5_Intron)$p.value,fisher.test(UTR5_Promoter)$p.value),
                                 c(fisher.test(UTR3_Distal)$p.value,fisher.test(UTR5_Distal)$p.value,NA,fisher.test(Distal_Downstream)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Distal_Promoter)$p.value),
                                 c(fisher.test(UTR3_Downstream)$p.value,fisher.test(UTR5_Downstream)$p.value,fisher.test(Distal_Downstream)$p.value,NA,fisher.test(Downstream_Exon)$p.value,fisher.test(Downstream_Intron)$p.value,fisher.test(Downstream_Promoter)$p.value),
                                 c(fisher.test(UTR3_Exon)$p.value,fisher.test(UTR5_Exon)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Downstream_Exon)$p.value,NA,fisher.test(Exon_Intron)$p.value,fisher.test(Exon_Promoter)$p.value),
                                 c(fisher.test(UTR3_Intron)$p.value,fisher.test(UTR5_Intron)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Downstream_Intron)$p.value, fisher.test(Exon_Intron)$p.value,NA,fisher.test(Intron_Promoter)$p.value),
                                 c(fisher.test(UTR3_Promoter)$p.value,fisher.test(UTR5_Promoter)$p.value,fisher.test(Distal_Promoter)$p.value,fisher.test(Downstream_Promoter)$p.value,fisher.test(Exon_Promoter)$p.value,fisher.test(Intron_Promoter)$p.value,NA))
#format to e number so the table is smaller in size
Stage25_CephNonSyn_fisher<-as.data.frame(round(Stage25_CephNonSyn_fisher, digits=4))
colnames(Stage25_CephNonSyn_fisher)<-c("3'UTR","5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
rownames(Stage25_CephNonSyn_fisher)<-c("3'UTR","5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
Stage25_CephNonSyn_fisher

#Stage 29 
#generates a wide formate by annotation and frequency with the first 3 columns 
counts_Stage29_wide<-spread(counts_Stage29[,1:3],annotation,frequency)
#replace NA with 0 
counts_Stage29_wide[is.na(counts_Stage29_wide)]<-0
#generates a data frame from column 2 to eight 
counts_Stage29_wide2<-as.data.frame(counts_Stage29_wide[,2:8])
row.names(counts_Stage29_wide2)<-counts_Stage29_wide$Synteny
counts_Stage29_wide2

#Stage 29: Meta to Ceph 
UTR3_UTR5<-counts_Stage29_wide2[1:2,c(1,2)]
UTR3_Distal<-counts_Stage29_wide2[1:2,c(1,3)]
UTR3_Downstream<-counts_Stage29_wide2[1:2,c(1,4)]
UTR3_Exon<-counts_Stage29_wide2[1:2,c(1,5)]
UTR3_Intron<-counts_Stage29_wide2[1:2,c(1,6)]
UTR3_Promoter<-counts_Stage29_wide2[1:2,c(1,7)]
UTR5_Distal<-counts_Stage29_wide2[1:2,c(2,3)]
UTR5_Downstream<-counts_Stage29_wide2[1:2,c(2,4)]
UTR5_Exon<-counts_Stage29_wide2[1:2,c(2,5)]
UTR5_Intron<-counts_Stage29_wide2[1:2,c(2,6)]
UTR5_Promoter<-counts_Stage29_wide2[1:2,c(2,7)]
Distal_Downstream<-counts_Stage29_wide2[1:2,c(3,4)]
Distal_Exon<-counts_Stage29_wide2[1:2,c(3,5)]
Distal_Intron<-counts_Stage29_wide2[1:2,c(3,6)]
Distal_Promoter<-counts_Stage29_wide2[1:2,c(3,7)]
Downstream_Exon<-counts_Stage29_wide2[1:2,c(4,5)]
Downstream_Intron<-counts_Stage29_wide2[1:2,c(4,6)]
Downstream_Promoter<-counts_Stage29_wide2[1:2,c(4,7)]
Exon_Intron<-counts_Stage29_wide2[1:2,c(5,6)]
Exon_Promoter<-counts_Stage29_wide2[1:2,c(5,7)]
Intron_Promoter<-counts_Stage29_wide2[1:2,c(6,7)]
#fisher test 
Stage29_CephMeta_fisher<-cbind(c(NA,fisher.test(UTR3_UTR5)$p.value,fisher.test(UTR3_Distal)$p.value,fisher.test(UTR3_Downstream)$p.value,fisher.test(UTR3_Exon)$p.value,fisher.test(UTR3_Intron)$p.value,fisher.test(UTR3_Promoter)$p.value),
                               c(fisher.test(UTR3_UTR5)$p.value,NA,fisher.test(UTR5_Distal)$p.value,fisher.test(UTR5_Downstream)$p.value,fisher.test(UTR5_Exon)$p.value,fisher.test(UTR5_Intron)$p.value,fisher.test(UTR5_Promoter)$p.value),
                               c(fisher.test(UTR3_Distal)$p.value,fisher.test(UTR5_Distal)$p.value,NA,fisher.test(Distal_Downstream)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Distal_Promoter)$p.value),
                               c(fisher.test(UTR3_Downstream)$p.value,fisher.test(UTR5_Downstream)$p.value,fisher.test(Distal_Downstream)$p.value,NA,fisher.test(Downstream_Exon)$p.value,fisher.test(Downstream_Intron)$p.value,fisher.test(Downstream_Promoter)$p.value),
                               c(fisher.test(UTR3_Exon)$p.value,fisher.test(UTR5_Exon)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Downstream_Exon)$p.value,NA,fisher.test(Exon_Intron)$p.value,fisher.test(Exon_Promoter)$p.value),
                               c(fisher.test(UTR3_Intron)$p.value,fisher.test(UTR5_Intron)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Downstream_Intron)$p.value, fisher.test(Exon_Intron)$p.value,NA,fisher.test(Intron_Promoter)$p.value),
                               c(fisher.test(UTR3_Promoter)$p.value,fisher.test(UTR5_Promoter)$p.value,fisher.test(Distal_Promoter)$p.value,fisher.test(Downstream_Promoter)$p.value,fisher.test(Exon_Promoter)$p.value,fisher.test(Intron_Promoter)$p.value,NA))
#format to e number so the table is smaller in size
Stage29_CephMeta_fisher<-as.data.frame(round(Stage29_CephMeta_fisher, digits=4))
colnames(Stage29_CephMeta_fisher)<-c("3'UTR","5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
rownames(Stage29_CephMeta_fisher)<-c("3'UTR","5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
Stage29_CephMeta_fisher
#Stage 25: Meta and non-syntenic
counts_Stage29_wide2
UTR3_UTR5<-counts_Stage29_wide2[2:3,c(1,2)]
UTR3_Distal<-counts_Stage29_wide2[2:3,c(1,3)]
UTR3_Downstream<-counts_Stage29_wide2[2:3,c(1,4)]
UTR3_Exon<-counts_Stage29_wide2[2:3,c(1,5)]
UTR3_Intron<-counts_Stage29_wide2[2:3,c(1,6)]
UTR3_Promoter<-counts_Stage29_wide2[2:3,c(1,7)]
UTR5_Distal<-counts_Stage29_wide2[2:3,c(2,3)]
UTR5_Downstream<-counts_Stage29_wide2[2:3,c(2,4)]
UTR5_Exon<-counts_Stage29_wide2[2:3,c(2,5)]
UTR5_Intron<-counts_Stage29_wide2[2:3,c(2,6)]
UTR5_Promoter<-counts_Stage29_wide2[2:3,c(2,7)]
Distal_Downstream<-counts_Stage29_wide2[2:3,c(3,4)]
Distal_Exon<-counts_Stage29_wide2[2:3,c(3,5)]
Distal_Intron<-counts_Stage29_wide2[2:3,c(3,6)]
Distal_Promoter<-counts_Stage29_wide2[2:3,c(3,7)]
Downstream_Exon<-counts_Stage29_wide2[2:3,c(4,5)]
Downstream_Intron<-counts_Stage29_wide2[2:3,c(4,6)]
Downstream_Promoter<-counts_Stage29_wide2[2:3,c(4,7)]
Exon_Intron<-counts_Stage29_wide2[2:3,c(5,6)]
Exon_Promoter<-counts_Stage29_wide2[2:3,c(5,7)]
Intron_Promoter<-counts_Stage29_wide2[2:3,c(6,7)]

Stage29_MetaNonSyn_fisher<-cbind(c(NA,fisher.test(UTR3_UTR5)$p.value,fisher.test(UTR3_Distal)$p.value,fisher.test(UTR3_Downstream)$p.value,fisher.test(UTR3_Exon)$p.value,fisher.test(UTR3_Intron)$p.value,fisher.test(UTR3_Promoter)$p.value),
                                 c(fisher.test(UTR3_UTR5)$p.value,NA,fisher.test(UTR5_Distal)$p.value,fisher.test(UTR5_Downstream)$p.value,fisher.test(UTR5_Exon)$p.value,fisher.test(UTR5_Intron)$p.value,fisher.test(UTR5_Promoter)$p.value),
                                 c(fisher.test(UTR3_Distal)$p.value,fisher.test(UTR5_Distal)$p.value,NA,fisher.test(Distal_Downstream)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Distal_Promoter)$p.value),
                                 c(fisher.test(UTR3_Downstream)$p.value,fisher.test(UTR5_Downstream)$p.value,fisher.test(Distal_Downstream)$p.value,NA,fisher.test(Downstream_Exon)$p.value,fisher.test(Downstream_Intron)$p.value,fisher.test(Downstream_Promoter)$p.value),
                                 c(fisher.test(UTR3_Exon)$p.value,fisher.test(UTR5_Exon)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Downstream_Exon)$p.value,NA,fisher.test(Exon_Intron)$p.value,fisher.test(Exon_Promoter)$p.value),
                                 c(fisher.test(UTR3_Intron)$p.value,fisher.test(UTR5_Intron)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Downstream_Intron)$p.value, fisher.test(Exon_Intron)$p.value,NA,fisher.test(Intron_Promoter)$p.value),
                                 c(fisher.test(UTR3_Promoter)$p.value,fisher.test(UTR5_Promoter)$p.value,fisher.test(Distal_Promoter)$p.value,fisher.test(Downstream_Promoter)$p.value,fisher.test(Exon_Promoter)$p.value,fisher.test(Intron_Promoter)$p.value,NA))
#format to e number so the table is smaller in size
Stage29_MetaNonSyn_fisher<-as.data.frame(round(Stage29_MetaNonSyn_fisher,digits=4))
colnames(Stage29_MetaNonSyn_fisher)<-c("3'UTR","5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
rownames(Stage29_MetaNonSyn_fisher)<-c("3'UTR","5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
Stage29_MetaNonSyn_fisher

#Stage29: Ceph to non-syntenic
counts_Stage29_wide2
UTR3_UTR5<-counts_Stage29_wide2[c(1,3),c(1,2)]
UTR3_Distal<-counts_Stage29_wide2[c(1,3),c(1,3)]
UTR3_Downstream<-counts_Stage29_wide2[c(1,3),c(1,4)]
UTR3_Exon<-counts_Stage29_wide2[c(1,3),c(1,5)]
UTR3_Intron<-counts_Stage29_wide2[c(1,3),c(1,6)]
UTR3_Promoter<-counts_Stage29_wide2[c(1,3),c(1,7)]
UTR5_Distal<-counts_Stage29_wide2[c(1,3),c(2,3)]
UTR5_Downstream<-counts_Stage29_wide2[c(1,3),c(2,4)]
UTR5_Exon<-counts_Stage29_wide2[c(1,3),c(2,5)]
UTR5_Intron<-counts_Stage29_wide2[c(1,3),c(2,6)]
UTR5_Promoter<-counts_Stage29_wide2[c(1,3),c(2,7)]
Distal_Downstream<-counts_Stage29_wide2[c(1,3),c(3,4)]
Distal_Exon<-counts_Stage29_wide2[c(1,3),c(3,5)]
Distal_Intron<-counts_Stage29_wide2[c(1,3),c(3,6)]
Distal_Promoter<-counts_Stage29_wide2[c(1,3),c(3,7)]
Downstream_Exon<-counts_Stage29_wide2[c(1,3),c(4,5)]
Downstream_Intron<-counts_Stage29_wide2[c(1,3),c(4,6)]
Downstream_Promoter<-counts_Stage29_wide2[c(1,3),c(4,7)]
Exon_Intron<-counts_Stage29_wide2[c(1,3),c(5,6)]
Exon_Promoter<-counts_Stage29_wide2[c(1,3),c(5,7)]
Intron_Promoter<-counts_Stage29_wide2[c(1,3),c(6,7)]


Stage29_CephNonSyn_fisher<-cbind(c(NA,fisher.test(UTR3_UTR5)$p.value,fisher.test(UTR3_Distal)$p.value,fisher.test(UTR3_Downstream)$p.value,fisher.test(UTR3_Exon)$p.value,fisher.test(UTR3_Intron)$p.value,fisher.test(UTR3_Promoter)$p.value),
                                 c(fisher.test(UTR3_UTR5)$p.value,NA,fisher.test(UTR5_Distal)$p.value,fisher.test(UTR5_Downstream)$p.value,fisher.test(UTR5_Exon)$p.value,fisher.test(UTR5_Intron)$p.value,fisher.test(UTR5_Promoter)$p.value),
                                 c(fisher.test(UTR3_Distal)$p.value,fisher.test(UTR5_Distal)$p.value,NA,fisher.test(Distal_Downstream)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Distal_Promoter)$p.value),
                                 c(fisher.test(UTR3_Downstream)$p.value,fisher.test(UTR5_Downstream)$p.value,fisher.test(Distal_Downstream)$p.value,NA,fisher.test(Downstream_Exon)$p.value,fisher.test(Downstream_Intron)$p.value,fisher.test(Downstream_Promoter)$p.value),
                                 c(fisher.test(UTR3_Exon)$p.value,fisher.test(UTR5_Exon)$p.value,fisher.test(Distal_Exon)$p.value,fisher.test(Downstream_Exon)$p.value,NA,fisher.test(Exon_Intron)$p.value,fisher.test(Exon_Promoter)$p.value),
                                 c(fisher.test(UTR3_Intron)$p.value,fisher.test(UTR5_Intron)$p.value,fisher.test(Distal_Intron)$p.value,fisher.test(Downstream_Intron)$p.value, fisher.test(Exon_Intron)$p.value,NA,fisher.test(Intron_Promoter)$p.value),
                                 c(fisher.test(UTR3_Promoter)$p.value,fisher.test(UTR5_Promoter)$p.value,fisher.test(Distal_Promoter)$p.value,fisher.test(Downstream_Promoter)$p.value,fisher.test(Exon_Promoter)$p.value,fisher.test(Intron_Promoter)$p.value,NA))
#format to e number so the table is smaller in size
Stage29_CephNonSyn_fisher<-as.data.frame(round(Stage29_CephNonSyn_fisher, digits=4))
colnames(Stage29_CephNonSyn_fisher)<-c("3'UTR","5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")
rownames(Stage29_CephNonSyn_fisher)<-c("3'UTR","5'UTR","Distal Intergenic","Downstream","Exon","Intron","Promoter")

datatable(Data, rownames = FALSE) %>%
  formatStyle(columns = numbofStars, 
              background = styleInterval(c(0.05,0.000000000000000001), c("white","magenta","white")))%>%
  formatStyle(columns = test, 
              background = styleInterval(c(0.05,0.000000000000000001), c("white","magenta","white")))

----
#display significance in red
Stage29_new<- formattable(Stage29_CephNonSyn_fisher, list(
  "3'UTR" = formatter("span", style = x ~ ifelse(x <= 0.05, style(color = "red", font.weight = "bold"), NA)),
  "5'UTR" = formatter("span", style = x ~ ifelse(x <= 0.05, style(color = "red", font.weight = "bold"), NA)),
  "Distal Intergenic" = formatter("span", style = x ~ ifelse(x <= 0.05, style(color = "red", font.weight = "bold"), NA)),
  "Downstream" = formatter("span", style = x ~ ifelse(x <= 0.05, style(color = "red", font.weight = "bold"), NA)),
  "Exon" = formatter("span", style = x ~ ifelse(x <= 0.05, style(color = "red", font.weight = "bold"), NA)),
  "Intron" = formatter("span", style = x ~ ifelse(x <= 0.05, style(color = "red", font.weight = "bold"), NA)),
  "Promoter" = formatter("span", style = x ~ ifelse(x <= 0.05, style(color = "red", font.weight = "bold"), NA))
))

library(gridBase)
library(gridExtra)
library(ggplotify)
library(cowplot)

#put tables into one figure:
th<-ttheme_minimal()

#Fishers exact test in a table
Stage20_fishers <- arrangeGrob(tableGrob(Stage29_CephMeta_fisher,theme = th),
                               tableGrob(Stage20_CephNonSyn_fisher,theme = th),
                               tableGrob(Stage20_MetaNonSyn_fisher,theme=th),
                               layout_matrix = cbind(c(1,2,3),c(1,2,3)))
# Add labels to the arranged plots
# transform to a ggplot
Fishers_table_20<-as.ggplot(Stage20_fishers) + 
  draw_plot_label(label = c("Stage 20",
                            "A Fisher's exact test: Cephalopod to Metazoan synteny", 
                            "B Fisher's exact test: Cephalopod synteny to non-syntenic",
                            "C Fisher's exact test: Metazoan synteny to non-syntenic"), 
                  size = 12,
                  # Add labels with the coordinates 
                  x = c(0,0.08, 0.08, 0.08), y = c(1, 1, 0.66, 0.33)) 
Fishers_table_20

#STage 25
Stage25_fishers <- arrangeGrob(tableGrob(Stage25_CephMeta_fisher,theme=th),
                               tableGrob(Stage25_CephNonSyn_fisher,theme = th),
                               tableGrob(Stage25_MetaNonSyn_fisher,theme=th),
                               layout_matrix = cbind(c(1,2,3), c(1,2,3)))
# Add labels to the arranged plots
# transform to a ggplot
Fishers_table_25<- as.ggplot(Stage25_fishers) +                               
  draw_plot_label(label = c("Stage 25", 
                            "A Fisher's exact test: Cephalopod to Metazoan synteny", 
                            "B Fisher's exact test: Cephalopod synteny to non-syntenic",
                            "C Fisher's exact test: Metazoan synteny to non-syntenic"), 
                  size = 12,
                  # Add labels with the coordinates 
                  x = c(0,0.08, 0.08, 0.08), y = c(1, 1, 0.66, 0.33)) 
Fishers_table_25
#Stage29:
Stage29_fishers <- arrangeGrob(tableGrob(Stage29_CephMeta_fisher,theme=th),
                               tableGrob(Stage29_CephNonSyn_fisher,theme = th),
                               tableGrob(Stage29_MetaNonSyn_fisher,theme=th),
                               layout_matrix = cbind(c(1,2,3), c(1,2,3)))
# Add labels to the arranged plots
# transform to a ggplot
Fishers_table_29<- as.ggplot(Stage29_fishers) +                               
  draw_plot_label(label = c("Stage 29", 
                            "A Fisher's exact test: Cephalopod to Metazoan synteny", 
                            "B Fisher's exact test: Cephalopod synteny to non-syntenic",
                            "C Fisher's exact test: Metazoan synteny to non-syntenic"), 
                  size = 12,
                  # Add labels with the coordinates 
                  x = c(0,0.08, 0.08, 0.08), y = c(1, 1, 0.66, 0.33)) 
Fishers_table_29
------
  #option 1: how to add titles  
  #save the plots individually then add title using ggplot
  #converts table to a graph with minimal theme
  a<-tableGrob(Stage20_CephMeta_fisher,theme=th,)
b<-tableGrob(Stage20_CephNonSyn_fisher,theme = th)
#adding the labels
p<-as.ggplot(a)+
  draw_plot_label(label = c("A"), size = 10)
p1<-as.ggplot(b)+
  draw_plot_label(label = c("A"), size = 10)



#using the minimal theme black and white
th<-ttheme_minimal()

#option 2: To add the label for all the graphs
#to add the table titles at once
Stage20_grobs <- arrangeGrob(counts_Stage20_plot,
                             tableGrob(Stage25_CephMeta_fisher,theme=th),
                             tableGrob(Stage20_CephNonSyn_fisher,theme = th),
                             tableGrob(Stage20_MetaNonSyn_fisher,theme=th),
                             layout_matrix = cbind(c(1,1,1), c(2,3,4)))
# Add labels to the arranged plots
# transform to a ggplot
Plot<-as.ggplot(Stage20_grobs) +                               
  draw_plot_label(label = c("A", 
                            "B Fisher's exact test: Cephalopod to Metazoan synteny", 
                            "C Fisher's exact test: Cephalopod synteny to non-syntenic",
                            "D Fisher's exact test: Metazoan synteny to non-syntenic"), 
                  size = 12,
                  # Add labels with the coordinates 
                  x = c(0, 0.4, 0.4, 0.4), y = c(1, 1, 0.66, 0.33)) 
Plot
#STage 25
Stage25_grobs <- arrangeGrob(counts_Stage25_plot,
                             tableGrob(Stage25_CephMeta_fisher,theme=th),
                             tableGrob(Stage25_CephNonSyn_fisher,theme = th),
                             tableGrob(Stage25_MetaNonSyn_fisher,theme=th),
                             layout_matrix = cbind(c(1,1,1), c(2,3,4)))
# Add labels to the arranged plots
# transform to a ggplot
Plot<- as.ggplot(Stage25_grobs) +                               
  draw_plot_label(label = c("A", 
                            "B Fisher's exact test: Cephalopod to Metazoan synteny", 
                            "C Fisher's exact test: Cephalopod synteny to non-syntenic",
                            "D Fisher's exact test: Metazoan synteny to non-syntenic"), 
                  size = 12,
                  # Add labels with the coordinates 
                  x = c(0, 0.4, 0.4, 0.4), y = c(1, 1, 0.66, 0.33)) 
Plot
#Stage29:
Stage29_grobs <- arrangeGrob(counts_Stage29_plot,
                             tableGrob(Stage29_CephMeta_fisher,theme=th),
                             tableGrob(Stage29_CephNonSyn_fisher,theme = th),
                             tableGrob(Stage29_MetaNonSyn_fisher,theme=th),
                             layout_matrix = cbind(c(1,1,1), c(2,3,4)))
# Add labels to the arranged plots
# transform to a ggplot
Plot<- as.ggplot(Stage29_grobs) +                               
  draw_plot_label(label = c("A", 
                            "B Fisher's exact test: Cephalopod to Metazoan synteny", 
                            "C Fisher's exact test: Cephalopod synteny to non-syntenic",
                            "D Fisher's exact test: Metazoan synteny to non-syntenic"), 
                  size = 12,
                  # Add labels with the coordinates 
                  x = c(0, 0.4, 0.4, 0.4), y = c(1, 1, 0.66, 0.33)) 
Plot
