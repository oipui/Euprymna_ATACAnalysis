#need:
#list of gene identifiers,gene scores, list of differentially expressed genes/criterion for selecting genes based on their scores; -> could make one list with ceph score = 1, all others 0 and the same for meta
#gene GO annotations

#for the first plots I made I ran everything with "classic"; now I am repeating everything using weighted01
library(topGO)
library(ALL)
library(affy)
library(hgu95av2.db)
library(ggplot2)
#use python script gene2go, then remove ' and [] from the file in text wrangler
geneID2GO <- readMappings("gene2go.txt",sep = "\t")
#now we need a list of interesting genes, e.g. all genes in ceph synteny vs. all genes in meta synteny
#use script "get_gene_list.py" to get a list of either ceph or metazoan specific genes
#also need a list of all genes in the genome, made this also with get_gene_list.py from the expression matrix
#read list of all genes
allgenes <- scan("/Users/pui/Documents/Lab_Vienna/ATAC/GO-Analysis/allgenes_list_GO.txt", what=character(1), sep=",")
#read list of cephsynteny genes
cephsynteny<-scan("/Users/pui/Documents/Lab_Vienna/ATAC/Cluster/ceph_list.txt", what=character(1), sep=",")
#same for metazoan synteny
metasynteny<-scan("meta_list.txt", what=character(1), sep=",")
#now we want to sample ceph_synteny -> this gives us a factor with 2 levels: 0 if gene is not of interest, 1 if it is of interest
cephgeneList <- factor(as.integer(allgenes %in% cephsynteny))
names(cephgeneList) <- allgenes
str(cephgeneList)
#now we can make a topgodata object
GOdata_ceph <- new("topGOdata", ontology = "MF", allGenes = cephgeneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
#test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
#resultFisher <- getSigGroups(GOdata, test.stat)#this is the same as the run test results
#resultFis <- runTest(GOdata_ceph, algorithm = "classic", statistic = "fisher")
#change this to weight01 and see if the results are better; used classic in first pltos
resultFis <- runTest(GOdata_ceph, algorithm = "weight01", statistic = "fisher") 

pvalFis <- score(resultFis)
hist(pvalFis, 50, xlab = "p-values")
#changed this too for the weighted01 analysis
#allRes_ceph <- GenTable(GOdata_ceph, classic = resultFis, orderBy = "classic",topNodes=127 )
#allRes_ceph$classic<- as.numeric(allRes_ceph$classic)

#use topNodes if you only want the significant genes! check the table to chose them :)
allRes_ceph <- GenTable(GOdata_ceph, pvalues = resultFis, orderBy = "pvalues",topNodes=29 )
allRes_ceph$pvalues<- as.numeric(allRes_ceph$pvalues)

#to visualize (forget about this pot and scroll down)
ggplot(allRes_ceph, aes(x=Term, y=-log10(pvalues))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Molecular function") +
  ylab("Enrichment") +
  ggtitle("Title") +
  scale_y_continuous(breaks = round(seq(0, max(-log10(allRes_ceph$pvalues)), by = 2), 1)) +
  theme_bw(base_size=24) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=18, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=18, face="bold", vjust=0.5),
    axis.title=element_text(size=24, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()


#do the same thing for metazoan genes

#now we want to sample ceph_synteny -> this gives us a factor with 2 levels: 0 if gene is not of interest, 1 if it is of interest
metageneList <- factor(as.integer(allgenes %in% metasynteny))
names(metageneList) <- allgenes
str(metageneList)
#now we can make a topgodata object
GOdata_meta <- new("topGOdata", ontology = "MF", allGenes = metageneList ,annot = annFUN.gene2GO, gene2GO = geneID2GO)
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")


#resultFisher <- getSigGroups(GOdata, test.stat)#this is the same as the run test results
#resultFis <- runTest(GOdata_meta, algorithm = "classic", statistic = "fisher") # this was the first go analysis I did
#try this also: use weight01 instead of classic, which takes the hierarchy of GOs into account
resultFis <- runTest(GOdata_meta, algorithm = "weight01", statistic = "fisher") 
pvalFis <- score(resultFis)
hist(pvalFis, 50, xlab = "p-values")
#allRes_meta <- GenTable(GOdata_meta, classic = resultFis, orderBy = "classic",  topNodes = 43 )
#allRes_meta$classic<- as.numeric(allRes_meta$classic)
# changed to before: the name of the column previously called classic is now called pvalues, which is also used to order the table
# topNodes shows you the first X number of GO terms with their corresponding pvalues, independently if they are significant (<0.05) or not
allRes_meta <- GenTable(GOdata_meta, pvalues = resultFis, orderBy = "pvalues", topNodes= 24)
allRes_meta$pvalues<- as.numeric(allRes_meta$pvalues)



allRes_meta$type<-"metazoan"
allRes_ceph$type<-"cephalopod"

allRes_meta$Term <- factor(allRes_meta$Term, levels = allRes_meta$Term)
allRes_ceph$Term <- factor(allRes_ceph$Term, levels = allRes_ceph$Term)

# y=-(-log10(classic)))) need to change this from "classic" to pvalues 
ceph_col=c("#007D6D")
ceph_fill=c("#65C7C2")
ggplot(allRes_ceph, aes(x=Term, y=-(-log10(pvalues)))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge",col=ceph_col,fill=ceph_fill) +
  xlab("Molecular function") +
  ylab("Enrichment -log10 p-value") +
  ggtitle("GO enrichment in cephalopod synteny") +
  scale_y_continuous(breaks = round(seq(0, -(max(-log10(allRes_ceph$pvalues))), by = -2), 1)) +
  scale_x_discrete(name = "", position = "top")+
  theme_bw(base_size=24) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position='none',
    legend.background=element_rect(),
    axis.text.x=element_text(angle=0, size=12, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=6, face="bold", vjust=0.5),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) + coord_flip()+geom_hline(aes(yintercept=-(-log10(0.05)), linetype="p-value 0.05"))



meta_col=c("#892888")
meta_fill=c("#CE8BBC")

ggplot(allRes_meta, aes(x=Term, y=(-log10(pvalues)))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge",col=meta_col,fill=meta_fill) +
  xlab("Molecular function") +
  ylab("Enrichment -log10 p-value") +
  ggtitle("GO enrichment in metazoan synteny") +
  scale_y_continuous(breaks = round(seq(0, (max(-log10(allRes_meta$pvalues))), by = 0.5), 1)) +
  theme_bw(base_size=24) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position='none',
    legend.background=element_rect(),
    axis.text.x=element_text(angle=0, size=12, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=8, face="bold", vjust=0.5),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) + coord_flip()+ 
  geom_hline(aes(yintercept=-log10(0.05), linetype="p-value 0.05"))# #adds line

#do the same for random synteny
random1 <- scan("random1_list.txt", what=character(1), sep=",")

#now we want to sample ceph_synteny -> this gives us a factor with 2 levels: 0 if gene is not of interest, 1 if it is of interest
random1List <- factor(as.integer(allgenes %in% random1))
names(random1List) <- allgenes
str(random1List)
#now we can make a topgodata object
GOdata_random1<- new("topGOdata", ontology = "BP", allGenes = random1List,annot = annFUN.gene2GO, gene2GO = geneID2GO)
#resultFisher <- getSigGroups(GOdata, test.stat)#this is the same as the run test results
resultFis <- runTest(GOdata_random1, algorithm = "classic", statistic = "fisher")
pvalFis <- score(resultFis)
hist(pvalFis, 50, xlab = "p-values")
allRes_random1 <- GenTable(GOdata_random1, classic = resultFis, orderBy = "classic",topNodes=50 )
allRes_random1$classic<- as.numeric(allRes_random1$classic)
allRes_random1$Term <- factor(allRes_random1$Term, levels =allRes_random1$Term)
#to visualize
ggplot(allRes_random1, aes(x=Term, y=(-log10(classic)))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge",fill=c("darkgreen")) +
  xlab("Biological process") +
  ylab("Enrichment -log10 p-value") +
  ggtitle("GO enrichment in random synteny") +
  scale_y_continuous(breaks = round(seq(0, (max(-log10(allRes_random1$classic))), by = 0.5), 1)) +
  theme_bw(base_size=24) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position='none',
    legend.background=element_rect(),
    axis.text.x=element_text(angle=0, size=12, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=8, face="bold", vjust=0.5),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=18),  #Text size
    title=element_text(size=18)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) + coord_flip()+ 
  geom_hline(aes(yintercept=(-log10(0.05)), linetype="p-value 0.05"))# #adds line






