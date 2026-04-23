
#####                      loading packages                           ####

library(phyloseq)
library(ape)
library(vegan)
library(microbiome)
library(microbiomeSeq)
library(factoextra)
library(FactoMineR)
library(ggplot2)
library(plyr)
library(ggpubr)
library(knitr)
library(dplyr)
library(labdsv)                                     
library(fso)
library(picante)


#####                      importing data                           ####
abund_table<-read.csv("input/feature-table_OTU97.csv",row.names=1,check.names=FALSE)
# abund_table<-t(abund_table)#Transpose the data to have sample names on rows
meta_table<-read.csv("input/AUS_ENV_DATA_analisi_LD .csv",row.names=1, check.names=FALSE)
#OTU_tree <- read.tree("tree.nwk")
OTU_taxonomy<-read.csv("input/otu97taxonomy.csv",row.names=1,check.names=FALSE)

#Convert the data to phyloseq format
OTU = otu_table(as.matrix(abund_table), taxa_are_rows = T)
TAX = tax_table(as.matrix(OTU_taxonomy))
SAM = sample_data(meta_table)
physeq<-merge_phyloseq(phyloseq(OTU, TAX, SAM))

sample_names(physeq)
#####rimuovo A1-3 perché sono TOCA e non GLTO
#e A22-25 perchè hanno 1 sola OTU
#A7 perchè  è un outlier

Samples_toRemove<-c("A1", "A2","A3","A7","A22" , "A23" ,"A25")
physeq<-subset_samples(physeq, !(Station %in% Samples_toRemove))

sample_names(physeq)





#####                      environmental variables                           ####

#boxplot per ogni variabile !!!

# screening con   ggpubr, se è significativo nel caso si fa il post hoc test
p <- ggboxplot(poll3_Al, x = "site", y = "V",
               color = "black",  palette = c("#66ccff", "#66ff65", "#ff66ff"),fill ="site"  )# ,fill ="site",add = "jitter", palette = c("#00AFBB", "#E7B800", "#FC4E07") per colorarep + stat_compare_means()
p + stat_compare_means()  #label.y = 35000

kruskal.test(As~ site,  data = poll3_Al)


alfaKdiv<-read.csv("alfa_diversity.csv",row.names=1,check.names=FALSE)

kruskal.test(pielou~port,alfaKdiv)



# Per quali variabili ambientali i siti son diversi?

# OC
# CPE
# S_W
# DIP_W
# Naph
# SUM_16
# LPAH
# HPAH
# IdAnt
# IdBaA
# IdInp
# Cr 
# Cu
# Fe
# Ni
# Mn
# Pb
# V 

# #correlazioni pearson standardizzate
# 
# colnames(poll3_Al)
# 
# ENV_table_scaled<-scale((as.matrix(poll3_Al[,-c(1,2,29:30,35:65)])))
# 
# library(corrplot)
# library(psych)
# library(RColorBrewer)
# 
# M <- corr.test(ENV_table_scaled, method = "pearson", adjust = "none")
# 
# r<- M$r
# p <- M$p
# 
# 
# corrplot(M$r, order = "hclust", diag = F, insig = "pch",
#          method = "color",
#          p.mat  = M$p,  sig.level = .05,
#          pch.cex = 0.5,
#          tl.col = "black",
#          is.corr=T,
#          hclust.method = "complete",
#          col=brewer.pal(n=11, name="Spectral"),
#          tl.cex= 0.8,
#          mar=c(0,0,0,0),
#          addgrid.col = "gray")
# 
# corrplot(M$r, order = "hclust", diag = F, insig = "blank",
#          method = "color",
#          p.mat  = M$p,  sig.level = .05,
#          pch.cex = 0.5,
#          tl.col = "black",
#          is.corr=T,
#          hclust.method = "complete",
#          col=brewer.pal(n=11, name="Spectral"),
#          tl.cex= 0.8,
#          mar=c(0,0,0,0),
#          addgrid.col = "gray")
# 
# 
# #plot correlazione
# 
# 
# a <- ggscatter(poll3_Al, x = "Cu", y = "IdInp", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "pearson",
#           xlab = "Cu", ylab = "IdInp")
# 
# b <- ggscatter(poll3_Al, x = "Cu", y = "IdBaA", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "pearson",
#           xlab = "Cu", ylab = "IdBaA")
# 
# c <- ggscatter(poll3_Al, x = "Cu", y = "IdFlu", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "pearson",
#           xlab = "Cu", ylab = "IdFlu")
# 
# d <- ggscatter(poll3_Al, x = "Cu", y = "IdAnt", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "pearson",
#           xlab = "Cu", ylab = "IdAnt")
# 
# a <- ggscatter(poll3_Al, x = "Cu", y = "T_S", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "pearson",
#           xlab = "Cu", ylab = "T_S")



## con variabili non normalizzate ad AL, rimuovo MN e MO

colnames(meta_table)

ENV_table_scaled_noAl <-scale((as.matrix(meta_table[,-c(1,2,30:31,39:69)])))

library(corrplot)
library(psych)
library(RColorBrewer)

## 14/01/2020: come da diapositive, rimuovo Naph (10), DIN (8),DIP (9),S_D (33),S_YE (32),S_PAH (34)
colnames(ENV_table_scaled_noAl)

# è da fare coi dati non scalati ad Al.
## File di partenza è meta_table che contiene i dati non normalizzati, quindi ok

ENV_table_scaled_noAl_selected <- ENV_table_scaled_noAl[,-c(8,9,10,32,33,34)]
colnames(ENV_table_scaled_noAl_selected)

M <- corr.test(ENV_table_scaled_noAl_selected, method = "pearson", adjust = "none")

r<- M$r
p <- M$p


corrplot(M$r, order = "hclust", diag = F, insig = "pch",
         method = "color",
         p.mat  = M$p,  sig.level = .05,
         pch.cex = 0.5,
         tl.col = "black",
         is.corr=T,
         hclust.method = "complete",
         col=brewer.pal(n=11, name="Spectral"),
         tl.cex= 0.8,
         mar=c(0,0,0,0),
         addgrid.col = "gray")

corrplot(M$r, order = "hclust", diag = F, insig = "blank",
         method = "color",
         p.mat  = M$p,  sig.level = .05,
         pch.cex = 0.5,
         tl.col = "black",
         is.corr=T,
         hclust.method = "complete",
         col=brewer.pal(n=11, name="Spectral"),
         tl.cex= 0.8,
         mar=c(0,0,0,0),
         addgrid.col = "gray")



#                                             PCA                               #####

colnames(ENV_table_scaled_noAl)
metalli <- ENV_table_scaled_noAl[,21:31]

res.pca <- prcomp(metalli, scale = TRUE)

fviz_eig(res.pca) #Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component
fviz_eig(res.pca,choice = "variance",addlabels = T)
ind<-get_pca_ind(res.pca)
var<-get_pca_var(res.pca)
var

#contributions of variables to PC1
fviz_contrib(res.pca, choice = "var",axes = 1)

fviz_contrib(res.pca, choice = "var",axes = 1, top = 10)

#contributions of variables to PC2
fviz_contrib(res.pca, choice = "var",axes =2)

fviz_contrib(res.pca, choice = "var",axes =2 , top = 10)

#contributions of variables to PC1&2
fviz_contrib(res.pca, choice = "var",axes =1:2 )

fviz_contrib(res.pca, choice = "var",axes =1:2 , top = 10)



fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality on the factor map
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_biplot(res.pca,addEllipses = TRUE,col.var =  "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
fviz_pca_biplot(res.pca,col.var =  "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)




####COLORARE I siteI IN BASE AL LORO COLORE!!!!!!!!
fviz_pca(res.pca,col.var =  "contrib",  
         gradient.cols = c("#f0f0f0", "#969696", "#525252"), repel = TRUE) + 
  theme_bw() + 
  ggtitle("") +
  geom_point(size=3, shape = 21,aes( fill = poll3_Al$site))


# #pca ggplot
# library(ggfortify)
# 
# metalli<-as.data.frame(poll3_Al[,c(23:34)])
# 
# autoplot(prcomp(metalli,scale. = T), data = poll3_Al, colour = 'site',
#          loadings = TRUE, loadings.colour = 'black',
#          loadings.label = TRUE, loadings.label.size = 5)+scale_color_manual(values=allGroupsColors)+ geom_text(mapping = aes(label = station), size = 4, vjust = 1,hjust=1.5) +stat_ellipse(type = "t")


### 14/01/2020: PCA su tutte le variabili con metalli non normalizzati e cerchio
# 
# ## Questa sui normalizzati
# colnames(ENV_table_scaled_noAl)
# metalli <- ENV_table_scaled_noAl[,c(12:15,22:31)]
# 
# res.pca <- prcomp(metalli, scale = T)
# 
# fviz_eig(res.pca,choice = "variance",addlabels = T)
# ind<-get_pca_ind(res.pca)
# var<-get_pca_var(res.pca)
# var
# 
# fviz_pca(res.pca,
#          addEllipses = TRUE,
#          fill.ind = poll3_Al$site,
#          gradient.cols = c("#f0f0f0", "#969696", "#525252"), 
#          col.var = "gray25",
#          repel = TRUE) + 
#   theme_bw() + 
#   ggtitle("") +
#   geom_point(size=3, shape = 21,aes( fill = poll3_Al$site)) 
#   xlab("PC1 (35.0%)") + 
#   ylab("PC2 (27.2%)")
# 
# #contributions of variables to PC1
# fviz_contrib(res.pca, choice = "var",axes = 1)
# 
# #contributions of variables to PC2
# fviz_contrib(res.pca, choice = "var",axes =2)
# 
# #contributions of variables to PC3
# fviz_contrib(res.pca, choice = "var",axes =3)


## Questa sui NON normalizzati

# 22/01/20 : figura 2 definitiva

unnormalized <- read.csv("variabili_ambientali_terza_stagione_140619.csv",row.names=1,check.names=FALSE)
ENV_table_scaled_unnorm <-scale((as.matrix(unnormalized[,-c(1,2,30:31,39:69)])))

colnames(ENV_table_scaled_noAl) == colnames(ENV_table_scaled_unnorm)

colnames(env_corr)
colnames(ENV_table_scaled_unnorm)

metalli <- ENV_table_scaled_unnorm[,c(12:31)]
metalli <- metalli[,c(1:4,10:20)] ## 24/01/20 versione definitiva
colnames(metalli)

res.pca <- prcomp(metalli, scale = F)

fviz_eig(res.pca,choice = "variance",addlabels = T)
ind<-get_pca_ind(res.pca)
var<-get_pca_var(res.pca)
var

library(cowplot)
library(ggpubr)

palette(allGroupsColors)
fviz_pca(res.pca,
         axes = c(1,2),
         addEllipses = F,
         fill.ind = poll3_Al$site,
         #gradient.cols = c("#f0f0f0", "#969696", "#525252"), 
         col.var = "grey40",
         repel = TRUE) + 
  theme_bw() + 
  ggtitle("") +
  geom_point(size=3, shape = 21,aes( fill = poll3_Al$site)) + 
  scale_fill_manual(values = allGroupsColors) + 
  xlab("PC1 (37.8%)") + 
  ylab("PC2 (26%)") +
  theme(legend.title = element_blank())

fviz_pca(res.pca,
         axes = c(1,3),
         addEllipses = F,
         fill.ind = poll3_Al$site,
         #gradient.cols = c("#f0f0f0", "#969696", "#525252"), 
         col.var = "grey40",
         repel = TRUE) + 
  theme_bw() + 
  ggtitle("") +
  geom_point(size=3, shape = 21,aes( fill = poll3_Al$site)) + 
  scale_fill_manual(values = allGroupsColors) + 
  xlab("PC1 (37.8%)") + 
  ylab("PC3 (13.4%)") +
  theme(legend.title = element_blank())

#contributions of variables to PC1
fviz_contrib(res.pca, choice = "var",axes = 1)

#contributions of variables to PC2
fviz_contrib(res.pca, choice = "var",axes =2)

#contributions of variables to PC3
fviz_contrib(res.pca, choice = "var",axes =3 )


###################
######## stessa immagine ma includendo anche

# T_S
# Eh 
# OC%
# CPE
# SC%
# S_W

ENV_table_scaled_unnorm <-scale((as.matrix(unnormalized[,-c(1,2,30:31,39:69)])))

colnames(ENV_table_scaled_noAl) == colnames(ENV_table_scaled_unnorm)

colnames(env_corr)
colnames(ENV_table_scaled_unnorm)

metalli_amb <- ENV_table_scaled_unnorm[,c(1:4,6,7,12:15, 21:31)]
colnames(metalli_amb)

res.pca <- prcomp(metalli_amb, scale = F)

fviz_eig(res.pca,choice = "variance",addlabels = T)
ind<-get_pca_ind(res.pca)
var<-get_pca_var(res.pca)
var

library(cowplot)
library(ggpubr)

palette(allGroupsColors)
fviz_pca(res.pca,
         axes = c(1,2),
         addEllipses = F,
         fill.ind = poll3_Al$site,
         #gradient.cols = c("#f0f0f0", "#969696", "#525252"), 
         col.var = "grey40",
         repel = TRUE) + 
  theme_bw() + 
  ggtitle("") +
  geom_point(size=3, shape = 21,aes( fill = poll3_Al$site)) + 
  scale_fill_manual(values = allGroupsColors) + 
  xlab("PC1 (37.6%)") + 
  ylab("PC2 (24%)") +
  theme(legend.title = element_blank())

fviz_pca(res.pca,
         axes = c(1,3),
         addEllipses = F,
         fill.ind = poll3_Al$site,
         #gradient.cols = c("#f0f0f0", "#969696", "#525252"), 
         col.var = "grey40",
         repel = TRUE) + 
  theme_bw() + 
  ggtitle("") +
  geom_point(size=3, shape = 21,aes( fill = poll3_Al$site)) + 
  scale_fill_manual(values = allGroupsColors) + 
  xlab("PC1 (37.6%)") + 
  ylab("PC3 (13.5%)") +
  theme(legend.title = element_blank())

#contributions of variables to PC1
fviz_contrib(res.pca, choice = "var",axes = 1)

#contributions of variables to PC2
fviz_contrib(res.pca, choice = "var",axes =2)

#contributions of variables to PC3
fviz_contrib(res.pca, choice = "var",axes =3 )



### a sto punto anche con solo var ambientali

amb <- ENV_table_scaled_unnorm[,c(1,2,3,4,6,7)]
colnames(amb)

res.pca <- prcomp(amb, scale = F)

fviz_eig(res.pca,choice = "variance",addlabels = T)
ind<-get_pca_ind(res.pca)
var<-get_pca_var(res.pca)
var

library(cowplot)
library(ggpubr)

palette(allGroupsColors)
fviz_pca(res.pca,
         axes = c(1,2),
         addEllipses = F,
         fill.ind = poll3_Al$site,
         #gradient.cols = c("#f0f0f0", "#969696", "#525252"), 
         col.var = "grey40",
         repel = TRUE) + 
  theme_bw() + 
  ggtitle("") +
  geom_point(size=3, shape = 21,aes( fill = poll3_Al$site)) + 
  scale_fill_manual(values = allGroupsColors) + 
  xlab("PC1 (53.3%)") + 
  ylab("PC2 (20.6%)") +
  theme(legend.title = element_blank())

fviz_pca(res.pca,
         axes = c(1,3),
         addEllipses = F,
         fill.ind = poll3_Al$site,
         #gradient.cols = c("#f0f0f0", "#969696", "#525252"), 
         col.var = "grey40",
         repel = TRUE) + 
  theme_bw() + 
  ggtitle("") +
  geom_point(size=3, shape = 21,aes( fill = poll3_Al$site)) + 
  scale_fill_manual(values = allGroupsColors) + 
  xlab("PC1 (53.3%)") + 
  ylab("PC3 (16.3%)") +
  theme(legend.title = element_blank())

#contributions of variables to PC1
fviz_contrib(res.pca, choice = "var",axes = 1)

#contributions of variables to PC2
fviz_contrib(res.pca, choice = "var",axes =2)

#contributions of variables to PC3
fviz_contrib(res.pca, choice = "var",axes =3) 


###                      ANALISI COMUNITA MICROBICA                           ####


####                                 normalizzazione                 ###############
library(metagenomeSeq)
phyloseq_obj<-physeq

# bokulich filtration

minTotRelAbun = 5e-5
x = taxa_sums(phyloseq_obj)
keepTaxa = (x / sum(x)) > minTotRelAbun
prunedSet = prune_taxa(keepTaxa, phyloseq_obj)
otu_table(prunedSet)
prunedSet


#filtering singletons
doubleton <- genefilter_sample(prunedSet, filterfun_sample(function(x) x > 1), A=1)
doubleton <- prune_taxa(doubleton, prunedSet) 
doubleton



## Controllo mitocondrio e cloroplasto

grep(pattern = "Mitochondria", tax_table(doubleton)) 
grep(pattern = "Chloroplast", tax_table(doubleton)) 

#rimuovo mitochondria 
mito <- taxa_names(subset_taxa(doubleton, Family  != " f__Mitochondria")) #2 mitok 
doubleton <- prune_taxa(mito, doubleton)
#rimuovo chloroplats 
subset_taxa(doubleton, Order != " o__Chloroplast") 
chloro <- taxa_names(subset_taxa(doubleton, Order  != " o__Chloroplast"))
doubleton <- prune_taxa(chloro, doubleton)
doubleton

# transforming 
data.metagenomeSeq = phyloseq_to_metagenomeSeq(doubleton)
p = cumNormStat(data.metagenomeSeq) #default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
#data.cumnorm
data.CSS = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
head(data.CSS)
dim(data.CSS)  # make sure the data are in a correct formal: number of samples in rows
phyloseq_obj_css <- phyloseq_obj
otu_table(phyloseq_obj_css) <- otu_table(data.CSS, taxa_are_rows = T)

phyloseq_obj_css

sample_names(phyloseq_obj_css)



plot_heatmap(tax_glom(phyloseq_obj_css, taxrank = "Genus"),sample.order = c("A4","A5","A6","A7","A8","A9","A10","A11","A12","A13","A14","A15","A16","A17","A18","A20","A21"), title="Heatmap at genus level")+
  theme(plot.title = element_text(hjust = 0.5))


# newSTorder = c("A4","A5","A6","A7","A8","A9","A10","A11","A12","A13","A14","A15","A16","A17","A18","A20","A21")



#salvare l'otu table 
# Extract abundance matrix from the phyloseq object
# OTU_matrix = as(otu_table(phyloseq_obj_css), "matrix")
# OTUdf = as.data.frame(OTU_matrix)
# head(OTUdf)
# write.csv(OTUdf,"OTU_table_bacteria_normalized.csv")

#export dist-matrix
B_C_COMUNITA_BATTERICA<-vegdist(t(OTUdf),method = "bray")
#write.csv(as.matrix(B_C_COMUNIT?_BATTERICA), "dist_matrix_BC.csv")

sample_names(phyloseq_obj_css)

write.csv(as.data.frame(tax_table(phyloseq_obj_css)),"tax_normalized.csv")
#rifare physeq con numero otu assolute e non scalate CSS in modo che siano interi
# e che vadano bene per il calcolo indici alpha div dopo

otu.absol<-round(otu_table(phyloseq_obj_css))
head(otu.absol)
OTU_absol = otu_table(as.matrix(otu.absol), taxa_are_rows = T)
physeq_normalized<-merge_phyloseq(phyloseq(OTU_absol, TAX),SAM)
# sample_data(phyloseq_obj_css)<-sample_data(physeq_normalized)

#questo lo usavo per bioenv
#df_otu table rounded 
batt1<-as.data.frame(otu_table(physeq_normalized))
head(batt1)

#                                     Alpha diversity                                  ####
head(otu_table(physeq_normalized))
sample_data(physeq_normalized)


grep(pattern = "Propionibacteriaceae", tax_table(physeq_normalized)) 
physeq_normalized <- subset_taxa(physeq_normalized, Genus != "Propionibacteriaceae")


p<-plot_richness(physeq_normalized,x="Tow" ,measures=c("Observed", "Shannon", "Simpson"),color="Tow")+ 
  scale_color_manual(values = c("#007FFF","#1aff00","#ff0008","#ff7003"))
p
write.csv(p$data, "alfadiv.csv")
p$data
# newSTorder = c("Before", "During", "After")

p$data$Period<- as.character(p$data$Period)
p$data$Period <- factor(p$data$Period, levels=newSTorder)


alpfadiv<-p+geom_boxplot(alpha = 0) + theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(legend.title.align=0.5) +theme(axis.text.x=element_text(colour="black"),axis.text.y =element_text(colour="black") )
alpfadiv

#statistica per l'alfa diversity
alf_div<-estimate_richness(physeq_normalized,measures=c("Observed", "Shannon", "Simpson"))
write.csv(alf_div,"alfa_diversity1.csv")
a<-as.data.frame(sample_data(physeq_normalized))


alf_div$Period<-a$Period


kruskal.test(Observed~ Period,  data = alf_div)
kruskal.test(Shannon~ Period,  data = alf_div)

pairwise.wilcox.test(alf_div$Shannon, 
                          alf_div$Period, 
                          p.adjust.method="none")


kruskal.test(Simpson~ Period,  data = alf_div)





#tutti gli indici di alfa diversita con pacchetto microbiome:

#fare l'ggetto PSEQ!!!!!!!!!!!!!!!
pseq_tutt<-physeq_normalized


#meta_tutt<-meta(SRB_physeq_silva_rounded)
#otu_tutt<-abundances(SRB_physeq_silva_rounded)


alpha(pseq_tutt,)
tab <- alpha(pseq_tutt, index = c("observed","diversities_shannon","evenness_pielou"))
tab$diversities_shannon<-log(tab$diversities_shannon)

tab$station<-c("C1","C2","C3","C4","C5","E1","E2","E3","H1","H3","H4","H5")
tab$site<-c(rep("C",5),rep("E",3),rep("H",4))
colnames(tab)
#write.csv(tab,"ALFADIVROUNDED_SRB_H.csv")
#alfa_div_HHH<-read.csv("ALFADIVROUNDED_SRB_H.csv",row.names=1,check.names=FALSE)

#kruskal.test(observed ~site,alfa_div_HHH)
#kruskal.test(diversities_shannon ~site,alfa_div_HHH)
#kruskal.test(evenness_pielou ~site,alfa_div_HHH)

#pairwise.wilcox.test(alfa_div_HHH$observed, 
#                    alfa_div_HHH$site, 
#                     p.adjust.method="bonf")
#pairwise.wilcox.test(alfa_div_HHH$diversities_shannon, 
#                    alfa_div_HHH$site, 
#                    p.adjust.method="bonf")
#pairwise.wilcox.test(alfa_div_HHH$evenness_pielou, 
#                    alfa_div_HHH$site, 
#                    p.adjust.method="bonf")

colnames(tab)[colnames(tab)=="observed"] <- "OTU Richness"
colnames(tab)[colnames(tab)=="diversities_shannon"] <- "Shannon H"
colnames(tab)[colnames(tab)=="evenness_pielou"] <- "Pielou J'"
colnames(tab)
tab$number<-c("1","2","3","4","5","1","2","3","1","3","4","5")


library(reshape2)

#per fare unico grafico con tutto
df.c<-melt(tab, c("site","station","number"))

library(ggrepel)

#boxplot porti
p1<-ggplot(df.c, aes(site,value ,color= site )) +
  geom_boxplot(alpha = 0) + 
  facet_wrap(~variable, scales = "free")+
  scale_color_manual(values=allGroupsColors)+
  theme_bw()+
  geom_point()
# + geom_text(aes(label=number),hjust=0, vjust=0,size=3)
p1+ geom_text_repel(aes(label=number),size=2.4)







#####                      BETA DIVERSITY                           ####

#PCoA on Bray-Curtis Dissimilarity
phyloseq_obj_css

sample_data(phyloseq_obj_css)




  
  #PCoA on Bray-Curtis Dissimilarity

phyloseq_obj_css
grep(pattern = "Propionibacteriaceae", tax_table(phyloseq_obj_css)) 
phyloseq_obj_css <- subset_taxa(phyloseq_obj_css, Genus != "Propionibacteriaceae")

  otu.ord <- ordinate(physeq = phyloseq_obj_css, "PCoA")
  

  ##beta diversity  Period

  #asse 1-2  
  a<-plot_ordination(physeq = phyloseq_obj_css, otu.ord,color = "Tow",
                  axes =c(1,2))+
    theme_bw()+ geom_point(size = 2)+
    geom_text(aes(label=Station), size = 3, vjust = 0,hjust=0)+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_color_manual(values = c("#007FFF","#1aff00","#ff0008","#ff7003"))
  
  #asse 1-3 
  b<-plot_ordination(physeq = phyloseq_obj_css, otu.ord,color = "Tow",
                  axes =c(1,3))+
    theme_bw()+ geom_point(size = 2)+
    geom_text(aes(label=Station), size = 3, vjust = 0,hjust=0)+
    theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("#007FFF","#1aff00","#ff0008","#ff7003"))
  
  title ="PCoA on Bray-Curtis Dissimilarity axes 1&3 \n Permanova p = 0.001" 
  c<-ggarrange(a, b , 
            labels = c("A", "B"),
            ncol = 2, nrow = 1,
            legend = "bottom",
            common.legend=T	)
  c
  betadivimg<-annotate_figure(c,top = text_grob("PCoA on Bray-Curtis Dissimilarity \n Permanova p = 0.001", color = "Black", face = "bold", size = 14))
  
  tiff("betadiv.tiff", units="in", width=9, height=5, res=300)
  betadivimg
  dev.off()
  
  
  
  
  
  ##beta diversity  STAGIONALITà

  
  #asse 1-2  
  plot_ordination(physeq = phyloseq_obj_css, otu.ord,color = "Season",
                  axes =c(1,2),title ="PCoA on Bray-Curtis Dissimilarity axes 1&2 \n Permanova p = 0.001" )+
    theme_bw()+ geom_point(size = 2)+
    geom_text(aes(label=Station), size = 3, vjust = 0,hjust=0)+
    theme(plot.title = element_text(hjust = 0.5))#)+
    scale_color_manual(values = c("#007FFF","#1aff00","#ff0008"))
  
  #asse 1-3 
  plot_ordination(physeq = phyloseq_obj_css, otu.ord,color = "Season",
                  axes =c(1,3),title ="PCoA on Bray-Curtis Dissimilarity axes 1&3 \n Permanova p = 0.001")+
    theme_bw()+ geom_point(size = 2)+
    geom_text(aes(label=Station), size = 3, vjust = 0,hjust=0)+
    theme(plot.title = element_text(hjust = 0.5))# +
    scale_color_manual(values = c("#007FFF","#1aff00","#ff0008"))
  
  
  ############ extraction vibrio for biplot

    VIBRIO_extr = subset_taxa(phyloseq_obj_css, Genus==" g__Vibrio")
    
  
    
    otu.ord_vibrio <- ordinate(physeq = VIBRIO_extr, "PCoA")
    
    sample_data(phyloseq_obj_css)$Period = factor(sample_data(phyloseq_obj_css)$Period, levels=c("Before", "During", "After"))
    
    
    ##beta diversity  Period
    
    #asse 1-2  
    a<-plot_ordination(physeq = VIBRIO_extr, otu.ord_vibrio, type="biplot",color = "Tow",
                       axes =c(1,2))+
      theme_bw()+ geom_point(size = 2)+
      geom_text(aes(label=Station), size = 3, vjust = 0,hjust=0)+
      theme(plot.title = element_text(hjust = 0.5))+
      scale_color_manual(values = c("#007FFF","#1aff00","#ff0008","#ff7003","#E52B50"))
    
    a
    
    
    
    
    
  
###### test difference between ports!!!!!!

OTU_matrix = as(otu_table(phyloseq_obj_css), "matrix")
OTUdf = as.data.frame(OTU_matrix)
braycurtis_bact<-vegdist(t(OTUdf),method = "bray")
env_data<-as.data.frame(sample_data(phyloseq_obj_css))


#permanova:
#Number of permutations: 999  

set.seed(123)
adonis(braycurtis_bact~env_data$Tow) #p: 0.001
set.seed(123)
adonis(braycurtis_bact~env_data$Season) # p:0.001


##plot taxa most abund

genus.sum = tapply(taxa_sums(phyloseq_obj_css_barplot), tax_table(phyloseq_obj_css_barplot)[, "Genus"], sum, na.rm=TRUE)
top5phyla = names(sort(genus.sum, TRUE))[1:15]
GP1 = prune_taxa((tax_table(phyloseq_obj_css_barplot)[, "Genus"] %in% top5phyla), phyloseq_obj_css_barplot)

GP2 <- subset_taxa(GP1, Genus != "Unassigned")
GP2 <- subset_taxa(GP2, Genus != "uncultured")
GP2 <- subset_taxa(GP2, Genus != "Propionibacteriaceae")
grep(pattern = "Propionibacteriaceae", tax_table(GP2)) 

GP.ord <- ordinate(GP2, "PCoA", "bray")
plot_ordination(GP2, GP.ord, type="taxa", color="Phylum", title="taxa")+ facet_wrap(~Genus, 3)
p3 = plot_ordination(GP2, GP.ord, type="biplot", color="Genus",  title="biplot")
p3 #  + facet_wrap(~Genus, 3)

otu_table(GP2)

tax_table(GP2)
p<-plot_bar(GP2, "Genus", fill="Genus", facet_grid=~Tow) #Tow 

# newSTorder = c("A4","A5","A6","A7","A8","A9","A10","A11","A12","A13","A14","A15","A16","A17","A18","A20","A21")

p$data$Sample<- as.character(p$data$Sample)
p$data$Sample <- factor(p$data$Sample, levels=newSTorder)

p+scale_fill_manual(values =c( "#708090","#7FFF00","#00A86B","#e6a910","#7B1B02",
                               "#aedce8","#884DA7","#E52B50","#293133","#7FFFD4",
                               "#708090","#1025e6","#fc0317"))


p + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")


write.csv(as.data.frame(p$data), "data_barplot_mostabundant.csv")






## togliere linee morte 


my_plot_bar = function (GP2, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                        facet_grid = NULL) {
  mdf = psmelt(GP2)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}



p<-my_plot_bar(GP2, "Genus", fill = "Genus", facet_grid=~Tow)+
  scale_fill_manual(values =c( "#7FFF00","#94755c","#fbff00","#7B1B02",
                               "#aedce8","#00A86B","#ff9d4d","#4B0082","#293133","#7FFFD4",
                               "#708090","#1025e6","#fc0317","#FFBF00","#ff0303","#4B0082"))

a<-p+theme_bw()+ theme(axis.text.x = element_text(angle = -90, hjust = 0,colour = "black"))
a<-p+theme_bw()+theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())
a+theme(axis.text=element_text(size=12,colour = "black"),
        axis.title=element_text(size=14,face="bold"))




#884DA7


########  phyloseq object con solo VIBRIO

Phylo_vibrio = subset_taxa(GP2, Genus=="Vibrio")
Phylo_vibrio_genus<-tax_glom(Phylo_vibrio, "Genus")
sample_data(Phylo_vibrio_genus) 
Phylo_vibrio_genusDF<-t(as.data.frame(otu_table(Phylo_vibrio_genus)))

write.csv(Phylo_vibrio_genusDF,"Phylo_vibrio_genusDF.csv")
mediavai_sumgenus<-read.csv("input/mediavai&sumgenus.csv",sep = ";")


##corr vai e otu

cor(mediavai_sumgenus$VAI_average_TOW,mediavai_sumgenus$VIBRIO_GENUS_SUM_TOW, method = "pearson")
cor.test(mediavai_sumgenus$VAI_average_TOW,mediavai_sumgenus$VIBRIO_GENUS_SUM_TOW, method = "pearson")

sp<-ggplot(mediavai_sumgenus, aes(x=VAI_average_TOW, y=VIBRIO_GENUS_SUM_TOW)) + 
  geom_point()+
  geom_smooth(method=lm)
sp + stat_cor(method = "pearson",r.digits = 3)




############ lefse ################## 
library(microbiomeMarker)

otu_css_names <- taxa_names(GP2) 
gp2_physeq <- prune_taxa(otu_css_names,physeq)
write.csv(as.data.frame(otu_table(gp2_physeq)), "otu_12_most_abund_gp2_physeq.csv")
taxa_names(GP2) == taxa_names(gp2_physeq)

colnames(tax_table(GP2)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus","Species")
write.csv(as.data.frame(otu_table(GP2)), "otu_12_most_abund.csv")
normalize()



mm <- lefse(
  GP2, 
  norm = "CSS", 
  class = "Tow", 
  multicls_strat = T)

head(marker_table(mm))
otu_table(mm)
#write.csv(as.data.frame(otu_table(mm)), "otu_12_most_abund_css_lefse.csv")

  

plot_ef_bar(mm, label_level = 7) +
  scale_fill_manual(values = c("GLTO20141129" = "blue", "GLTO20150414" = "green", "GLTO20160416"="red", "GLTO20160816"="orange"))

plot_cladogram(mm, color = c(GLTO20141129 = "blue", GLTO20150414 = "green", GLTO20160416="red", GLTO20160816="orange"),
                  group_legend_param = ggplot2::guides(fill=guide_legend(nrow=4)),
                  marker_legend_param = ggplot2::guides(fill=guide_legend(nrow=8))) + 
  theme(legend.position = "top")  
 


plot_cladogram(
  mm,only_marker = T,clade_label_font_size = 5,
  color=c(GLTO20141129 = "blue", GLTO20150414 = "green", GLTO20160416="red", GLTO20160816="orange"),
  branch_size = 1,
  alpha = 0.2,
  node_size_scale = 1,
  node_size_offset = 1,
  clade_label_level = 4,
  annotation_shape = 22,
  annotation_shape_size = 5, #per la legenda 
  group_legend_param = list(),
  marker_legend_param = list()) 



geom_text(aes(), size = 3)
geom_tiplab(size =2)













library(cowplot)
my_legend <- get_legend(p)
library(ggpubr)
as_ggplot(my_legend)





#####         ANALISI LINK TRA COMUNITà MICROBICA E VARIABILI AMBIENTALI   ####
#plot temperatura:

data_T<-read.csv("input/temperature.csv",row.names=1, check.names=FALSE)

data_T$elnino <- factor(data_T$elnino, levels=c("Before", "During", "After"))

p<-ggplot(data = data_T, aes(x = Date, y = Daily_SST, group=1, color= elnino))   


  geom_line()

  Torder = c("29/11/2014",
             "30/11/2014",
             "15/04/2015",
             "16/04/2015",
             "16/08/2015",
             "17/08/2015",
             "18/10/2015",
             "19/10/2015",
             "16/04/2016",
             "17/04/2016",
             "17/08/2016")
  p$data$Date<- as.character(p$data$Date)
  p$data$Date <- factor(p$data$Date, levels=Torder)
  a <- ifelse(data_T$SEQ == "si", "red", "black")
  
  #plot grafico temperatura  
  p+ labs(y="Weekly SST (°C)")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))+
    geom_path()+
    geom_point()+scale_linetype_discrete("El Niño") +
    scale_shape_discrete("El Niño") +
    scale_colour_discrete("El Niño")+
    scale_color_manual(values = c("#4595b0","#fcba03","#5e0219"))
 
  

  
  

  
  
 z<- env_data[,c(2,4,5,8,9,11,12,16)]
 z<-env_data[,4]
set.seed(123)

adonis(braycurtis_bact~env_data$Daily_SST,perm = 9999)



adonis(braycurtis_bact~env_data$Daily_SST,perm = 9999)
adonis(braycurtis_bact~env_data$Weekly_SST,perm = 9999)
adonis(braycurtis_bact~env_data$PCI_idx,perm = 9999)

#https://fromthebottomoftheheap.net/slides/advanced-vegan-webinar-2020/advanced-vegan#57

by = "terms"
by = "margin"
adonis(braycurtis_bact~Daily_SST+Weekly_SST+PCI_idx,data=env_data,perm = 9999)

adonis()
adonis(braycurtis_bact~env_data$Daily_SST*env_data$Weekly_SST,perm = 9999)
adonis(braycurtis_bact~env_data$Weekly_SST*env_data$Daily_SST*env_data$PCI_idx,perm = 9999)





############### VAI index e temperatura

data_T1<-read.csv("input/luigi environmental_data_GLTO.csv",row.names=1, check.names=FALSE)

dati_t1_VAI<- merge(as.data.frame(data_T1),as.data.frame(VAI), by="Trip_code")
dati_t1_VAI_WO_NA<-dati_t1_VAI[-c(3,6,8,14,20),]
as.list(dati_t1_VAI_WO_NA$VAI)
##correlazioni temperatura VAI

cor.test( dati_t1_VAI_WO_NA$VAI,dati_t1_VAI_WO_NA$Daily_SST, method = c("pearson"))

cor.test(dati_t1_VAI_WO_NA$VAI,dati_t1_VAI_WO_NA$Weekly_SST, method = c("pearson"))


cor.test(dati_t1_VAI_WO_NA$VAI,dati_t1_VAI_WO_NA$Weekly_SST, method = c("spearman"))

cor.test(dati_t1_VAI_WO_NA$VAI,dati_t1_VAI_WO_NA$Daily_SST,  method = c("spearman"))

ggplot(data = dati_t1_VAI_WO_NA, aes(x = Event, y = VAI, label=Trip_code))+ geom_boxplot()

library("ggpubr")
ggscatter(dati_t1_VAI_WO_NA, x = "VAI", y = "Daily_SST", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman")#  ,
        #  xlab = "Daily SST °C", ylab = "VAI")
# Shapiro-Wilk normality test for mpg
shapiro.test(my_data$mpg) # => p = 0.1229

# p-values are greater than the significance level 0.05 implying that the distribution
#  of the data are not significantly different from normal distribution.
# In other words, we can assume the normality

VAI<-read.csv("input/DATI_VAI_SOLO_STAZIONI.csv",row.names=1, check.names=FALSE)


 VAI$Station <- factor(VAI$Station, levels =c("A4","A5","A6","A7","A8","A9","A10","A11","A12","A13","A14","A15","A16","A17","A18","A20","A21"))
                  
p<-ggplot(data=VAI, aes(x=Station, y=VAI)) +
  geom_bar(stat="identity")
  
p+facet_grid(~Tow, scale='free_x')+labs(y= "VAI %", x = "Sample") + theme_bw()
#  coord_cartesian( ylim = c(0, 100))

library(dplyr)
gd <- VAI %>% 
  group_by(Tow) %>% 
  summarise(VAI = mean(VAI))
gd


p<-ggplot(data=gd, aes(x=Tow, y=VAI)) +
  geom_bar(stat="identity", fill = "#cf0000")+facet_grid(~Tow, scale='free_x')

b<-p+facet_grid(~Tow, scale='free_x')+labs(y= "VAI %", x = "Tow") +theme_bw()# oord_cartesian( ylim = c(0, 120))

b
b















###############################   MAPPA   ##################### 
library(dplyr) # A staple for modern data management in R
library(lubridate) # Useful functions for dealing with dates
library(ggplot2) # The preferred library for data visualisation
library(tidync) # For easily dealing with NetCDF data
library(rerddap) # For easily downloading subsets of data
library(doParallel) # For parallel processing
library(ggOceanMaps)


# MAPdata<-as.data.frame(meta_table[-c(1:3),c(1,3,4,10,11)])
# write.csv(MAPdata,"mapdata.csv")

#  con valori di lat e long fissi!!!!!--->MAPdata1
MAPdata1<-read.csv("mapdata.csv",row.names=1,check.names=FALSE)









###mappa base GLTO
### COSTA GLTO
basemap(limits = c(143, 155, -25, -15),grid.col = NA,land.col = "#783d01")


##mappa BASE queensland
basemap(limits = c(137.99465, 153.611604, -29.17927, -9.088012),bathymetry = T,grid.col = NA,land.col = "#eeeac4")

##mappa BASE AUS CON BATIMETRIA
basemap(limits = c(110.1, 155.1, -40.17927, -8.088012),bathymetry = T,grid.col = NA,land.col = "#eeeac4")+
  theme(legend.position="left") +labs(x = "Longitude", y = "Latitude")


## queensland, stazioni per latitudine e tow
basemap(limits = c(137.99465, 153.611604, -29.17927, -9.088012)
        ,bathymetry = F, grid.col = NA,land.col = "#eeeac4")+
  labs(x = "Longitude", y = "Latitude")+
  geom_point(data = MAPdata1, aes(x = Longitude, y = Lat,shape=Tow, color=Latitude))

## GLTO, stazioni per latitudine
basemap(limits = c(143, 155, -25, -15)
        ,bathymetry = F, grid.col = NA,land.col = "#783d01",land.border.col="#000000",land.size=0.4)+
  labs(x = "Longitude", y = "Latitude") # +
  geom_point(data = MAPdata1, aes(x = Long, y = Lat, color=Latitude))


## GLTO, stazioni per latitudine e tow
basemap(limits = c(143, 155, -25, -15)
        ,bathymetry = F, grid.col = NA,land.col = "#eeeac4")+
  labs(x = "Longitude", y = "Latitude")+
  geom_point(data = MAPdata1, aes(x = Longitude, y = Lat,shape=Tow,color=Latitude,size=2))+
  scale_shape_manual(values=c(0, 1, 2,6)) + scale_size(guide="none")



## GLTO, stazioni per latitudine e tow--------> quella per la mappa
mappaIniziale<-basemap(limits = c(143, 155, -25, -15)
        ,bathymetry = F, grid.col = NA,land.col = "#eeeac4")+
  labs(x = "Longitude", y = "Latitude")+
  geom_point(data = MAPdata1, aes(x = Longitude, y = Lat,shape=Tow, color=Tow), size = 1.5, stroke = 1)+
  scale_shape_manual(values=c(0, 1, 2,6)) + scale_size(guide="none")+
  geom_point(data = citta, aes(x = Long, y = Lat), size = 3, shape=1.5, stroke = 1)+
  geom_text(data = citta, aes(x = Long, y = Lat, label= Cityc, hjust= 1.2))+
  scale_color_manual(values = c("#007FFF","#1aff00","#ff0008","#ff7003"))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y =element_text(colour="black"))+
  theme(legend.title.align=0.5) +theme(
                                       axis.text.y=element_blank(),
                                       axis.title.y=element_blank(),axis.ticks.y=element_blank()
                                       )
mappaIniziale<-mappaIniziale+ scale_y_continuous(expand = c(0, 0))+scale_x_continuous(expand = c(0, 0))
mappaIniziale



##mappa BASE AUS CON BATIMETRIA
mappaAUS<-basemap(limits = c(110.1, 155.1, -40.17927, -8.088012),bathymetry = T,grid.col = NA,land.col = "#eeeac4")+
  theme(legend.position="none") +labs(x = "Longitude", y = "Latitude")+
  geom_point(data = citta[1,], aes(x = Long, y = Lat), size = 15, shape=1.5, stroke = 1,colour="red")+
  theme(axis.text.x=element_text(colour="black"),axis.text.y =element_text(colour="black"))+ scale_y_continuous(expand = c(0, 0))
 
patchwork<-mappaAUS+mappaIniziale
patchwork  + plot_annotation(tag_levels = 'A')

mappaAUS + plot_spacer() + mappaIniziale + plot_layout(widths = c(5.2, -0.6 ,4.5),guides = "collect")& theme(legend.position = "none")
tiff("mappa2.tiff", units="in", width=5, height=5, res=300)
patchwork
dev.off()


## GLTO, stazioni per latitudine e nome (una schifezza)
basemap(limits = c(143, 155, -25, -15)
        ,bathymetry = F, grid.col = NA,land.col = "#eeeac4")+
  labs(x = "Longitude", y = "Latitude")+
  geom_point(data = MAPdata1, aes(x = Longitude, y = Lat,color=Latitude,))+
  geom_text(data = MAPdata1, aes(x = Long, y = Lat,label=Station))

#dayly T 





### scaricare dati temperatura:
#quello usato      ncdcOisst21Agg_LonPM180
###CORDINATE GLTO 143, 155, -25, -15

##funzione per impostare le COORDINATE

OISST_sub_dl <- function(time_df){
  OISST_dat <- griddap(x = "NOAA_DHW", 
                       url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                       time = c(time_df$start, time_df$end), 
                       zlev = c(0, 0),
                       latitude = c(-15, -25),
                       longitude = c(143, 155),
                       fields = "CRW_SST")$data %>% 
    mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>% 
    dplyr::rename(t = time, temp = CRW_SST) %>% 
    select(lon, lat, t, temp) %>% 
    na.omit()}


# per impostare le DATE
dl_years <- data.frame(date_index = 1:1,
                       start = as.Date("2016-08-16"),
                       end = as.Date("2016-08-16"))

#2014-11-29= OISST_data1
#2015-04-15   OISST_data2
#2016-04-16 OISST_data3
#2016-08-16 OISST_data4



# Download all of the data with one nested request

system.time(
  OISST_data <- dl_years %>% 
    group_by(date_index) %>% 
    group_modify(~OISST_sub_dl(.x)) %>% 
    ungroup() %>% 
    select(lon, lat, t, temp))


### MULTIMAPPA MAPPA SCELTA!!!
# COSTA GLTO segnalini neri senza latitudine

#novembre 2014

I2014<-basemap(limits = c(143, 155, -25, -15),grid.col = NA,land.col = "#783d01")+
  geom_spatial_polygon(data = as.data.frame(OISST_data), aes(x = lon, y = lat))+
  geom_tile(data = as.data.frame(OISST_data),aes(x = lon, y = lat,fill = temp))+
  scale_fill_viridis_c(limits = c(18,30)) +
  geom_point(data = MAPdata1, aes(x = Long, y = Lat)) +
  scale_size(guide="none")+
  scale_color_manual(values = c("#000000","##000000","##000000","##000000","##000000","##000000"))
  theme(legend.position = "bottom")


  #aprile 2015

  A2015<-basemap(limits = c(143, 155, -25, -15),grid.col = NA,land.col = "#783d01")+
    geom_spatial_polygon(data = as.data.frame(OISST_data), aes(x = lon, y = lat))+
    geom_tile(data = as.data.frame(OISST_data),aes(x = lon, y = lat,fill = temp))+
    scale_fill_viridis_c(limits = c(18,30)) +
    geom_point(data = MAPdata1, aes(x = Long, y = Lat)) +
    scale_size(guide="none")+
    scale_color_manual(values = c("#000000","##000000","##000000","##000000","##000000","##000000"))

  #aprile 2016
  A2016<-basemap(limits = c(143, 155, -25, -15),grid.col = NA,land.col = "#783d01")+
    geom_spatial_polygon(data = as.data.frame(OISST_data), aes(x = lon, y = lat))+
    geom_tile(data = as.data.frame(OISST_data),aes(x = lon, y = lat,fill = temp))+
    scale_fill_viridis_c(limits = c(18,30)) +
    geom_point(data = MAPdata1, aes(x = Long, y = Lat)) +
    scale_size(guide="none")+
    scale_color_manual(values = c("#000000","##000000","##000000","##000000","##000000","##000000"))

  #agosto 2016
  AG2016<-basemap(limits = c(143, 155, -25, -15),grid.col = NA,land.col = "#783d01")+
    geom_spatial_polygon(data = as.data.frame(OISST_data), aes(x = lon, y = lat))+
    geom_tile(data = as.data.frame(OISST_data),aes(x = lon, y = lat,fill = temp))+
    scale_fill_viridis_c(limits = c(18,30)) +
    geom_point(data = MAPdata1, aes(x = Long, y = Lat)) +
    scale_size(guide="none")+
    scale_color_manual(values = c("#000000","##000000","##000000","##000000","##000000","##000000"))

  
  
  ## unire le 4 mappe

  c<-ggarrange(I2014, A2015 , A2016,AG2016,
               labels = c("A", "B","C","D"),
               ncol = 2, nrow = 2,
               legend = "bottom",
               common.legend=T	)
  
  annotate_figure(c,top = text_grob("PCoA on Bray-Curtis Dissimilarity \n Permanova p = 0.001", color = "Black", face = "bold", size = 14))
  
  
  
  
  
  
  
  

  ###  MAPPA SCELTA!!!
  # COSTA GLTO segnalini neri senza latitudine
  basemap(limits = c(143, 155, -25, -15),grid.col = NA,land.col = "#783d01")+
    geom_spatial_polygon(data = as.data.frame(OISST_data), aes(x = lon, y = lat))+
    geom_tile(data = as.data.frame(OISST_data),aes(x = lon, y = lat,fill = temp))+
    scale_fill_viridis_c(limits = c(18,30)) +
    geom_point(data = MAPdata1, aes(x = Longitude, y = Lat,size=2,shape=Tow)) +
    scale_shape_manual(values=c(0, 1, 2,6))+scale_size(guide="none")+
    scale_color_manual(values = c("#000000","##000000","##000000","##000000","##000000","##000000"))
  theme(legend.position = "bottom")


####mappa con temperatura
### COSTA queensland
basemap(limits = c(141, 160, -29.17927, -9.088012),grid.col = NA,land.col = "#783d01")+
  geom_spatial_polygon(data = as.data.frame(OISST_data), aes(x = lon, y = lat))+
  geom_tile(data = as.data.frame(OISST_data),aes(x = lon, y = lat,fill = temp))+
  scale_fill_viridis_c(limits = c(18,30)) +
  geom_point(data = MAPdata1, aes(x = Long, y = Lat,color=Latitude)) + scale_size(guide="none")




## COSTA GLTO
basemap(limits = c(143, 155, -25, -15),grid.col = NA,land.col = "#783d01")+
  geom_spatial_polygon(data = as.data.frame(OISST_data), aes(x = lon, y = lat))+
  geom_tile(data = as.data.frame(OISST_data),aes(x = lon, y = lat,fill = temp))+
  scale_fill_viridis_c(limits = c(18,30)) +
  geom_point(data = MAPdata1, aes(x = Long, y = Lat,color=Latitude,)) + scale_size(guide="none")






###MAPPA CON TEMPERATURA E STAZIONI
#QUEENSLAND
c(137.99465, 153.611604
  -29.17927,  -9.088012)



### COSTA queensland#
#  basemap(limits = c(137.99465, 155, -29.17927, -9.088012),grid.col = NA,land.col = "#783d01")+
  geom_spatial_polygon(data = as.data.frame(OISST_data), aes(x = lon, y = lat))+
  geom_tile(data = as.data.frame(OISST_data),aes(x = lon, y = lat,fill = temp))+
  scale_fill_viridis_c(limits = c(20,30)) +
  geom_point(data = MAPdata1, aes(x = Long, y = Lat,color=Latitude,)) + scale_size(guide="none")













####################    metagenome enrichemnt##########
  
  #####                      importing data                           ####
  abund_table_enrich<-read.csv("input/target_enrichment_taxa_table_only_vibrio.csv",row.names=1,check.names=FALSE)
  OTU_taxonomy_enrich<-read.csv("input/target_enrichment_taxonomy_only_vibrio.csv",row.names=1,check.names=FALSE)
  
  #Convert the data to phyloseq format
  OTU_enrich = otu_table(as.matrix(abund_table_enrich), taxa_are_rows = T)
  TAX_enrich = tax_table(as.matrix(OTU_taxonomy_enrich))
  physeq_enrich<-merge_phyloseq(phyloseq(OTU_enrich, TAX_enrich))
  
 p<- plot_bar(physeq_enrich, "Species", fill="Species")

 p + guides(fill=guide_legend(title="Vibrio species")) +ggtitle("Vibriome") +theme(plot.title = element_text(hjust = 0.5))


















#####BIOENV


variabili_bioenv<-env_data[,c(3:9,14:22,23:28,31:34)]
as.data.frame(colnames(variabili_bioenv))
#write.csv(as.data.frame(colnames(variabili_bioenv)), "variabili_BIOENV.csv")
















#bioenv comunit? batterica CON MET NORMALIZZATI





bioenv_comu_batt<-bioenv(t(OTU_DF_BIOENV),variabili_bioenv, method = "spearman", index = "bray")

bat.dist_bact<-vegdist(x = t(OTU_DF_BIOENV),method = "bray")
env.dist_bact<-dist(x =variabili_bioenv[,c(3,16,20)],method = "euclidean")
set.seed(1234)
mantel(xdis = bat.dist_bact,ydis = env.dist_bact,method = "spearman",permutations = 9999)
#p=0.0311










##############  RDA:                            #####
physeq_RDA<-subset_samples(phyloseq_obj_css, station != "5" | site != "H")


#!!!qui ? scomparsa la traformazione hellinger, la rifaccio !!!!
OTU_TABLE_hel_trasf = decostand(otu_table(phyloseq_obj_css), "hel")
physeq_RDA<-phyloseq_obj_css

otu_table(physeq_RDA) <- otu_table(OTU_TABLE_hel_trasf, taxa_are_rows = T)
otu_table(physeq_RDA)==otu_table(phyloseq_obj_css)


#RDA model, con variabili che vengono furoi da BIOENV/relate

df <- otu_table(physeq_RDA)
head(df)
# usiamo qui il file scalato su alluminio, senza Mo, Mn, DIP e DIN

RDA.model <- rda(t(df[,-12])~ OC*IdInp*Cu, poll3_Al[-12,], scale = T)


#test dei modelli tramite ordistep:

#non capisco bene come funziona!!!!!!

ordistep(object = RDA.model,direction = "forward",permutations = 9999)






palette(allGroupsColors)
allGroupsColors
plot(RDA.model, type = "n")
points(RDA.model, pch=21, col="black", bg=poll3_Al[-12,]$site, cex=2.5)
text(RDA.model, dis="cn")
palette(allGroupsColors)
text(RDA.model, "sites", col="black", cex=0.8)

#anova per constrained models 
set.seed(1234)
anova.cca(RDA.model)
anova.cca(RDA.model,by="margin")
anova.cca(RDA.model,by="term")


#VARIATION PARTITION:
var.part<-varpart(Y =t(df[,-12]), X = ~ OC,~IdInp,~Cu,data = poll3_Al[-12,])
plot(var.part, digit=2)
legend("topleft", legend=c("X1=OC", "X2=IdInp","X3=Cu"))#, lty=1:2, cex=0.8)

# # significance of partitions
anova.cca(rda(t(df[,-12]), poll3_Al[-12,]$OC, step=1000))
anova.cca(rda(t(df[,-12]), poll3_Al[-12,]$IdInp, step=1000))
anova.cca(rda(t(df[,-12]), poll3_Al[-12,]$Cu, step=1000))
summary(var.part)



# 
ordination.rda<-ordinate(physeq_RDA_SRB, formula = ~ T_S*IdInp*Cu, poll3_Al, method = "RDA")
# 
# 
# ## overall test
# anova.cca(ordination.rda)
# ## Test for axes
# anova.cca(ordination.rda, by="axis")
# ## Sequential test for terms
# anova.cca(ordination.rda, by="terms")
# ## Marginal effects
# anova.cca(ordination.rda, by="margin")
# 
# 
# #plottare RDA
p0=plot_ordination(physeq_RDA_SRB, ordination.rda, color = "site" )+ scale_color_manual(values = allGroupsColors)+  geom_point( size = 4)+ geom_text(mapping = aes(label = station), size = 4, hjust = 2) #,geom_text(label = "station", size = 5)
# #
#???p0=plot_ordination(physeq_normalized, ordination.cca, color = "site",label = "station")
# 
# 
# 
# # Now add the environmental variables as arrows
arrowmat = vegan::scores(ordination.rda, display = "bp")
# # Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# # Define the arrow aesthetic mapping
arrow_map = aes(xend = RDA1, yend = RDA2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 1.2 * RDA1, y = 1.2 * RDA2, shape = NULL, color = NULL, 
                label = labels)
# # Make a new graphic
arrowhead = arrow(length = unit(0.05, "npc"))
p1 =p0 + geom_segment(arrow_map, size = 1, data = arrowdf, color = "black", 
                      arrow = arrowhead)  +theme_bw()#+ geom_text(label_map, size = 3, data = arrowdf, vjust = 0,hjust=-0.5)

p1 + geom_text_repel(label_map, size = 2.5, data = arrowdf, vjust = 0,hjust=-0,force = 1)
#+ scale_x_reverse()
# 
# 
# #per aumentare la grandezza delle frecce 
# arrowhead = arrow(length = unit(0.05, "npc"))
# p1 = p0 + geom_segment(arrow_map, size =1, data = arrowdf, color = "black", 
#                        arrow = arrowhead) + geom_text(label_map, size = 3, data = arrowdf)
# p1+ scale_x_reverse()
# 
# 
# 
# 
# 
# 
# 
# #-> rda at phylum level
# 
# 
# p0=plot_ordination(physeq_RDA, ordination.rda, color = "phylum",type = "taxa")
# 
# 
# p0=plot_ordination(physeq_RDA, ordination.rda, color = "site",type = "biplot",label = "station")#+ scale_color_manual(values = allGroupsColors)
# 
# #per ingrandire le freccie
# # Define the arrow aesthetic mapping
# arrow_map = aes(xend = RDA1, yend = RDA2, x = 0, y = 0, shape = NULL, color = NULL, 
#                 label = labels)
# label_map = aes(x = 1.2 * RDA1, y = 1.2 * RDA2, shape = NULL, color = NULL, 
#                 label = labels)#+ scale_color_manual(values = allGroupsColors)
# # Make a new graphic
# arrowhead = arrow(length = unit(0.05, "npc"))
# p1 = p0 + geom_segment(arrow_map, size = 1, data = arrowdf, color = "black", 
#                        arrow = arrowhead) + geom_text(label_map, size = 5, data = arrowdf)
# p1



###########    RDA SRB
SRB_physeq_silva_woH5 
SRB_physeq_silva



OTU_TABLE_hel_trasf_SRB = decostand(otu_table(SRB_physeq_silva_woH5), "hel")
physeq_RDA_SRB<-SRB_physeq_silva_woH5

otu_table(physeq_RDA_SRB) <- otu_table(OTU_TABLE_hel_trasf_SRB, taxa_are_rows = T)
otu_table(physeq_RDA_SRB)==otu_table(SRB_physeq_silva_woH5)


#RDA model, con variabili che vengono furoi da BIOENV/relate

df_SRB <- otu_table(physeq_RDA_SRB)
head(df)
# usiamo qui il file scalato su alluminio, senza Mo, Mn, DIP e DIN

#RDA.model_SRB <- rda(t(df_SRB[,-12])~ OC*IdInp*Cu, poll3_Al[-12,], scale = T)


RDA.model_SRB <- rda(t(df_SRB)~ T_S*IdInp*Cu, poll3_Al[-12,], scale = T)



# allGroupsColors
# plot(RDA.model_SRB, type = "n")
# text(RDA.model_SRB, dis="cn")
# palette(allGroupsColors)
# points(RDA.model_SRB, pch=21, col="black", bg=poll3_Al$site, cex=2.5)
# text(RDA.model_SRB, "sites", col="blue", cex=0.8)

palette(allGroupsColors)
allGroupsColors
plot(RDA.model_SRB, type = "n")
points(RDA.model_SRB, pch=21, col="black", bg=poll3_Al[-12,]$site, cex=2.5)
text(RDA.model_SRB, dis="cn")
palette(allGroupsColors)
text(RDA.model_SRB, "sites", col="black", cex=0.8)



#anova per constrained models 
set.seed(1234)
anova.cca(RDA.model_SRB)
anova.cca(RDA.model_SRB,by="margin")
anova.cca(RDA.model_SRB,by="term")


#VARIATION PARTITION:
var.part<-varpart(Y =t(df[,-12]), X = ~ OC,~IdInp,~Cu,data = poll3_Al[-12,])
plot(var.part, digit=2)
legend("topleft", legend=c("X1=OC", "X2=IdInp","X3=Cu"))#, lty=1:2, cex=0.8)

# # significance of partitions
anova.cca(rda(t(df[,-12]), poll3_Al[-12,]$OC, step=1000))
anova.cca(rda(t(df[,-12]), poll3_Al[-12,]$IdInp, step=1000))
anova.cca(rda(t(df[,-12]), poll3_Al[-12,]$Cu, step=1000))
summary(var.part)







####  MANTEL E PARTIAL MANTER SRB CON SILVA!!!!!!

SRB_physeq_silva = subset_taxa(phyloseq_obj_css, family=="D_4__Desulfobacteraceae" |family=="D_4__Desulfarculaceae" |family=="D_4__Syntrophaceae" |family=="D_4__Desulfobulbaceae" |family=="D_4__Syntrophobacteraceae")

# Extract abundance matrix from the phyloseq object
OTU_matrix_SRB_silva = as(otu_table(subset_samples(SRB_physeq_silva, station != "5" | site != "H")), "matrix")
OTU_SRB_df = as.data.frame(OTU_matrix_SRB_silva)
head(OTU_matrix_SRB_silva)



#Costruisco le dist-matrix per analisi mantel e partial mantel:
B_C_DIST_SRB_silva<-vegdist(x = t(OTU_matrix_SRB_silva),method = "bray")
ind.PAH.dist_1<-dist(scale(x =ind.sources.emission_1[-12,]),method = "euclidean")
#ind.PAH.dist_2<-dist(scale(x =ind.sources.emission_2),method = "euclidean")
coord_dist<-dist(x =coordinates[-12,],method = "euclidean")

rownames(as.matrix(B_C_DIST_SRB_silva)) == rownames(as.matrix(ind.PAH.dist_1))
rownames(as.matrix(B_C_DIST_SRB_silva)) == rownames(as.matrix(coord_dist))

#mantel fatta con vegan SRB_silva contro indici
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_SRB_silva,
              ydis = ind.PAH.dist_1,
              method = "pearson",
              permutations = 9999)

# R = 0.5087
# p = 0.0018

#mantel fatta con vegan SRB_silva contro coordinate
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_SRB_silva,
              ydis = coord_dist,
              method = "pearson",
              permutations = 9999)
# R = 0.345  
# p = 0.0262

#partial.mantel fatta con vegan SRB_silva contro indici controllando coordinate

set.seed(1234)
vegan::mantel.partial(xdis = B_C_DIST_SRB_silva,
                      ydis = ind.PAH.dist_1,
                      zdis =coord_dist,
                      method = "pearson",
                      permutations = 9999 )
# R = 0.4608 
# p = 0.0045










########## MANTEL E PARTIAL MANTEL CON INTERA COMUNITA      ####


### MANTEL E PARTIAL MANTEL CON INTERA COMUNITA

# Extract abundance matrix from the phyloseq object
OTU_matrix.ALL = as(otu_table(phyloseq_obj_css), "matrix")
OTU_ALL_df = as.data.frame(OTU_matrix.ALL)
head(OTU_ALL_df)

#Costruisco le dist-matrix per analisi mantel e partial mantel:
B_C_DIST_ALL<-vegdist(x = t(OTU_ALL_df[,-12]),method = "bray")
ind.PAH.dist_1<-dist(scale(x =ind.sources.emission_1[-12,]),method = "euclidean")
#ind.PAH.dist_2<-dist(scale(x =ind.sources.emission_2),method = "euclidean")
coord_dist<-dist(x =coordinates[-12,],method = "euclidean")

rownames(as.matrix(B_C_DIST_ALL)) == rownames(as.matrix(ind.PAH.dist_1))
rownames(as.matrix(B_C_DIST_ALL)) == rownames(as.matrix(coord_dist))

#mantel fatta con vegan TUTTI contro indici
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_ALL,
              ydis = ind.PAH.dist_1,
              method = "pearson",
              permutations = 9999)

# R = 0.4045 
# p = 0.0072

# NB rispetto a T-RFLP mantel è piu "sensibile" e non tornava significativo senza rimuovere H5

#mantel fatta con vegan TUTTI contro coordinate
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_ALL,
              ydis = coord_dist,
              method = "pearson",
              permutations = 9999)

# R = 0.4794 
# p = 0.0057

#partial.mantel fatta con vegan TUTTI contro indici controllando coordinate

set.seed(1234)
vegan::mantel.partial(xdis = B_C_DIST_ALL,
                      ydis = ind.PAH.dist_1,
                      zdis =coord_dist,
                      method = "pearson",
                      permutations = 9999 )

# R = 0.3273 
# p = 0.0267





####  NB un "errore" che ci vedo, è che l'analisi con tutti comprende ANCHE gli SRB
#### per questo secondo me con lo stesso metodo e logica sarebbe da provare
#### tutta la comunità eccetto SRB

### MANTEL E PARTIAL MANTEL CON INTERA COMUNITA - SRB

# Extract abundance matrix from the phyloseq object

OTU_matrix.ALL = as(otu_table(subset_samples(phyloseq_obj_css)), "matrix")#, station != "5" | site != "H" )), "matrix")
dim(OTU_matrix.ALL) # 2482 OTU

dim(OTU_matrix_SRB_silva) # 321 OTU

dim(OTU_matrix.ALL) - dim(OTU_matrix_SRB_silva) # 2161 OTU

OTU_matrix.noSRB <- subset(OTU_matrix.ALL, !(rownames(OTU_matrix.ALL) %in% rownames(OTU_matrix_SRB_silva)))
dim(OTU_matrix.noSRB) # 2161 OTU


#Costruisco le dist-matrix per analisi mantel e partial mantel:
B_C_DIST_noSRB<-vegdist(x = t(OTU_matrix.noSRB),method = "bray")
ind.PAH.dist_1<-dist(scale(x =ind.sources.emission_1[-12,]),method = "euclidean")
#ind.PAH.dist_2<-dist(scale(x =ind.sources.emission_2),method = "euclidean")
coord_dist<-dist(x =coordinates[-12,],method = "euclidean")

rownames(as.matrix(B_C_DIST_noSRB)) == rownames(as.matrix(ind.PAH.dist_1))
rownames(as.matrix(B_C_DIST_noSRB)) == rownames(as.matrix(coord_dist))

#mantel fatta con vegan TUTTI contro indici
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_noSRB,
              ydis = ind.PAH.dist_1,
              method = "pearson",
              permutations = 9999)

# R = 0.3794 
# p = 0.0127

# NB rispetto a T-RFLP mantel è piu "sensibile" e non tornava significativo senza rimuovere H5

#mantel fatta con vegan TUTTI contro coordinate
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_noSRB,
              ydis = coord_dist,
              method = "pearson",
              permutations = 9999)
# R = 0.5081 
# p = 0.0057 

#partial.mantel fatta con vegan TUTTI contro indici controllando coordinate

set.seed(1234)
vegan::mantel.partial(xdis = B_C_DIST_noSRB,
                      ydis = ind.PAH.dist_1,
                      zdis =coord_dist,
                      method = "pearson",
                      permutations = 9999 )

# R = 0.2942 
# p = 0.0361 





#MANTEL E PARTIAL MANTEL TRFLP terza stagione 


#import TRFLP Bacteria
abund_table_trflp_bacteria<-read.csv("TRFLP_BACTERIA_TERZA_STAGIONE.csv",row.names=1,check.names=FALSE)
abund_table_trflp_bacteria<-t(abund_table_trflp_bacteria)


#mantel trflp terza stagione BACTERIA:


#costruisco le dist-matrix per analisi mantel e partial mantel:

B_C_DIST_TRFLP_BACTERIA<-vegdist(x = (as.data.frame(abund_table_trflp_bacteria[-11,])),method = "bray")
ind.PAH.dist_1<-dist(scale(x =ind.sources.emission_1[-c(7,12),]),method = "euclidean")
coord_dist<-dist(x =coordinates[-c(7,12),],method = "euclidean")



rownames(as.matrix(B_C_DIST_TRFLP_BACTERIA)) == rownames(as.matrix(ind.PAH.dist_1))
rownames(as.matrix(B_C_DIST_TRFLP_BACTERIA)) == rownames(as.matrix(coord_dist))

#mantel fatta con vegan trflp BACTERIA contro indici
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_TRFLP_BACTERIA,
              ydis = ind.PAH.dist_1,
              method = "pearson",
              permutations = 9999)

# R = 0.1871
# p = 0.1526

#mantel fatta con vegan trflp BACTERIA contro coordinate
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_TRFLP_BACTERIA,
              ydis = coord_dist,
              method = "pearson",
              permutations = 9999)
# R = 0.4212
# p = 0.0226

#partial.mantel fatta con vegan trflp BACTERIA contro indici controllando coordinate

set.seed(1234)
vegan::mantel.partial(xdis = B_C_DIST_TRFLP_BACTERIA,
                      ydis = ind.PAH.dist_1,
                      zdis =coord_dist,
                      method = "pearson",
                      permutations = 9999 )
# R = 0.05106 
# p = 0.373 




#import TRFLP dsrAB
abund_table_trflp_dsrab<-read.csv("TRFLP_dsrAB_TERZA_STAGIONE.csv",row.names=1,check.names=FALSE)
abund_table_trflp_dsrab<-t(abund_table_trflp_dsrab)


#costruisco le dist-matrix per analisi mantel e partial mantel:

B_C_DIST_TRFLP_DSRAB<-vegdist(x = (as.data.frame(abund_table_trflp_dsrab[-12,])),method = "bray")
ind.PAH.dist_1<-dist(scale(x =ind.sources.emission_1[-12,]),method = "euclidean")
coord_dist<-dist(x =coordinates[-12,],method = "euclidean")



rownames(as.matrix(B_C_DIST_TRFLP_DSRAB)) == rownames(as.matrix(ind.PAH.dist_1))
rownames(as.matrix(B_C_DIST_TRFLP_DSRAB)) == rownames(as.matrix(coord_dist))

#mantel fatta con vegan trflp DSRAB contro indici
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_TRFLP_DSRAB,
              ydis = ind.PAH.dist_1,
              method = "pearson",
              permutations = 9999)

# R = 0.4303
# p = 0.0077

#mantel fatta con vegan trflp DSRAB contro coordinate
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_TRFLP_DSRAB,
              ydis = coord_dist,
              method = "pearson",
              permutations = 9999)
# R = 0.6631
# p = 1e-04

#partial.mantel fatta con vegan trflp DSRAB contro indici controllando coordinate

set.seed(1234)
vegan::mantel.partial(xdis = B_C_DIST_TRFLP_DSRAB,
                      ydis = ind.PAH.dist_1,
                      zdis =coord_dist,
                      method = "pearson",
                      permutations = 9999 )
# R = 0.3518 
# p = 0.0183 



#import TRFLP ARCHEA
abund_table_trflp_archaea<-read.csv("TRFLP_ARCHAEA_TERZA_STAGIONE.csv",row.names=1,check.names=FALSE)
abund_table_trflp_archaea<-t(abund_table_trflp_archaea)

#costruisco le dist-matrix per analisi mantel e partial mantel:

B_C_DIST_TRFLP_ARCHEA<-vegdist(x = (as.data.frame(abund_table_trflp_archaea[-11,])),method = "bray")
ind.PAH.dist_1<-dist(scale(x =ind.sources.emission_1[c(-7,-12),]),method = "euclidean")
coord_dist<-dist(x =coordinates[c(-7,-12),],method = "euclidean")



rownames(as.matrix(B_C_DIST_TRFLP_ARCHEA)) == rownames(as.matrix(ind.PAH.dist_1))
rownames(as.matrix(B_C_DIST_TRFLP_ARCHEA)) == rownames(as.matrix(coord_dist))

#mantel fatta con vegan trflp ARCHEA contro indici
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_TRFLP_ARCHEA,
              ydis = ind.PAH.dist_1,
              method = "pearson",
              permutations = 9999)

# R = 0.1055
# p = 0.2321

#mantel fatta con vegan trflp ARCHEA contro coordinate
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_TRFLP_ARCHEA,
              ydis = coord_dist,
              method = "pearson",
              permutations = 9999)
# R = 0.4955
# p = 0.0096

#partial.mantel fatta con vegan trflp ARCHEA contro indici controllando coordinate

set.seed(1234)
vegan::mantel.partial(xdis = B_C_DIST_TRFLP_ARCHEA,
                      ydis = ind.PAH.dist_1,
                      zdis =coord_dist,
                      method = "pearson",
                      permutations = 9999 )
# R = -0.07771
# p = 0.6704 





##mantel e partial mantel   bioenv   NORMALIZZATO ####

# Extract abundance matrix from the phyloseq object
OTU_matrix.ALL = as(otu_table(phyloseq_obj_css), "matrix")
OTU_ALL_df = as.data.frame(OTU_matrix.ALL)
head(OTU_ALL_df)

#Costruisco le dist-matrix per analisi mantel e partial mantel:
B_C_DIST_ALL_bioenv<-vegdist(x = t(OTU_ALL_df[,-12]),method = "bray")
bioenv_bact_normalizzati<-dist(scale(x =variabili_bioenv[,c(3,16,20)]),method = "euclidean")


coord_dist<-dist(x =coordinates[-12,],method = "euclidean")

rownames(as.matrix(B_C_DIST_ALL_bioenv)) == rownames(as.matrix(bioenv_bact_normalizzati))
rownames(as.matrix(B_C_DIST_ALL_bioenv)) == rownames(as.matrix(coord_dist))

#mantel fatta con vegan TUTTI contro indici
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_ALL_bioenv,
              ydis = bioenv_bact_normalizzati,
              method = "pearson",
              permutations = 9999)

# R =  0.8463 
# p = 1e-04


#mantel fatta con vegan TUTTI contro coordinate
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_ALL_bioenv,
              ydis = coord_dist,
              method = "pearson",
              permutations = 9999)

# R =  0.4794 
# p = 0.0053

#partial.mantel fatta con vegan TUTTI contro bioenv controllando coordinate

set.seed(1234)
vegan::mantel.partial(xdis = B_C_DIST_ALL_bioenv,
                      ydis = bioenv_bact_normalizzati,
                      zdis =coord_dist,
                      method = "pearson",
                      permutations = 9999 )

# R = 0.8038
# p = 1e-04






#BACTNON NON NORM

bioenv_bact_NO_normalizzati<-dist(scale(x =variabili_bioenv_met_non_norm[,c(1,3,16,20,21)]),method = "euclidean")

rownames(as.matrix(B_C_DIST_ALL_bioenv)) == rownames(as.matrix(bioenv_bact_NO_normalizzati))
rownames(as.matrix(B_C_DIST_ALL_bioenv)) == rownames(as.matrix(coord_dist))

#mantel fatta con vegan TUTTI contro indici
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_ALL_bioenv,
              ydis = bioenv_bact_NO_normalizzati,
              method = "pearson",
              permutations = 9999)

# R =   0.8615 
# p = 1e-04


#mantel fatta con vegan TUTTI contro coordinate
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_ALL_bioenv,
              ydis = coord_dist,
              method = "pearson",
              permutations = 9999)

# R =  0.4794
# p = 0.0053

#partial.mantel fatta con vegan TUTTI contro bioenv controllando coordinate

set.seed(1234)
vegan::mantel.partial(xdis = B_C_DIST_ALL_bioenv,
                      ydis = bioenv_bact_NO_normalizzati,
                      zdis =coord_dist,
                      method = "pearson",
                      permutations = 9999 )

# R = 0.8183
# p = 1e-04




###
####  MANTEL E PARTIAL MANTER SRB CON SILVA!!!!!!

SRB_physeq_silva = subset_taxa(phyloseq_obj_css, family=="D_4__Desulfobacteraceae" |family=="D_4__Desulfarculaceae" |family=="D_4__Syntrophaceae" |family=="D_4__Desulfobulbaceae" |family=="D_4__Syntrophobacteraceae")

# Extract abundance matrix from the phyloseq object
OTU_matrix_SRB_silva = as(otu_table(subset_samples(SRB_physeq_silva, station != "5" | site != "H")), "matrix")
OTU_SRB_df = as.data.frame(OTU_matrix_SRB_silva)
head(OTU_matrix_SRB_silva)



#Costruisco le dist-matrix per analisi mantel e partial mantel:
B_C_DIST_SRB_silva<-vegdist(x = t(OTU_matrix_SRB_silva),method = "bray")


BIOENV_SRB_NORMALIZZATO<-dist(scale(x =variabili_bioenv[,c(1,16,20)]),method = "euclidean")



rownames(as.matrix(B_C_DIST_SRB_silva)) == rownames(as.matrix(BIOENV_SRB_NORMALIZZATO))
rownames(as.matrix(B_C_DIST_SRB_silva)) == rownames(as.matrix(coord_dist))

#mantel fatta con vegan SRB_silva contro indici
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_SRB_silva,
              ydis = BIOENV_SRB_NORMALIZZATO,
              method = "pearson",
              permutations = 9999)

# R = 0.8438
# p = 2e-04

#mantel fatta con vegan SRB_silva contro coordinate
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_SRB_silva,
              ydis = coord_dist,
              method = "pearson",
              permutations = 9999)
# R =  0.3450


# p = 0.0262

#partial.mantel fatta con vegan SRB_silva contro indici controllando coordinate

set.seed(1234)
vegan::mantel.partial(xdis = B_C_DIST_SRB_silva,
                      ydis = BIOENV_SRB_NORMALIZZATO,
                      zdis =coord_dist,
                      method = "pearson",
                      permutations = 9999 )
# R = 0.8378
# p = 2e-04



#srb NON NORM
env.dist_srb<-dist(x =variabili_bioenv_met_non_norm[,c(1,16,20:22)],method = "euclidean")


BIOENV_SRB_NO_NORMALIZZATO<-dist(scale(x =variabili_bioenv_met_non_norm[,c(1,16,20:22)]),method = "euclidean")



rownames(as.matrix(B_C_DIST_SRB_silva)) == rownames(as.matrix(BIOENV_SRB_NORMALIZZATO))
rownames(as.matrix(B_C_DIST_SRB_silva)) == rownames(as.matrix(coord_dist))

#mantel fatta con vegan SRB_silva contro indici
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_SRB_silva,
              ydis = BIOENV_SRB_NO_NORMALIZZATO,
              method = "pearson",
              permutations = 9999)

# R = 0.8706
# p = 2e-04 

#mantel fatta con vegan SRB_silva contro coordinate
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_SRB_silva,
              ydis = coord_dist,
              method = "pearson",
              permutations = 9999)
# R =  0.3450


# p = 0.0262

#partial.mantel fatta con vegan SRB_silva contro indici controllando coordinate

set.seed(1234)
vegan::mantel.partial(xdis = B_C_DIST_SRB_silva,
                      ydis = BIOENV_SRB_NO_NORMALIZZATO,
                      zdis =coord_dist,
                      method = "pearson",
                      permutations = 9999 )
# R = 0.8531 
# p =  2e-04 





variabili_bioenv_met_non_norm[,c(1,3,16,20,21)]


#Costruisco le dist-matrix per analisi mantel e partial mantel:
B_C_DIST_noSRB<-vegdist(x = t(OTU_matrix.noSRB),method = "bray")


BIOENV_batt_NO_SRB_NO_NORMALIZZATO<-dist(scale(x =variabili_bioenv_met_non_norm[,c(1,3,16,20,21)]),method = "euclidean")




#mantel fatta con vegan SRB_silva contro indici
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_noSRB,
              ydis = BIOENV_batt_NO_SRB_NO_NORMALIZZATO,
              method = "pearson",
              permutations = 9999)

# R = 0.8534 
# p = 1e-04 

#mantel fatta con vegan contro coordinate
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_noSRB,
              ydis = coord_dist,
              method = "pearson",
              permutations = 9999)
# R =  0.4994


# p = 0.0066

#partial.mantel fatta con vegan SRB_silva contro indici controllando coordinate

set.seed(1234)
vegan::mantel.partial(xdis = B_C_DIST_noSRB,
                      ydis = BIOENV_batt_NO_SRB_NO_NORMALIZZATO,
                      zdis =coord_dist,
                      method = "pearson",
                      permutations = 9999 )
# R = 0.7993 
# p = 1e-04





## batt No SRB MET NORMALIZZATI
BIOENV_batt_NO_SRB_NORMALIZZATO<-dist(scale(x =variabili_bioenv[,c(3,16,20)]),method = "euclidean")



#mantel fatta con vegan SRB_silva contro indici
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_noSRB,
              ydis = BIOENV_batt_NO_SRB_NORMALIZZATO,
              method = "pearson",
              permutations = 9999)

# R = 0.836 
# p = 1e-04 

#mantel fatta con vegan contro coordinate
set.seed(1234)
vegan::mantel(xdis = B_C_DIST_noSRB,
              ydis = coord_dist,
              method = "pearson",
              permutations = 9999)
# R =  0.4994


# p = 0.0066

#partial.mantel fatta con vegan SRB_silva contro indici controllando coordinate

set.seed(1234)
vegan::mantel.partial(xdis = B_C_DIST_noSRB,
                      ydis = BIOENV_batt_NO_SRB_NORMALIZZATO,
                      zdis =coord_dist,
                      method = "pearson",
                      permutations = 9999 )
# R = 0.789 
# p = 1e-04


## 14/01/20: correlazione tra variabili uscite da bioenv e lat/long NORMALIZZATE ad al

# input di variabili
colnames(ENV_table_scaled_noAl) == colnames(ENV_table_scaled_unnorm)

colnames(ENV_table_scaled_noAl)
env_var_correl <- ENV_table_scaled_noAl[-12, c(3,20,25)]

env_met_correl <- ENV_table_scaled_unnorm[-12, 21:31]

# input di coordinate
coordinate <- read.csv("COORDINATE.csv", row.names=1,check.names=FALSE)
coordinate <- coordinate[-12,]

## unisco

correl_df <- cbind(coordinate, env_var_correl)

library("PerformanceAnalytics")
chart.Correlation(correl_df[,-6], histogram=TRUE, method = "spearman")
chart.Correlation(correl_df[,-6], histogram=TRUE, method = "pearson")


pairs(correl_df)

correl_df$Site <- factor(c(rep("C",5), rep("E",3), rep("H",3)))
p1 <- ggplot(correl_df, aes(x=latitude, y= OC)) + geom_point(aes(col = Site)) +  geom_smooth(method = lm, linetype = 2, se = F)
cor.test(correl_df$latitude,correl_df$OC, method = "spearman")
p3 <- ggplot(correl_df, aes(x=latitude, y= IdInp )) + geom_point(aes(col = Site)) +  geom_smooth(method = lm, linetype = 2, se = F)
cor.test(correl_df$latitude,correl_df$IdInp, method = "spearman")
p5 <- ggplot(correl_df, aes(x=latitude, y= Cu)) + geom_point(aes(col = Site)) + geom_smooth(method = lm, linetype = 2, se = F)
cor.test(correl_df$latitude,correl_df$Cu, method = "spearman")
p2 <- ggplot(correl_df, aes(x=longitude, y= OC)) + geom_point(aes(col = Site)) +  geom_smooth(method = lm, linetype = 2, se = F)
cor.test(correl_df$latitude,correl_df$OC, method = "spearman")
p4 <- ggplot(correl_df, aes(x=longitude, y= IdInp )) + geom_point(aes(col = Site)) +  geom_smooth(method = lm, linetype = 2, se = F)
cor.test(correl_df$latitude,correl_df$IdInp, method = "spearman")
p6 <- ggplot(correl_df, aes(x=longitude, y= Cu)) + geom_point(aes(col = Site)) +  geom_smooth(method = lm, linetype = 2, se = F)
cor.test(correl_df$latitude,correl_df$Cu, method = "spearman")

library(cowplot)
plot_grid(p1,p2,p3,p4,p5,p6, ncol = 2)


## con i metalli

correl_df <- cbind(coordinate, env_met_correl)

library("PerformanceAnalytics")
chart.Correlation(correl_df[,-6], histogram=F, method = "spearman")
chart.Correlation(correl_df[,-6], histogram=F, method = "pearson")



# regioni <- c("Italy", "Greece", "Tunisia","Crete","Albania","Croatia", "Egypt",
#              "Egitto","Israel")
# usa <- map_data("world", region = regioni) # we already did this, but we can do it again
# ggplot() + 
#   geom_polygon(data = usa, aes(x = long, y = lat, group = group)) + 
#   coord_quickmap()


### modello MVABUND ####

## MVabund # problema è che questo pacchetto mi va in conflitto con un altro
## pacchetto che è installato dal sistema, tipo rcpp. 

# dopo fa un test univariato per ogni OTU, per questo è lento. 
# Forse sarebbe meglio da lavorare sui generi anzi che le OTU?

# library(mvabund)
# tabella_mvabb <- otu_table(phyloseq_obj_css, taxa_are_rows = F)
# tabella_mvabb<-tabella_mvabb[,-12]
# max(tabella_mvabb)
# tabella_mvabb <- as.data.frame(tabella_mvabb@.Data)
# #rimuovo taxa a 0
# tabella_mvabb <- tabella_mvabb[rowSums(tabella_mvabb) != 0,]
# 
# #rimuovo taxa che hanno 0 al 20% dei campioni. Da verificare se soglia è ok
# min <- round(ncol(tabella_mvabb )*0.2)
# tabella_mvabb <- tabella_mvabb[rowSums(tabella_mvabb != 0) > min, ]
# 
# # rimodello i livelli tassonomici
# taxa_mvabb <- as.data.frame(tax_table(phyloseq_obj_css))
# taxa_mvabb <- taxa_mvabb[row.names(taxa_mvabb) %in% row.names(tabella_mvabb), ]
# row.names(taxa_mvabb) == row.names(tabella_mvabb)
# nrow(taxa_mvabb)
# tax_vector <- paste(row.names(taxa_mvabb),sep="|",paste("f_", sep = "", taxa_mvabb$family))
# tax_vector <- paste(tax_vector, sep = "|", paste("g_",taxa_mvabb$genus, sep=""))
# row.names(tabella_mvabb) <- tax_vector
# tabella_mvabb <- as.data.frame(t(tabella_mvabb))
# 
# var_mvabb_all <- sample_data(phyloseq_obj_css)
# var_mvabb<-var_mvabb_all[-12,]
# ## modelliamo OC*IdInp*Cu
# mvabun1 <- mvabund(tabella_mvabb)
# manyglm_1 <- manyglm(as.matrix(tabella_mvabb) ~ var_mvabb$OC *var_mvabb$IdInp*var_mvabb$Cu, 
#                      family="negative.binomial")
# plot(manyglm_1)
# 
# # SLOW!!!!!!!!
# summary.model_unadj <- summary(manyglm_1, 
#                                test="LR", 
#                                p.uni = "unadjusted", 
#                                show.cor = T)
# 
# # see https://rdrr.io/rforge/mvabund/man/summary.manylm.html
# summary.model_unadj$uni.p
# summary.model_unadj$uni.test
# test.dev.T1 <- as.data.frame(summary.model_unadj$uni.test)
# test.p.T1 <- as.data.frame(summary.model_unadj$uni.p)
# 
# ## output 
# out_mvabund_univar <- cbind(summary.model_unadj$uni.test, summary.model_unadj$uni.p, summary.model_unadj$uni.p)
# write.table(x = out_mvabund_univar, file = "./mvabund_unviaraite_res.csv")
# 
# ## SLOW!!!!!!!
# anova.model.noadjust <- anova(manyglm_1, p.uni = "unadjusted", nCores = 7)
# 
# # qui fa un test per ogni OTU, per questo è lento. Forse sarebbe da lavorare sui generi anzi che le OTU?


######## CORRELATIONS

## idea è di fare correlazioni spearman tra tutte le otu e tra le otu e le variabili selezionate,
## fare poi un filtro su R e p-value e poi trasformare la matrice in modo da ottenere le coppie

# uso come input phyloseq_obj_css come da linea 384

phylo_correl <- phyloseq_obj_css
colnames(tax_table(phylo_correl))
phylo_correl <- tax_glom(phylo_correl, "genus")
phylo_correl <- subset_samples(phylo_correl, station != "5" | site != "H")
sample_names(phylo_correl)
abb_table_corr <- t(otu_table(phylo_correl))

# serve per dopo, cosi so che posso sotituire otu_names alla tax table, dato che l'ordine è lo stesso
rownames(tax_table(phylo_correl)) == colnames(abb_table_corr)

rownames(abb_table_corr)
head(tax_table(phylo_correl)[,6])
#colnames(abb_table_corr) <- tax_table(phylo_correl)[,6]
# cambio tattica, faccio il buon vecchio otu1 - otuX e poi esporto la tassonomia e la importo come node attribute
length(colnames(abb_table_corr))
otu_names <- paste("OTU",seq(1:length(colnames(abb_table_corr))), sep = "")
colnames(abb_table_corr) <- otu_names

env_corr <-scale(poll3_Al[-12,c("OC","IdInp","Cu")])

corr_df <- cbind(abb_table_corr,env_corr)

rownames(corr_df)
colnames(corr_df)

# esporto quindi la tabella tassonomia

tax_correlation <- as.data.frame(tax_table(phylo_correl)@.Data)
rownames(tax_correlation) <- otu_names

library(psych)

correl_all <- corr.test(x = corr_df, use = "pairwise",method="spearman",
                         adjust="fdr", alpha=.05,ci=F,minlength=5)

r<- correl_all$r
p <- correl_all$p

library(reshape2)
melt_r <- melt(replace(correl_all$r, lower.tri(correl_all$r, TRUE), NA), na.rm = TRUE)
colnames(melt_r)[3] <- "r_Value"
melt_p <- melt(replace(correl_all$p, lower.tri(correl_all$p, TRUE), NA), na.rm = TRUE)
colnames(melt_p)[3] <- "p_Value"

head(melt_p)
head(melt_r)

melt_r$p_Value <- melt_p$p_Value
head(melt_r)

## ora filtriamo
# 1: seleziono solo le correlazioni significative
melt_sign <- melt_r[melt_r$p_Value < 0.05,] #1303 correlazioni

# 2: selezioniamo solo le correlazioni medio-forti
melt_sign_strong <- melt_sign[abs(melt_sign$r_Value) > 0.8,] #1303 correlazioni
# le significative son tutte strong in sostanza

# output tables
write.table(melt_sign_strong, "./correlations_spearman_ALL.csv", sep = ";")


# per i taxa, direi inutile portarsi dietro tutto
names <- unique(c(as.character(melt_sign_strong$Var1), as.character(melt_sign_strong$Var2))) 
# considera che una stessa OTU ha piu connesisoni, quindi compare piu di una volta
# in teoria in names ci sono tutti i nodi della network 

tax_correlation_subset <- subset(tax_correlation, rownames(tax_correlation) %in% names)

write.table(tax_correlation, "./node_attr_correlation.csv", sep = ";")


#### CORRELATION CON TUTTE VARIABILI
colnames(poll3_Al)
env_corr <-scale(poll3_Al[-12,-c(1,2,10:13,29,30,35:65)])

colnames(env_corr)

corr_df <- cbind(abb_table_corr,env_corr)

rownames(corr_df)
colnames(corr_df)

# esporto quindi la tabella tassonomia

tax_correlation <- as.data.frame(tax_table(phylo_correl)@.Data)
rownames(tax_correlation) <- otu_names

library(psych)

correl_all <- corr.test(x = corr_df, use = "pairwise",method="spearman",
                        adjust="fdr", alpha=.05,ci=F,minlength=5)

r<- correl_all$r
p <- correl_all$p

library(reshape2)
melt_r <- melt(replace(correl_all$r, lower.tri(correl_all$r, TRUE), NA), na.rm = TRUE)
colnames(melt_r)[3] <- "r_Value"
melt_p <- melt(replace(correl_all$p, lower.tri(correl_all$p, TRUE), NA), na.rm = TRUE)
colnames(melt_p)[3] <- "p_Value"

head(melt_p)
head(melt_r)

melt_r$p_Value <- melt_p$p_Value
head(melt_r)

## ora filtriamo
# 1: seleziono solo le correlazioni significative
melt_sign <- melt_r[melt_r$p_Value < 0.05,] #1303 correlazioni

# 2: selezioniamo solo le correlazioni medio-forti
melt_sign_strong <- melt_sign[abs(melt_sign$r_Value) > 0.8,] #1303 correlazioni
# le significative son tutte strong in sostanza

# output tables
write.table(melt_sign_strong, "./correlations_spearman_ALLVar.csv", sep = ";")


# per i taxa, direi inutile portarsi dietro tutto
names <- unique(c(as.character(melt_sign_strong$Var1), as.character(melt_sign_strong$Var2))) 
# considera che una stessa OTU ha piu connesisoni, quindi compare piu di una volta
# in teoria in names ci sono tutti i nodi della network 

tax_correlation_subset <- subset(tax_correlation, rownames(tax_correlation) %in% names)

write.table(tax_correlation, "./node_attr_correlation_ALLvar.csv", sep = ";")

## calcolo CV% per tutte le OTU

ixs <- as.data.frame(t(abb_table_corr)@.Data)
tax_correlation$SD <- apply(ixs,1,sd)

co.var <- function(x) {(sd(x)/mean(x))*100}

tax_correlation$CV <- apply(ixs,1,co.var)

write.table(tax_correlation, "./node_attr_correlation_ALLvar.csv", sep = ";")


#### visualmente sulla network non si vede nessuna distribuzione strana rispetto a CV o SD ######


#otu_cluster_6 <- read.delim("clipboard")

keep <- otu_cluster_6$shared.name[-c(72,74,76,78,79,80)]

keep <- as.character(keep)

pippo <- phylo_correl
otu_names <- paste("OTU",seq(1:length(colnames(abb_table_corr))), sep = "")

taxa_names(pippo) <- otu_names

pippo_subset <- subset(otu_table(pippo), rownames(otu_table(pippo)) %in% keep)

inp <- t(otu_table(pippo_subset))@.Data

inp[inp == 0] <- NA

library(pheatmap)
pheatmap(inp,
         gaps_row = c(5,8),
         color = brewer.pal(n = 9, "YlGn"),
         cluster_cols=T, 
         clustering_distance_cols = "euclidean",
         cluster_rows = F, 
         show_rownames = T,
         border_color = "grey27",
         na_col = "black",
         scale = "column")


out <- pheatmap(inp,
                gaps_row = c(5,8),
                color = brewer.pal(n = 9, "YlGn"),
                cluster_cols=T, 
                clustering_distance_cols = "euclidean",
                cluster_rows = F, 
                show_rownames = T,
                border_color = "grey27",
                na_col = "black",
                scale = "column")

colnames(inp[,out$tree_col[["order"]]])

cagliari_heatmap <- colnames(inp[,out$tree_col[["order"]]])[1:44]
Elk_her_heatmap <- colnames(inp[,out$tree_col[["order"]]])[45:77]

subset_cluster6_heatmap1 <- subset_taxa(pippo, taxa_names(pippo)%in%cagliari_heatmap)

tax_table(subset_cluster6_heatmap1)

library(pals)
paletta <- colorRampPalette(as.character(alphabet(26)))(30)
otu_table_collapsed <- merge_samples(subset_cluster6_heatmap1, group = "site")
otu_table_collapsed <- tax_glom(otu_table_collapsed, taxrank="family")
#qiime_file_proportional <- transform_sample_counts(otu_table_collapsed, function(x) 100 * x/sum(x))
#qiime_file_proportional_oneperc <-  filter_taxa(qiime_file_proportional, function(x) max(x) > 1, TRUE)
p_phylum <- plot_bar(otu_table_collapsed, fill= "family")
p_phylum <- p_phylum + theme_linedraw() + scale_fill_manual(values= paletta) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="right") + ylab("Abundance (%)") + xlab("")
p_phylum + guides(fill = guide_legend(nrow = 15)) 



subset_cluster6_heatmap2 <- subset_taxa(pippo, taxa_names(pippo)%in%Elk_her_heatmap)

tax_table(subset_cluster6_heatmap2)

library(pals)
paletta <- colorRampPalette(as.character(alphabet(26)))(30)
otu_table_collapsed <- merge_samples(subset_cluster6_heatmap2, group = "site")
otu_table_collapsed <- tax_glom(otu_table_collapsed, taxrank="family")
#qiime_file_proportional <- transform_sample_counts(otu_table_collapsed, function(x) 100 * x/sum(x))
#qiime_file_proportional_oneperc <-  filter_taxa(qiime_file_proportional, function(x) max(x) > 1, TRUE)
p_phylum <- plot_bar(otu_table_collapsed, fill= "family")
p_phylum <- p_phylum + theme_linedraw() + scale_fill_manual(values= paletta) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="right") + ylab("Abundance (%)") + xlab("")
p_phylum + guides(fill = guide_legend(nrow = 15)) 


## faccio anche su cluster OC

#otu_cluster_2 <- read.delim("clipboard")

keep <- otu_cluster_2$shared.name[-c(29)]

keep <- as.character(keep)

pippo <- phylo_correl
otu_names <- paste("OTU",seq(1:length(colnames(abb_table_corr))), sep = "")

taxa_names(pippo) <- otu_names

pippo_subset <- subset(otu_table(pippo), rownames(otu_table(pippo)) %in% keep)

inp <- t(otu_table(pippo_subset))@.Data

inp[inp == 0] <- NA

library(pheatmap)
pheatmap(inp,
         gaps_row = c(5,8),
         color = brewer.pal(n = 9, "YlGn"),
         cluster_cols=T, 
         clustering_distance_cols = "euclidean",
         cluster_rows = F, 
         show_rownames = T,
         border_color = "grey27",
         na_col = "black",
         scale = "column")


subset_cluster2 <- subset_taxa(pippo, taxa_names(pippo)%in%keep)

plot_heatmap(subset_cluster2) # torna forse un pelo meglio, ma la visualizzazione pheatmap è piu bella

sample_data(subset_cluster2)
tax_table(subset_cluster2)
plot_bar(subset_cluster2, fill = "family")

library(pals)
paletta <- colorRampPalette(as.character(alphabet(26)))(30)
otu_table_collapsed <- merge_samples(subset_cluster2, group = "site")
otu_table_collapsed <- tax_glom(otu_table_collapsed, taxrank="phylum")
qiime_file_proportional <- transform_sample_counts(otu_table_collapsed, function(x) 100 * x/sum(x))
qiime_file_proportional_oneperc <-  filter_taxa(qiime_file_proportional, function(x) max(x) > 1, TRUE)
p_phylum <- plot_bar(qiime_file_proportional_oneperc, fill= "phylum")
p_phylum <- p_phylum + theme_linedraw() + scale_fill_manual(values= paletta) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="right") + ylab("Abundance (%)") + xlab("")
p_phylum + guides(fill = guide_legend(nrow = 15)) 



## faccio anche su cluster 1 solo OTU

otu_cluster_1 <- read.delim("clipboard")

keep <- otu_cluster_1$shared.name

keep <- as.character(keep)

pippo <- phylo_correl
otu_names <- paste("OTU",seq(1:length(colnames(abb_table_corr))), sep = "")

taxa_names(pippo) <- otu_names

pippo_subset <- subset(otu_table(pippo), rownames(otu_table(pippo)) %in% keep)

inp <- t(otu_table(pippo_subset))@.Data

inp[inp == 0] <- NA

library(pheatmap)
pheatmap(inp,
         gaps_row = c(5,8),
         color = brewer.pal(n = 9, "YlGn"),
         cluster_cols=T, 
         clustering_distance_cols = "euclidean",
         cluster_rows = F, 
         show_rownames = T,
         border_color = "grey27",
         na_col = "black",
         scale = "column")


subset_cluster1 <- subset_taxa(pippo, taxa_names(pippo)%in%keep)

plot_heatmap(subset_cluster1) # torna forse un pelo meglio, ma la visualizzazione pheatmap è piu bella

sample_data(subset_cluster1)
tax_table(subset_cluster1)
#plot_bar(subset_cluster1, fill = "family")

library(pals)
paletta <- colorRampPalette(as.character(alphabet(26)))(30)
otu_table_collapsed <- merge_samples(subset_cluster1, group = "site")
otu_table_collapsed <- tax_glom(otu_table_collapsed, taxrank="family")
qiime_file_proportional <- transform_sample_counts(otu_table_collapsed, function(x) 100 * x/sum(x))
qiime_file_proportional_oneperc <-  filter_taxa(qiime_file_proportional, function(x) max(x) > 1, TRUE)
p_phylum <- plot_bar(qiime_file_proportional_oneperc, fill= "family")
p_phylum <- p_phylum + theme_linedraw() + scale_fill_manual(values= paletta) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="right") + ylab("Abundance (%)") + xlab("")
p_phylum + guides(fill = guide_legend(nrow = 15)) 


## CORRELATION SOLO SRB

phylo_correl_SRB <- SRB_physeq_silva
colnames(tax_table(phylo_correl_SRB))
#phylo_correl_SRB <- tax_glom(phylo_correl_SRB, "genus")
phylo_correl_SRB <- subset_samples(phylo_correl_SRB, station != "5" | site != "H")
sample_names(phylo_correl_SRB)
abb_table_corr_SRB <- t(otu_table(phylo_correl_SRB))

# serve per dopo, cosi so che posso sotituire otu_names alla tax table, dato che l'ordine è lo stesso
rownames(tax_table(phylo_correl_SRB)) == colnames(abb_table_corr_SRB)

rownames(abb_table_corr)
head(tax_table(phylo_correl_SRB)[,6])
#colnames(abb_table_corr) <- tax_table(phylo_correl)[,6]
# cambio tattica, faccio il buon vecchio otu1 - otuX e poi esporto la tassonomia e la importo come node attribute
length(colnames(abb_table_corr_SRB))
otu_names <- paste("OTU",seq(1:length(colnames(abb_table_corr_SRB))), sep = "")
colnames(abb_table_corr_SRB) <- otu_names

env_corr <-scale(poll3_Al[-12,c("OC","IdInp","Cu")])

corr_df_SRB <- cbind(abb_table_corr_SRB,env_corr)

rownames(corr_df_SRB)
colnames(corr_df_SRB)

# esporto quindi la tabella tassonomia

tax_correlation_SRB <- as.data.frame(tax_table(phylo_correl_SRB)@.Data)
rownames(tax_correlation_SRB) <- otu_names

library(psych)

correl_SRB <- corr.test(x = corr_df_SRB, use = "pairwise",method="spearman",
                        adjust="fdr", alpha=.05,ci=F,minlength=5)

r<- correl_SRB$r
p <- correl_SRB$p

library(reshape2)
melt_r <- melt(replace(correl_SRB$r, lower.tri(correl_all$r, TRUE), NA), na.rm = TRUE)
colnames(melt_r)[3] <- "r_Value"
melt_p <- melt(replace(correl_SRB$p, lower.tri(correl_all$p, TRUE), NA), na.rm = TRUE)
colnames(melt_p)[3] <- "p_Value"

head(melt_p)
head(melt_r)

melt_r$p_Value <- melt_p$p_Value
head(melt_r)

## ora filtriamo
# 1: seleziono solo le correlazioni significative
melt_sign <- melt_r[melt_r$p_Value < 0.05,] #1303 correlazioni

# 2: selezioniamo solo le correlazioni medio-forti
melt_sign_strong <- melt_sign[abs(melt_sign$r_Value) > 0.8,] #1303 correlazioni
# le significative son tutte strong in sostanza

# output tables
write.table(melt_sign_strong, "./correlations_spearman_SRB.csv", sep = ";")


# per i taxa, direi inutile portarsi dietro tutto
names <- unique(c(as.character(melt_sign_strong$Var1), as.character(melt_sign_strong$Var2))) 
# considera che una stessa OTU ha piu connesisoni, quindi compare piu di una volta
# in teoria in names ci sono tutti i nodi della network 

tax_correlation_subset <- subset(tax_correlation_SRB, rownames(tax_correlation_SRB) %in% names)

write.table(tax_correlation_subset, "./node_attr_correlation_SRB.csv", sep = ";")



#####################################################
#### provare anche correlationTest in metagenomeseq  --> non mi piace perchè il metodo non aggiusta, e viene fuori una palla di pelo
# 
# correlations_metagenomeseq <- correlationTest(obj = t(corr_df), method = "spearman", cores = 7)
# 
# # elaboro per essere input di correlazione
# 
# coppie <- row.names(correlations_metagenomeseq)
# 
# metagenomeseq_corr_df<- data.frame(t(sapply(strsplit(coppie,split = "-"), `[`)))
# colnames(metagenomeseq_corr_df) <- c("SOURCE","TARGET")
# metagenomeseq_corr_df$INTERACTION <- rep("pp", nrow(metagenomeseq_corr_df))
# metagenomeseq_corr_df$R <- correlations_metagenomeseq[,1]
# metagenomeseq_corr_df$p <- correlations_metagenomeseq[,2]
# 
# # filtro
# dim(metagenomeseq_corr_df)
# metagenomeseq_corr_df_sign <- metagenomeseq_corr_df[metagenomeseq_corr_df$p < 0.05,] #1303 correlazioni
# dim(metagenomeseq_corr_df_sign)
# 
# metagenomeseq_corr_df_sign_strong <- metagenomeseq_corr_df_sign[abs(metagenomeseq_corr_df_sign$R) > 0.5,] #1303 correlazioni
# dim(metagenomeseq_corr_df_sign_strong)
# 
# # output tables
# write.table(metagenomeseq_corr_df_sign_strong, "./correlations_spearman_ALL_metagenomeseq.csv", sep = ";")
# 
# 
# # per i taxa, direi inutile portarsi dietro tutto
# names <- unique(c(as.character(metagenomeseq_corr_df_sign_strong$SOURCE), as.character(metagenomeseq_corr_df_sign_strong$TARGET))) 
# # considera che una stessa OTU ha piu connesisoni, quindi compare piu di una volta
# # in teoria in names ci sono tutti i nodi della network 
# 
# tax_correlation_subset2 <- subset(tax_correlation, rownames(tax_correlation) %in% names)
# 
# write.table(tax_correlation_subset2, "./node_attr_correlation.csv", sep = ";")

