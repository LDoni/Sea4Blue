### analisi atlatintico 16S
library(phyloseq)
library(ggside)
library(vegan)
library(ggplot2)
library(readr)
library(plyr)
library(dplyr)
library(hrbrthemes)
library(ggpubr)
library(viridis)
library(ggforce)
library(concaveman)
library(ggside)
library(ggdist)
library("foreach")
library("doParallel")
library(metagenomeSeq)
library(dplyr)
library(microbiome)
library(vegan)
library(scales)
library(microbiomeMarker)

#https://github.com/joey711/phyloseq/issues/1163

#load data PS
#ps <- readRDS("ps.rds")


#phyloseq_ill<-ps



phyloseq_ill <- phyloseq(otu_table(ps), tax_table(ps), sample_data(ps))




ps_nanopore <- readRDS("physeq_nanopore_spaghetti.rds")
head(otu_table(ps_nanopore))

taxa_names(ps_nanopore) <- paste0("ps_nan_", taxa_names(ps_nanopore))

taxa_names(phyloseq_ill) <- paste0("ps_ill_", taxa_names(phyloseq_ill))


phyloseq_obj_ALL<-merge_phyloseq(ps_nanopore,phyloseq_ill)


phyloseq_obj_ALL<- tax_glom(phyloseq_obj_ALL, taxrank = 'Genus')





### l'ho modificato 
#write.csv(sample_data(phyloseq_obj_ALL),"metadata_sea4blue.csv")
sample_data(phyloseq_obj_ALL)<-sample_data(read.csv("metadata_sea4blue.csv",row.names=1,check.names=FALSE,sep = ";"))



#saveRDS(phyloseq_obj_ALL, "phyloseq_obj_ALL.rds")



phyloseq_obj_ALL <- readRDS("phyloseq_obj_ALL.rds")









###                      ANALISI COMUNITA MICROBICA                           ####


####                                 normalizzazione                 ###############


### RIMOZIONE CHL E MITOK

phyloseq_obj<-phyloseq_obj_ALL
grep(pattern = "Mitochondria", tax_table(phyloseq_obj)) 
grep(pattern = "Chloroplast", tax_table(phyloseq_obj)) 
phyloseq_obj <- phyloseq_obj %>% subset_taxa( Family!= "Mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class) ) 

phyloseq_obj <- subset_taxa(phyloseq_obj, (tax_table(phyloseq_obj)[,"Order"]!="Chloroplast") | is.na(tax_table(phyloseq_obj)[,"Order"]))

phyloseq_obj

####                                 normalizzazione                 ###############

#filtering singletons
doubleton <- genefilter_sample(phyloseq_obj, filterfun_sample(function(x) x > 1), A=1)
doubleton <- prune_taxa(doubleton, phyloseq_obj) 



## rimuovo i pos e neg

sample_variables(doubleton)
(sample_data(doubleton)$samplenames)

#rimuovo i positivi e negativi
PosNeg = subset_samples(doubleton, samplenames == "Positive" | samplenames=="Negative")
head(sample_data(PosNeg))
Pos = subset_samples(PosNeg, samplenames == "Lapospositive" )

#rimuovo i bluehole
doubleton = subset_samples(doubleton,!( samplenames=="Positive" | samplenames=="Negative"))
tail(sample_data(doubleton))

sample_names(doubleton)


# transforming 
data.metagenomeSeq = phyloseq_to_metagenomeSeq(doubleton)

p = cumNormStat(data.metagenomeSeq) #default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
#data.cumnorm
data.CSS = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
head(data.CSS)
dim(data.CSS)  # make sure the data are in a correct formal: number of samples in rows
phyloseq_obj_css_comparison <- phyloseq_obj
otu_table(phyloseq_obj_css) <- otu_table(data.CSS, taxa_are_rows = T)

phyloseq_obj_css

#cambio i nomi
#sample_names(phyloseq_obj_css_comparison)<-sample_data(phyloseq_obj_css_comparison)$Nominativo.campione.






phyloseq_obj_css_comparison = subset_samples(phyloseq_obj_css, Giorno.di.navigazione. == "1" |
                                    Giorno.di.navigazione.=="11"|
                                    Giorno.di.navigazione.=="17"|
                                    Giorno.di.navigazione.=="25"|
                                    Giorno.di.navigazione.=="26")


sample_data(phyloseq_obj_css_comparison)

#plot heatmap  
p<-plot_heatmap(tax_glom(phyloseq_obj_css_comparison, taxrank = "Genus"),taxa.label = "Genus",first.taxa = "Vibrio",sample.order = c("eDNA1","eDNA3","eDNA4","eDNA6","eDNA7","eDNA9","eDNA15","eDNA16","eDNA17","eDNA18","eDNA19","eDNA25","eDNA26","eDNA27","eDNA28"))

p+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+theme_bw(base_size = 9,"black")
p<-p+theme(axis.text.y = element_text(color = "black", size = 6, angle = 0, hjust = 1, vjust = 0, face = "plain"),axis.text.x = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"))
p




#rifare physeq con numero otu assolute e non scalate CSS in modo che siano interi
# e che vadano bene per il calcolo indici alpha div dopo

otu.absol<-round(otu_table(phyloseq_obj_css_comparison))
head(otu.absol)
physeq_normalized <- phyloseq_obj_css_comparison
otu_table(physeq_normalized) <- otu_table(as.matrix(otu.absol), taxa_are_rows = T)


#                                     Alpha diversity                                  ####
head(otu_table(physeq_normalized))

sample_variables(physeq_normalized)

plot_richness(physeq_normalized,x="Longitude" ,measures=c("Observed", "Shannon", "Simpson"),color = "Zona")+
  geom_boxplot(aes(alpha = 1/10),show.legend = FALSE)+theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(values = c("#C37F7B","#814e07", "#9F914B" ,"#58A069" , "#7B90C4","#00A1A4"))



plot_richness(physeq_normalized,x="sequencing" ,measures=c("Observed", "Shannon", "Simpson"),color = "DNA")+
  geom_boxplot(aes(alpha = 1/10),show.legend = FALSE)+theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(values = c("#C37F7B","#814e07", "#9F914B" ,"#58A069" , "#7B90C4","#00A1A4"))




################################  bellissime
###plot GGPLOT con geom ponit alpha 
p1<-plot_richness(physeq_normalized,x="Latitude",color = "Zona",measures=c("Observed"))

#Longitudine 
newSTorder = c( "ANW","ANC","ANE")
p1$data$Zona<- as.character(p1$data$Zona)
p1$data$Zona <- factor(p1$data$Zona, levels=newSTorder) 


min(p1$data$Longitude)
max(p1$data$Longitude)
#facet_grid(~Depth)+
ggplot(p1$data,aes(Longitude,value))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour = Zona),alpha =0.45)+
  scale_color_manual(values=c("#4280fc", "#ffb452","#f7170a","#000000"))+
  geom_ysideboxplot(aes(x=Zona,y=value,colour = Zona), orientation = "x") +
  theme(        ggside.panel.scale.y = .4)+scale_ysidex_discrete()+
  geom_smooth(aes(colour = variable ),size=1.5)+ ylab('Richness') 
# scale_x_continuous(breaks = seq(-65, 80, by = 20))#, expand = c(0, 0)    theme(axis.text.x = element_text(angle = 45, hjust = 1))+


ggboxplot(p1$data, x = "Zona", y = "value",
          color = "Zona", palette = "jco",
          add = "jitter")+ stat_compare_means()  #### Kruskal-Wallis p=0.08












scale_color_hue(l=60, c=48)
scale_color_brewer(palette='Spectral')


show_col(hue_pal(l=70, c=90,direction = -1)(50))





show_col(hue_pal(l=60, c=48)(6))
pseq_tutt<-physeq_normalized










#####                      BETA DIVERSITY                           ####

#PCoA on Bray-Curtis Dissimilarity
phyloseq_obj_css_comparison
metad

library(ggpubr)
library(ggrepel)
head(sample_data(phyloseq_obj_css_comparison))


### bray-curtis

otu.ord <- ordinate(physeq = phyloseq_obj_css_comparison, "PCoA", distance = "bray")


#plot ordination

A<-plot_ordination(phyloseq_obj_css_comparison, otu.ord,  color="Zona",axes =c(1,2))+
  scale_color_manual(values = c("#4280fc", "#ffb452","#f7170a"))+ geom_point(size=3)+ 
  geom_text_repel(aes(label=samples ),max.overlaps = Inf, show.legend = FALSE)+
 theme_bw() + 
  geom_xsidedensity(aes(y=stat(density),fill=Zona), alpha = 0.5, show.legend = FALSE) +
  geom_ysidedensity(aes(x=stat(density),fill=Zona), alpha = 0.5, show.legend = FALSE) +
  scale_xsidey_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
  scale_ysidex_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +scale_ysidex_discrete()+
  ggside::theme_ggside_void()  +
  scale_fill_manual(values = c("#4280fc", "#ffb452","#f7170a"))+ggtitle("bray-Curtis")
A

# geom_mark_ellipse(aes(color = Zona), show.legend = FALSE)+ theme_void()+


B<-plot_ordination(phyloseq_obj_css_comparison, otu.ord,  color="Zona",axes =c(1,3))+
  scale_color_manual(values = c("#4280fc", "#ffb452","#f7170a"))+ geom_point(size=3)+ 
  geom_text_repel(aes(label=samples),max.overlaps = Inf, show.legend = FALSE)+
  geom_mark_ellipse(aes(color = Zona), show.legend = FALSE)+ theme_void()+theme_bw() + 
  geom_xsidedensity(aes(y=stat(density),fill=Zona), alpha = 0.5, show.legend = FALSE) +
  geom_ysidedensity(aes(x=stat(density),fill=Zona), alpha = 0.5, show.legend = FALSE) +
  scale_xsidey_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
  scale_ysidex_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +scale_ysidex_discrete()+
  ggside::theme_ggside_void()  +
  scale_fill_manual(values = c("#4280fc", "#ffb452","#f7170a"))+ggtitle("bray-Curtis")
B
library(patchwork)

A+B
# colori "#C37F7B","#814e07", "#9F914B" ,"#58A069" , "#7B90C4","#00A1A4"


##### permanova

set.seed(1)

# Calculate bray curtis distance matrix
bray <- phyloseq::distance(phyloseq_obj_css_comparison, method = "bray")
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(phyloseq_obj_css_comparison))

# Adonis test
adonis2(bray ~ Zona, data = sampledf) #0.001 ***

adonis2(bray ~ Zona*Giorno.di.navigazione., data = sampledf)

#Zona                         2  1.17766 0.50695 8.4339  0.001 ***
#Giorno.di.navigazione.       1  0.18168 0.07821 2.6022  0.038 *  
#  Zona:Giorno.di.navigazione.  2  0.33533 0.14435 2.4015  0.006 ** 

adonis2(bray ~ Zona+Giorno.di.navigazione.+Latitude+Longitude, data = sampledf)

#Zona                    2  1.17766 0.50695 8.7154  0.001 ***
# Giorno.di.navigazione.  1  0.18168 0.07821 2.6891  0.031 *  
# Latitude                1  0.14997 0.06456 2.2198  0.058 .  
#Longitude               1  0.20565 0.08853 3.0438  0.021 * 



################# unifracchete
wunifrac

otu.ord <- ordinate(physeq = phyloseq_obj_css, "PCoA", distance = "wunifrac")


#plot ordination

plot_ordination(phyloseq_obj_css, otu.ord,  color="Zona",axes =c(1,2))+
  scale_color_manual(values = c("#4280fc", "#ffb452","#f7170a"))+ geom_point(size=3)+ 
  geom_text_repel(aes(label=Giorno.di.navigazione.),max.overlaps = Inf, show.legend = FALSE)+
  geom_mark_ellipse(aes(color = Zona), show.legend = FALSE)+ theme_void()+theme_bw() + 
  geom_xsidedensity(aes(y=stat(density),fill=Zona), alpha = 0.5, show.legend = FALSE) +
  geom_ysidedensity(aes(x=stat(density),fill=Zona), alpha = 0.5, show.legend = FALSE) +
  scale_xsidey_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
  scale_ysidex_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +scale_ysidex_discrete()+
  ggside::theme_ggside_void()  +
  scale_fill_manual(values = c("#4280fc", "#ffb452","#f7170a"))+ggtitle("wunifrac")






#### barplot###



genus.sum = tapply(taxa_sums(phyloseq_obj_css_comparison), tax_table(phyloseq_obj_css_comparison)[, "Genus"], sum, na.rm=TRUE)
top5phyla = names(sort(genus.sum, TRUE))[1:30]
GP1 = prune_taxa((tax_table(phyloseq_obj_css_comparison)[, "Genus"] %in% top5phyla), phyloseq_obj_css_comparison)


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



p<-my_plot_bar(GP1, "Genus", fill = "Genus", facet_grid=~Zona)+
  scale_fill_manual(values =c( "#7FFF00","#94755c","#fbff00","#7B1B02",
                               "#aedce8","#00A86B","#ff9d4d","#4B0082","#293133","#7FFFD4",
                               "#708090","#1025e6","#fc0317","#FFBF00","#ff0303","#4B0082"))


newSTorder = c( "ANW","ANC","ANE")
p$data$Zona<- as.character(p$data$Zona)
p$data$Zona <- factor(p$data$Zona, levels=newSTorder)

p
a<-p+theme_bw()+ theme(axis.text.x = element_text(angle = -90, hjust = 0,colour = "black"))
a<-p+theme_bw()+theme(axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())
a+theme(axis.text=element_text(size=12,colour = "black"),
        axis.title=element_text(size=14,face="bold"))




##### lefse

?run_lefse
mm_lefse <- run_lefse(
  phyloseq_obj_css_comparison,
  taxa_rank="Genus",
  wilcoxon_cutoff = 0.05,
  group = "Zona",
  kw_cutoff = 0.05,
  multigrp_strat = T,
  lda_cutoff = 4)
rm(mm_lefse)
length(rownames(marker_table(mm_lefse)))


mm_lefse <- run_lefse(
  phyloseq_obj_css_comparison,
  wilcoxon_cutoff = 0.05,
  group = "Zona",
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 4,taxa_rank="Genus"
)




plot_ef_bar(mm_lefse)
plot_ef_dot(mm_lefse)

plot_cladogram(mm_lefse, color = c(ANC = "blue", ANW = "red",ANE="orange"),only_marker = T) +
  theme(plot.margin = margin(0, 0, 0, 0))

plot_cladogram(mm_lefse, color = c(ANC = "blue", ANW = "red",ANE="orange")) +
  theme(plot.margin = margin(0, 0, 0, 0))
