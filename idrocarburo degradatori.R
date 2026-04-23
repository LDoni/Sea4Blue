### idrocarburo degradatori

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



#load data PS
#ps <- readRDS("ps.rds")




phyloseq_obj<-ps
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



doubleton = subset_samples(doubleton,!( Nominativo.campione.=="Positive" | Nominativo.campione.=="Negative"))
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
phyloseq_obj_css <- phyloseq_obj
otu_table(phyloseq_obj_css) <- otu_table(data.CSS, taxa_are_rows = T)

phyloseq_obj_css


#taxa_da_estrarre :::: Alcanivorax Oleibacter Methylophaga Cycloclasticus Marinobacter


IDROCARB_phyloseq = subset_taxa(phyloseq_obj_css, Genus=="Alcanivorax" | Genus=="Oleibacter"| Genus=="Methylophaga"| Genus=="Cycloclasticus"| Genus=="Marinobacter" )

#refseq(IDROCARB_phyloseq)



## salvare fasta
IDROCARB_phyloseq %>%
  refseq() %>%
  Biostrings::writeXStringSet("IDROCARB_phyloseq.fasta", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")




sample_data(IDROCARB_phyloseq)

IDROCARB_phyloseq


## per alfa div
otu.absol<-round(otu_table(IDROCARB_phyloseq))
head(otu.absol)
physeq_normalized <- IDROCARB_phyloseq
otu_table(physeq_normalized) <- otu_table(as.matrix(otu.absol), taxa_are_rows = T)


################################  bellissime
###plot GGPLOT con geom ponit alpha 
p1<-plot_richness(physeq_normalized,x="Longitude",color = "Zona",measures=c("Observed"))
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










############## barplot su mappa










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



IDROCARB_phyloseq

IDROCARB_phyloseq_bar = subset_samples(IDROCARB_phyloseq,!( Nominativo.campione.=="eDNA19" | Nominativo.campione.=="eDNA26"))


p<-my_plot_bar(IDROCARB_phyloseq_bar, "Genus", fill = "Genus")+facet_grid(~ Giorno.di.navigazione., scales = "free")
p+ggtitle("giorni di navigazione")


p$data <- p$data %>%
  filter(Giorno.di.navigazione. != "19")

 


ggplot(p$data, aes(x = Giorno.di.navigazione., y = Abundance, group = Genus, color = Genus)) +
  geom_smooth(method = "loess", se = FALSE) #geom_line() 

  

  

head(p$data)

newSTorder = c( "ANW","ANC","ANE")
p$data$Zona<- as.character(p$data$Zona)
p$data$Zona <- factor(p$data$Zona, levels=newSTorder)
p




