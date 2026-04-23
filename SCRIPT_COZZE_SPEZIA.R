
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

abund_table<-read.csv("feature-table.tsv.csv",row.names=1,check.names=FALSE)
#  abund_table<-t(abund_table)#Transpose the data to have sample names on rows
#  head(abund_table)
meta_table<-read.csv("metadata.csv",row.names=1, check.names=FALSE)
# manca questo file nella cartella dropbox
# OTU_tree <- read.tree("tree.nwk")
OTU_taxonomy<-read.csv("taxonomy.tsv.csv",row.names=1,check.names=FALSE)
# sort(colnames(abund_table) )==sort(rownames(OTU_taxonomy) )
#Convert the data to phyloseq format
OTU = otu_table(as.matrix(abund_table), taxa_are_rows = T)
TAX = tax_table(as.matrix(OTU_taxonomy))
SAM = sample_data(meta_table)
physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM)

###                      ANALISI COMUNITA MICROBICA                           ####







####                                 normalizzazione                 ###############
library(metagenomeSeq)
phyloseq_obj<-physeq
tax_table(physeq)
#filtering singletons
doubleton <- genefilter_sample(phyloseq_obj, filterfun_sample(function(x) x > 1), A=1)
doubleton <- prune_taxa(doubleton, phyloseq_obj) 

## Controllo mitocondrio e cloroplasto


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

#plot heatmap 
 p<-plot_heatmap(tax_glom(phyloseq_obj_css, taxrank = "Genus"),taxa.label = "Genus",sample.order = c("FEB","MAY","JUL","SEP","NOV"),first.taxa = "Vibrio")

p+theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold"))+theme_bw(base_size = 9,"black")
p<-p+theme(axis.text.y = element_text(color = "black", size = 6, angle = 0, hjust = 1, vjust = 0, face = "plain"),axis.text.x = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"))




tiff("heatmap.tiff", units="in", width=5, height=5, res=700)
p
dev.off()



#salvare l'otu table 
# Extract abundance matrix from the phyloseq object
OTU_matrix = as(otu_table(tax_glom(phyloseq_obj_css, taxrank = "Genus")), "matrix")
head(OTU_matrix)
# transpose if necessary
#if(taxa_are_rows(OTU_matrix)){OTU_matrix <- t(OTU_matrix)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU_matrix)
head(OTUdf)
write.csv(OTUdf,"OTU_table_bacteria_normalized_GENUS.csv")
OTU_taxonomy
head(OTU_taxonomy)


abund_table_normalized<-read.csv("OTU_table_bacteria_normalized_GENUS.csv",check.names=FALSE)
head(abund_table_normalized)
mergedTAX_OTU<-merge(abund_table_normalized,OTU_taxonomy,by = "OTU")
head(mergedTAX_OTU)

write.csv(as.data.frame(mergedTAX_OTU),"TAX_OTU_table_bacteria_normalized_GENUS.csv")

#export dist-matrix
B_C_COMUNITA_BATTERICA<-vegdist(t(OTUdf),method = "bray")
#write.csv(as.matrix(B_C_COMUNIT?_BATTERICA), "dist_matrix_BC.csv")



#rifare physeq con numero otu assolute e non scalate CSS in modo che siano interi
# e che vadano bene per il calcolo indici alpha div dopo

otu.absol<-round(otu_table(phyloseq_obj_css))
head(otu.absol)
OTU_absol = otu_table(as.matrix(otu.absol), taxa_are_rows = T)
physeq_normalized<-merge_phyloseq(phyloseq(OTU_absol, TAX),SAM,OTU_tree)


#questo lo usavo per bioenv
#df_otu table rounded 
batt1<-as.data.frame(otu_table(physeq_normalized))
head(batt1)

#                                     Alpha diversity                                  ####
head(otu_table(physeq_normalized))

sample_variables(physeq_normalized)
p<-plot_richness(physeq_normalized,x="Mese" ,measures=c("Observed", "Shannon", "Simpson"),color = "Mese")
library(microbiome)



pseq_tutt<-physeq_normalized


#meta_tutt<-meta(SRB_physeq_silva_rounded)
#otu_tutt<-abundances(SRB_physeq_silva_rounded)



tab <- alpha(pseq_tutt, index = c("observed","diversity_shannon","evenness_pielou"))

colnames(tab)[colnames(tab)=="observed"] <- "Richness"
colnames(tab)[colnames(tab)=="diversity_shannon"] <- "Shannon Diversity"
colnames(tab)[colnames(tab)=="evenness_pielou"] <- "Pielou Evenness"


tab$Mese<-c("FEB","JUL","MAY","NOV","SEP")

newSTorder = c( "FEB","MAY","JUL" ,"SEP", "NOV")

 

#per fare unico grafico con tutto
df.c<-melt(tab,  "Mese")

df.c$Mese<- as.character(df.c$Mese)
df.c$Mese <- factor(df.c$Mese, levels=newSTorder) 



library(ggrepel)

#boxplot porti
p<-ggplot(df.c, aes(Mese,value ,color= Mese )) +
  geom_point() +
  theme_bw()+ 
  facet_wrap(~variable, scales = "free")+ylab("Value")+xlab("Month")


getwd()
p







#####                      BETA DIVERSITY                           ####

#PCoA on Bray-Curtis Dissimilarity
phyloseq_obj_css

library(ggpubr)

otu.ord <- ordinate(physeq = phyloseq_obj_css, "PCoA", distance = "bray")

a<-plot_ordination(physeq = phyloseq_obj_css, otu.ord,title ="PCoA on Bray-Curtis Dissimilarity"  ,axes =c(1,2))+
  geom_point(aes(fill = Temperatura),size = 4, pch = 21)+ 
  geom_text(mapping = aes(label = Mese), size = 4, vjust = -0.5,hjust=0.5)+
  theme_bw() 
a+ scale_fill_viridis()
library(vegan)

OTU_matrix = as(otu_table(phyloseq_obj_css), "matrix")
OTUdf = as.data.frame(OTU_matrix)
braycurtis_bact<-vegdist(t(OTUdf),method = "bray")
env_data<-as.data.frame(sample_data(phyloseq_obj_css))


#permanova:
#Number of permutations: 999  

set.seed(123)
adonis2(braycurtis_bact~env_data$Temperatura, permutations=9999, method="bray") #p: 0.01667 *









##### test difference between ports!!!!!!


set.seed(123)
# adonis2(braycurtis_bact~env_data$Mese, permutations=9999, method="bray") # p:0.001














##plot taxa most abund
phyloseq_obj_css_barplot<-phyloseq_obj_css
genus.sum = tapply(taxa_sums(phyloseq_obj_css_barplot), tax_table(phyloseq_obj_css_barplot)[, "Genus"], sum, na.rm=TRUE)
top5phyla = names(sort(genus.sum, TRUE))[1:10]
GP1 = prune_taxa((tax_table(phyloseq_obj_css_barplot)[, "Genus"] %in% top5phyla), phyloseq_obj_css_barplot)

GP2 <- subset_taxa(GP1, Genus != "Unassigned")
GP2 <- subset_taxa(GP2, Genus != "uncultured")
GP2 <- subset_taxa(GP2, Genus != "")


GP.ord <- ordinate(GP2, "PCoA", "bray")
plot_ordination(GP2, GP.ord, type="taxa", color="Genus", title="taxa")+ facet_wrap(~Genus, 3)
p3 = plot_ordination(GP2, GP.ord, type="biplot", color="Genus")
p3 #  + facet_wrap(~Genus, 3)
p3+ geom_text(aes(label=Mese), size = 3, vjust = 0,hjust=0, color="#000000")



otu_table(GP2)

tax_table(GP2)
p<-plot_bar(GP2, "Genus", fill="Genus", facet_grid=~Mese) #Tow 

newSTorder = c( "FEB","MAY","JUL", "SEP", "NOV")

p$data$Mese<- as.character(p$data$Mese)
p$data$Mese <- factor(p$data$Mese, levels=newSTorder)  

p



# write.csv(as.data.frame(p$data), "data_barplot_mostabundant.csv")






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



p<-my_plot_bar(GP2, "Genus", fill = "Genus", facet_grid=~Mese)+
  scale_fill_manual(values =c( "#7FFF00","#94755c","#fbff00","#7B1B02",
                               "#aedce8","#00A86B","#ff9d4d","#4B0082","#293133","#7FFFD4",
                               "#708090","#1025e6","#fc0317","#FFBF00","#ff0303","#4B0082"))

a<-p+theme_bw()+ theme(axis.text.x = element_text(angle = -90, hjust = 0,colour = "black"))
a<-p+theme_bw()+theme(axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())
a+theme(axis.text=element_text(size=12,colour = "black"),
        axis.title=element_text(size=14,face="bold"))

newSTorder = c( "FEB","MAY","JUL", "SEP", "NOV")

p$data$Mese<- as.character(p$data$Mese)
p$data$Mese <- factor(p$data$Mese, levels=newSTorder)  

p










##### core microbiome


phyloseq_obj_css

# Calculate compositional version of the data
# (relative abundances)

x1 <- tax_glom(phyloseq_obj_css, taxrank="Genus")
pseq <- x1
pseq.rel <- microbiome::transform(pseq, "compositional")

# core.taxa.standard <- core_members(pseq.rel, detection = 0.0001, prevalence = 50/100)

# core.taxa.standard

# Use the microbiome function add_besthit to get taxonomic identities of ASVs.
ps.m3.rel.f <- microbiome::add_besthit(pseq.rel)

# Check 
taxa_names(ps.m3.rel.f)[1:10]


#tabella con i core taxa 

pseq.core <- core(ps.m3.rel.f, detection = 0.0001, prevalence = 0.9)
core.taxa <- taxa(pseq.core)
class(core.taxa)
# get the taxonomy data
tax.mat <- tax_table(pseq.core)
tax.df <- as.data.frame(tax.mat)

# add the OTus to last column
tax.df$OTU <- rownames(tax.df)

# select taxonomy of only 
# those OTUs that are core memebers based on the thresholds that were used.
core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
knitr::kable((core.taxa.class))


nrow(core.taxa.class)




### plot core taxa
# Core with compositionals:
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-2), log10(.2), length = 10), 3)

#Deletes "ASV" from taxa_names, e.g. ASV1 --> 1
#taxa_names(ps.m3.rel) = taxa_names(ps.m3.rel) %>% str_replace("ASV", "")
# Also define gray color palette







gray <- gray(seq(0,1,length=5))

p1 <- plot_core(ps.m3.rel.f,
                plot.type = "heatmap",
                colours = gray,
                prevalences = prevalences,
                detections = detections, min.prevalence = .8) +
  xlab("Detection Threshold (Relative Abundance (%))")

p1 <- p1 + theme_bw() + ylab("ASVs")
p1

library(viridis)
print()


p1 + scale_fill_viridis()




mycol <- c("white", "darkgrey")




heatmap_manon<-p1+scale_fill_gradientn(colours = mycol)





#cambio nomi 
heatmap_manon$data$Taxa<- gsub(".*:",
                               "", heatmap_manon$data$Taxa)

heatmap_manon+ theme_bw() +
  theme(
    axis.text.x = element_text(color="black"),
    axis.ticks = element_line(color = "black"))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y =element_text(colour="black") )+
 ggtitle("Core Microbiome") +
  geom_tile(
    aes(
      width = 0.5,
      height = 0.5))+theme(plot.title = element_text(hjust = 0.5))+ ylab("Genera")


    
theme(axis.text=element_text(size=14),axis.text.x = element_text(size=10),
      axis.title=element_text(size=14))




######



#correlazione vibrio con T

library(ggpubr)
cor(meta_table$Temperatura,meta_table$Richness_vibrio, method = "pearson")
cor.test(meta_table$Temperatura,meta_table$Richness_vibrio, method = "pearson")

sp<-ggplot(meta_table, aes(x=Richness_vibrio, y=Temperatura)) + 
  geom_point(aes(color=Mese),size=3)+
  geom_smooth(method=lm)
p<-sp + stat_cor(method = "pearson",r.digits = 3)+
  xlab("Vibrio") + ylab("Temperature (°C)")
p<-p+scale_color_manual(values = c("#D95F02", "#7570B3" ,"#66A61E" ,"#E6AB02" ,"#A6761D")) + labs(color='Samples')

newSTorder = c("FEB","MAY","JUL","SEP","NOV")
p$data$Mese<- as.character(p$data$Mese)
p$data$Mese <- factor(p$data$Mese, levels=newSTorder)
p







#multi scatter plot  

cor(psdf$Temperatura,psdf$Vibrio, method = "pearson")
cor.test(psdf$Temperatura,psdf$Vibrio, method = "pearson")






newSTorder = c("Amphritea", "Colwellia", "Photobacterium", "Pseudoalteromonas", 
               "Uncultured Arcobacteraceae", "Shewanella", 
               "Vibrio", "Tenacibaculum")
psdf$Genus<- as.character(psdf$Genus)
psdf$Genus <- factor(psdf$Genus, levels=newSTorder)




ggplot(psdf, aes(x=Temperatura, y=Abundance, color=Genus)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method=lm, se=FALSE,show.legend = FALSE)+ 
  stat_cor(method = "pearson",r.digits = 3,label.y = c(40,30,75,85,80,40,350,90),show.legend = FALSE)+
  facet_wrap(vars(Genus), scales = "free_y", ncol = 2, strip.position = "top")+ scale_color_hue(l=70, c=35)+
  theme(
    axis.text.x = element_text(color="black"),
    axis.ticks = element_line(color = "black"))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y =element_text(colour="black") )+ xlab("Temperature")



### color=c("#77DD77","#836953","#89cff0","#99c5c4","#9adedb","#aa9499","#aaf0d1","#b2fba5")









#microviz



library(dplyr)
library(phyloseq)
library(microViz)



x1



p<-GP2_X %>% 
  tax_fix(unknowns = c("uncultured"))%>%

    cor_heatmap(vars = c("Temperatura"))




x1 %>% 
  tax_fix(unknowns = c("uncultured"))%>%
  cor_heatmap(
    taxa = tax_top(x1, 15, by = max, rank = "Genus"),
   cor = "spearman"
  )





# set up the data with numerical variables and filter to top taxa
psq <- GP2 %>%
  tax_transform("identity", rank = "Genus")
#> Proportional min_prevalence given: 0.1 --> min 23/222 samples.

# randomly select 30 taxa from the 50 most abundant taxa (just for an example)
set.seed(123)
taxa <- sample(tax_top(psq))
# actually draw the heatmap
cor_heatmap(
  data = GP2, taxa = taxa)




otu_table(GP1)

GP2_X<-tax_glom(GP2, taxrank = "Genus")
psdf <- psmelt(GP2_X)

# Plot
psdf %>%
  ggplot(aes(x = Abundance,
                    y = Temperatura)) +
  geom_point() + 
  geom_smooth(method=lm)















