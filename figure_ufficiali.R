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
library(ggpubr)
library(ggside)
library(ggdist)
#####                      importing data                           ####






devtools::install_github("jbgb13/peRReo") 

devtools::install_github("tylermorganwall/rayshader")



abund_table<-read.csv("input/OTUtable_ALL_ENTERO_derep.csv",row.names=1, check.names=FALSE,sep = ";")
abund_table<-t(abund_table)

#TAXONOMY table
OTU_taxonomy<-read.csv("input/braken_all_ENTERO_prokEprot_merged_fract_INPUT_R_OTU_TAX.csv",row.names=1,check.names=FALSE,sep = ";")
nrow(OTU_taxonomy)


meta_table<-read.csv("input/braken_all_REFSEQ_prokEprot_merged_fract_METAWcooRd.csv",row.names=1,check.names=FALSE,sep = ";")

#Convert the data to phyloseq format
OTU = otu_table(as.matrix(abund_table), taxa_are_rows = T)
TAX = tax_table(as.matrix(OTU_taxonomy))
SAM = sample_data(meta_table)
physeq<-merge_phyloseq(phyloseq(OTU, TAX, SAM))

#####################################################################################################################

#####   alfadiversity

#####################################################################################################################

###############   OTUTABLE CEILING

OTU_ceiling = otu_table(as.matrix(ceiling(abund_table)), taxa_are_rows = T)

physeqCEILING<-merge_phyloseq(phyloseq(OTU_ceiling, TAX, SAM))

physeqCEILING<-subset_samples(physeqCEILING, !(  Fraction=="NA" ))
# 



#####################################################################################################################
###FRACTION 
#####################################################################################################################
#utilizzo i dati nel plot di phyloseq per plottarli in ggpubrr 

p<-plot_richness(physeqCEILING,x="Fraction3",color = "Fraction3",measures=c("Observed"),title = "Vibrio Richness")


my_comparisons <- list( c("PAB", "FLB"))

#tutto insieme 
ggboxplot(p$data, x = "Fraction3", y = "value",
          color = "Fraction3", add = "jitter",palette = c("#FF3355", "#9C9596"),outlier.shape = NA)+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y =250, label.x =0.5)     + 
  ggtitle("Vibrio richness x fraction")+ theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(colour="black"))+xlab("Fraction") + ylab('Richness')
  



#per tutte le zone: 

ggboxplot(p$data, x = "Fraction3", y = "value",
          color = "Fraction3", add = "jitter",palette = c("#FF3355", "#9C9596"),outlier.shape = NA)+
  stat_compare_means(comparisons = my_comparisons,label.y = 200)+ # Add pairwise comparisons p-value
  ggtitle("Vibrio richness x fraction")+ theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(colour="black"))+
  facet_wrap( ~Zone)+  xlab("Fraction") + ylab('Richness%')








#####################################################################################################################
###    DEPTH 
#####################################################################################################################
physeqCEILING<-merge_phyloseq(phyloseq(OTU_ceiling, TAX, SAM))

physeqCEILING<-subset_samples(physeqCEILING, !( ( Depth=="MIX" | Depth=="ZZZ"| Depth=="NA")))

#utilizzo i dati nel plot di phyloseq per plottarli in ggpubrr 
p<-plot_richness(physeqCEILING,x="Fraction3",color = "Depth",measures=c("Observed"),title = "Vibrio Richness")

my_comparisons1 <- list( c("SRF", "DCM"), c("SRF", "MES"), c("DCM", "MES") )



ggboxplot(p$data, x = "Depth", y = "value",
          color = "Depth",add = "jitter",  palette = c("#66aaee", "#76d05a","#032d83"),outlier.shape = NA)+
  stat_compare_means(comparisons = my_comparisons1)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y =325)     + 
  ggtitle("Vibrio richness x Depth")+ theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(colour="black"))+  xlab("Depth") + ylab('Richness')



#all zones
##
ggboxplot(p$data, x = "Depth", y = "value",
          color = "Depth",add = "jitter",  palette = c("#66aaee", "#76d05a","#032d83"),outlier.shape = NA)+
  stat_compare_means(comparisons = my_comparisons1, label.y = c(180,190,210))+ # Add pairwise comparisons p-value
  stat_compare_means(label.y =50,)     + 
  ggtitle("Vibrio richness x Depth")+ theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(colour="black"))+facet_wrap( ~Zone)+  xlab("Depth") + ylab('Richness')




# RICHNESS ALL OCEANS per latutude and long 
newSTorder =  c("SRF",
                "DCM","MES")
p1<-plot_richness(physeqCEILING,x="Latitude",color = "Depth",measures=c("Observed"),title = "Vibrio Species Richness")+facet_grid(~Fraction3)+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+geom_point(size = 2,alpha=1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_smooth()
p1$data$Depth<- as.character(p1$data$Depth)
p1$data$Depth <- factor(p1$data$Depth, levels=newSTorder)
p1+scale_color_manual(values=c("#66aaee", "#76d05a","#032d83"))


p1<-plot_richness(physeqCEILING,x="Longitude",color = "Depth",measures=c("Observed"),title = "Vibrio Species Richness")+facet_grid(~Fraction3)+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+geom_point(size = 2,alpha=1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_smooth()+scale_color_manual(values=c("#66aaee", "#76d05a","#032d83"))



p1$data$Depth<- as.character(p1$data$Depth)
p1$data$Depth <- factor(p1$data$Depth, levels=newSTorder)
p1+scale_color_manual(values=c("#66aaee", "#76d05a","#032d83"))


# RICHNESS ALL OCEANS per latutude and long  COME ISME!!!!!

p1<-plot_richness(physeqCEILING,x="Latitude",color = "Fraction3",measures=c("Observed"),title = "Vibrio Species Richness")+facet_grid(~Depth)+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+geom_point(size = 2,alpha=1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_smooth()+scale_color_manual(values=c("#FF3355", "#9C9596"))
p1$data$Depth<- as.character(p1$data$Depth)
p1$data$Depth <- factor(p1$data$Depth, levels=newSTorder)
#con ggside 
p1+scale_color_manual(values=c("#FF3355", "#9C9596"))+geom_ysideboxplot(aes(y =value), orientation = "x") +
  theme(        ggside.panel.scale.y = .4)


#longitudine 
p1<-plot_richness(physeqCEILING,x="Longitude",color = "Fraction3",measures=c("Observed"),title = "Vibrio Species Richness")+facet_grid(~Depth)+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+geom_point(size = 2,alpha=1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_smooth(size=2)
p1$data$Depth<- as.character(p1$data$Depth)
p1$data$Depth <- factor(p1$data$Depth, levels=newSTorder)
#con ggside 
p1+scale_color_manual(values=c("#FF3355", "#9C9596"))+geom_ysideboxplot(aes(y =value), orientation = "x") +
  theme(        ggside.panel.scale.y = .4)+  scale_xsidey_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
  scale_ysidex_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) 




################################  bellissime
###plot GGPLOT con geom ponit alpha 
p1<-plot_richness(physeqCEILING,x="Latitude",color = "Fraction3",measures=c("Observed"),title = "Vibrio Species Richness")

#latitudine 

min(p1$data$Latitude)
max(p1$data$Latitude)

ggplot(p1$data,aes(Latitude,value))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+  geom_point(aes(colour = Fraction3),alpha =0.45)+
  scale_color_manual(values=c("#FF3355", "#9C9596"))+facet_grid(~Depth)+
  geom_ysideboxplot(aes(x=Fraction3,y=value,colour = Fraction3), orientation = "x") +
  theme(        ggside.panel.scale.y = .4)+scale_ysidex_discrete()+
  geom_smooth(aes(colour = Fraction3),size=1.5)+ ylab('Richness') + labs(color='Fraction')  +
  scale_x_continuous(breaks = seq(-65, 80, by = 20))#, expand = c(0, 0)  




#longitudine 

min(p1$data$Longitude)
max(p1$data$Longitude)
ggplot(p1$data,aes(Longitude,value))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+  geom_point(aes(colour = Fraction3),alpha =0.45)+
  scale_color_manual(values=c("#FF3355", "#9C9596"))+facet_grid(~Depth)+
  geom_ysideboxplot(aes(x=Fraction3,y=value,colour = Fraction3), orientation = "x") +
  theme(        ggside.panel.scale.y = .4)+scale_ysidex_discrete()+
  geom_smooth(aes(colour = Fraction3),size=1.5)+ ylab('Richness')+ labs(color='Fraction')+
  scale_x_continuous(breaks = seq(-170, 175, by = 45))#, expand = c(0, 0) 





#####################################################################################################################
###    BETA DIV 
#####################################################################################################################

############merge average per zona


variable1 = as.character(get_variable(physeq, "Zone"))
variable2 = as.character(get_variable(physeq, "Fraction3"))
# variable3 = as.character(get_variable(physeq, "Distance1"))

sample_data(physeq)$NewPastedVar <- mapply(paste, variable1, variable2  
                                           , sep = "_")
sample_data(physeq)
# write.csv(sample_data(physeq),"sample_data_physeq.csv")#ho messo i nomi come quella dei kmers
# sample_data(physeq)<-sample_data(read.csv("sample_data_physeq.csv",row.names=1, check.names=FALSE,sep = ";"))


physeq_MERGED_zone_FRACTION<-merge_samples_mean(physeq, "NewPastedVar")#cambiato er fare la mantel, per fare b div ricalcolarlo DHN 







sample_variables(physeq_MERGED_zone_FRACTION)
sample_data(physeq_MERGED_zone_FRACTION)
sample_names(physeq_MERGED_zone_FRACTION)


#per far quello che facevo su excel
rm(df)
RR<-nrow(as.data.frame(sample_names(physeq_MERGED_zone_FRACTION)))
df <- data.frame(matrix(ncol = 1, nrow = RR))
df$names<-(as.data.frame(sample_names(physeq_MERGED_zone_FRACTION)))

df<-df[,2]
colnames(df)[1] <- "Samples"

foo <- data.frame(do.call('rbind', strsplit(as.character(df$Samples),'_',fixed=TRUE)))
df$Zone <-foo$X1
df$Fraction <-foo$X2
colnames(df)

rownames(df) <- df[,1]

# df$Ocean<-c(rep("Atlantic",9),rep("Polar",5),rep("Atlantic",8),rep("Indian",9),rep("Sea",5),rep("Pacific",12),
#            rep("Sea",4),rep("Polar",4))
df$Ocean<-c(rep("Atlantic",4),rep("Polar",2),rep("Atlantic",4),rep("Indian",4),rep("Sea",2),rep("Pacific",6),
            rep("Sea",2),rep("Polar",2))

physeq_MERGED_zone_FRACTION
#cambio il metadata




sample_data(physeq_MERGED_zone_FRACTION)<-sample_data(df)
physeq_MERGED_zone_FRACTION

sample_variables(physeq_MERGED_zone_FRACTION)
sample_data(physeq_MERGED_zone_FRACTION)



physeq_MERGED_zone_FRACTION1<-subset_samples(physeq_MERGED_zone_FRACTION, !(  Fraction=="NA" ))
# physeq_MERGED_zone_FRACTION1<-subset_samples(physeq_MERGED_zone_FRACTION, !(  Samples=="IOS_20-180" ))



##beta diversity  

otu.ord <- ordinate(physeq = physeq_MERGED_zone_FRACTION1, "PCoA")

#asse 1-2  
library(ggforce)
a<-plot_ordination(physeq = physeq_MERGED_zone_FRACTION1, otu.ord,color = "Fraction",shape = "Ocean",
                   axes =c(1,2))+
  geom_text(aes(label=Zone), Fraction = 3, vjust = 0,hjust=0, show.legend = FALSE)+
  theme(plot.title = element_text(hjust = 0.0))+ geom_point(Fraction = 2)+ theme_void()+theme_bw()
a1<-a+  geom_mark_ellipse(aes(color = Fraction), show.legend = FALSE)+ theme_void()+theme_bw() + geom_xsidedensity(aes(y=stat(density),fill=Fraction), alpha = 0.5, show.legend = FALSE) +
  geom_ysidedensity(aes(x=stat(density),fill=Fraction), alpha = 0.5, show.legend = FALSE) +
  theme_bw() +  scale_xsidey_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
  scale_ysidex_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +scale_ysidex_discrete()+
  ggside::theme_ggside_void() 


a1

#####  dovrei correlare le due matridi di distanza fatte con kmers e taxonomia 



OTU_matrix = as(otu_table(physeq_MERGED_zone_FRACTION), "matrix") #FRACTION no FRACTION1  
# transpose if necessary
#if(taxa_are_rows(OTU_matrix)){OTU_matrix <- t(OTU_matrix)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU_matrix)
head(OTUdf)
# write.csv(OTUdf,"OTU_table_bacteria_normalized.csv")
library(vegan)
#export dist-matrix
BC_TAX<-vegdist(t(OTUdf),method = "bray")
head(BC_TAX)

BC_kmers=as.matrix(read.table("mat_abundance_braycurtis.csv",sep=";", header=TRUE, row.names=1))

BC_kmers[upper.tri(BC_kmers)] <- 0
BC_kmersDist <- as.dist(BC_kmers, diag = TRUE)



set.seed(1234)
mantel(xdis =BC_TAX ,ydis = BC_kmersDist,method = "pearson",permutations = 9999)

Mantel statistic r: 0.7614 
Significance: 1e-04 


### plottare kmers betadiv con ggplot





physeq_MERGED_zone_FRACTION
sample_data(physeq_MERGED_zone_FRACTION1)

iMDS  <- ordinate(physeq_MERGED_zone_FRACTION, "PCoA", distance=BC_kmersDist) #metto la matrice di distanza  KMERS con l'altro phyloseq  



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ho messo i metadati come il phyloseq merged1 = tax quello fatto con la tax


# write.csv(sample_data(physeq_MERGED_zone_FRACTION),"sample_data_physeq_MERGED_zone_FRACTION.csv")#ho messo i nomi come quella dei kmers
# sample_data(physeq_MERGED_zone_FRACTION)<-sample_data(read.csv("sample_data_physeq_MERGED_zone_FRACTION.csv",row.names=1, check.names=FALSE,sep = ";"))


#lo plotto

p<-plot_ordination(physeq_MERGED_zone_FRACTION, iMDS,color = "Fraction",shape = "Ocean",
                   axes =c(1,2))+
  geom_text(aes(label=Zone), Fraction = 3, vjust = 0,hjust=0, show.legend = FALSE)+
  theme(plot.title = element_text(hjust = 0.0))+ geom_point(Fraction = 2)+ theme_void()+theme_bw()
p1<-p+  geom_mark_ellipse(aes(color = Fraction), show.legend = FALSE)+ theme_void()+theme_bw() + geom_xsidedensity(aes(y=stat(density),fill=Fraction), alpha = 0.5, show.legend = FALSE) +
  geom_ysidedensity(aes(x=stat(density),fill=Fraction), alpha = 0.5, show.legend = FALSE) +
  theme_bw() +  scale_xsidey_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
  scale_ysidex_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +scale_ysidex_discrete()+
  ggside::theme_ggside_void() 





#betadiversity metaphlan----> no buona!!!!!!!!

BC_metaphlan=as.matrix(read.csv2("metaphlan4_merged_abundance_table_bray-curtis.1tsv.csv", sep=";",header=TRUE, row.names=1,check.names = F))
head(BC_metaphlan)


BC_metaphlan[upper.tri(BC_metaphlan)] <- 0
BC_metaphlanDist <- as.dist(BC_metaphlan, diag = TRUE)

#tolgo quelli non classificati da metaphlan
physeq_MERGED_zone_FRACTION_metaphlan = subset_samples(physeq_MERGED_zone_FRACTION, Samples != "SOC_Protist Fraction" & Samples != "SOC_Prokaryotes Fractions" &  Samples !="ASE_Prokaryotes Fractions")

sample_names(physeq_MERGED_zone_FRACTION_metaphlan)




iMDS_metaphlan  <- ordinate(physeq_MERGED_zone_FRACTION_metaphlan, "PCoA", distance=BC_metaphlanDist) #metto la matrice di distanza  KMERS con l'altro phyloseq  

physeq_MERGED_zone_FRACTION

plot_ordination(physeq_MERGED_zone_FRACTION_metaphlan, iMDS_metaphlan,color = "Size",shape = "Ocean",
                axes =c(1,2))+
  geom_text(aes(label=Zone), Size = 3, vjust = 0,hjust=0, show.legend = FALSE)+
  theme(plot.title = element_text(hjust = 0.0))+ geom_point(Size = 2)+ theme_void()+theme_bw()












########################################################################################################################
#################################  CORRELAZIONI CON VAR AMB ###############################################  


  https://david-barnett.github.io/microViz/reference/cor_heatmap.html




library(ggplot2)
library(microViz)
library(dplyr)
#per vedere online il phyloseq obj 
ord_explore(physeq)



#barplot 
physeq %>%
  comp_barplot(
    tax_level = "Species", n_taxa = 15, other_name = "Other",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    palette = distinct_palette(n = 15, add = "grey90"),
    merge_other = FALSE, bar_outline_colour = "darkgrey"
  ) +
  coord_flip() +
  facet_wrap("Fraction3", nrow = 1, scales = "free") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())



library(viridis)


htmp <- physeq %>%
  ps_mutate(Fraction3 = as.character(Fraction3)) %>%
  tax_transform("log2", add = 1, chain = TRUE) %>%
  comp_heatmap(
    taxa = tax_top(physeq, n = 30), grid_col = NA, name = "Log2p",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    colors = heat_palette(palette = viridis::turbo(11)),
    row_names_side = "left", row_dend_side = "right", sample_side = "bottom",
    sample_anno = sampleAnnotation(
      Fraction = anno_sample_cat(
        var = "Fraction3", col = c( FLB= "grey35", PAB = "grey85"),
        box_col = NA, legend_title = "Fraction", size = grid::unit(4, "mm")
      )
    )
  )



ComplexHeatmap::draw(
  object = htmp, annotation_legend_list = attr(htmp, "AnnoLegends"),
  merge_legends = TRUE
)







library(dplyr)

####### correlazioni

physeq_1<-subset_samples(physeq, !( ( Depth=="MIX" | Depth=="ZZZ"| Depth=="NA")))




# set up the data with numerical variables and filter to top taxa
psq <- physeq_1 %>%
  ps_mutate(
    Depth = recode(Depth, SRF = 3, DCM = 2, MES = 1),
    Fraction = if_else(Fraction3 == "FLB", true = 1, false = 0),
  ) %>%
  tax_filter(
    tax_level = "Species", min_prevalence = 1 / 10, min_sample_abundance = 1 / 10
  )  %>%
  tax_transform("identity", rank = "Species")
#> Proportional min_prevalence given: 0.1 --> min 23/222 samples.

# randomly select 30 taxa from the 50 most abundant taxa (just for an example)
set.seed(123)
taxa <- sample(tax_top(psq, n = 15), size = 15)
# actually draw the heatmap
cor_heatmap(
  data = psq, taxa = taxa,
  taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
  tax_anno = taxAnnotation(
    Prev. = anno_tax_prev(undetected = 15),
    Log2 = anno_tax_box(undetected = 15, trans = "log2", zero_replace = 1)
  )
)


####################################  RDA




physeq %>%
  ps_mutate(
    Depth = as.numeric(Depth == "SRF"),
    Fraction = as.numeric(Fraction3 == "FLB"),
  ) %>%
  tax_transform("clr", rank = "Species") %>%
  ord_calc(
    constraints = c("Depth", "Fraction"),
    # method = "RDA", # Note: you can specify RDA explicitly, and it is good practice to do so, but microViz can guess automatically that you want an RDA here (helpful if you don't remember the name?)
    scale_cc = FALSE # doesn't make a difference
  ) %>%
  ord_plot(
    colour = "Depth", size = 2, alpha = 0.5, shape = "Fraction",
    plot_taxa = 1:8)




#quello buono 
physeq %>%
  ps_mutate(
    Depth = recode(Depth, SRF = 3, DCM = 2, MES = 1),
    Fraction = if_else(Fraction3 == "FLB", true = 1, false = 0)
    ) %>%
  tax_transform("clr", rank = "Species") %>%
  ord_calc(
    constraints = c("Depth", "Fraction"),
    # method = "RDA", # Note: you can specify RDA explicitly, and it is good practice to do so, but microViz can guess automatically that you want an RDA here (helpful if you don't remember the name?)
    scale_cc = FALSE # doesn't make a difference
  ) %>%
  ord_plot(
    colour = "Fraction", size = 2, alpha = 0.5, shape = "Depth",
    plot_taxa = 1:30)














physeq %>%
  ps_calc_dominant(
    rank = "Species", other = "Other", none = "Not dominated",
    threshold = 0.4, n_max = 3
  ) %>%
  tax_transform(rank = "Species", trans = "clr") %>%
  ord_calc("PCA") %>%
  ord_plot(colour = "Fraction3", size = 3, alpha = 0.6) 











#Contour colour only
ggplot(CPF, aes(x = day, y = depth, z = temp)) +
  geom_tile(aes(fill = temp)) + 
  stat_contour() +
  scale_y_reverse(name="depth (m)") +
  theme_bw()

#Contour lines only
contourplot(temp ~ day * depth, data = CPF)








#cca 


plot_cca(physeq = physeq, grouping_column = "Fraction3", pvalueCutoff = 0.01, 
         env.variables = NULL, num.env.variables = NULL, exclude.variables = NULL, 
         draw_species = F)


sample_variables(physeq)

class(sample_data(physeq)$Fraction3)

library(phylosmith)
abundance_heatmap(physeq, classification = 'Species',  treatment = 'Fraction3', transformation = 'log10')

variable_correlation_heatmap(physeq, treatment = 'Fraction3',
                             method = 'spearman')
sample_data(physeq)
