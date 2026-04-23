
#####                      loading packages                           ####

library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(readr)
library(plyr)
library(dplyr)
library(metagMisc)

#####                      importing data                           ####



### load data in the system:
#frequencies table

#######NULL per importare rows duplicate dhn bastardo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 abund_table<-read.csv("input/braken_all_REFSEQ_prokEprot_merged_fract_OTU.csv",row.names=NULL, check.names=FALSE,sep = ";")

## aggregare colonne ripetute facendo la media ( per dereplicare le repliche dei campioni)

 write.csv(aggregate(.~Genus,abund_table,mean), "input/OTU_REFSEQ_ALL_derep.csv")


 
 
 
abund_table<-read.csv("input/OTU_REFSEQ_ALL_derep.csv",row.names=1, check.names=FALSE,sep = ";")
abund_table<-t(abund_table)

head(abund_table)
ncol(abund_table)
nrow(abund_table)
#TAXONOMY table
OTU_taxonomy<-read.csv("input/braken_all_REFSEQ_prokEprot_merged_fract_TAX.csv",row.names=1,check.names=FALSE,sep = ";")
nrow(OTU_taxonomy)

##  Metatable daereplicato su excel, filtro avanzato, non mostrare duplicati!
meta_table<-read.csv("input/braken_all_REFSEQ_prokEprot_merged_fract_META.csv",row.names=1,check.names=FALSE,sep = ";")
colnames(meta_table)
nrow(meta_table)
ncol(meta_table)




#Convert the data to phyloseq format
OTU = otu_table(as.matrix(abund_table), taxa_are_rows = T)
TAX = tax_table(as.matrix(OTU_taxonomy))
SAM = sample_data(meta_table)
physeq<-merge_phyloseq(phyloseq(OTU, TAX, SAM))



# physeq<-subset_samples(physeq, (  Fraction=="0.2-5" ))




genus.sum = tapply(taxa_sums(physeq), tax_table(physeq)[, "Genus"], sum, na.rm=TRUE)
top5phyla = names(sort(genus.sum, TRUE))[1:7]
GP1 = prune_taxa((tax_table(physeq)[, "Genus"] %in% top5phyla), physeq)



############ FUNZIONE PER UNIRE I CAMPIONI FACENDONE LA MEDIA DELLE OTU



# https://github.com/joey711/phyloseq/issues/465


merge_samples_mean <- function(physeq, group){
  group_sums <- as.matrix(table(sample_data(physeq)[ ,group]))[,1]
  merged <- merge_samples(physeq, group)
  x <- as.matrix(otu_table(merged))
  if(taxa_are_rows(merged)){ x<-t(x) }
  out <- t(x/group_sums)
  out <- otu_table(out, taxa_are_rows = TRUE)
  otu_table(merged) <- out
  return(merged)
}
  





sample_data(GP1)
GP1

#  BUBLE PLOT TAXA PIU ABBONDANTI 
phyloseq_obj_css_Zone<- merge_samples_mean(GP1,"Zone")
sample_data(phyloseq_obj_css_Zone)




# bubble_plot(x = as.data.frame(otu_table(transform_sample_counts(phyloseq_obj_css_Zone,function(x) 100 * x))))
bubble_plot(x = as.data.frame(otu_table(phyloseq_obj_css_Zone)))   

BUBLE_plot<-bubble_plot(x = as.data.frame(otu_table(phyloseq_obj_css_Zone)))+scale_x_discrete(position = "top")

BUBLE_plot_5_2000<-BUBLE_plot
BUBLE_plot_0.2_5<-BUBLE_plot


ggarrange(BUBLE_plot_0.2_5,BUBLE_plot_5_2000,
          ncol = 2, nrow = 1,
          legend = "right",
          common.legend=F, align = "v"	)






#############  con faced grid  FRACTION




GP3<-subset_samples(GP1, !(  Fraction=="NA" ))

data_glom<- psmelt(GP3) # create dataframe from phyloseq object
data_glom$Genus <- as.character(data_glom$Genus)



 FIG1S<-ggplot(data_glom, aes(x = Zone, y = Genus)) + geom_point(aes(size = Abundance, 
      colour = Abundance), shape = 19, alpha = 0.9) + scale_size_continuous(name = "Counts ", 
         range = c(0, 16)) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = NULL, y = NULL)+  facet_wrap( ~factor(Fraction, levels=c("0.2-5",
                                       "5-20",
                                       "20-180","180-2000")),1)+
   scale_y_discrete(limits=rev)


 FIG1S+ theme(axis.text=element_text(colour="black"))
 
 
 
 #############  con faced grid  FRACTION
 
 sample_data(GP1)
 GP4<-subset_samples(GP1, !( ( Depth=="MIX" | Depth=="ZZZ"| Depth=="NA")))
 
 data_glom<- psmelt(GP4) # create dataframe from phyloseq object
 data_glom$Genus <- as.character(data_glom$Genus)
 
 
 
 FIG1S_2<-ggplot(data_glom, aes(x = Zone, y = Genus)) + geom_point(aes(size = Abundance, 
  colour = Abundance), shape = 19, alpha = 0.9) + scale_size_continuous(name = "Counts ", 
 range = c(0, 16)) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + 
   labs(x = NULL, y = NULL)+  facet_wrap( ~factor(Depth, levels=c(
     "SRF",
     "DCM",
     "MES")),1)+
   scale_y_discrete(limits=rev)
 
 FIG1S_2+ theme(axis.text=element_text(colour="black"))
 
 
 
 
 
 
 
 
 
 
 



newSTorder = rep(c( "Prochlorococcus","Pseudoalteromonas","Candidatus Pelagibacter",
                    "Polaribacter","Synechococcus" , "Alteromonas","Pseudomonas","Vibrio" ),13)
newSTorder =c( "Prochlorococcus","Pseudoalteromonas","Candidatus Pelagibacter",
   "Polaribacter","Synechococcus" , "Alteromonas","Pseudomonas","Vibrio" )
BUBLE_plot$data$Spec<- as.character(BUBLE_plot$data$Spec)
BUBLE_plot$data$Spec <- factor(BUBLE_plot$data$Spec, levels=newSTorder)  



##### subset VIBRIPPR HEATMAPS

VIBRIO <- subset_taxa(GP1, Genus=="Vibrio" )





#https://github.com/joey711/phyloseq/issues/293


#############              ZONA & DEPTH
VIBRIO <- subset_taxa(GP1, Genus=="Vibrio" )
sample_variables(VIBRIO)
sample_data(VIBRIO)


variable1 = as.character(get_variable(VIBRIO, "Zone"))
variable2 = as.character(get_variable(VIBRIO, "Water_Layer"))
# variable3 = as.character(get_variable(VIBRIO, "Distance1"))

sample_data(VIBRIO)$NewPastedVar <- mapply(paste, variable1, variable2  
                                           , sep = "_")
sample_data(VIBRIO)

VIBRIO_MERGED_zone_depth<-merge_samples_mean(VIBRIO, "NewPastedVar")

sample_variables(VIBRIO_MERGED_zone_depth)
sample_data(VIBRIO_MERGED_zone_depth)
sample_names(VIBRIO_MERGED_zone_depth)



#dati->testo in colonna così ho le nuove variabli

# write.csv(sample_names(VIBRIO_MERGED_zone_depth),"names_Zone_Water_layer.csv" )
#le importo
meta_zone_depth<-read.csv("file_conversione_nomi/names_Zone_Water_layer.csv",row.names=1, check.names=FALSE,sep = ";")

VIBRIO_MERGED_zone_depth
#cambio il metadata

sample_data(VIBRIO_MERGED_zone_depth)<-sample_data(meta_zone_depth)
VIBRIO_MERGED_zone_depth

sample_variables(VIBRIO_MERGED_zone_depth)
sample_data(VIBRIO_MERGED_zone_depth)



#   rimuovo campioni MIX e ZZZ perché non sono rappresentati



Samples_toRemove<-c("ANE_MIX", "PSW_MIX","ION_ZZZ","ARC_ZZZ","ARC_MIX" , "ANE_ZZZ" )
VIBRIO_MERGED_zone_depth_NOMIX_ZZZ<-subset_samples(VIBRIO_MERGED_zone_depth, !(  sample_names(VIBRIO_MERGED_zone_depth) %in% Samples_toRemove))


####################   forse dovrei trafsormalo in percentuale del totale ????

# physeq_perc=transform_sample_counts(VIBRIO_MERGED_zone_depth,function(x) 100 * x)

data_glom<- psmelt(VIBRIO_MERGED_zone_depth_NOMIX_ZZZ) # create dataframe from phyloseq object
data_glom$Genus <- as.character(data_glom$Genus) #convert to character
# data_glom$depth1 <- as.factor(data_glom$depth1 ,levels = c("5-10M","10-50M","50-100M","100-200M","200-500M","500-1000M")) #convert to character



p<-ggplot(data_glom, aes(Zone, Water_Layer, fill= Abundance)) + 
  geom_tile()+ theme(axis.text.x = element_text(angle = 90))+scale_fill_viridis_c() 
         p                     

       
         
         
#limits = c(0,6)
                  
p$data
  
  p$data
newSTorder = c("SRF","DCM","MES")  #  ,"100-200M","200-500M","500-1000M")c("100-1000M","10-100M","5-10M") 

p$data$Water_Layer<- as.character(p$data$Water_Layer)
p$data$Water_Layer <- factor(p$data$Water_Layer, levels=  rev(newSTorder))  

  p  + theme(panel.background = element_rect(fill = 'black'),panel.grid.major = element_blank(),
             panel.grid.minor = element_blank()) 

 DEPTH<- p  + theme(panel.background = element_rect(fill = 'black'),panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) 
  
    

ggplot(data_glom, aes(zone, depth, fill= Abundance)) + 
  geom_tile()+
  theme_ipsum() +scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 90))





#############              ZONA & Distance1
VIBRIO <- subset_taxa(GP1, Genus=="Vibrio" )
sample_variables(VIBRIO)
sample_data(VIBRIO)


variable1 = as.character(get_variable(VIBRIO, "Zone"))
variable2 = as.character(get_variable(VIBRIO, "Distance3"))
# variable3 = as.character(get_variable(VIBRIO, "Distance1"))

sample_data(VIBRIO)$NewPastedVar1 <- mapply(paste, variable1, variable2  
                                            , sep = "_")
VIBRIO_MERGED_zone_Distance1<-merge_samples_mean(VIBRIO, "NewPastedVar1")

sample_variables(VIBRIO_MERGED_zone_Distance1)
sample_data(VIBRIO_MERGED_zone_Distance1)
sample_names(VIBRIO_MERGED_zone_Distance1)



#dati->testo in colonna così ho le nuove variabli

write.csv(sample_names(VIBRIO_MERGED_zone_Distance1),"names_Zone_Distance.csv" )
#le importo
meta_zone_depth<-read.csv("file_conversione_nomi/names_Zone_Distance.csv",row.names=1, check.names=FALSE,sep = ";")

VIBRIO_MERGED_zone_Distance1
#cambio il metadata

sample_data(VIBRIO_MERGED_zone_Distance1)<-sample_data(meta_zone_depth)
VIBRIO_MERGED_zone_Distance1

sample_variables(VIBRIO_MERGED_zone_Distance1)
sample_data(VIBRIO_MERGED_zone_Distance1)




####################   forse dovrei trafsormalo in percentuale del totale ????

physeq_perc=transform_sample_counts(VIBRIO_MERGED_zone_Distance1,function(x) 100 * x)

data_glom<- psmelt(VIBRIO_MERGED_zone_Distance1) # create dataframe from phyloseq object
data_glom$Genus <- as.character(data_glom$Genus) #convert to character
# data_glom$depth1 <- as.factor(data_glom$depth1 ,levels = c("5-10M","10-50M","50-100M","100-200M","200-500M","500-1000M")) #convert to character

library(dplyr)
library(tidyr)

p<-ggplot(data_glom, aes(x=zone, y=distance, fill= Abundance)) + 
  geom_tile()+ theme(axis.text.x = element_text(angle = 90))+scale_fill_gradient( na.value="black")


p

p$data
newSTorder = c(
  "5-100Km",
"100-400Km",
"400-1400Km")



p$data$distance<- as.character(p$data$distance)
p$data$distance <- factor(p$data$distance, levels=rev(newSTorder) )

p  + theme(panel.background = element_rect(fill = 'black'),panel.grid.major = element_blank(),
           panel.grid.minor = element_blank()) +scale_fill_viridis_c() 
  


DIST<-p  + theme(panel.background = element_rect(fill = 'black'),panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank()) +scale_fill_viridis_c() 

######  salinità

VIBRIO <- subset_taxa(GP1, Genus=="Vibrio" )
VIBRIO_SRF_DCM<-subset_samples(VIBRIO, Env_feature == "DCM" | Env_feature == "SRF")


#############              ZONA & salinità
sample_variables(VIBRIO_SRF_DCM)
sample_data(VIBRIO_SRF_DCM)


variable1 = as.character(get_variable(VIBRIO_SRF_DCM, "Zone"))
variable2 = as.character(get_variable(VIBRIO_SRF_DCM, "S_1"))
# variable3 = as.character(get_variable(VIBRIO, "Distance1"))

sample_data(VIBRIO_SRF_DCM)$NewPastedVar <-  mapply(paste, variable1, variable2, sep = "_")
sample_data(VIBRIO_SRF_DCM)

VIBRIO_SRF_DCM_MERGED_zone_Sal<-merge_samples_mean(VIBRIO_SRF_DCM, "NewPastedVar")

sample_variables(VIBRIO_SRF_DCM_MERGED_zone_Sal)
sample_data(VIBRIO_SRF_DCM_MERGED_zone_Sal)
sample_names(VIBRIO_SRF_DCM_MERGED_zone_Sal)



#dati->testo in colonna così ho le nuove variabli

write.csv(sample_names(VIBRIO_SRF_DCM_MERGED_zone_Sal),"names_Zone_Sal.csv" )
#le importo
meta_zone_Sal<-read.csv("file_conversione_nomi/names_Zone_Sal.csv",row.names=1, check.names=FALSE,sep = ";")

VIBRIO_SRF_DCM_MERGED_zone_Sal
#cambio il metadata

sample_data(VIBRIO_SRF_DCM_MERGED_zone_Sal)<-sample_data(meta_zone_Sal)
VIBRIO_SRF_DCM_MERGED_zone_Sal

sample_variables(VIBRIO_SRF_DCM_MERGED_zone_Sal)
sample_data(VIBRIO_SRF_DCM_MERGED_zone_Sal)




####################   forse dovrei trafsormalo in percentuale del totale ????

# physeq_perc=transform_sample_counts(VIBRIO_MERGED_zone_depth,function(x) 100 * x)

data_glom<- psmelt(VIBRIO_SRF_DCM_MERGED_zone_Sal) # create dataframe from phyloseq object
data_glom$Genus <- as.character(data_glom$Genus) #convert to character
# data_glom$depth1 <- as.factor(data_glom$depth1 ,levels = c("5-10M","10-50M","50-100M","100-200M","200-500M","500-1000M")) #convert to character



p<-ggplot(data_glom, aes(Zone, Sal, fill= Abundance)) + 
  geom_tile()+ theme(axis.text.x = element_text(angle = 90))+scale_fill_viridis_c() 

#limits = c(0,6)

p$data

p$data
newSTorder = c("20-30","30-35","35-40")  

p$data$Sal<- as.character(p$data$Sal)
p$data$Sal <- factor(p$data$Sal, levels=rev(newSTorder))  

p  + theme(panel.background = element_rect(fill = 'black'),panel.grid.major = element_blank(),
           panel.grid.minor = element_blank()) 

SALY<- p  + theme(panel.background = element_rect(fill = 'black'),panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank()) 



ggplot(data_glom, aes(Zone, Sal, fill= Abundance)) + 
  geom_tile()+
  theme_ipsum() +scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 90))
















#### GGARRANGE 4 plots!!!!!!!!!!!!!!!


P1<-BUBLE_plot+theme(axis.line=element_blank(),
                 axis.text.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks=element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank())                

P2<-DEPTH+theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())
P3<-DIST+theme(axis.line=element_blank(),
           axis.text.x=element_blank(),
           axis.text.y=element_blank(),
           axis.ticks=element_blank(),
           axis.title.x=element_blank(),
           axis.title.y=element_blank())
#temp l'ho fatto più in basso!!!

P4<-TEMP+theme(axis.line=element_blank(),
               axis.text.x=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks=element_blank(),
               axis.title.x=element_blank(),
               axis.title.y=element_blank())

ggarrange(BUBLE_plot, DEPTH , DIST,TEMP,
          ncol = 1, nrow = 4,
          legend = "none",
          common.legend=F,heights = c(1.3, 0.3, 0.3, 0.3), align = "v"	)

library(ggpubr)
figuraAbund_genera<-ggarrange(P1, P2 , P3,P4,
             ncol = 1, nrow = 4,
             legend = "none",
             common.legend=F,heights = c(1.3, 0.3, 0.3, 0.3)	)

annotate_figure(c,top = text_grob("sea surface temperature", color = "Black", face = "bold", size = 14))



BUBLE_plot1<-BUBLE_plot+theme(axis.line=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks=element_blank(),
                     axis.title.x=element_blank() ) +
  theme(text=element_text(color="black"),axis.text=element_text(color="black"))           

DEPTH1<-DEPTH+theme(axis.line=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank())+
                  theme(text=element_text(color="black"),axis.text=element_text(color="black"))
                
DIST1<-DIST+theme(axis.line=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks=element_blank(),
               axis.title.x=element_blank())  +
  theme(text=element_text(color="black"),axis.text=element_text(color="black"))

TEMP1<-TEMP+theme(axis.line=element_blank(),
                  axis.text.x=element_blank(),
               axis.ticks=element_blank(),
               axis.title.x=element_blank()
)        
TEMP1$data
newSTorder = c(
  "-2-10",
  "10-20",
  "20-30")

TEMP1$data$Temp<- as.character(TEMP1$data$Temp)
TEMP1$data$Temp <- factor(TEMP1$data$Temp, levels=rev(newSTorder) )
TEMP1<-TEMP1+  theme(text=element_text(color="black"),axis.text=element_text(color="black"))

SALY1<-SALY+theme(axis.line=element_blank(),
                  axis.ticks=element_blank(),
                  axis.title.x=element_blank())+
  theme(text=element_text(color="black"),axis.text=element_text(color="black",size = 16))
SALY1




ALL=ggarrange(BUBLE_plot1, DEPTH1,DIST1,TEMP1,SALY1,
          ncol = 1, nrow = 5,
          legend = "none",
          common.legend=F,heights = c(1.3, 0.3, 0.3,0.3,0.42), align = "v"	)


ALL

ALL






























#############              ZONA & Temperature
sample_variables(VIBRIO)
sample_data(VIBRIO)
#                                   SRF 


VIBRIO_SRF<-subset_samples(VIBRIO, Env_feature == "SRF" )
sample_variables(VIBRIO_SRF)
sample_data(VIBRIO_SRF)


variable1 = as.character(get_variable(VIBRIO_SRF, "Zone"))
variable2 = as.character(get_variable(VIBRIO_SRF, "T_AVG"))
#variable3 = as.character(get_variable(VIBRIO_SRF, "T_AVG"))

sample_data(VIBRIO_SRF)$NewPastedVar <- mapply(paste, variable1, variable2, sep = "_")
sample_data(VIBRIO_SRF)

VIBRIO_SRF_MERGED_T<-merge_samples_mean(VIBRIO_SRF, "NewPastedVar")

sample_variables(VIBRIO_SRF_MERGED_T)
sample_data(VIBRIO_SRF_MERGED_T)
sample_names(VIBRIO_SRF_MERGED_T)


#dati->testo in colonna così ho le nuove variabli

###     fprse dovrei salvare questo             sample_data

write.csv(sample_names(VIBRIO_SRF_MERGED_T),"SRF_names_Zone_T.csv" )
#le importo
meta_zone_T<-read.csv("SRF_names_Zone_T.csv",row.names=1, check.names=FALSE,sep = ";")


z<-as.data.frame(sample_data(VIBRIO_SRF_MERGED_T))
z$T_AVG==meta_zone_T$T_SRF


VIBRIO_SRF_MERGED_T
#cambio il metadata

sample_data(VIBRIO_SRF_MERGED_T)<-sample_data(meta_zone_T)
VIBRIO_SRF_MERGED_T

sample_variables(VIBRIO_SRF_MERGED_T)
sample_data(VIBRIO_SRF_MERGED_T)




####################   forse dovrei trafsormalo in percentuale del totale ????

# physeq_perc=transform_sample_counts(VIBRIO_SRF_MERGED_T,function(x) 100 * x)

data_glom<- psmelt(VIBRIO_SRF_MERGED_T) # create dataframe from phyloseq object
data_glom$Genus <- as.character(data_glom$Genus) #convert to character
# data_glom$depth1 <- as.factor(data_glom$depth1 ,levels = c("5-10M","10-50M","50-100M","100-200M","200-500M","500-1000M")) #convert to character


############ buble plot SRF TEMP

  
 B<- ggplot(data_glom, aes(Zone, T_SRF)) + 
    geom_point(aes(size = Abundance,color=Abundance)) +
    scale_size_continuous(range = c(1, 10))
  
  
  
  
  
  
########    DCM
  
  
  
  
  #                                   DCM 
  
  
  VIBRIO_DCM<-subset_samples(VIBRIO, Env_feature == "DCM" )
  sample_variables(VIBRIO_DCM)
  sample_data(VIBRIO_DCM)
  
  
  variable1 = as.character(get_variable(VIBRIO_DCM, "Zone"))
  variable2 = as.character(get_variable(VIBRIO_DCM, "T_AVG"))
  #variable3 = as.character(get_variable(VIBRIO_DCM, "T_AVG"))
  
  sample_data(VIBRIO_DCM)$NewPastedVar <- mapply(paste, variable1, variable2, sep = "_")
  sample_data(VIBRIO_DCM)
  
  VIBRIO_DCM_MERGED_T<-merge_samples_mean(VIBRIO_DCM, "NewPastedVar")
  
  sample_variables(VIBRIO_DCM_MERGED_T)
  sample_data(VIBRIO_DCM_MERGED_T)
  sample_names(VIBRIO_DCM_MERGED_T)
  
  
  #dati->testo in colonna così ho le nuove variabli
  
  ###     fprse dovrei salvare questo             sample_data
  
  write.csv(sample_names(VIBRIO_DCM_MERGED_T),"DCM_names_Zone_T.csv" )
  #le importo
  meta_zone_T<-read.csv("DCM_names_Zone_T.csv",row.names=1, check.names=FALSE,sep = ";")
  
  
  z<-as.data.frame(sample_data(VIBRIO_DCM_MERGED_T))
  z$T_AVG==meta_zone_T$T_DCM
  
  
  VIBRIO_DCM_MERGED_T
  #cambio il metadata
  
  sample_data(VIBRIO_DCM_MERGED_T)<-sample_data(meta_zone_T)
  VIBRIO_DCM_MERGED_T
  
  sample_variables(VIBRIO_DCM_MERGED_T)
  sample_data(VIBRIO_DCM_MERGED_T)
  
  
  
  
  ####################   forse dovrei trafsormalo in percentuale del totale ????
  
  # physeq_perc=transform_sample_counts(VIBRIO_DCM_MERGED_T,function(x) 100 * x)
  
  data_glom<- psmelt(VIBRIO_DCM_MERGED_T) # create dataframe from phyloseq object
  data_glom$Genus <- as.character(data_glom$Genus) #convert to character
  # data_glom$depth1 <- as.factor(data_glom$depth1 ,levels = c("5-10M","10-50M","50-100M","100-200M","200-500M","500-1000M")) #convert to character
  
  
  ############ buble plot DCM TEMP
  
  
  A<-ggplot(data_glom, aes(Zone, T_DCM)) + 
    geom_point(aes(size = Abundance,color=Abundance)) +
    scale_size_continuous(range = c(1, 10))
  
  
A1<-A+  ggtitle("Vibrio frequencies in DCM") +
    xlab("Zone") + ylab("T (°C)  \n (from 17 to 188 m)")+ theme(plot.title = element_text(hjust = 0.5))

B1<-B+ggtitle("Vibrio frequencies in SRF") +
  xlab("Zone") + ylab("T (°C) \n (at a depth of 10 m) ") +theme(plot.title = element_text(hjust = 0.5))
theme(axis.line=element_blank(),
                                           axis.text.x=element_blank(),axis.ticks.x = element_blank(),
                                           axis.title.x=element_blank())
                                           







  library(ggpubr)
  figuraAbund_genera<-ggarrange(B1, A1,
                                ncol = 1, nrow = 2,
                                legend = "right",
                                common.legend=F	)
  
  annotate_figure(c,top = text_grob("sea surface temperature", color = "Black", face = "bold", size = 14))
  




####################    prova senzasubset SRF e DCM


  
  
  VIBRIO<-subset_samples(VIBRIO, Env_feature == "DCM" | Env_feature == "SRF")
  sample_variables(VIBRIO)
  sample_data(VIBRIO)
  
  
  variable1 = as.character(get_variable(VIBRIO, "Zone"))
  variable2 = as.character(get_variable(VIBRIO, "T_AVG"))
  variable3 = as.character(get_variable(VIBRIO, "depth3"))
  
  sample_data(VIBRIO)$NewPastedVar <- mapply(paste, variable1, variable2,variable3, sep = "_")
  sample_data(VIBRIO)
  
  VIBRIO_MERGED_T<-merge_samples_mean(VIBRIO, "NewPastedVar")
  
  sample_variables(VIBRIO_MERGED_T)
  sample_data(VIBRIO_MERGED_T)
  sample_names(VIBRIO_MERGED_T)
  
  
  #dati->testo in colonna così ho le nuove variabli
  
  ###     fprse dovrei salvare questo             sample_data
  
  write.csv(sample_names(VIBRIO_MERGED_T),"names_Zone_T_depth3.csv" )
  #le importo
  meta_zone_T<-read.csv("names_Zone_T_depth3.csv",row.names=1, check.names=FALSE,sep = ";")
  
  
  z<-as.data.frame(sample_data(VIBRIO_MERGED_T))
  z$T_AVG==meta_zone_T$Temp
  
  
  VIBRIO_MERGED_T
  #cambio il metadata
  
  sample_data(VIBRIO_MERGED_T)<-sample_data(meta_zone_T)
  VIBRIO_MERGED_T
  
  sample_variables(VIBRIO_MERGED_T)
  sample_data(VIBRIO_MERGED_T)
  
  
  
  
  ####################   forse dovrei trafsormalo in percentuale del totale ????
  genus.sum = tapply(taxa_sums(VIBRIO_MERGED_T), tax_table(VIBRIO_MERGED_T)[, "Genus"], sum, na.rm=TRUE)
  taxa_sums(VIBRIO_MERGED_T)

  physeq_perc=transform_sample_counts(VIBRIO_MERGED_T,function(x)  100 * x/ taxa_sums(VIBRIO_MERGED_T))
  otu_table(VIBRIO_MERGED_T)
  data_glom<- psmelt(VIBRIO_MERGED_T) # create dataframe from phyloseq object
  data_glom$Genus <- as.character(data_glom$Genus) #convert to character
  # data_glom$depth1 <- as.factor(data_glom$depth1 ,levels = c("5-10M","10-50M","50-100M","100-200M","200-500M","500-1000M")) #convert to character
  
  
  ############ buble plot DCM TEMP
  
  
MIXed<-ggplot(data_glom, aes(Zone, Temp)) + 
    geom_point(aes(size = Abundance,color=Abundance)) +
    scale_size_continuous(range = c(1, 10))

MIXed



diplox<-ggplot(data_glom, aes(Zone, Temp)) + 
  geom_point(aes(size = Abundance,color=Env)) +
  scale_size_continuous(range = c(1, 10))




Depth

dip_depth<-ggplot(data_glom, aes(Zone, Temp)) + 
  geom_point(aes(size = Abundance,color=Depth)) +
  scale_size_continuous(range = c(0, 10))




write.csv(  otu_table(physeq_perc),"OTU_TAB_temp%.csv" )



meta_table

ggscatter(meta_table, x = "T_AVG", y = "Depth", 
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = "T_AVG", ylab = "Depth")

ggplot(meta_table, aes(x=sort(T_AVG), y=(Depth))) + geom_point()

colnames(meta_table)

duplicated(colnames(meta_table))





###############     CON RANGES DI TEMPERATURA E SALINIT°



variable1 = as.character(get_variable(VIBRIO, "Zone"))
variable2 = as.character(get_variable(VIBRIO, "T_2"))


sample_data(VIBRIO)$NewPastedVar <- mapply(paste, variable1, variable2, sep = "_")
sample_data(VIBRIO)

VIBRIO_MERGED_T<-merge_samples_mean(VIBRIO, "NewPastedVar")

sample_variables(VIBRIO_MERGED_T)
sample_data(VIBRIO_MERGED_T)
sample_names(VIBRIO_MERGED_T)


#dati->testo in colonna così ho le nuove variabli

###     fprse dovrei salvare questo             sample_data

write.csv(sample_names(VIBRIO_MERGED_T),"names_Zone_T_1.csv" )
#le importo
meta_zone_T<-read.csv("names_Zone_T_1.csv",row.names=1, check.names=FALSE,sep = ";")


z<-as.data.frame(sample_data(VIBRIO_MERGED_T))
z$T_AVG==meta_zone_T$Temp


VIBRIO_MERGED_T
#cambio il metadata

sample_data(VIBRIO_MERGED_T)<-sample_data(meta_zone_T)
VIBRIO_MERGED_T

sample_variables(VIBRIO_MERGED_T)
sample_data(VIBRIO_MERGED_T)




####################   forse dovrei trafsormalo in percentuale del totale ????


data_glom<- psmelt(VIBRIO_MERGED_T) # create dataframe from phyloseq object
data_glom$Genus <- as.character(data_glom$Genus) #convert to character
# data_glom$depth1 <- as.factor(data_glom$depth1 ,levels = c("5-10M","10-50M","50-100M","100-200M","200-500M","500-1000M")) #convert to character


############ buble plot DCM TEMP

p<-ggplot(data_glom, aes(Zone, Temp, fill= Abundance)) + 
  geom_tile()+ theme(axis.text.x = element_text(angle = 90))+scale_fill_gradient( na.value="black")

TEMP<-p+ theme(panel.background = element_rect(fill = 'black'),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +scale_fill_viridis_c() 


ggplot(data_glom, aes(Zone, Temp)) + 
  geom_point(aes(size = Abundance,color=Abundance)) +
  scale_size_continuous(range = c(1, 10))














###                      ANALISI COMUNITA MICROBICA                           ####
#########  min taxa sample filter function




####                                 normalizzazione                 ###############

library(metagenomeSeq)
phyloseq_obj<-noempty_physeq


#filtering singletons #NON CI SONO! 
doubleton <- genefilter_sample(phyloseq_obj, filterfun_sample(function(x) x > 1), A=1)
doubleton <- prune_taxa(doubleton, phyloseq_obj) 
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

phyloseq_obj_css_contigs_spades
head(sample_names(phyloseq_obj_css))

####    add the metadata

meta_table_norm<-read.csv("input/metatable_table_ATL_normalized_sizeEdepth1_withtaradata.csv",row.names=1, check.names=FALSE,sep = ";")
sample_data(phyloseq_obj_css_contigs_spades)<-sample_data(meta_table_norm)
sample_variables(phyloseq_obj_css)



##salvare OTU_TABLE


#salvare l'otu table 
# Extract abundance matrix from the phyloseq object
#OTU_matrix = as(sample_data(phyloseq_obj_css), "matrix")
#OTUdf = as.data.frame(OTU_matrix)
# head(OTUdf)
#write.csv(sample_data(phyloseq_obj_css_contigs_spades),"output_csv/Contigs_spades_samplDATA_table_ATL_normalized.csv")













p<-plot_heatmap(phyloseq_obj_css, title="TARA SAMPLES  ordered by Depth",method = "PCoA")+
  theme(plot.title = element_text(hjust = 0.5))

p+scale_fill_gradientn(     colours=c("black","darkgreen","green","red","darkred"),
                            limits=c(0,15),
                            breaks = c(0,5,10,15))   +labs(y="Species")+
  theme (axis.text.y = element_text(size=7,colour = "black"))+
  theme(axis.text.x = element_text(size=7,colour = "black"))












#                                     Alpha diversity                                  ####



otu.absol<-round(otu_table(phyloseq_obj_css))
head(otu.absol)
# physeq_normalized<-merge_phyloseq(phyloseq(OTU_absol, TAX),SAM)
phyloseq_obj_css <-phyloseq_obj_css
otu_table(phyloseq_obj_css) <- otu_table(otu.absol)





head(otu_table(physeq_normalized))
sample_data(physeq_normalized)
sample_names(physeq_normalized)


### plottare alfa div per stazione con colore zona

plot_richness(physeq_normalized,x="Station",color = "Zone",measures=c("Observed", "Shannon", "Simpson"),title = "Alfa Diversity")+
  geom_boxplot()+scale_color_manual(values = c("#007FFF","#1aff00","#ff0008","#ff7003"))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))







plot_richness(physeq_normalized,measures=c("Observed", "Shannon", "Simpson"),title = "Alfa Diversity",color="Station1") 




#write.csv(p$data, "alfadiv.csv")



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

#####                      BETA DIVERSITY                           ####
sample_variables(phyloseq_obj_css)
#PCoA on Bray-Curtis Dissimilarity
phyloseq_obj_css
phyloseq_obj_css_Zone
phyloseq_obj_css_Station
head(sample_data(phyloseq_obj_css_Station))

write.csv(sample_names(phyloseq_obj_css_Station), "samples_station_filtered.csv")
sample_names(phyloseq_obj_css_Station)


meta_table_station<-read.csv("samples_station_filtered.csv",row.names=1, check.names=FALSE)
SAM_station = sample_data(meta_table_station)
sample_data(phyloseq_obj_css_Station) <- sample_data(SAM_station)
head(sample_data(phyloseq_obj_css_Station))


sample_variables(phyloseq_obj_css_Zone)
sample_data(phyloseq_obj_css_ANE_surface)




#PCoA on Bray-Curtis Dissimilarity
phyloseq_obj_css_surface = subset_samples(phyloseq_obj_css, depth1=="Surface")
phyloseq_obj_css_surface_fraction02.5 = subset_samples(phyloseq_obj_css, depth1=="Surface" | fraction1=="0.22-5")
phyloseq_obj_css_ANE_surface = subset_samples(phyloseq_obj_css_surface, Zone=="ANE")

# phyloseq_obj_css_surface_only5M<-subset_samples(phyloseq_obj_css_surface, Depth != "9")


otu.ord <- ordinate(physeq = phyloseq_obj_css_ANE_surface, "PCoA")



##beta diversity  

#asse 1-2  

plot_ordination(physeq = phyloseq_obj_css_ANE_surface, otu.ord,color = "Station",
                axes =c(1,2))+
  theme_bw()+ geom_point(size = 2)+
  geom_text(aes(label=Station), size = 3, vjust = 0,hjust=0)+
  theme(plot.title = element_text(hjust = 0.5))+ stat_ellipse(type = "t", linetype = 2) 
stat_ellipse(type = "t") 

scale_color_manual(values = c("#007FFF","#1aff00","#ff0008","#ff7003"))

#asse 1-3 
plot_ordination(physeq = phyloseq_obj_css_surface, otu.ord,color = "Zone",
                axes =c(1,3))+
  theme_bw()+ geom_point(size = 2)+
  geom_text(aes(label=Station), size = 3, vjust = 0,hjust=0)+
  theme(plot.title = element_text(hjust = 0.5))
scale_color_manual(values = c("#007FFF","#1aff00","#ff0008","#ff7003"))

title ="PCoA on Bray-Curtis Dissimilarity axes 1&3 \n Permanova p = 0.001" 
c<-ggarrange(a, b , 
             labels = c("A", "B"),
             ncol = 2, nrow = 1,
             legend = "bottom",
             common.legend=T	)

annotate_figure(c,top = text_grob("PCoA on Bray-Curtis Dissimilarity \n Permanova p = 0.001", color = "Black", face = "bold", size = 14))
##beta diversity  



GP.ord <- ordinate(phyloseq_obj_css, "PCoA", "bray")
plot_ordination(phyloseq_obj_css, GP.ord, type="taxa", color="Species", title="taxa")



plot_ordination(phyloseq_obj_css, GP.ord, type="split", color="Species", shape="Zone", label="Station", title="split") 






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

genus.sum = tapply(taxa_sums(physeq_normalized_contigs_spades), tax_table(physeq_normalized_contigs_spades)[, "Species"], sum, na.rm=TRUE)
top5phyla = names(sort(genus.sum, TRUE))[1:20]
GP1 = prune_taxa((tax_table(physeq_normalized_contigs_spades)[, "Species"] %in% top5phyla), physeq_normalized_contigs_spades)

GP2 <- subset_taxa(GP1, Species != "Unassigned")
GP2 <- subset_taxa(GP2, Species != "uncultured")



p<-plot_bar(GP2, "Species", fill="Species", facet_grid=~Zone) #Tow , facet_grid=~Sample
sample_data(GP2)
p + geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



p+scale_fill_manual(values =c( "#542750","#7FFF00","#00A86B","#e6a910","#7B1B02",
                               "#aedce8","#884DA7","#E52B50","#293133","#7FFFD4",
                               "#708090","#1025e6","#fc0317","#568b8c","#f59de5",
                               "#960018","#000000","#007FFF","#ff0303","#66562b",
                               "#F4C430","#FFBF00"
))

newSTorder =  c("5-160M",
                "5M","9M","25M","30M",
                "35M","40M","50M",
                "80M","85M","100M",
                "120M","125M","150M",
                "250M","300M","390M",
                "590M","640M","700M",
                "740M","800M")

p$data$Sample<- as.character(p$data$Sample)
p$data$Sample <- factor(p$data$Sample, levels=newSTorder)




write.csv(as.data.frame(p$data), "data_barplot_mostabundant.csv")
#


##############rimovere linee orizzontali nere noiose

my_plot_bar = function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                        facet_grid = NULL) {
  mdf = psmelt(physeq)
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



p<-my_plot_bar(GP2, "Species", fill = "Species")
scale_fill_manual(values =c( "#7FFF00","#00A86B","#e6a910","#7B1B02",
                             "#aedce8","#884DA7","#E52B50","#293133","#7FFFD4",
                             "#708090","#1025e6","#fc0317","#FFBF00","#ff0303","#4B0082"))


plot_bar(GP2, "Species", fill="Species") #Tow , facet_grid=~Sample













############ lefse ################## 

otu_css_names <- taxa_names(GP2) 
gp2_physeq <- prune_taxa(otu_css_names,physeq)
write.csv(as.data.frame(otu_table(gp2_physeq)), "otu_12_most_abund_gp2_physeq.csv")
taxa_names(GP2) == taxa_names(gp2_physeq)

colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus","Species")
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



plot_ef_bar(mm, label_level = 6) +
  scale_fill_manual(values = c("GLTO20141129" = "blue", "GLTO20150414" = "green", "GLTO20160416"="red", "GLTO20160816"="orange"))

plot_cladogram(mm, color = c(GLTO20141129 = "blue", GLTO20150414 = "green", GLTO20160416="red", GLTO20160816="orange")) + 
  theme(legend.position = "none")





