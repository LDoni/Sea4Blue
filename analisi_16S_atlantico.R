#da provare https://biorgeo.github.io/bioregion/index.html









https://docs.oceanparcels.org/en/latest/
https://github.com/andrew-s28/nccs-transport
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')

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
library(sf)
library("rnaturalearth")
 library(raster)
library("RColorBrewer")
 



#https://github.com/joey711/phyloseq/issues/1163

#load data PS
ps <- readRDS("ps.rds")

head(t(otu_table(ps)))


#phyloseq_ill<-ps


#####---------------------------------------------------->>>>> fatto per il confronto, ma poi l'ho fatto nel nuovo file
#phyloseq_ill <- phyloseq(otu_table(ps), tax_table(ps), sample_data(ps))
#ps_nanopore <- readRDS("physeq_nanopore_spaghetti.rds")
#head(otu_table(ps_nanopore))

#taxa_names(ps_nanopore) <- paste0("ps_nan_", taxa_names(ps_nanopore))

#taxa_names(phyloseq_ill) <- paste0("ps_ill_", taxa_names(phyloseq_ill))


#phyloseq_obj_ALL<-merge_phyloseq(ps_nanopore,phyloseq_ill)

#phyloseq_obj_ALL<- tax_glom(phyloseq_obj_ALL, taxrank = 'Genus')

### l'ho modificato 
#write.csv(sample_data(phyloseq_obj_ALL),"metadata_sea4blue.csv")
#sample_data(phyloseq_obj_ALL)<-sample_data(read.csv("metadata_sea4blue.csv",row.names=1,check.names=FALSE,sep = ";"))

#saveRDS(phyloseq_obj_ALL, "phyloseq_obj_ALL.rds")

#phyloseq_obj_ALL <- readRDS("phyloseq_obj_ALL.rds")




###        ----   >>>>   mappa longrust

sample_data(ps)
sample_names(ps)
 
   

coord.df <- ps %>%
  subset_samples(!(sample.description %in% c("Laposnegative", "Lapospositive"))) %>%
  meta() %>%
  dplyr::select(sample.description, Latitude, Longitude)


ps %>%
  subset_samples(!(sample.description %in% c("Laposnegative", "Lapospositive"))) %>%
  meta() %>%
  dplyr::select(Nominativo.campione., Latitude, Longitude,data)




colnames(coord.df)<-c("id","lat", "lon")
head(coord.df)



points_sf <- st_as_sf(coord.df, coords = c("lon", "lat"), crs = 4326)

longhurst <- st_read("longhurst_v4_2010/Longhurst_world_v4_2010.shp")
longhurst <- st_make_valid(longhurst)
longhurst <- st_transform(longhurst, 4326)
# Spatial join
result <- st_join(points_sf, longhurst["ProvCode"], join = st_intersects)

head(as.data.frame(result))

world <- ne_countries(scale = "medium", returnclass = "sf")
crs_used <- 4326  # EPSG:4326
province_usate <- unique(result$ProvCode)
longhurst_subset <- longhurst %>% filter(ProvCode %in% province_usate)

ggplot() +
  geom_sf(data = world, fill = "#D9D9D9", color = "gray70", size = 0.3) +
  geom_sf(data = longhurst_subset, aes(fill = ProvCode), color = NA, alpha = 0.4) +
  geom_sf(data = result, color = "black", size = 2) +
  geom_sf_text(data = result, aes(label = id), size = 3, nudge_y = 1, check_overlap = TRUE) +
  scale_fill_viridis_d(name = "Longhurst Province") +
  coord_sf(crs = crs_used, xlim = c(-90, 10), ylim = c(25, 45)) +
  theme_minimal() +
  labs(title = "Longhurst Provinces and Sample Points",
       subtitle = "Proiezione WGS84 (EPSG:4326)",
       x = "Longitude", y = "Latitude")


sample_data(ps)$longhurst_provinces <- as.data.frame(result)$ProvCode[match(sample_data(ps)$sample.description, as.data.frame(result)$id)]

#scarico sst annuale da https://oceandata.sci.gsfc.nasa.gov/file_search/ 
wget --auth-no-challenge --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies "https://oceandata.sci.gsfc.nasa.gov/ob/getfile/AQUA_MODIS.20220101_20221231.L3m.YR.SST.sst.4km.nc"

 






crop.vals_Atl <- c(lat = c(-90,10), 
                   lon = c(25, 45))

SST_Atl <- ncdf4::nc_open("AQUA_MODIS.20220101_20221231.L3m.YR.SST.sst.4km.nc") %>%
  oceanmap::nc2raster(., "sst") %>%
  raster::flip(., "y") %>%
  raster::crop(., raster::extent(crop.vals_Atl))

library(oceanmap)

oceanmap::v(SST_Atl, cbpos = "b", pal = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(300)),
            zlim = c(0,35), 
            cb.xlab = expression("Annual SST (°C)"),
            bwd = 0.01, grid = F, replace.na = F, border = "#504f4f",
            cex.ticks = 1, axeslabels = F, figdim = c(4,5), show.colorbar = T)

SST_df <- as.data.frame(SST_Atl, xy = TRUE)
names(SST_df)[3] <- "sst" 
points_sf$longhurst_provinces <- as.data.frame(result)$ProvCode


ggplot() +
  # Sfondo SST
  geom_raster(data = SST_df, aes(x = x, y = y, fill = sst)) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")),
                       name = "Annual SST (°C)",
                       limits = c(0, 35),
                       na.value = "transparent") +
  
  # Mappa dei continenti
  geom_sf(data = world, fill = "#D9D9D9", color = "gray70", size = 0.3) +
  
  # Sovrapposizione poligoni Longhurst senza fill
  geom_sf(data = longhurst_subset, color = "gray60", fill = NA, size = 0.4) +
  
  # Punti colorati per provincia
  geom_sf(data = points_sf, aes(color = longhurst_provinces), size = 3) +
  scale_color_viridis_d(name = "Longhurst Province") +
  
  # Etichette
  geom_sf_text(data = points_sf, aes(label = id), size = 3, nudge_y = 1, check_overlap = TRUE) +
  
  # Coordinate e stile
  coord_sf(crs = 4326, xlim = c(-90, 10), ylim = c(25, 45), expand = FALSE) +
  theme_minimal() +
  labs(title = "Sea Surface Temperature and Sample Points",
       subtitle = "SST from MODIS + Longhurst Provinces",
       x = "Longitude", y = "Latitude")




###                      ANALISI COMUNITA MICROBICA                           ####


####                                 normalizzazione                 ###############
library(dplyr)
library(umap)
library(phyloseq)
library(microbiome)
library(MiscMetabar)
library(tidyverse)
library(vegan)
library(picante)
library(Biostrings)
library(ape)
#importoape#importo i dati scaricati in reticulate_var_amb_giuste.R
all_data<-read.csv("dati_scaricati_copernicus.csv")
ps <- readRDS("ps.rds")
otu_table(ps) <- otu_table(t(otu_table(ps)), taxa_are_rows = TRUE)
head(tax_table(ps))
phyloseq_obj<-ps
### RIMOZIONE CHL E MITOK
grep(pattern = "Mitochondria", tax_table(phyloseq_obj)) 
grep(pattern = "Chloroplast", tax_table(phyloseq_obj)) 
phyloseq_obj <- phyloseq_obj %>% subset_taxa( Family!= "Mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class) ) 
phyloseq_obj <- subset_taxa(phyloseq_obj, (tax_table(phyloseq_obj)[,"Order"]!="Chloroplast") | is.na(tax_table(phyloseq_obj)[,"Order"]))
phyloseq_obj = subset_samples(phyloseq_obj,!( Nominativo.campione.=="Positive" | Nominativo.campione.=="Negative"))
sample_names(phyloseq_obj)<-sample_data(phyloseq_obj)$Nominativo.campione. 
# Unisci i metadati esistenti con all_data
merged_metadata <- merge(sample_data(phyloseq_obj)%>% meta(), all_data, 
                         by = "Nominativo.campione.", 
                         all.x = TRUE)
rownames(merged_metadata) <- merged_metadata$Nominativo.campione.
# Aggiorna i metadati nell'oggetto phyloseq
sample_data(phyloseq_obj) <- merged_metadata

# Verifica che i nuovi dati siano stati aggiunti
sample_data(phyloseq_obj)
phyloseq_obj<-clean_pq(phyloseq_obj)




### da questo phyloseq_obj estraggo i dati che poi verrranno usati per l'analisi dei meccanismi ecologici di comunità
dir.create("Data", recursive = TRUE)
dir.create("Data/Meta_Data/Prok", recursive = TRUE)
dir.create("Data/Count_Data/Processed/Prok", recursive = TRUE)
dir.create("Data/Count_Data/Fasta/Prok", recursive = TRUE)
dir.create("Data/Trees/Prok", recursive = TRUE)
dir.create("Data/Processed/Prok", recursive = TRUE)
dir.create("Data/Taxonomy", recursive = TRUE)

#metadata
sample_data(phyloseq_obj) %>%
   as_tibble() %>%
  rename(Sample_ID = Nominativo.campione.) %>%
  mutate(Sample_ID = as.character(Sample_ID)) %>%
  relocate(Sample_ID, .before = 1)%>%
  write_tsv("Data/Meta_Data/Prok/Meta_Data.tsv")

# OTU TABLE

otu<-otu_table(phyloseq_obj) %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID") 
  write_tsv(otu,"Data/Count_Data/Processed/Prok/Full_Prok_Count.tsv")  
  write_tsv(otu,"Data/Processed/Prok/Full_Prok_Count.tsv")
 
#   Sequenze FASTA
writeXStringSet(refseq(phyloseq_obj), "Data/Count_Data/Fasta/Prok/Full_Prok_Sequences.fasta")

# 4.4 Albero filogenetico
write.tree(phy_tree(phyloseq_obj), "Data/Trees/Prok/Prok_Combined.tree")

# 4.5 Tassonomia
tax_table(phyloseq_obj) %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID") %>%
  write_tsv("Data/Taxonomy/Prok_Taxonomy.tsv")







####                                 normalizzazione                 ###############
 
#filtering singletons
doubleton <- genefilter_sample(phyloseq_obj, filterfun_sample(function(x) x > 1), A=1)
doubleton <- prune_taxa(doubleton, phyloseq_obj) 



## rimuovo i pos e neg
sample_data(doubleton)
sample_variables(doubleton)
(sample_data(doubleton)$Nominativo.campione.)

 #rimuovo i positivi e negativi
PosNeg = subset_samples(doubleton, Nominativo.campione. == "Positive" | Nominativo.campione.=="Negative")
head(sample_data(PosNeg))
Pos = subset_samples(PosNeg, sample.description == "Lapospositive" )

#rimuovo i bluehole
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

#cambio i nomi
#sample_names(phyloseq_obj_css)<-sample_data(phyloseq_obj_css)$Nominativo.campione.






#phyloseq_obj_css = subset_samples(phyloseq_obj_css, Giorno.di.navigazione. == "1" |
#                                                 Giorno.di.navigazione.=="11"|
 #                                                Giorno.di.navigazione.=="17"|
 #                                                Giorno.di.navigazione.=="25"|
 #                                                Giorno.di.navigazione.=="26")
sample_names(phyloseq_obj_css)<-sample_data(phyloseq_obj_css)$Nominativo.campione.

#plot heatmap  
p<-plot_heatmap(tax_glom(phyloseq_obj_css, taxrank = "Genus"),taxa.label = "Genus",first.taxa = "Vibrio",sample.order = c("eDNA1","eDNA3","eDNA4","eDNA6","eDNA7","eDNA9","eDNA15","eDNA16","eDNA17","eDNA18","eDNA19","eDNA25","eDNA26","eDNA27","eDNA28"))

p+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+theme_bw(base_size = 9,"black")
p<-p+theme(axis.text.y = element_text(color = "black", size = 6, angle = 0, hjust = 1, vjust = 0, face = "plain"),axis.text.x = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"))
p




#rifare physeq con numero otu assolute e non scalate CSS in modo che siano interi
# e che vadano bene per il calcolo indici alpha div dopo

otu.absol<-round(otu_table(phyloseq_obj_css))
head(otu.absol)
physeq_normalized <- phyloseq_obj_css
otu_table(physeq_normalized) <- otu_table(as.matrix(otu.absol), taxa_are_rows = T)


#                                     Alpha diversity                                  ####

# Trasforma l'otu table in presenza/assenza (1=presente, 0=assente)
physeq_pa <- transform_sample_counts(phyloseq_obj_css, function(x) ifelse(x > 0, 1, 0))
head(otu_table(physeq_pa))
sample_variables(physeq_pa)
sample_data(physeq_pa)



###plot GGPLOT con geom ponit alpha 
p1<-plot_richness(physeq_pa,x="Longitude",color = "Zona",measures=c("Observed"))
#Longitudine 
newSTorder = c( "ANW","ANC","ANE")
p1$data$Zona<- as.character(p1$data$Zona)
p1$data$Zona <- factor(p1$data$Zona, levels=newSTorder) 


 
#facet_grid(~Depth)+
ggplot(p1$data,aes(Longitude,value))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour = Zona))+
  scale_color_manual(values=c("#000000","#4280fc", "#ffb452","#f7170a"))+
  geom_ysideboxplot(aes(x=Zona,y=value,colour = Zona), orientation = "x") +
  theme(        ggside.panel.scale.y = .4)+scale_ysidex_discrete()+
  geom_smooth(aes(colour = variable ),size=1.5)+ ylab('Richness') 
 # scale_x_continuous(breaks = seq(-65, 80, by = 20))#, expand = c(0, 0)    theme(axis.text.x = element_text(angle = 45, hjust = 1))+


ggboxplot(p1$data, x = "Zona", y = "value",
               color = "Zona", palette = "jco",
               add = "jitter")+ stat_compare_means()  #### Kruskal-Wallis p=0.06








scale_color_hue(l=60, c=48)
scale_color_brewer(palette='Spectral')


show_col(hue_pal(l=70, c=90,direction = -1)(50))





show_col(hue_pal(l=60, c=48)(6))
pseq_tutt<-physeq_normalized










#####                      BETA DIVERSITY                           ####

#PCoA on Bray-Curtis Dissimilarity
phyloseq_obj_css
metad
library(ggforce)

library(ggpubr)
library(ggrepel)
head(sample_data(phyloseq_obj_css))


### bray-curtis

otu.ord <- ordinate(physeq = phyloseq_obj_css, "PCoA", distance = "bray")


#plot ordination

A<-plot_ordination(phyloseq_obj_css, otu.ord,  color="Zona",axes =c(1,2))+
  scale_color_manual(values = c("#4280fc", "#ffb452","#f7170a"))+ geom_point(size=3)+ 
  geom_text_repel(aes(label=Giorno.di.navigazione.),max.overlaps = Inf, show.legend = FALSE)+
  geom_mark_ellipse(aes(color = Zona), show.legend = FALSE)+ theme_void()+theme_bw() + 
  geom_xsidedensity(aes(y=stat(density),fill=Zona), alpha = 0.5, show.legend = FALSE) +
  geom_ysidedensity(aes(x=stat(density),fill=Zona), alpha = 0.5, show.legend = FALSE) +
 scale_xsidey_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
  scale_ysidex_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +scale_ysidex_discrete()+
  ggside::theme_ggside_void()  +
  scale_fill_manual(values = c("#4280fc", "#ffb452","#f7170a"))+ggtitle("bray-Curtis")
A
B<-plot_ordination(phyloseq_obj_css, otu.ord,  color="Zona",axes =c(1,3))+
  scale_color_manual(values = c("#4280fc", "#ffb452","#f7170a"))+ geom_point(size=3)+ 
  geom_text_repel(aes(label=Giorno.di.navigazione.),max.overlaps = Inf, show.legend = FALSE)+
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
bray <- phyloseq::distance(phyloseq_obj_css, method = "bray")
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(phyloseq_obj_css))

colonne_da_selezionare <- c("Nominativo.campione.", "data.x", "Latitude.x", "Longitude.x")




write.csv(sampledf %>% 
            dplyr::select(all_of(colonne_da_selezionare)), "parcels/stations.csv")




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

sample_data(phyloseq_obj_css)$Zona <-
  factor(sample_data(phyloseq_obj_css)$Zona)
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



set.seed(1)

# Calculate bray curtis distance matrix
uniF <- phyloseq::distance(phyloseq_obj_css, method = "wunifrac")
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(phyloseq_obj_css))

adonis2(uniF ~ Zona, data = sampledf) #0.001 ***

adonis2(uniF ~ Zona+Giorno.di.navigazione.+thetao +Longitude.x, data = sampledf,by = "margin")












library(phyloseq)
library(Matrix)
library(foreach)
library(doParallel)
library(ape)
library(vegan)
library(phyloseq)
library(dplyr)


###distance decay

library(phyloseq)
library(dplyr)

# Zona come fattore
sample_data(phyloseq_obj_css)$Zona <-
  factor(sample_data(phyloseq_obj_css)$Zona)

# Coordinate e temperatura come numeriche
sample_data(phyloseq_obj_css)$Latitude.y  <- as.numeric(sample_data(phyloseq_obj_css)$Latitude.y)
sample_data(phyloseq_obj_css)$Longitude.y <- as.numeric(sample_data(phyloseq_obj_css)$Longitude.y)
sample_data(phyloseq_obj_css)$thetao      <- as.numeric(sample_data(phyloseq_obj_css)$thetao)

bray_dist <- phyloseq::distance(
  phyloseq_obj_css,
  method = "bray"
)
unifrac_dist <- phyloseq::distance(
  phyloseq_obj_css,
  method = "unifrac",
  weighted = FALSE
)

library(geosphere)
library(tibble)
#distanza geografica
coords <- data.frame(
  lon = sample_data(phyloseq_obj_css)$Longitude.y,
  lat = sample_data(phyloseq_obj_css)$Latitude.y
)

geo_dist <- distm(coords, fun = distHaversine) / 1000
geo_dist <- as.dist(geo_dist)


#distanza ambientale T
temp_dist <- dist(sample_data(phyloseq_obj_css)$thetao)

library(tibble)

dd_df <- tibble(
  geographic_km = as.vector(geo_dist),
  delta_temp    = as.vector(temp_dist),
  bray          = as.vector(bray_dist),
  unifrac       = as.vector(unifrac_dist)
) %>%
  drop_na()

#similaruty
dd_df <- dd_df %>%
  mutate(
    bray_sim    = 1 - bray,
    unifrac_sim = 1 - unifrac
  )

#Bray–Curtis vs distanza geografica

ggplot(dd_df, aes(geographic_km, bray_sim)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() 



#UniFrac vs distanza geografica
ggplot(dd_df, aes(geographic_km, unifrac_sim)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw()


#Bray–Curtis vs differenza di temperatura
ggplot(dd_df, aes(delta_temp, bray_sim)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() 



#-->(più le condizioni di temperatura sono simili più le cumunità sono simili)


#UniFrac vs differenza di temperatura
ggplot(dd_df, aes(delta_temp, unifrac_sim)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() 



ggplot(dd_df, aes(delta_temp, geographic_km)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() 


dd_long <- dd_df %>%
  dplyr::select(
    geographic_km,
    delta_temp,
    bray_sim,
    unifrac_sim
  ) %>%
  tidyr::pivot_longer(
    cols = c(bray_sim, unifrac_sim),
    names_to = "metric",
    values_to = "similarity"
  ) %>%
  tidyr::pivot_longer(
    cols = c(geographic_km, delta_temp),
    names_to = "distance_type",
    values_to = "distance"
  ) %>%
  dplyr::mutate(
    relation = paste(metric, distance_type, sep = " vs ")
  ) %>%
  dplyr::group_by(relation) %>%
  dplyr::mutate(
    distance_z = scale(distance)[,1],
    similarity_z = scale(similarity)[,1]
  ) %>%
  ungroup()



ggplot(dd_long, aes(distance_z, similarity_z)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ relation, scales = "free") +
  labs(
    x = "Normalized distance (z-score)",
    y = "Normalized community similarity (z-score)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )



dd_plot <- dd_df %>%
  dplyr::select(
    geographic_km,
    delta_temp,
    bray_sim,
    unifrac_sim
  ) %>%
  tidyr::pivot_longer(
    cols = c(bray_sim, unifrac_sim),
    names_to = "metric",
    values_to = "similarity"
  ) %>%
  tidyr::pivot_longer(
    cols = c(geographic_km, delta_temp),
    names_to = "driver",
    values_to = "distance"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(metric,
                           bray_sim = "Bray–Curtis",
                           unifrac_sim = "UniFrac"),
    driver = dplyr::recode(driver,
                           geographic_km = "Geographic distance",
                           delta_temp   = "Temperature difference"),
    group_id = paste(metric, driver, sep = "_")
  ) %>%
  dplyr::group_by(group_id) %>%
  dplyr::mutate(
    distance_z   = scale(distance)[,1],
    similarity_z = scale(similarity)[,1]
  ) %>%
  dplyr::ungroup()
ggplot(dd_plot,
       aes(x = distance_z,
           y = similarity_z,
           color = metric,
           linetype = driver)) +
  
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  
  labs(
    x = "Normalized distance (z-score)",
    y = "Normalized community similarity (z-score)",
    color = "Similarity metric",
    linetype = "Driver"
  ) +
  
  scale_linetype_manual(
    values = c(
      "Geographic distance" = "solid",
      "Temperature difference" = "dashed"
    )
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    panel.grid.minor = element_blank()
  )


#statistics

summary(lm(bray ~ geographic_km, data = dd_df))
summary(lm(unifrac ~ geographic_km, data = dd_df))

library(vegan)

mantel(bray_dist, geo_dist, permutations = 9999)
mantel(unifrac_dist, geo_dist, permutations = 9999)

# Controllando per temperatura
mantel.partial(bray_dist, geo_dist, temp_dist,permutations = 9999)
mantel.partial(unifrac_dist, geo_dist, temp_dist,permutations = 9999)






head(meta(phyloseq_obj_css))




##environmental drivers


# Variabili ambientali (escluso spazio)
env_data <- sample_data(phyloseq_obj_css)[, c(
  "thetao", "so", "chl", "o2", "no3", "po4", "si"
)]

env_dist <- dist(scale(env_data))

# Mantel ambientale
mantel(bray_dist, env_dist,permutations = 9999)
mantel(unifrac_dist, env_dist,permutations = 9999)





###oceanography

# Velocità risultante delle correnti
curr_speed <- with(sample_data(phyloseq_obj_css),
                   sqrt(uo^2 + vo^2))

curr_dist <- dist(curr_speed)

# Test: correnti vs comunità
mantel(bray_dist, curr_dist,permutations = 9999)
mantel(unifrac_dist, curr_dist,permutations = 9999)






# Subset tassonomico--> trovare il genus più abbondante

ps_genus <- tax_glom(phyloseq_obj_css, taxrank = "Genus")
abund_mat <- as(otu_table(ps_genus), "matrix")

# Se taxa sono righe
if (!taxa_are_rows(ps_genus)) {
  abund_mat <- t(abund_mat)
}

total_abundance <- rowSums(abund_mat)
most_abundant_taxon <- names(which.max(total_abundance))
most_abundant_taxon


total_abundance[most_abundant_taxon]

tax_table(ps_genus)[most_abundant_taxon, ]

pattern_list <- c("Prochlorococcus", "Synechococcus", "Clade Ia", "Clade Ib")

tax <- as.data.frame(tax_table(ps_genus))
otu_mat <- as(otu_table(ps_genus), "matrix")
if (!taxa_are_rows(ps_genus)) {
  otu_mat <- t(otu_mat)
}
keep_idx <- grepl(paste(pattern_list, collapse="|"),
                  tax$Genus,
                  ignore.case = TRUE)
keep_taxa <- rownames(tax)[keep_idx]

# --- 5. Calcola abbondanza totale dei generi selezionati ---
selected_abund <- rowSums(otu_mat[keep_taxa, , drop=FALSE])
total_community <- sum(otu_mat)
relative_abund <- selected_abund / total_community * 100
plot_df <- data.frame(
  Genus = tax$Genus[keep_idx],         # usa il nome del genus
  Total_Abundance = selected_abund,
  Relative_Abundance = relative_abund,
  row.names = NULL                     # opzionale per rimuovere ASV IDs
)

print(plot_df)
ggplot(plot_df, aes(x = reorder(Genus, -Relative_Abundance), y = Relative_Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  ylab("Relative Abundance (%) of Entire Community") +
  xlab("Genus / Clade") +
  ggtitle("Relative Abundance of Selected Genera / Clades vs Entire Community") +
  theme_minimal() +
  theme(legend.position = "none")



##
####script per le analisi/tutto su proclorococcus_vs_clade1A.R
##

library(vegan)

# Distance matrix della community
library(vegan)

# 1. OTU table come matrice
otu_mat <- as(otu_table(phyloseq_obj_css), "matrix")
if (!taxa_are_rows(phyloseq_obj_css)) {
  otu_mat <- t(otu_mat)
}

# 2. Variabili ambientali
meta_df <- as.data.frame(sample_data(phyloseq_obj_css))

temp <- data.frame(thetao = as.numeric(meta_df$thetao), row.names = rownames(meta_df))
nutrients <- data.frame(
  chl = as.numeric(meta_df$chl),
  o2 = as.numeric(meta_df$o2),
  no3 = as.numeric(meta_df$no3),
  po4 = as.numeric(meta_df$po4),
  si = as.numeric(meta_df$si),
  row.names = rownames(meta_df)
)

# Controllo righe
all(rownames(temp) == rownames(otu_mat))   # deve essere TRUE
all(rownames(nutrients) == rownames(otu_mat))  # deve essere TRUE
# 3. Variance partitioning
varpart_res <- varpart(t(otu_mat), temp, nutrients)

# 4. Plot
plot(varpart_res)





#### barplot###



genus.sum = tapply(taxa_sums(phyloseq_obj_css), tax_table(phyloseq_obj_css)[, "Genus"], sum, na.rm=TRUE)
top5phyla = names(sort(genus.sum, TRUE))[1:30]
GP1 = prune_taxa((tax_table(phyloseq_obj_css)[, "Genus"] %in% top5phyla), phyloseq_obj_css)


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
  phyloseq_obj_css,
  taxa_rank="Genus",
  wilcoxon_cutoff = 0.05,
  group = "Zona",
  kw_cutoff = 0.05,
  multigrp_strat = T,
  lda_cutoff = 4)
rm(mm_lefse)
length(rownames(marker_table(mm_lefse)))


mm_lefse <- run_lefse(
  phyloseq_obj_css,
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
