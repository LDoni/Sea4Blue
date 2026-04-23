

library(ggplot2)
library(phyloseq)
library(tidyverse)
library(randomcoloR)


#   OTU al 97 %
#load data
OTU_taxonomy_barplot<-read.csv("input/tax_normalized97_unass.csv",row.names=1,check.names=FALSE)
TAX_BP = tax_table(as.matrix(OTU_taxonomy_barplot))
dim(TAX_BP)


#FARE IL PHYLOSEQ
phyloseq_obj_css_barplot<-merge_phyloseq(otu_table(phyloseq_obj_css),sample_data(phyloseq_obj_css) ,TAX_BP)





# MERGING PHYLOSEQ
phyloseq_obj_css_barplot_merged_period<- merge_samples(phyloseq_obj_css_barplot,"Tow")




#trasformazione in%

#GUARDARE CHE PHYLOSEQ STO TRSFORMANDO!!!!

physeq_perc=transform_sample_counts(phyloseq_obj_css,function(x) 100 * x/sum(x))


 (names(sort(taxa_sums(glom), TRUE))[1:15])
 
#Turn all OTUs into phylum (or phylum or order level) counts

#GUARDARE CHE PHYLOSEQ STO UTILIZZANDO

glom <- tax_glom(physeq_perc, taxrank = 'Species')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Genus <- as.character(data_glom$Genus) #convert to character
# data_glom$Sample <- factor(data_glom$Sample, levels =c("Before","During","After"))

data_glom$Sample <- factor(data_glom$Sample, levels=c("A4","A5","A6","A8","A9","A10","A11","A12","A13","A14","A15","A16","A17","A18","A20","A21") )
# data_glom$Sample <- factor(data_glom$Sample, levels =c("GLTO20141129" ,"GLTO20150414", "GLTO20160416","GLTO20160816"))

write.csv(data_glom, "dataglom_merged_tow1.csv")


Count = length(unique(data_glom$Genus))
Count
sort(unique(data_glom$Genus))



data_glom$Genus <- factor(data_glom$Genus, levels =c("0319-6G20",
                                                     "4572-13",
                                                     "Acidovorax",
                                                     "Acinetobacter",
                                                     "Actinomyces",
                                                     "ADurb.Bin063-1",
                                                     "AEGEAN-169_marine_group",
                                                     "Aeromonas",
                                                     "Alicycliphilus",
                                                     "Aliiglaciecola",
                                                     "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
                                                     "Alteromonas",
                                                     "Anoxybacillus",
                                                     "AqS1",
                                                     "Aquabacterium",
                                                     "Aquitalea",
                                                     "Arctic97B-4_marine_group",
                                                     "Arenimonas",
                                                     "AT-s2-59",
                                                     "Aurantimicrobium",
                                                     "B2M28",
                                                     "Bacillus",
                                                     "Balneola",
                                                     "BD2-11_terrestrial_group",
                                                     "BD7-11",
                                                     "BD7-8",
                                                     "Bdellovibrio",
                                                     "BIrii41",
                                                     "Blastopirellula",
                                                     "Bradymonadaceae",
                                                     "Bradymonadales",
                                                     "Bradyrhizobium",
                                                     "Brevibacillus",
                                                     "Brevundimonas",
                                                     "Burkholderia-Caballeronia-Paraburkholderia",
                                                     "C39",
                                                     "Candidatus_Accumulibacter",
                                                     "Candidatus_Flaviluna",
                                                     "Candidatus_Fritschea",
                                                     "Candidatus_Hepatoplasma",
                                                     "Candidatus_Kaiserbacteria",
                                                     "Candidatus_Methylopumilus",
                                                     "Candidatus_Obscuribacter",
                                                     "Candidatus_Omnitrophus",
                                                     "Candidatus_Planktoluna",
                                                     "Candidatus_Planktophila",
                                                     "Candidatus_Rhabdochlamydia",
                                                     "Candidatus_Tenderia",
                                                     "Cellvibrio",
                                                     "Cetobacterium",
                                                     "Chthoniobacter",
                                                     "CL500-29_marine_group",
                                                     "CL500-3",
                                                     "Clade_Ib",
                                                     "Clade_IV",
                                                     "Clostridiisalibacter",
                                                     "Cobetia",
                                                     "Collinsella",
                                                     "Coraliomargarita",
                                                     "Corynebacterium",
                                                     "Coxiella",
                                                     "Crocosphaera",
                                                     "Cupriavidus",
                                                     "Curvibacter",
                                                     "Cutibacterium",
                                                     "cvE6",
                                                     "Cyanobium_PCC-6307",
                                                     "Cyanothece_ATCC_51142_(UCYN-C)",
                                                     "Dechloromonas",
                                                     "Defluviitaleaceae_UCG-011",
                                                     "Delftia",
                                                     "DEV007",
                                                     "Diaphorobacter",
                                                     "DSSD61",
                                                     "Duganella",
                                                     "Eel-36e1D6",
                                                     "Ellin6055",
                                                     "Endozoicomonas",
                                                     "Enhydrobacter",
                                                     "Enterococcus",
                                                     "Enterovibrio",
                                                     "Epulopiscium",
                                                     "Erysipelotrichaceae_UCG-003",
                                                     "Escherichia-Shigella",
                                                     "Ethanoligenens",
                                                     "EUB33-2",
                                                     "Exiguobacterium",
                                                     "FCPU426",
                                                     "Ferrimonas",
                                                     "Finegoldia",
                                                     "Friedmanniella",
                                                     "Fusobacterium",
                                                     "Galbitalea",
                                                     "Gemella",
                                                     "Gemmatimonas",
                                                     "Geobacillus",
                                                     "Georgenia",
                                                     "Grimontia",
                                                     "Haemophilus",
                                                     "Halioglobus",
                                                     "Halodesulfovibrio",
                                                     "Halomonas",
                                                     "Hathewaya",
                                                     "Herbaspirillum",
                                                     "hgcI_clade",
                                                     "Hirschia",
                                                     "Holophaga",
                                                     "Hydrocarboniphaga",
                                                     "Hyphomicrobium",
                                                     "Ilumatobacter",
                                                     "IMCC26207",
                                                     "Inquilinus",
                                                     "Janthinobacterium",
                                                     "JTB23",
                                                     "Kocuria",
                                                     "Lactobacillus",
                                                     "Lacunisphaera",
                                                     "Latescibacterota",
                                                     "Lautropia",
                                                     "Lawsonella",
                                                     "Leeia",
                                                     "Legionella",
                                                     "Lentisphaera",
                                                     "Leucobacter",
                                                     "Limnohabitans",
                                                     "Lineage_IIb",
                                                     "Longivirga",
                                                     "Luminiphilus",
                                                     "Luteibacter",
                                                     "LWQ8",
                                                     "Lysinibacillus",
                                                     "Marine_Methylotrophic_Group_3",
                                                     "Marinimicrobia_(SAR406_clade)",
                                                     "Marinomonas",
                                                     "Massilia",
                                                     "Mesorhizobium",
                                                     "Methylobacterium-Methylorubrum",
                                                     "Methylophilus",
                                                     "Methylotenera",
                                                     "Microbacteriaceae",
                                                     "Microbacterium",
                                                     "Microbulbifer",
                                                     "Micrococcus",
                                                     "MWH-UniP1_aquatic_group",
                                                     "Mycobacterium",
                                                     "Nannocystis",
                                                     "NB1-j",
                                                     "Neisseria",
                                                     "Neorickettsia",
                                                     "Neptuniibacter",
                                                     "Nocardioides",
                                                     "Nodularia_PCC-9350",
                                                     "Novosphingobium",
                                                     "Oceanirhabdus",
                                                     "Oceanospirillum",
                                                     "Ochrobactrum",
                                                     "Oligoflexus",
                                                     "OM182_clade",
                                                     "OM190",
                                                     "OM27_clade",
                                                     "OM43_clade",
                                                     "OM60(NOR5)_clade",
                                                     "Opitutus",
                                                     "P3OB-42",
                                                     "Paenibacillus",
                                                     "Paraclostridium",
                                                     "Paracoccus",
                                                     "Parahaliea",
                                                     "PB19",
                                                     "Pectobacterium",
                                                     "Pelagicoccus",
                                                     "Pelomonas",
                                                     "PeM15",
                                                     "Peredibacter",
                                                     "Photobacterium",
                                                     "Phreatobacter",
                                                     "Phycisphaera",
                                                     "pItb-vmat-80",
                                                     "Pla3_lineage",
                                                     "Polynucleobacter",
                                                     "Porphyrobacter",
                                                     "possible_genus_04",
                                                     "Prochlorococcus_MIT9313",
                                                     "Propionibacterium",
                                                     "Propionigenium",
                                                     "Prosthecobacter",
                                                     "Proteus",
                                                     "Pseudarcicella",
                                                     "Pseudoalteromonas",
                                                     "Pseudobacteriovorax",
                                                     "Pseudomonas",
                                                     "Pseudorhodobacter",
                                                     "Psychrilyobacter",
                                                     "Psychrobacter",
                                                     "R76-B128",
                                                     "R7C24",
                                                     "Ralstonia",
                                                     "Reyranella",
                                                     "Rheinheimera",
                                                     "Rhizobacter",
                                                     "Rhodococcus",
                                                     "Rhodoferax",
                                                     "Rhodoluna",
                                                     "Richelia_HH01",
                                                     "Rickettsiella",
                                                     "Romboutsia",
                                                     "Rothia",
                                                     "Rubritalea",
                                                     "Sandaracinus",
                                                     "SAR116_clade",
                                                     "SAR324_clade(Marine_group_B)",
                                                     "SAR86_clade",
                                                     "SCGC_AAA164-E04",
                                                     "Sericytochromatia",
                                                     "Shewanella",
                                                     "SM1A02",
                                                     "Sorangium",
                                                     "Sphaerotilus",
                                                     "Sphingobium",
                                                     "Sphingopyxis",
                                                     "Sphingorhabdus",
                                                     "Sporanaerobacter",
                                                     "Sporichthyaceae",
                                                     "Sporosalibacterium",
                                                     "Staphylococcus",
                                                     "Stenotrophobacter",
                                                     "Stenotrophomonas",
                                                     "Streptococcus",
                                                     "Subgroup_17",
                                                     "Subgroup_22",
                                                     "Sva0081_sediment_group",
                                                     "Sva1033",
                                                     "Synechococcus_CC9902",
                                                     "Tepidibacter",
                                                     "Tepidimonas",
                                                     "Tepidiphilus",
                                                     "Terrimicrobium",
                                                     "Thermicanus",
                                                     "Timonella",
                                                     "Trichodesmium_IMS101",
                                                     "UBA10353_marine_group",
                                                     "Uliginosibacterium",
                                                     "Unassigned",
                                                     "uncultured",
                                                     "Undibacterium",
                                                     "Urania-1B-19_marine_sediment_group",
                                                     "Vallitalea",
                                                     "Variovorax",
                                                     "Vibrio",
                                                     "Vicinamibacteraceae",
                                                     "Weissella",
                                                     "WN-HWB-116",
                                                     "Woeseia",
                                                     "Xanthobacteraceae",
                                                     "Xanthomonas","Sphingomonas"))



spatial_plot <- ggplot(data=data_glom, aes(x=Sample, y=Abundance, fill=Species)) + facet_grid(~Tow, scales = "free")



#RANDOM COLOR PER OGNI TAXA
spatial_plot + geom_bar(aes(), stat="identity", position="stack",width = 0.9) +
  scale_fill_manual(values=distinctColorPalette(length(unique(data_glom$Genus))), breaks=c("Acinetobacter",
                                                                                           "AEGEAN-169_marine_group",
                                                                                           "Bacillus",
                                                                                           "Blastopirellula",
                                                                                           "Burkholderia-Caballeronia-Paraburkholderia",
                                                                                           "CL500-3",
                                                                                           "Corynebacterium",
                                                                                           "Cutibacterium",
                                                                                           "Duganella",
                                                                                           "Methylotenera",
                                                                                           "Nannocystis",
                                                                                           "OM190",
                                                                                           "Photobacterium",
                                                                                           "Pseudomonas",
                                                                                           "SAR86_clade",
                                                                                           "Sphingomonas",
                                                                                           "Streptococcus",
                                                                                           "Synechococcus_CC9902",
                                                                                           "Unassigned",
                                                                                           "Vibrio","< 0.2% abund."))+
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))

#




#  scale_fill_manual(values\ = c("#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#646B63","#00A86B","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#646B63","#00A86B","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#646B63","#00A86B","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#646B63","#00A86B","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#646B63","#00A86B","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#646B63","#00A86B","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#646B63","#00A86B","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#646B63","#00A86B","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#646B63","#00A86B","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#646B63","#00A86B","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#646B63","#00A86B","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#646B63","#00A86B")) +


library(cowplot)
my_legend <- get_legend(p)
library(ggpubr)
as_ggplot(my_legend)




#####
#simple way to rename phyla with < 1% abundance
data_glom$Genus[data_glom$Abundance < 0.2] <- "< 0.2% abund."

#Count # phyla to set color palette
Count = length(unique(data_glom$Genus))
Count
sort(unique(data_glom$Genus))




##plot with condensed phyla into "unknown" category
spatial_plot <- ggplot(data=data_glom, aes(x=Sample, y=Abundance, fill=Genus)) + facet_grid(~Lat, scales = "free")

spatial_plot + geom_bar(aes(), stat="identity", position="stack",width = 0.7) +
  scale_fill_manual(values = c( "#708090","#7FFF00","#00A86B","#960018","#7B1B02",
                               "#4B0082","#884DA7","#E52B50","#293133","#7FFFD4",
                               "#708090","#960018","#7FFF00","#00A86B","#7B1B02",
                               "#4B0082","#000000","#B20000","#00FF00","#0F52BA",
                               "#FF9933","#ff42ad","#FFD700","#646B63",
                               "#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF",
                               "#FFBF00","#884DA7","#E52B50","#7FFFD4",
                               "#4B0082",
                               "#ffa526","#B20000","#00FF00","#0F52BA","#FF9933","#177245",
                               "#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430",
                               "#708090","#960018","#7FFF00","#00A86B","#ff0303","#4B0082",
                               "#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000",
                               "#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082",
                               "#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4",
                               "#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082",
                               "#045e08","#B20000","#00FF00","#0F52BA","#FF9933","#ff42ad",
                               "#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430",
                               "#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133",
                               "#7FFFD4","#708090","#750004","#960018","#7FFF00","#00A86B","#FFBF00",
                               "#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018",
                               "#7FFF00","#00A86B","#7B1B02","#4B0082","#884DA7","#E52B50",
                               "#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B",
                               "#7B1B02","#4B0082","#ffa526","#f2d6ff","#ff242b",
                               "#B20000","#00FF00","#0F52BA","#FF9933","#ff42ad","#FF0000",
                               "#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082",
                               "#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4",
                               "#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082",
                               "#750004","#B20000","#000000","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#ff42ad","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#00FF00","#0F52BA","#FF9933","#ff42ad","#FF0000","#FFD700","#646B63","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933","#ff42ad","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000","#646B63","#00A8)" ), breaks=c("Acinetobacter",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "AEGEAN-169_marine_group",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "Bacillus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "Blastopirellula",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "Burkholderia-Caballeronia-Paraburkholderia",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "CL500-3",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "Corynebacterium",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "Cutibacterium",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "Duganella",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "Methylotenera",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "Nannocystis",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "OM190",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "Photobacterium",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "Pseudomonas",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "SAR86_clade",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "Sphingomonas",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "Streptococcus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "Synechococcus_CC9902",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "Unassigned",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "uncultured",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "Vibrio")) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))

#,legend.key.size = unit(0.5,"line")

#  per genus

spatial_plot + geom_bar(aes(), stat="identity", position="stack",width = 0.9) +
  scale_fill_manual(values = c("#7FFFD4","#708090","#7B1B02",
                               "#4B0082","#00FF00","#0F52BA","#FF9933",
                               "#ff42ad","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F",
                               "#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7",
                               "#E52B50","#293133","#7FFFD4","#708090", "#4B0082","#884DA7","#E52B50",
                               "#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#d4cce0",
                               "#B20000","#00FF00","#0F52BA","#FF9933","#ff42ad","#FF0000",
                               "#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082",
                               "#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4",
                               "#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082",
                               "#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245",
                               "#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430",
                               "#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133",
                               "#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02",
                               "#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933",
                               "#ff42ad","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F",
                               "#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50",
                               "#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B",
                               "#7B1B02","#4B0082","#000000","#000000","#B20000","#00FF00",
                               "#0F52BA","#FF9933","#ff42ad","#FF0000","#FFD700","#646B63",
                               "#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00",
                               "#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018",
                               "#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000",
                               "#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700",
                               "#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF",
                               "#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090",
                               "#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000",
                               "#B20000","#00FF00","#0F52BA","#FF9933","#ff42ad","#FF0000",
                               "#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082",
                               "#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4",
                               "#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082",
                               "#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B",
                               "#7B1B02","#4B0082","#000000","#B20000","#00FF00","#0F52BA",
                               "#FF9933","#177245","#FF0000","#FFD700","#646B63","#00A86B",
                               "#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7",
                               "#E52B50","#293133","#7FFFD4","#708090","#960018","#7FFF00",
                               "#00A86B","#7B1B02","#4B0082","#000000","#B20000","#00FF00",
                               "#0F52BA","#FF9933","#ff42ad","#FF0000","#FFD700","#646B63",
                               "#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00",
                               "#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018",
                               "#7FFF00","#00A86B","#FFBF00","#884DA7","#E52B50","#293133",
                               "#7FFFD4","#708090","#960018","#7FFF00","#00A86B",
                               "#4B0082","#884DA7","#E52B50","#293133","#7FFFD4","#708090",
                               "#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000",
                               "#B20000","#00FF00","#0F52BA","#FF9933","#ff42ad","#FF0000",
                               "#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082",
                               "#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4",
                               "#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082",
                               "#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245",
                               "#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430",
                               "#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133",
                               "#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02",
                               "#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933",
                               "#ff42ad","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F",
                               "#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50",
                               "#4B0082","#884DA7","#E52B50","#293133","#7FFFD4","#708090",
                               "#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000",
                               "#B20000","#00FF00","#0F52BA","#FF9933","#ff42ad","#FF0000",
                               "#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082",
                               "#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4",
                               "#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082",
                               "#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245",
                               "#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430",
                               "#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133",
                               "#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02",
                               "#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933",
                               "#ff42ad","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F",
                               "#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#003333",
                               "#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B",
                               "#7B1B02","#4B0082","#000000","#000000","#B20000","#00FF00",
                               "#0F52BA","#FF9933","#ff42ad","#FF0000","#FFD700","#646B63",
                               "#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00",
                               "#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018",
                               "#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000",
                               "#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700",
                               "#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF",
                               "#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090",
                               "#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000",
                               "#B20000","#00FF00","#0F52BA","#FF9933","#ff42ad","#FF0000",
                               "#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082",
                               "#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4",
                               "#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082",
                               "#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B",
                               "#7B1B02","#4B0082","#000000","#000000","#B20000","#00FF00",
                               "#0F52BA","#FF9933","#ff42ad","#FF0000","#FFD700","#646B63",
                               "#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00",
                               "#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018",
                               "#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000",
                               "#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700",
                               "#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF",
                               "#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090",
                               "#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000",
                               "#4B0082","#884DA7","#E52B50","#293133","#7FFFD4","#708090",
                               "#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000",
                               "#B20000","#00FF00","#0F52BA","#FF9933","#ff42ad","#FF0000",
                               "#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082",
                               "#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4",
                               "#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082",
                               "#000000","#B20000","#00FF00","#0F52BA","#FF9933","#177245",
                               "#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430",
                               "#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50","#293133",
                               "#7FFFD4","#708090","#960018","#7FFF00","#00A86B","#7B1B02",
                               "#4B0082","#000000","#B20000","#00FF00","#0F52BA","#FF9933",
                               "#ff42ad","#FF0000","#FFD700","#646B63","#00A86B","#3D2B1F",
                               "#F4C430","#4B0082","#007FFF","#FFBF00","#884DA7","#E52B50",
                               "#293133","#7FFFD4","#708090","#960018","#7FFF00","#00A86B",
                               "#7B1B02","#4B0082","#000000","#000000","#B20000","#00FF00",
                               "#0F52BA","#FF9933","#ff42ad","#FF0000","#FFD700","#646B63",
                               "#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF","#FFBF00",
                               "#884DA7","#E52B50","#293133","#7FFFD4","#708090","#960018",
                               "#7FFF00","#00A86B","#7B1B02","#4B0082","#000000","#B20000",
                               "#00FF00","#0F52BA","#FF9933","#177245","#FF0000","#FFD700",
                               "#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#007FFF",
                               "#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4","#708090",
                               "#960018","#7FFF00","#00A86B","#7B1B02","#4B0082","#000000",
                               "#B20000","#00FF00","#0F52BA","#FF9933","#ff42ad","#FF0000",
                               "#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082",
                               "#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4",
                               "#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082",
                               "#B20000","#00FF00","#0F52BA","#FF9933","#ff42ad","#FF0000",
                               "#FFD700","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082",
                               "#007FFF","#FFBF00","#884DA7","#E52B50","#293133","#7FFFD4",
                               "#708090","#960018","#7FFF00","#00A86B","#7B1B02","#4B0082",
                               "#000000","#B20000","#646B63","#00A86B"), breaks=c("Acinetobacter",
                                                                                  "AEGEAN-169_marine_group",
                                                                                  "Bacillus",
                                                                                  "Blastopirellula",
                                                                                  "Burkholderia-Caballeronia-Paraburkholderia",
                                                                                  "CL500-3",
                                                                                  "Corynebacterium",
                                                                                  "Cutibacterium",
                                                                                  "Duganella",
                                                                                  "Methylotenera",
                                                                                  "Nannocystis",
                                                                                  "OM190",
                                                                                  "Photobacterium",
                                                                                  "Pseudomonas",
                                                                                  "SAR86_clade",
                                                                                  "Sphingomonas",
                                                                                  "Streptococcus",
                                                                                  "Synechococcus_CC9902",
                                                                                  "Unassigned",
                                                                                  "Vibrio")) +
  theme(legend.position="bottom") +  guides(fill=guide_legend(nrow=5)) 

  
library(cowplot)
library(ggpubr)
my_legend <- get_legend(p)

as_ggplot(my_legend)



# Using the cowplot package
legend <- cowplot::get_legend(p)

grid.newpage()
grid.draw(legend)

## genus più abbondanti da plottare
c("AEGEAN-169_marine_group",
"Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
"Anaerococcus",
"Bacillus",
"BD7-11",
"Blastopirellula",
"Bradyrhizobium",
"Candidatus_Omnitrophus",
"CL500-3",
"Lineage_IIb",
"OM190",
"OM27_clade",
"P3OB-42",
"Polynucleobacter",
"Propionigenium",
"Ralstonia",
"SAR116_clade",
"SM1A02",
"SM2D12",
"Sphingobium",
"Sporanaerobacter",
"Unassigned",
"uncultured",
"Urania-1B-19_marine_sediment_group",
"vadinHA49",
"Xanthobacteraceae")




##OTU97


c("Acinetobacter",
"AEGEAN-169_marine_group",
"Bacillus",
"Blastopirellula",
"Burkholderia-Caballeronia-Paraburkholderia",
"CL500-3",
"Corynebacterium",
"Cutibacterium",
"Duganella",
"Methylotenera",
"Nannocystis",
"OM190",
"Photobacterium",
"Pseudomonas",
"SAR86_clade",
"Sphingomonas",
"Streptococcus",
"Synechococcus_CC9902",
"Unassigned",
"Vibrio")









##simple way to rename phyla with < 2% abundance
data_glom$phylum[data_glom$Abundance < 0.02] <- "< 2% abund"

#Count # phyla to set color palette
Count = length(unique(data_glom$phylum))
Count
unique(data_glom$phylum)
#data_glom$phylum <- factor(data_glom$phylum, levels = c(" p__Proteobacteria"," p__Spirochaetes" ," p__Acidobacteria"," p__Actinobacteria"," p__Firmicutes" , " p__WS3"," p__Chloroflexi"," p__Planctomycetes"," p__Bacteroidetes"," p__Gemmatimonadetes"," p__GN04"," p__Caldithrix"," p__OP1"," p__Verrucomicrobia"," p__" , "< 2% abund", "Unassigned"))
##plot with condensed phyla into "unknown" category
spatial_plot <- ggplot(data=data_glom, aes(x=Sample, y=Abundance, fill=phylum)) #+ facet_grid(~site, scales = "free")

spatial_plot + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = colors_bar_plot2)#+
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
#  scale_fill_manual(values = c( "#1900ff","#ff0000", "#ffb300", "#ff00f2", "#40ff00", "#00d9ff", "#121112", "#9e344d", "#5d731d", "#925abf", "#a9f28a", "#381938", "#007311", "#613314", "#708090", "#1477b5")) +


#"#007FFF","#FFBF00","#884DA7","#E52B50","#7FFFD4","#708090","#960018","#ff0000","#ff5500","#7B1B02","#4B0082","#000000","#B20000","#FF9933","#00FF00","#0F52BA","#177245","#FF0000","#646B63","#00A86B","#3D2B1F","#F4C430","#4B0082","#293133"










