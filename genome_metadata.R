#Gather metadata on published gastropod and intertidal species genomes
# install new packages
install.packages("taxize")
install.packages("readr")
install.packages("geiger")
install.packages("myTAI")
install.packages("usethis")
install.packages("plyr")
install.packages("rredlist")
install.packages("forcats")

##### loading dependencies #####
library(forcats)
library(ape)
library(dplyr)
library(taxize)
library(ggplot2)
library(readr)
library(geiger)
library(myTAI)
library(usethis)
library(plyr)
library(rredlist)

getwd()
setwd("/Users/emily/Dropbox/School/Thesis/Genomics-Ch1/genome_metadata")
ENTREZ_KEY='5ac727e2b17036171bfb2d56d5113649fb08'

##### Read in data generated using ncbi dataset genome search tool #####
mol_nome <- read_tsv("mollusca_genomes_02.tsv")
head(mol_nome)
colnames(mol_nome) <- c("Assembly_Accession", "Assembly_Name", "Annotation_Name", "Release_Date", "Organism_Name",
                         "BUSCO_C", "BUSCO_D", "BUSCO_F", "BUSCO_Lineage", "Prot_coding_count", "Gene_Count", "Level",
                         "Date","Seq_Tech","GC","Scaf_L50","Scaf_N50","Chrm_No")
##### Remove entries with genome level equal to contig
mollusca1<-mol_nome[!(mol_nome$Level=="Contig"),]
head(mollusca1)
summary(mollusca1)

##### Subset the dataframe in batches of 10 species b/c otherwise taxize crashes... not sure why??!!! #####
n <- 10
split_df<-split(mollusca1, rep(1:(nrow(mollusca1)/n+nrow(mollusca1)%%n),each=n))
for (i in 1:length(split_df)) {
  assign(paste0("split_df", i), as.data.frame(split_df[[i]]))
}
#You now have 28 dataframes of no more than 10 species each. They are called "split_df#"
#######Retrieve taxonomic information for the species listed in mollusca1$Organism_Name ####
df_list_split1 <- list()
for (i in 1:nrow(split_df1)){
  tax  <- taxonomy(organism = split_df1[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split1[[i]] <- df
}
df_list_split2 <- list()
for (i in 1:nrow(split_df2)){
  tax  <- taxonomy(organism = split_df2[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split2[[i]] <- df
}
df_list_split3 <- list()
for (i in 1:nrow(split_df3)){
  tax  <- taxonomy(organism = split_df3[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split3[[i]] <- df
}
df_list_split4 <- list()
for (i in 1:nrow(split_df4)){
  tax  <- taxonomy(organism = split_df4[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split4[[i]] <- df
}
df_list_split5 <- list()
for (i in 1:nrow(split_df5)){
  tax  <- taxonomy(organism = split_df5[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split5[[i]] <- df
}
df_list_split6 <- list()
for (i in 1:nrow(split_df6)){
  tax  <- taxonomy(organism = split_df6[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split6[[i]] <- df
}
df_list_split7 <- list()
for (i in 1:nrow(split_df7)){
  tax  <- taxonomy(organism = split_df7[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split7[[i]] <- df
}
df_list_split8 <- list()
for (i in 1:nrow(split_df8)){
  tax  <- taxonomy(organism = split_df8[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split8[[i]] <- df
}
df_list_split9 <- list()
for (i in 1:nrow(split_df9)){
  tax  <- taxonomy(organism = split_df9[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split9[[i]] <- df
}
df_list_split10 <- list()
for (i in 1:nrow(split_df10)){
  tax  <- taxonomy(organism = split_df10[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split10[[i]] <- df
}
df_list_split11 <- list()
for (i in 1:nrow(split_df11)){
  tax  <- taxonomy(organism = split_df11[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split11[[i]] <- df
}
df_list_split12 <- list()
for (i in 1:nrow(split_df12)){
  tax  <- taxonomy(organism = split_df12[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split12[[i]] <- df
}
df_list_split13 <- list()
for (i in 1:nrow(split_df13)){
  tax  <- taxonomy(organism = split_df13[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split13[[i]] <- df
}
df_list_split14 <- list()
for (i in 1:nrow(split_df14)){
  tax  <- taxonomy(organism = split_df14[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split14[[i]] <- df
}
df_list_split15 <- list()
for (i in 1:nrow(split_df15)){
  tax  <- taxonomy(organism = split_df15[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split15[[i]] <- df
}
df_list_split16 <- list()
for (i in 1:nrow(split_df16)){
  tax  <- taxonomy(organism = split_df16[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split16[[i]] <- df
}
df_list_split17 <- list()
for (i in 1:nrow(split_df17)){
  tax  <- taxonomy(organism = split_df17[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split17[[i]] <- df
}
df_list_split18 <- list()
for (i in 1:nrow(split_df18)){
  tax  <- taxonomy(organism = split_df18[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split18[[i]] <- df
}
df_list_split19 <- list()
for (i in 1:nrow(split_df19)){
  tax  <- taxonomy(organism = split_df19[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split19[[i]] <- df
}
df_list_split20 <- list()
for (i in 1:nrow(split_df20)){
  tax  <- taxonomy(organism = split_df20[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split20[[i]] <- df
}
df_list_split21 <- list()
for (i in 1:nrow(split_df21)){
  tax  <- taxonomy(organism = split_df21[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split21[[i]] <- df
}
df_list_split22 <- list()
for (i in 1:nrow(split_df22)){
  tax  <- taxonomy(organism = split_df22[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split22[[i]] <- df
}
df_list_split23 <- list()
for (i in 1:nrow(split_df23)){
  tax  <- taxonomy(organism = split_df23[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split23[[i]] <- df
}
df_list_split24 <- list()
for (i in 1:nrow(split_df24)){
  tax  <- taxonomy(organism = split_df24[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split24[[i]] <- df
}
df_list_split25 <- list()
for (i in 1:nrow(split_df25)){
  tax  <- taxonomy(organism = split_df25[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split25[[i]] <- df
}
df_list_split26 <- list()
for (i in 1:nrow(split_df26)){
  tax  <- taxonomy(organism = split_df26[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split26[[i]] <- df
}
df_list_split27 <- list()
for (i in 1:nrow(split_df27)){
  tax  <- taxonomy(organism = split_df27[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split27[[i]] <- df
}
df_list_split28 <- list()
for (i in 1:nrow(split_df28)){
  tax  <- taxonomy(organism = split_df28[i,]$Organism_Name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list_split28[[i]] <- df
}

####Make results into dataframes "combined_df_split#" #####
combined_df_split1 <- do.call(rbind.fill, df_list_split1)
combined_df_split2 <- do.call(rbind.fill, df_list_split2)
combined_df_split3 <- do.call(rbind.fill, df_list_split3)
combined_df_split4 <- do.call(rbind.fill, df_list_split4)
combined_df_split5 <- do.call(rbind.fill, df_list_split5)
combined_df_split6 <- do.call(rbind.fill, df_list_split6)
combined_df_split7 <- do.call(rbind.fill, df_list_split7)
combined_df_split8 <- do.call(rbind.fill, df_list_split8)
combined_df_split9 <- do.call(rbind.fill, df_list_split9)
combined_df_split10 <- do.call(rbind.fill, df_list_split10)
combined_df_split11 <- do.call(rbind.fill, df_list_split11)
combined_df_split12 <- do.call(rbind.fill, df_list_split12)
combined_df_split13 <- do.call(rbind.fill, df_list_split13)
combined_df_split14 <- do.call(rbind.fill, df_list_split14)
combined_df_split15 <- do.call(rbind.fill, df_list_split15)
combined_df_split16 <- do.call(rbind.fill, df_list_split16)
combined_df_split17 <- do.call(rbind.fill, df_list_split17)
combined_df_split18 <- do.call(rbind.fill, df_list_split18)
combined_df_split19 <- do.call(rbind.fill, df_list_split19)
combined_df_split20 <- do.call(rbind.fill, df_list_split20)
combined_df_split21 <- do.call(rbind.fill, df_list_split21)
combined_df_split22 <- do.call(rbind.fill, df_list_split22)
combined_df_split23 <- do.call(rbind.fill, df_list_split23)
combined_df_split24 <- do.call(rbind.fill, df_list_split24)
combined_df_split25 <- do.call(rbind.fill, df_list_split25)
combined_df_split26 <- do.call(rbind.fill, df_list_split26)
combined_df_split27 <- do.call(rbind.fill, df_list_split27)
combined_df_split28 <- do.call(rbind.fill, df_list_split28)

#####MERGE the dataframes created in the last step ####
merged_df_tax <- dplyr::bind_rows(combined_df_split1, 
                       combined_df_split2, 
                       combined_df_split3,
                       combined_df_split4,
                       combined_df_split5,
                       combined_df_split6,
                       combined_df_split7,
                       combined_df_split8,
                       combined_df_split9,
                       combined_df_split10,
                       combined_df_split11,
                       combined_df_split12,
                       combined_df_split13,
                       combined_df_split14,
                       combined_df_split15,
                       combined_df_split16,
                       combined_df_split17,
                       combined_df_split18,
                       combined_df_split19,
                       combined_df_split20,
                       combined_df_split21,
                       combined_df_split22,
                       combined_df_split23,
                       combined_df_split24,
                       combined_df_split25,
                       combined_df_split26,
                       combined_df_split27,
                       combined_df_split28)
head(merged_df_tax, n=11)

names(merged_df_tax)[names(merged_df_tax) == 'species'] <- 'Organism_Name'

#####MERGE that with mollusca1 #####
head(mollusca1)
merged_gnomes_tax <- dplyr::full_join(mollusca1, merged_df_tax, by="Organism_Name")
head(merged_gnomes_tax)
merged_gnomes_tax_uniq <- merged_gnomes_tax %>% distinct()
head(merged_gnomes_tax_uniq)

merged_gnomes_tax_uniq<-merged_gnomes_tax_uniq[!(merged_gnomes_tax_uniq$Organism_Name=="Argopecten irradians concentricus"|merged_gnomes_tax_uniq$Organism_Name=="Argopecten irradians irradians"),]

##### Add habitat ####
## automatic habitat retrieval does not work, doing this step by hand via the literature!!!
getwd()
habitat_small_list <- read_tsv("habitat.tsv")
head(habitat_small_list)

merged_gnomes_tax_uniq_hab <- dplyr::left_join(merged_gnomes_tax_uniq, habitat_small_list, by="Organism_Name")
head(merged_gnomes_tax_uniq_hab)
merged_gnomes_tax_uniq_hab <- merged_gnomes_tax_uniq_hab %>% distinct()

merged_gnomes_tax_uniq_hab_partial <- read_tsv("merged_gnomes_tax_uniq_hab.tsv")

##### Write table #######
write.table(merged_gnomes_tax_uniq_hab, file='/Users/emily/Dropbox/School/Thesis/Genomics-Ch1/genome_metadata/merged_gnomes_tax_uniq_hab.txt', row.names=FALSE, sep='\t', quote=FALSE)



##### Metadata ready, Analysis, Graphs ####
df1 <- read_tsv("/Users/emily/Dropbox/School/Thesis/Genomics-Ch1/genome_metadata/merged_gnomes_tax_uniq_hab_noscurria.tsv") 
head(df1)
str(df1)
summary(df1)

df1_scurrias <- read_tsv("/Users/emily/Dropbox/School/Thesis/Genomics-Ch1/genome_metadata/merged_gnomes_tax_uniq_hab.tsv") 

df1_reducedbysp <- read_tsv("/Users/emily/Dropbox/School/Thesis/Genomics-Ch1/genome_metadata/merged_gnomes_tax_uniq_hab_noscurrias_reducedbysp.tsv") 

#Count number of species
n_distinct(df1$Organism_Name)
n_distinct(df1_reducedbysp$Organism_Name)
#195

#Count number of assemblies with annotations
x<- df1 %>% group_by(Organism_Name, class, subclass) %>%
  summarise(total_non_na = sum(!is.na(Annotation_Name)))
y<- x[x$total_non_na != 0, ]
#40 assemblies with annotations

x<- df1_reducedbysp %>% group_by(Organism_Name, class, subclass) %>%
  summarise(total_non_na = sum(!is.na(Annotation_Name)))
y<- x[x$total_non_na != 0, ]
#30 assemblies with annotations

#Back to original DF
df2 <- df1 %>% 
  dplyr::count(Assembly_Accession, Annotation_Name,class,subclass,order,family,genus, Organism_Name, Habitat_1,Habitat_2,Habitat_3)
head(df2)
summary(df2)
df2$class
str(df2)

df2_reducedbysp <- df1_reducedbysp %>% 
  dplyr::count(Assembly_Accession, Annotation_Name,class,subclass,order,family,genus, Organism_Name, Habitat_1,Habitat_2,Habitat_3)
head(df2_reducedbysp)

# Remove duplicates based on Organism_Name 
df2_species<- df2[!duplicated(df2$Organism_Name), ]
str(df2_species)

##### Plot data #####
dev.off()
str(df1_scurrias)
df1_scurrias %>%
  ggplot(aes(class, Chrm_No, label=ifelse(df1_scurrias$subclass=='Patellogastropoda',genus," "))) + 
  geom_boxplot(outline=FALSE)+
  geom_jitter(alpha=0.6, size=4, color=df1_scurrias$color, position = position_jitter(seed = 4)) +
  geom_text(position = position_jitter(seed = 4), size=2, angle = 45, hjust=.75) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Chromosome Number") + 
  theme(panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10))

df1_scurrias$color <- ifelse(df1_scurrias$subclass=="Caenogastropoda", '#008B00',
                    ifelse(df1_scurrias$subclass=="Heterobranchia", '#8B8B00',
                           ifelse(df1_scurrias$subclass=='Neomphalina', '#008B8B',
                                  ifelse(df1_scurrias$subclass=='Patellogastropoda','hotpink',
                                         ifelse(df1_scurrias$subclass=='Vetigastropoda', '#8B008B', 'grey')))))
str(df1_scurrias)
df1_scurrias %>%
  ggplot(aes(class, Chrm_No, label=ifelse(df1_scurrias$subclass=='Patellogastropoda',genus," "))) + 
  geom_boxplot(outline=FALSE)+
  geom_jitter(alpha=0.6, size=4, color=df1_scurrias$color, position = position_jitter(seed = 4)) +
  geom_text(position = position_jitter(seed = 4), size=2, angle = 45, hjust=.75) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Chromosome Number") + 
  theme(panel.background=element_blank(),
    axis.line=element_line(colour="black"),
    axis.text=element_text(size=10),
    axis.title=element_text(size=10))
  

df1 %>%
ggplot(aes(class, Chrm_No, label=Organism_Name)) + 
  geom_point(alpha=0.5) +
  geom_point(data = df1 %>% filter(subclass == "Caenogastropoda"), color = "#008B00") +
  geom_point(data = df1 %>% filter(subclass == "Heterobranchia"), color = "#8B8B00") +
  geom_point(data = df1 %>% filter(subclass == "Neomphalina"), color = "#008B8B") +
  geom_point(data = df1 %>% filter(subclass == "Patellogastropoda"), color = 'hotpink') +
  geom_point(data = df1 %>% filter(subclass == "Vetigastropoda"), color = '#8B008B') +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Chromosome Number") + 
  theme(#legend.title=element_blank(), 
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        #axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=1),
        #axis.text.x=element_blank(),
        #legend.key.size = unit(0.5, 'cm'),
        legend.position="bottom",
        legend.key = element_blank())

ggplot(y2, aes(x = fct_relevel(Species, "Lottia gigantea", "Patella vulgata", 
                               "Haliotis rubra", "Haliotis rufescens", "Gigantopelta aegis","Pomacea canaliculata",
                               "Elysia marginata", "Elysia chlorotica","Plakobranchus ocellatus"), y=Annotations, fill=ifelse(Class=='Gastropoda', Subclass, 'other'))) + 
  geom_bar(stat="identity", width = 0.5, alpha=0.5, position="stack", color="black", show.legend = TRUE) +
  coord_flip() +
  scale_fill_manual(values=c('#008B00','#8B8B00','#008B8B','grey','hotpink','#8B008B')) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Number of Annotations in NCBI", limits=c(0, 2), breaks=c(0,1,2)) + 
  theme(legend.title=element_blank(), 
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=1),
        #axis.text.x=element_blank(),
        #legend.key.size = unit(0.5, 'cm'),
        legend.position="bottom",
        legend.key = element_blank())

ggplot(x, aes(x = Organism_Name, y=total_non_na, fill=subclass)) + 
  geom_bar(stat="identity", width = 0.5, alpha=0.5, position="stack", color="black", show.legend = TRUE) +
  coord_flip() +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Number of Assemblies in NCBI", limits=c(0, 2), breaks=c(0,1,2)) + 
  theme(legend.title=element_blank(), 
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        #axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=1),
        #axis.text.x=element_blank(),
        #legend.key.size = unit(0.5, 'cm'),
        legend.position="bottom",
        legend.key = element_blank())

ggplot(df2_reducedbysp, aes(x=reorder(Habitat_1, n), y=n, fill=(ifelse(Annotation_Name=='NA', "not_annotated", "annotated")))) + 
  geom_bar(stat="identity", width = 0.5, alpha=0.5, position="stack", color="black", show.legend = FALSE) +
  coord_flip() +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Species with Genome Assembly in NCBI", limits=c(0,60), breaks=c(0,20,40,60)) + 
  theme(legend.title=element_blank(), 
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=1),
        #axis.text.x=element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.position="bottom",
        legend.key = element_blank())

ggplot(df2_reducedbysp, aes(x=reorder(Habitat_3, n), y=n, fill=(ifelse(Annotation_Name=='NA', "not_annotated", "annotated")))) + 
  geom_bar(stat="identity", width = 0.5, alpha=0.5, position="stack", color="black", show.legend = FALSE) +
  coord_flip() +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Species with Genome Assembly in NCBI", limits=c(0,30), breaks=c(0,10,20,30)) + 
  theme(legend.title=element_blank(), 
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=1),
        #axis.text.x=element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.position="bottom",
        legend.key = element_blank())

ggplot(df2_species, aes(fill=(ifelse(class=='Gastropoda', subclass, 'other')), x = fct_relevel(class, "Gastropoda"), y=n)) + 
  geom_bar(stat="identity", width = 0.5, alpha=0.5, position="stack", color="black") +
  coord_flip() +
  scale_fill_manual(values=c('#008B00','#8B8B00','#008B8B','grey','hotpink','#8B008B')) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Species with Genome Assembly in NCBI", limits=c(1,100), breaks=c(0, 25, 50, 75, 100)) + 
  theme(legend.title=element_blank(), 
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=1),
        #axis.text.x=element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.position="bottom")


df2<-df1[!(df1$class=="Bivalvia"|df1$class=="Cephalopoda"|df1$class=="Polyplacophora"|df1$class=="Solenogastres"|df1$class=="Scaphopoda"),]
head(df2)

ggplot(df2, aes(fill=family, x=reorder(family, family), y=n)) + 
  geom_bar(stat="identity", width = 0.5, alpha=0.5, color="black") +
  coord_flip() +
  #scale_fill_manual(values=c('deeppink1','plum1','hotpink1','lightpink','maroon2')) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Number of Genome Assemblies in NCBI") + 
  theme(legend.title=element_blank(), 
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=1),
        #axis.text.x=element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.position="side")

head(merged_gnomes_tax_uniq_hab_partial)
ggplot(merged_gnomes_tax_uniq_hab_partial) +
  aes(x =class, y = (BUSCO_C)*100, color= (ifelse(genus=='Scurria', merged_gnomes_tax_uniq_hab_partial$Organism_Name...5, "other"))) +
  geom_point(show.legend = FALSE, aes(), size=5, alpha=0.5, pch=19) +
  scale_color_manual(values=c("darkgrey", "#FB61D7", "#A58AFF", "#00B6EB"))+
  #stat_summary(fun.data=mean_sdl, 
               #geom="pointrange", color="black", pch=16, size=.3) +
  #stat_summary(fun.y=median, geom="point", size=3, color="black", pch=17) +
  scale_x_discrete(name=" ") +
  scale_y_continuous(name="% Complete BUSCOs Mollusca_odb10") + 
  theme(legend.title=element_blank(),
        axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10))

head(merged_gnomes_tax_uniq_hab_partial)
ggplot(merged_gnomes_tax_uniq_hab_partial) +
  aes(x =class, y = (BUSCO_D)*100, color= (ifelse(genus=='Scurria', merged_gnomes_tax_uniq_hab_partial$Organism_Name...5, "other"))) +
  geom_point(show.legend = FALSE, aes(), size=5, alpha=0.5, pch=19) +
  scale_color_manual(values=c("darkgrey", "#FB61D7", "#A58AFF", "#00B6EB"))+
  #stat_summary(fun.data=mean_sdl, 
  #geom="pointrange", color="black", pch=16, size=.3) +
  #stat_summary(fun.y=median, geom="point", size=3, color="black", pch=17) +
  scale_x_discrete(name=" ") +
  scale_y_continuous(name="% Duplicated BUSCOs Mollusca_odb10") + 
  theme(legend.title=element_blank(),
        axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10))

head(merged_gnomes_tax_uniq_hab_partial)
ggplot(merged_gnomes_tax_uniq_hab_partial) +
  aes(x =class, y = (BUSCO_F)*100, color= (ifelse(genus=='Scurria', merged_gnomes_tax_uniq_hab_partial$Organism_Name...5, "other"))) +
  geom_point(show.legend = FALSE, aes(), size=5, alpha=0.5, pch=19) +
  scale_color_manual(values=c("darkgrey", "#FB61D7", "#A58AFF", "#00B6EB"))+
  #stat_summary(fun.data=mean_sdl, 
  #geom="pointrange", color="black", pch=16, size=.3) +
  #stat_summary(fun.y=median, geom="point", size=3, color="black", pch=17) +
  scale_x_discrete(name=" ") +
  scale_y_continuous(name="% Fragmented BUSCOs Mollusca_odb10") + 
  theme(legend.title=element_blank(),
        axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10))

merged_gnomes_tax_uniq_hab_partial2<-merged_gnomes_tax_uniq_hab_partial[!(merged_gnomes_tax_uniq_hab_partial$subclass=="NA"),]
ggplot(merged_gnomes_tax_uniq_hab_partial2) +
  aes(x =class, y = (Scaf_N50)/1000000, color= (ifelse(genus=='Scurria', merged_gnomes_tax_uniq_hab_partial$Organism_Name...5, "other"))) +
  geom_point(show.legend = FALSE, aes(), size=5, alpha=0.5, pch=19) +
  scale_color_manual(values=c("darkgrey", "#FB61D7", "#A58AFF", "#00B6EB"))+
  stat_summary(fun.data=mean_sdl, 
  geom="pointrange", color="black", pch=16, size=.3) +
  #stat_summary(fun.y=median, geom="point", size=3, color="black", pch=17) +
  scale_x_discrete(name=" ") +
  scale_y_continuous(name="Scaffold N50 (Mb)", limits=c(0, 170)) + 
  theme(legend.title=element_blank(),
        axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10))

ggplot(merged_gnomes_tax_uniq_hab_partial2) +
  aes(x =class, y = Scaf_L50, color= (ifelse(genus=='Scurria', merged_gnomes_tax_uniq_hab_partial$Organism_Name...5, "other"))) +
  geom_point(show.legend = FALSE, aes(), size=5, alpha=0.5, pch=19) +
  scale_color_manual(values=c("darkgrey", "#FB61D7", "#A58AFF", "#00B6EB"))+
  stat_summary(fun.data=mean_sdl, 
               geom="pointrange", color="black", pch=16, size=.3) +
  #stat_summary(fun.y=median, geom="point", size=3, color="black", pch=17) +
  scale_x_discrete(name=" ") +
  scale_y_continuous(name="Scaffold L50", limits=c(0, 1500)) + 
  theme(legend.title=element_blank(),
        axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10))
### Add tree -------------
phylo.Hedges <- read.tree("~/Dropbox/School/Thesis/Genomics-Ch1/comparative/Hedges_tree.txt")

tr.sup2 <- treedata(phylo.Hedges, 
                        setNames(sup2$sp, sup2$sp))[[1]]


