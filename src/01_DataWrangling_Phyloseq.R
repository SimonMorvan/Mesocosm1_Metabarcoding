# Simon Morvan
# July 2025
# R version 4.4.2 - Pile of Leaves
# 
#                          Title : Data wrangling and phyloseq object 
# Description : Script to create the phyloseq object that will serve for the other scripts 

source(here::here("src", "libraries.R"))




#### 16S solid samples  #####

data_path <- here("data/16s_sed/")

# Sample data
meta <- read.csv(file = paste0(data_path,"/mapping_file.csv"), 
                 dec = ".", header = T, row.names = 1, sep = ";", 
                 comment.char = "") #load 196 in 8 vars

# Remove some samples with undefined dates
meta <- subset(meta,meta$Sampling_date!="undefined") 

## Convert character vectors to factors
meta[sapply(meta, is.character)] <- lapply(meta[sapply(meta, is.character)], 
                                           as.factor) 

# Create group variable (Sample type + Water type)
meta$Group <- paste(meta$Compartment,meta$Sample_type,meta$Water_type,sep="_")
meta$Group <- factor(meta$Group, levels = c("Sediments_No_plant_Artificial_OSPW","Sediments_No_plant_OSPW",
                                            "Sediments_Carex_OSPW","Rhizosphere_Carex_OSPW","Roots_Carex_OSPW"))

# Re-order levels in Time variable
meta$Time <- factor(meta$Time, levels = c("D-0","D8","D28","D42","D84","BlankpcrCES"))
# 192 samples, 8 categories



# Abundance data
raw_ab <- read.csv(file = paste0(data_path,"/feature_table_filtered_clean_no-undefined-samples.csv"), dec = ".", sep = ",", header = T, row.names = 1, comment.char = "")

# Erase the "X" added by R 
names(raw_ab) <- gsub("X","",names(raw_ab)) 

# Add ASV in front of the feature identification number 
rownames(raw_ab) <- paste0("ASV",rownames(raw_ab)) 
#5149 ASvs


# Simplify sample names in both the abundance data and sample data

# Abundance data 

# First rename mispelled Mesocosm
colnames(raw_ab) <- gsub("Mesosocm", "Mesocosm", colnames(raw_ab), perl = TRUE)
colnames(raw_ab)<- gsub("Mesococm", "Mesocosm", colnames(raw_ab), perl = TRUE)

# Second remove the numbers before "Mesocosm"
colnames(raw_ab) <- gsub(".*(?=Mesocosm)", "", colnames(raw_ab), perl = TRUE)

# Third change D.0 to D0 (will help in further steps)
colnames(raw_ab) <- gsub("D.0", "D0", colnames(raw_ab), perl = TRUE)


# Sample data

meta$alias <- gsub("-",".",meta$alias)
meta$alias <- gsub("Mesosocm", "Mesocosm", meta$alias, perl = TRUE)
meta$alias<- gsub("Mesococm", "Mesocosm", meta$alias, perl = TRUE)

# Second remove the numbers before mesocosm
meta$alias <- gsub(".*(?=Mesocosm)", "", meta$alias, perl = TRUE)

# Third change D.0 to D0 
meta$alias <- gsub("D.0", "D0", meta$alias, perl = TRUE)
meta$alias <- gsub("TEST", "D84", meta$alias, perl = TRUE)
meta$alias <- gsub("Racines","Racine",meta$alias)


rownames(meta) <- meta$alias 

# Check if the samples found in the 
common_samples <- intersect(rownames(meta),colnames(raw_ab))
unique_to_meta <- rownames(meta)[!(rownames(meta)%in%common_samples)]
unique_to_ab <-  colnames(raw_ab)[!(colnames(raw_ab)%in%common_samples)]
# D84 and TEST need to match


colnames(raw_ab) <- gsub("TEST","D84",colnames(raw_ab))
common_samples <- intersect(rownames(meta),colnames(raw_ab))
unique_to_meta <- rownames(meta)[!(rownames(meta)%in%common_samples)]
unique_to_ab <-  colnames(raw_ab)[!(colnames(raw_ab)%in%common_samples)]
# Racine sometimes has an S in colnames(raw_ab)

colnames(raw_ab) <- gsub("Racines","Racine",colnames(raw_ab))
common_samples <- intersect(rownames(meta),colnames(raw_ab))
unique_to_meta <- rownames(meta)[!(rownames(meta)%in%common_samples)]
unique_to_ab <-  colnames(raw_ab)[!(colnames(raw_ab)%in%common_samples)]
# Now it's good.
# There's a taxonomy column in the abundance data, we'll put it in it's own
# data.frame, and split the individual taxonomy levels in individual columns


# Taxonomy 

# Last column of the abundance table contains the taxonomy
taxonomy <- as.data.frame(raw_ab[,ncol(raw_ab)]) 
rownames(taxonomy) <- rownames((raw_ab))
colnames(taxonomy) <- c("Taxonomy")

# Now we can remove the taxonomy column from raw_ab
raw_ab <- raw_ab[,-ncol(raw_ab)]

# Separate the taxonomic levels in individual columns
tax_clean <- separate(taxonomy, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep = ";")
tax_clean <- tax_clean[,-which(names(tax_clean) %in% c("Species", "Strain"))]
# Removes species and strains (No information there)

tax_clean <- as.data.frame(tax_clean)
levels(as.factor(tax_clean$Domain))    # 4   Domain 
levels(as.factor(tax_clean$Phylum))    # 50  Phyla 
levels(as.factor(tax_clean$Class))     # 129 Classes
levels(as.factor(tax_clean$Order))     # 273 Orders
levels(as.factor(tax_clean$Family))    # 402 Families
levels(as.factor(tax_clean$Genus))     # 710 Genera
str(tax_clean)

# Remove all the k__, g__ etc from taxonomy 
tax_clean[] <- lapply(tax_clean, function(x) substr(x, 4, nchar(x)))

# Replace empty and NA cells with "NA" 
tax_clean[]<- lapply(tax_clean[], function(x) ifelse(x == "" | is.na(x), "NA", x))

# Remove Chloroplasts and Mitochondria
tax_chloro <- subset(tax_clean,tax_clean$Order=="Chloroplast") # 79 Choloplast sequence
tax_mito <- subset(tax_clean,tax_clean$Family=="Mitochondria") # 0 Mitochondria sequence
summary(as.factor(tax_clean$Domain)) # Remove Eukaryotes and Viridiplantea
rm(tax_chloro,tax_mito)
tax_clean <- subset(tax_clean,
                    tax_clean$Order!="Chloroplast"&
                         tax_clean$Family!="Mitochondria"&
                         tax_clean$Domain!="Eukaryota"&
                         tax_clean$Domain!="Viridiplantae")

## Some taxa have been classified at a taxonomy level but with the previous
## taxonomy name for instance BurkholderialesOR at the genus level...
## I want to get rid of those and replace them by NA as their taxonomy was not 
## assigned at these lower taxonomy levels.

pattern <- "(KI|PH|OR|CL|FA)$"
# Replace matching values with "NA"
tax_clean <- tax_clean %>%
     mutate(across(everything(), ~ if_else(str_detect(.x, pattern), "NA", .x)))

# Similarly, a lot of the genera are labelled "uncultured-Familyname
# I want to get rid of this.
tax_clean$Genus[grepl("uncultured", tax_clean$Genus)] <- "NA"

# Noticed down the line that Hungateiclostridiaceae is a family not an order
# To remedy this, I'll move Hungateiclostridiaceae assignment to the family level.
# and replace Hungateiclostridiaceae by Oscillospirales at the order level 
tax_clean$Family[tax_clean$Order=="Hungateiclostridiaceae"]<- "Hungateiclostridiaceae"
tax_clean$Order[tax_clean$Order=="Hungateiclostridiaceae"]<- "Oscillospirales"

#Same thing for Desulfuromonadaceae
tax_clean$Family[tax_clean$Order=="Desulfuromonadaceae"]<- "Desulfuromonadaceae"
tax_clean$Order[tax_clean$Order=="Desulfuromonadaceae"]<- "Desulfuromonadales"

levels(as.factor(tax_clean$Genus))

#### ⚠ Pseudoreplicartes ####
# There are pseudoreplicates in the study,
# Mesocosms were sampled twice. 
# To deal with it we will sum the abundances of the pseudoreplicates

# What differs between the pseudoreplicates is the .1 .2 after the mesocosm id
names(raw_ab)

# Step 1: Group pseudoreplicates
# Extract core names to use as grouping factors
raw_ab_long <- raw_ab %>%
     rownames_to_column(var = "ASV") %>% # Create an ASV column based on rownames
     pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Abundance") %>% # Melt dataframe
     mutate(CoreName = gsub("(Mesocosm\\.\\d+\\.).*(D\\d+.*)", "\\1\\2", Sample))%>%
     group_by(ASV, CoreName) %>%
     summarize(Abundance = sum(Abundance), .groups = "drop")


# Step 3: Reshape back to wide format
raw_ab_wide <- raw_ab_long %>%
     pivot_wider(names_from = CoreName, values_from = Abundance)%>%
     as.data.frame()

# Check for example mesocoms 12, date d0, sediments
Meso12.d0.merged <- raw_ab_long[grepl("Mesocosm.12.D0.Sed", raw_ab_long$CoreName), ]

Meso12.d0.raw <- raw_ab[grepl("(Mesocosm.12).*(D0.Sed).*", names(raw_ab))]
Meso12.d0.raw$ASV <- row.names(Meso12.d0.raw)
Meso12.d0.raw$Abundance  <- rowSums(Meso12.d0.raw[ , !names(Meso12.d0.raw) %in% "ASV"])
Meso12.d0.raw <- Meso12.d0.raw[,c("ASV","Abundance")]

Meso12.d0.merged <- as.data.frame(Meso12.d0.merged[,c("ASV","Abundance")])
rownames(Meso12.d0.merged) <- Meso12.d0.merged$ASV

Meso12.d0.merged <- Meso12.d0.merged[order(rownames(Meso12.d0.merged)), ]
Meso12.d0.raw <- Meso12.d0.raw[order(rownames(Meso12.d0.raw)), ]
Meso12.d0.merged$Abundance <- as.numeric(Meso12.d0.merged$Abundance)
Meso12.d0.raw$Abundance <- as.numeric(Meso12.d0.raw$Abundance)
identical(Meso12.d0.raw,Meso12.d0.merged)


# Now let's make all these data.frames congruent
meta <- meta %>%
     mutate(Sample_name = gsub("(Mesocosm\\.\\d+\\.).*(D\\d+.*)", "\\1\\2", alias))%>%
     as.data.frame()
meta$Sample_name


meta_merged <- meta %>%
     group_by(Sample_name, Time,Compartment) %>%
     summarise(
          Water_type = first(Water_type),   # Keep the first value of Water_type
          Sample_type = first(Sample_type), # Keep the first value of Sample_type
          Mesocosm = first(Mesocosm),
          Group= first(Group),
          Sampling_date=first(Sampling_date),
          NA_concentration=first(NA_concentration),
          .groups = "drop"
     )%>%
     as.data.frame()


rownames(meta_merged) <- meta_merged$Sample_name

rownames(raw_ab_wide) <- raw_ab_wide$ASV
raw_ab_wide <- raw_ab_wide[,!names(raw_ab_wide)=="ASV"]


# Check congruence in sample names between sample data and Abundance data
identical(row.names(meta_merged),colnames(raw_ab_wide))

# Drop ASVs from Abundance table that we don't want based on taxonomy
raw_ab_clean <- raw_ab_wide[row.names(raw_ab_wide) %in% row.names(tax_clean),]



# ----------------------------------#
#### Phyloseq obj 
# ----------------------------------#

ASV_ab <- otu_table(raw_ab_clean,taxa_are_rows = T)

Tax <- tax_table(as.matrix(tax_clean))

Sdata <- sample_data(meta_merged)

ps <- phyloseq(ASV_ab,Tax,Sdata)
ps

save(ps, file = here(data_path,"ps_16S_Sed.RData"))


### Repeat for other datasets. 


##____####
#### ITS solid samples  #####

data_path <- here("data/16s_sed/")
data_path <- here("data/its_sed/")
list.files(path=data_path)

# Sample data
meta <- read.csv(file = paste0(data_path,"/mapping_file.csv"), 
                 dec = ".", header = T, row.names = 1, sep = ";", 
                 comment.char = "") #load 196 in 8 vars

str(meta)

meta <- subset(meta,meta$Sampling_date!="undefined") # Remove some samples with undefined dates.

meta$Group <- paste(meta$Compartment,meta$Sample_type,meta$Water_type,sep="_")
meta$Water_Stype <- paste(meta$Water_type,meta$Sample_type,sep="_")



meta[sapply(meta, is.character)] <- lapply(meta[sapply(meta, is.character)], 
                                           as.factor) ## Convert character vectors to factors


meta$Group <- factor(meta$Group, levels = c("Sediments_No_plant_Artificial_OSPW","Sediments_No_plant_OSPW",
                                            "Sediments_Carex_OSPW","Rhizosphere_Carex_OSPW","Roots_Carex_OSPW"))

meta$Time <- factor(meta$Time, levels = c("D-0","D8","D28","D42","D84"))
# 192 samples, 10 categories

# Abundance data
raw_ab <- read.csv(file = paste0(data_path,"/feature_table_filtered.tsv"), 
                   dec = ".", sep = "\t", header = T, row.names = 1,
                   comment.char = "")


names(raw_ab) <- gsub("X","",names(raw_ab)) #Erase the X added by R 

rownames(raw_ab) <- paste0("ASV",rownames(raw_ab))
# 3814 ASvs

colnames(raw_ab)
rownames(meta)



# Simplify sample names
# First rename mispelled Mesocosm
colnames(raw_ab) <- gsub("Mesosocm", "Mesocosm", colnames(raw_ab), perl = TRUE)
colnames(raw_ab)<- gsub("Mesococm", "Mesocosm", colnames(raw_ab), perl = TRUE)
# Second remove the numbers before mesocosm
colnames(raw_ab) <- gsub(".*(?=Mesocosm)", "", colnames(raw_ab), perl = TRUE)
# Third change D.0 to D0 
colnames(raw_ab) <- gsub("D.0", "D0", colnames(raw_ab), perl = TRUE)
colnames(raw_ab)<- gsub("TEST", "D84",colnames(raw_ab), perl = TRUE)
colnames(raw_ab) <- gsub("Racines","Racine",colnames(raw_ab))


# First have the original sample data names match the ASV abundance original sample names
meta$alias <- gsub("-",".",meta$alias)
meta$alias <- gsub("Mesosocm", "Mesocosm", meta$alias, perl = TRUE)
meta$alias<- gsub("Mesococm", "Mesocosm", meta$alias, perl = TRUE)
# Second remove the numbers before mesocosm
meta$alias <- gsub(".*(?=Mesocosm)", "", meta$alias, perl = TRUE)

# Third change D.0 to D0 and TEST to D84 
meta$alias <- gsub("D.0", "D0", meta$alias, perl = TRUE)
meta$alias <- gsub("TEST", "D84", meta$alias, perl = TRUE)
meta$alias <- gsub("Racines","Racine",meta$alias)

rownames(meta) <- meta$alias 


common_samples <- intersect(rownames(meta),colnames(raw_ab))
unique_to_meta <- rownames(meta)[!(rownames(meta)%in%common_samples)]
unique_to_ab <-  colnames(raw_ab)[!(colnames(raw_ab)%in%common_samples)]
# The samples unique to the abundance table are the ones that were removed 
# due to the undefined information. There's a also a taxonomy column in the abundance data, 
# we'll put it in it's owndata.frame, and split the individual taxonomy levels in individual columns

# Let's remove the undefined samples
raw_ab <- raw_ab[(colnames(raw_ab)%in%c(common_samples,"taxonomy"))]


taxonomy <- as.data.frame(raw_ab[,ncol(raw_ab)]) # taxonomy from last column as independent table
rownames(taxonomy) <- rownames((raw_ab))
colnames(taxonomy) <- c("Taxonomy")

# Now we can remove the taxonomy column from raw_ab
raw_ab <- raw_ab[,-ncol(raw_ab)]

# Separate the taxonomic levels in individual columns
tax_clean <- separate(taxonomy, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";")

tax_clean <- as.data.frame(tax_clean)
levels(as.factor(tax_clean$Domain)) # 8   Domains 
levels(as.factor(tax_clean$Phylum)) # 21  Phyla (one of which is empty) 
levels(as.factor(tax_clean$Class))  # 63  Classes
levels(as.factor(tax_clean$Order))  # 147 Orders
levels(as.factor(tax_clean$Family)) # 309 Families
levels(as.factor(tax_clean$Genus))  # 592 Genera
str(tax_clean)

# Remove all the k__, g__ etc from taxonomy 
tax_clean[] <- lapply(tax_clean, function(x) substr(x, 4, nchar(x)))

# Replace empty and NA cells with "NA" 
tax_clean[]<- lapply(tax_clean[], function(x) ifelse(x == "" | is.na(x), "NA", x))




# There are also pseudoreplicates in the study,
# Mesocosms were sampled twice. 
# We'll sum the abundances of this pseudoreplicates
names(raw_ab)

# The thing that differs is the .1 .2 after the mesocosm id

# Step 1: Group pseudoreplicates
# Extract core names to use as grouping factors
raw_ab_long <- raw_ab %>%
    rownames_to_column(var = "ASV") %>% # Create an ASV column based on rownames
    pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Abundance") %>% # Melt dataframe
    mutate(CoreName = gsub("(Mesocosm\\.\\d+\\.).*(D\\d+.*)", "\\1\\2", Sample))%>%
    group_by(ASV, CoreName) %>%
    summarize(Abundance = sum(Abundance), .groups = "drop")


# Step 3: Reshape back to wide format
raw_ab_wide <- raw_ab_long %>%
    pivot_wider(names_from = CoreName, values_from = Abundance)%>%
    as.data.frame()

# Check for example mesocoms 12, date d0, sediments
Meso12.d0.merged <- raw_ab_long[grepl("Mesocosm.12.D0.Sed", raw_ab_long$CoreName), ]

Meso12.d0.raw <- raw_ab[grepl("(Mesocosm.12).*(D0.Sed).*", names(raw_ab))]
Meso12.d0.raw$ASV <- row.names(Meso12.d0.raw)
Meso12.d0.raw$Abundance  <- rowSums(Meso12.d0.raw[ , !names(Meso12.d0.raw) %in% "ASV"])
Meso12.d0.raw <- Meso12.d0.raw[,c("ASV","Abundance")]

Meso12.d0.merged <- as.data.frame(Meso12.d0.merged[,c("ASV","Abundance")])
rownames(Meso12.d0.merged) <- Meso12.d0.merged$ASV

Meso12.d0.merged <- Meso12.d0.merged[order(rownames(Meso12.d0.merged)), ]
Meso12.d0.raw <- Meso12.d0.raw[order(rownames(Meso12.d0.raw)), ]
Meso12.d0.merged$Abundance <- as.numeric(Meso12.d0.merged$Abundance)
Meso12.d0.raw$Abundance <- as.numeric(Meso12.d0.raw$Abundance)
identical(Meso12.d0.raw,Meso12.d0.merged)


# Now let's make all these data.frames congruent
meta <- meta %>%
    mutate(Sample_name = gsub("(Mesocosm\\.\\d+\\.).*(D\\d+.*)", "\\1\\2", alias))%>%
    as.data.frame()
meta$Sample_name


meta_merged <- meta %>%
    group_by(Sample_name, Time,Compartment) %>%
    summarise(
        Water_type = first(Water_type),   # Keep the first value of Water_type
        Sample_type = first(Sample_type), # Keep the first value of Sample_type
        Mesocosm = first(Mesocosm), # etc
        Group= first(Group),
        Water_Stype=first(Water_Stype),
        Sampling_date=first(Sampling_date),
        NA_concentration=first(NA_concentration),
        .groups = "drop"
    )%>%
    as.data.frame()


rownames(meta_merged) <- meta_merged$Sample_name

rownames(raw_ab_wide) <- raw_ab_wide$ASV
raw_ab_wide <- raw_ab_wide[,!names(raw_ab_wide)=="ASV"]


# Check congruence in sample names between sample data and Abundance data
identical(row.names(meta_merged),colnames(raw_ab_wide))

# Drop ASVs from Abundance table that we don't want based on taxonomy
#raw_ab_clean <- raw_ab_wide[row.names(raw_ab_wide) %in% row.names(tax_clean),]


### Phyloseq object ###

ASV_ab <- otu_table(raw_ab_wide,taxa_are_rows = T)

Tax <- tax_table(as.matrix(tax_clean))

Sdata <- sample_data(meta_merged)

ps <- phyloseq(ASV_ab,Tax,Sdata)
save(ps, file = here(data_path,"ps_ITS.RData"))

ps <- subset_taxa(ps, Domain=="Fungi") # 2916 fungal ASVs
save(ps, file = here(data_path,"ps_fungi.RData"))



#### 18S water samples #####
data_path <- here("data/18s_wat/")

# Sample data
meta <- read.csv(file = paste0(data_path,"/mapping_file.tsv"), 
                 dec = ".", header = T, row.names = 1, sep = "\t", 
                 comment.char = "") #load 196 in 7 vars
str(meta)
meta[sapply(meta, is.character)] <- lapply(meta[sapply(meta, is.character)], 
                                           as.factor) ## Convert character vectors to factors
meta$Group <- paste(meta$Sample_type,meta$Water_type,sep="_")
meta$Group <- factor(meta$Group, levels = c("No_plant_Artificial_OSPW","Carex_OSPW","No_plant_OSPW"))

meta$Time <- factor(meta$Time, levels = c("D-0","D8","D28","D42","D84","BlankpcrCES"))

# 62 samples, 8 categories

# Abundance data
raw_ab <- read.csv(file = paste0(data_path,"/feature_table_filtered.tsv"), dec = ".", sep = "\t", header = T, row.names = 1, comment.char = "")
rownames(raw_ab) <- paste0("ASV",rownames(raw_ab))
# 60 samples, 484 ASvs

#Sample names differ slightly between Sample data and Abundance data
# so we'll fix it
rownames(meta) <- gsub("-",".",rownames(meta))

common_samples <- intersect(rownames(meta),colnames(raw_ab))
unique_to_meta <- rownames(meta)[!(rownames(meta)%in%common_samples)]
unique_to_ab <-  colnames(raw_ab)[!(colnames(raw_ab)%in%common_samples)]

# There's a taxonomy column in the abundance data, we'll put it in it's own
# data.frame, and split the individual taxonomy levels in individual columns

taxonomy <- as.data.frame(raw_ab[,ncol(raw_ab)]) # taxonomy from last column as independent table
rownames(taxonomy) <- rownames((raw_ab))
colnames(taxonomy) <- c("Taxonomy")

# Now we can remove the taxonomy column from raw_ab
raw_ab <- raw_ab[, !names(raw_ab) %in% "taxonomy"]

# Separate the taxonomic levels in individual columns
tax_clean <- separate(taxonomy, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
tax_clean <- tax_clean[,-which(names(tax_clean) %in% "Species")]
# Removes species (No information there)

tax_clean <- as.data.frame(tax_clean)
levels(as.factor(tax_clean$Domain)) # 2 Domain 
levels(as.factor(tax_clean$Phylum)) #  32 Phyla 
levels(as.factor(tax_clean$Class)) # 56 classes
levels(as.factor(tax_clean$Order)) # 80 orders
levels(as.factor(tax_clean$Family)) # 89 Families
levels(as.factor(tax_clean$Genus)) # 152 genus
str(tax_clean)

# Remove all the k__, g__ etc from taxonomy 
tax_clean[] <- lapply(tax_clean, function(x) substr(x, 4, nchar(x)))

# Replace empty and NA cells with "NA" 
tax_clean[]<- lapply(tax_clean[], function(x) ifelse(x == "" | is.na(x), "NA", x))

# Remove Chloroplasts and Mitochondria
tax_chloro <- subset(tax_clean,tax_clean$Order=="Chloroplast") # 0 Choloplast sequence
tax_mito <- subset(tax_clean,tax_clean$Family=="Mitochondria") # 1 Mitochondria sequence
rm(tax_chloro,tax_mito)
tax_clean <- subset(tax_clean,tax_clean$Order!="Chloroplast"&tax_clean$Family!="Mitochondria")


## Some taxa have been classified at a taxonomy level but with the previous
## taxonomy name for instance BacillariophyceaeCL at the family level...
## I want to get rid of those and replace them by NA as their taxonomy was not 
## assigned at these lower taxonomy levels.

pattern <- "(KI|PH|OR|CL|FA)$"
# Replace matching values with "NA"
tax_clean <- tax_clean %>%
    mutate(across(everything(), ~ if_else(str_detect(.x, pattern), "NA", .x)))

# Chlorophyta phylum as a space at the beginning
tax_clean$Phylum <- gsub(" Chlorophyta","Chlorophyta",tax_clean$Phylum)


# Some taxonomy assignment skip one taxonomy level, here is a way to remedy 
# some of them.
pattern="Trebouxiophyceae|Chlorodendrophyceae|Nephroselmidophyceae"
tax_clean <- tax_clean %>%
    mutate(Phylum = if_else(str_detect(Class, pattern), "Chlorophyta", Phylum))

pattern="Chrysophyceae|Eustigmatophyceae"
tax_clean <- tax_clean %>%
    mutate(Phylum = if_else(str_detect(Class, pattern), "Ochrophyta", Phylum))

pattern="Choanoflagellida"
tax_clean <- tax_clean %>%
    mutate(Phylum = if_else(str_detect(Class, pattern), "Choanozoa", Phylum))

pattern="Tubulinea"
tax_clean <- tax_clean %>%
    mutate(Phylum = if_else(str_detect(Class, pattern), "Amoebozoa", Phylum))

pattern="Ichthyosporea"
tax_clean <- tax_clean %>%
    mutate(Phylum = if_else(str_detect(Class, pattern), "Mesomycetozoa", Phylum))

# Update Class and Phylum based on order
pattern="Cryptomonadales"
tax_clean <- tax_clean %>%
    mutate(Phylum = if_else(str_detect(Order, pattern), "Cryptista", Phylum))%>%
    mutate(Class = if_else(str_detect(Order, pattern), "Cryptophyceae", Class))

# Update class 
pattern="Peronosporomycetes"
tax_clean <- tax_clean %>%
    mutate(Class = if_else(str_detect(Phylum, pattern), "Peronosporomycetes", Class))

pattern="Chytridiomycota"
tax_clean <- tax_clean %>%
    mutate(Class = if_else(str_detect(Phylum, pattern), "Chytridiomycetes", Class))

tax_clean$Class <- gsub("Incertae Sedis","NA",tax_clean$Class)

# Now let's make all these data.frames congruent

# Drop samples not present in abundance table
meta_clean <- meta[row.names(meta) %in% colnames(raw_ab),]

# Drop chloroplast and mitochondria ASVs from Abundance table
raw_ab_clean <- raw_ab[row.names(raw_ab) %in% row.names(tax_clean),]

### Phyloseq object ###

ASV_ab <- otu_table(raw_ab_clean,taxa_are_rows = T)

Tax <- tax_table(as.matrix(tax_clean))

Sdata <- sample_data(meta)

ps <- phyloseq(ASV_ab,Tax,Sdata)
ps

save(ps, file = here(data_path,"ps_18S_Wat.RData"))
# 59 samples, 483 taxa



#### 16S water samples #####

source(here::here("src", "libraries.R"))

data_path <- here("data/16s_wat")

R.version

# Sample data
meta <- read.csv(file = paste0(data_path,"/mapping_file.tsv"), 
                 dec = ".", header = T, row.names = 1, sep = "\t", 
                 comment.char = "") #load 196 in 7 vars
str(meta)
meta[sapply(meta, is.character)] <- lapply(meta[sapply(meta, is.character)], 
                                           as.factor) ## Convert character vectors to factors
meta$Group <- paste(meta$Sample_type,meta$Water_type,sep="_")
meta$Group <- factor(meta$Group, levels = c("No_plant_Artificial_OSPW","Carex_OSPW","No_plant_OSPW"))

meta$Time <- factor(meta$Time, levels = c("D-0","D8","D28","D42","D84","BlankpcrCES"))

# 62 samples, 8 categories

# Abundance data
raw_ab <- read.csv(file = paste0(data_path,"/feature_table_filtered.tsv"), dec = ".", sep = "\t", header = T, row.names = 1, comment.char = "")
rownames(raw_ab) <- paste0("ASV",rownames(raw_ab))
# 59 samples, 2957 ASvs

#Sample names differ slightly between Sample data and Abundance data
# so we'll fix it
rownames(meta) <- gsub("-",".",rownames(meta))

common_samples <- intersect(rownames(meta),colnames(raw_ab))
unique_to_meta <- rownames(meta)[!(rownames(meta)%in%common_samples)]
unique_to_ab <-  colnames(raw_ab)[!(colnames(raw_ab)%in%common_samples)]

# There's a taxonomy column in the abundance data, we'll put it in it's own
# data.frame, and split the individual taxonomy levels in individual columns

taxonomy <- as.data.frame(raw_ab[,ncol(raw_ab)]) # taxonomy from last column as independent table
rownames(taxonomy) <- rownames((raw_ab))
colnames(taxonomy) <- c("Taxonomy")

# Now we can remove the taxonomy column from raw_ab
raw_ab <- raw_ab[, !names(raw_ab) %in% "taxonomy"]

# Separate the taxonomic levels in individual columns
tax_clean <- separate(taxonomy, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep = ";")
tax_clean <- tax_clean[,-which(names(tax_clean) %in% c("Species", "Strain"))]
# Removes species and strains (No information there)

tax_clean <- as.data.frame(tax_clean)
levels(as.factor(tax_clean$Domain)) # 3 Domain 
levels(as.factor(tax_clean$Phylum)) #  39 Phyla 
levels(as.factor(tax_clean$Class)) # 87 classes
levels(as.factor(tax_clean$Order)) # 201 orders
levels(as.factor(tax_clean$Family)) # 321 Families
levels(as.factor(tax_clean$Genus)) # 617 genus
str(tax_clean)

# Remove all the k__, g__ etc from taxonomy 
tax_clean[] <- lapply(tax_clean, function(x) substr(x, 4, nchar(x)))

# Replace empty and NA cells with "NA" 
tax_clean[]<- lapply(tax_clean[], function(x) ifelse(x == "" | is.na(x) | x =="uncultured", "NA", x))

# Remove Chloroplasts, Mitochondria, Eukaryotes
tax_chloro <- subset(tax_clean,tax_clean$Order=="Chloroplast") # 64 Choloplast sequence
tax_mito <- subset(tax_clean,tax_clean$Family=="Mitochondria") # 25 Mitochondriq sequence
tax_euk <- subset(tax_clean,tax_clean$Domain=="Eukaryota") # 12 eukaryote sequence (Cryptomycota)
rm(tax_chloro,tax_mito,tax_euk)
tax_clean <- subset(tax_clean,tax_clean$Order!="Chloroplast"&
                        tax_clean$Family!="Mitochondria"&
                        tax_clean$Domain!="Eukaryota")

## Some taxa have been classified at a taxonomy level but with the previous
## taxonomy name for instance FimbriimonadaceaeFA at the genus level...
## I want to get rid of those and replace them by NA as their taxonomy was not 
## assigned at these lower taxonomy levels.
pattern <- "(KI|PH|OR|CL|FA)$"
# Replace matching values with "NA"
tax_clean <- tax_clean %>%
    mutate(across(everything(), ~ if_else(str_detect(.x, pattern), "NA", .x)))

# Similarly, a lot of the genera are labelled "uncultured-Familyname
# I want to get rid of this.
tax_clean$Genus[grepl("uncultured", tax_clean$Genus)] <- "NA"


# Now let's make all these data.frames congruent

# Drop samples not present in abundance table
meta_clean <- meta[row.names(meta) %in% colnames(raw_ab),]

# Drop chloroplast and mitochondria ASVs from Abundance table
raw_ab_clean <- raw_ab[row.names(raw_ab) %in% row.names(tax_clean),]



#### Phyloseq object 

ASV_ab <- otu_table(raw_ab_clean,taxa_are_rows = T)

Tax <- tax_table(as.matrix(tax_clean))

meta_clean$Sample_name <- paste(meta_clean$Mesocosm,meta_clean$Time,sep="_")
Sdata <- sample_data(meta_clean)

ps <- phyloseq(ASV_ab,Tax,Sdata)
ps
# 58 samples, 2868 taxa

names <- as.data.frame(sample_data(ps))
sample_names(ps) <- names$Sample_name

plot_bar(ps)
# Mesocosm_3-D-0 has almost no data, we can remove it. 

ps <- subset_samples(ps,sample_names(ps)!="Mesocosm_3_D-0")
ps<- prune_taxa(taxa_sums(ps)>0, ps)

save(ps, file = here(data_path,"ps_16S_Wat.RData"))


