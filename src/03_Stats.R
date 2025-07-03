# Simon Morvan
# July 2025
# R version 4.4.2 - Pile of Leaves
# 
#                           Title : Statistical analysis 
#                           
# Description : Script to generate the statistical analysis presented in the manuscript 

source(here::here("src", "libraries.R"))

##### Beta divesity -  Permanova ####


Permanova <- function(data_path, file_name) {
     
     load(here::here(data_path, file_name))
     
     # Data processing pipeline
     exclude_groups <-  c("Roots_Carex_OSPW",
                          "Sediments_No_plant_Artificial_OSPW",
                          "No_plant_Artificial_OSPW")
     
     # Groups that are not of interest for the analysis are removed 
     ps_sub <- prune_samples(!(sample_data(ps)$Group %in% exclude_groups), ps)
     ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
     
     # RCLR Distance
     ps_dist <- microbiome::transform(ps_sub, "rclr")
     
     rclr_dist_matrix <- phyloseq::distance(ps_dist, method = "euclidean") # Robust Aitchison distance
     metadata <- as(sample_data(ps_dist), "data.frame") # Metadata
     
     # Permanova
     Permanova_Sample<-adonis2(rclr_dist_matrix~Group,
                               data = metadata,
                               strata = metadata$Time,
                               by="terms",
                               permutations = 999)
     
     Permanova_Time<-adonis2(rclr_dist_matrix~Time,
                             data = metadata,
                             strata = metadata$Group,
                             by="terms",
                             permutations = 999)
     
     return(list(Permanova_Sample = Permanova_Sample, 
                 Permanova_Time = Permanova_Time))
          
     }

# Process all datasets
Bac_solid <- Permanova("data/16s_sed/", "ps_16S_Sed.RData")
Bac_water <- Permanova("data/16s_wat/", "ps_16S_Wat.RData")
Fun_solid <- Permanova("data/its_sed/", "ps_fungi.RData")
Euk_water <- Permanova("data/18s_wat/", "ps_18S_Wat.RData")



##### Pairwise Permanova ####

# Extracts p_value from pairwise permanova results
# to adjust for multiple testing
extract_p_values <- function(result) {
     if ("Model" %in% rownames(result)) {
          return(result["Model", "Pr(>F)"])
     } else {
          return(NA)
     }
}

# Pairwise permanova function
PW_Permanova <- function(data_path, file_name, factor){
     
     load(here::here(data_path, file_name))

     # Data processing pipeline
     exclude_groups <-  c("Roots_Carex_OSPW",
                          "Sediments_No_plant_Artificial_OSPW",
                          "No_plant_Artificial_OSPW")
     
     # Groups that are not of interest for the analysis are removed 
     ps_sub <- prune_samples(!(sample_data(ps)$Group %in% exclude_groups), ps)
     ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
     
     # RCLR Distance
     ps_dist <- microbiome::transform(ps_sub, "rclr")
     
     rclr_dist_matrix <- phyloseq::distance(ps_dist, method = "euclidean") # Robust Aitchison distance
     metadata <- as(sample_data(ps_dist), "data.frame") # Metadata

     # Pairwise testing
     pairwise_results <- pairwise.adonis2(as.formula(paste("rclr_dist_matrix ~", factor)), 
                                          data = metadata)

     # Need to adjust p_value for multiple testing
     # Extract p values
     p_values <- sapply(pairwise_results[-1], extract_p_values)
     
     # Apply Bonferroni correction
     p_adjusted <- p.adjust(p_values, method = "bonferroni")
     
     # Add adjusted p-values back to the pairwise_results obj
     for (i in seq_along(pairwise_results)[-1]) {
          pairwise_results[[i]]["Model", "P_adj_Bonf"] <- p_adjusted[i-1]
     }
     return(list(pairwise_results=pairwise_results))
     }

Bac_solid_time <- PW_Permanova("data/16s_sed/", "ps_16S_Sed.RData","Time")
Bac_solid_stype <- PW_Permanova("data/16s_sed/", "ps_16S_Sed.RData","Group")

Bac_water_time <- PW_Permanova("data/16s_wat/", "ps_16S_Wat.RData","Time")
Bac_water_stype <- PW_Permanova("data/16s_wat/", "ps_16S_Wat.RData","Group")

Fun_solid_time <- PW_Permanova("data/its_sed/", "ps_fungi.RData","Time")
Fun_solid_stype <- PW_Permanova("data/its_sed/", "ps_fungi.RData","Group")

Euk_water_time <- PW_Permanova("data/18s_wat/", "ps_18S_Wat.RData","Time")
Euk_water_group <- PW_Permanova("data/18s_wat/", "ps_18S_Wat.RData","Group")


#### Indicator species #### 

load(here(data_path,"ps_16S_Sed.RData"))

## Optional 
ps_sub <- subset_samples(ps,Group=="Rhizosphere_Carex_OSPW"|Group=="Sediments_No_plant_OSPW")
ps_sub <- prune_taxa(taxa_sums(ps_sub)>0, ps_sub)
sample_names(ps_sub)
ps <- ps_sub

# RCLR transformation
ps_rclr<- microbiome::transform(ps, "rclr") # Robust Aitchison transformation
ps <- ps_rclr

tax_lvl <- "Family"
ps_glom <- tax_glom(ps,taxrank=tax_lvl,NArm = T,bad_empty="NA")
#Removes Taxa with NA as tax_lvl

ps_glom <- prune_taxa(rowSums(otu_table(ps_glom) == 0) < ncol(otu_table(ps_glom)) * 0.75, ps_glom)
# This code removes (or "prunes") taxa that are absent (have zero abundance) in 75% 
# or more of the samples. 
# In other words, only taxa that appear in at least 25% of the samples 
# are retained in the phyloseq object.


tax <- as.data.frame(tax_table(ps_glom))
Abund <- as.data.frame((otu_table(ps_glom)))

merged_df <- merge(Abund, tax[, tax_lvl, drop = FALSE], by = "row.names", all = TRUE)
rownames(merged_df) <- merged_df[[tax_lvl]]
merged_df <- merged_df[,-c(1,ncol(merged_df))]


Sdata <- as.data.frame(sample_data(ps_glom))
merged_df <- as.data.frame(t(merged_df))

###### Sample Type ####
h <- how(nperm = 999, 
         blocks = Sdata$Time)

ASV_bac_indval<- multipatt(merged_df,Sdata$Group,
                           control=h,
                           #restcomb = c(1,2,3,4,5,25),
                           duleg =T,
                           print.perm =TRUE)

summary(ASV_bac_indval,alpha=1,indvalcomp=TRUE)
# Component ‘A’ (Specificity) is sample estimate of the probability that the surveyed site 
# belongs to the target site group given the fact that the species has been found. 
# This conditional probability is called the specificity or positive predictive 
# value of the species as indicator of the site group.
# So if A = 1; the specie is present only in this group 
# 
# Component ‘B’ is sample estimate of the probability of finding the species 
# in sites belonging to the site group. This second conditional probability 
# is called the fidelity or sensitivity of the species as indicator of the target site group.
# So if B = 1; the specie is present in all the sites of the group

P_adj <- ASV_bac_indval$sign
P_adj$p_val_adj_holm <- p.adjust(P_adj$p.value, "holm")
P_adj$p_val_adj_fdr <- p.adjust(P_adj$p.value, "fdr")
view(P_adj)

P_adj_signif <- subset(P_adj,P_adj$p_val_adj_fdr<0.05)



ASV_bac_indval<- multipatt(merged_df,Sdata$Group,
                           control=how(nperm=999),
                           duleg =F,
                           print.perm =TRUE)

summary(ASV_bac_indval,alpha=1,indvalcomp=TRUE)

P_adj <- ASV_bac_indval$sign
P_adj$p_val_adj_holm <- p.adjust(P_adj$p.value, "holm")
P_adj$p_val_adj_fdr <- p.adjust(P_adj$p.value, "fdr")
view(P_adj)

P_adj_signif <- subset(P_adj,P_adj$p_val_adj_fdr<0.1)







