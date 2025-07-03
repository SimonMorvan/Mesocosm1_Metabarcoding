# Simon Morvan
# July 2025
# R version 4.4.2 - Pile of Leaves
# 
#                           Title : Figures 
# Description : Script to generate the figures presented in the manuscript 

source(here::here("src", "libraries.R"))

#### Figure 1 - NA concentration #####
data_path <- here("data/its_sed/")

metadata <- read.csv(file = paste0(data_path,"/NA_predict.csv"), 
                     dec = ".", header = T, row.names = 1, sep = ";", 
                     comment.char = "") 

metadata_sub <- subset(metadata, Compartment=="Sediments")

metadata_sub$Group <- factor(metadata_sub$Group, levels = c("Sediments_No_plant_Artificial_OSPW","Sediments_No_plant_OSPW","Sediments_Carex_OSPW"), labels = c("Unplanted LPW","Unplanted OSPW","Carex OSPW "))
metadata_sub$Time <- factor(metadata_sub$Time, levels = c("D-0","D8","D28","D42","D84"))
metadata_sub <- subset(metadata_sub, Water_type=="OSPW")

NA_CONC <- ggplot(metadata_sub, aes(x=Time,y=NA_concentration,fill=Group))+
    theme_bw()+
    theme(#title = element_text(size=20),
        axis.text.x =element_blank(),
        axis.text.y =element_text(size=12),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=15),
        legend.text.align = 0,
        axis.ticks.x = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.background = element_rect(color = "black"),
        legend.position = "bottom",
        legend.spacing = unit(0.01,"cm"),
        strip.text = element_text(size=15),
        strip.background = element_rect(fill='white'))+
    geom_pwc(method="dunn_test",
             hide.ns=F,
             label="p.adj.format",
             p.adjust.method = "bonferroni",
             label.size = 5,
             #tip.length = 0, 
             #position="identity",
             step.increase = 0.05,
             linetype = 1,
             size = 1)+
    facet_grid(~Time, scales = "free")+
    geom_boxplot(outliers = F)+
    geom_point(position=position_dodge(width=0.75),aes(fill=Group),alpha=0.3,shape=21,size=2.5)+
    #geom_point(aes(x=Time,y=NA_concentration,fill=Group,stroke = 0.7),alpha=0.3,position=position_dodge(width=0.75,preserve = "total"),size=2.5)+
    scale_fill_manual(name="Sample type :",values=c("#687967","#9DEE74"))+
    labs(y = "NAFCs concentration (mg/L)")


ggsave(NA_CONC,device='svg',height = 6, width = 15,
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/NA_CONC.svg")


#  metadata_sub$Group <- factor(metadata_sub$Group, levels = c("Sediments_No_plant_Artificial_OSPW","Sediments_No_plant_OSPW","Sediments_Carex_OSPW"), labels = c("LPW","OSPW","Carex OSPW "))
#  
#  # If p-value > 0.05 variance is homogeneous 
#  bartlett.test(NA_concentration ~Group, data=metadata_sub) 
#  #p-value < 2.2e-16 --» Variance is not homogeneous
#  
#  Select_Time <- "D84"
#  metadata_sub_time <- subset(metadata_sub,metadata_sub$Time==Select_Time)
#  kruskal.test(NA_concentration ~Group,data=metadata_sub_time)
#  dunn.test(metadata_sub_time$NA_concentration,
#            g=metadata_sub_time$Group,
#            method="bonferroni")



###____####
#### Figure 2 - Beta diversity ####

source(here::here("src", "libraries.R"))

Beta_div_Ordination <- function(data_path, file_name) {

    # Load phyloseq object
    load(here::here(data_path, file_name))
   
    # Data processing pipeline
    exclude_groups <-  c("Roots_Carex_OSPW",
                         "Sediments_No_plant_Artificial_OSPW",
                         "No_plant_Artificial_OSPW")
    # Groups that are not of interest for the plot are removed 
    ps_sub <- prune_samples(!(sample_data(ps)$Group %in% exclude_groups), ps)
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
    # 

    # Ordination analysis
    ps_dist <- microbiome::transform(ps_sub, "rclr")
    ord_rclr <- phyloseq::ordinate(ps_dist, "PCoA", distance = "euclidean")
    
    # Calculate variance explained
    PC1 <- round(ord_rclr$values$Relative_eig[1] * 100, 1)
    PC2 <- round(ord_rclr$values$Relative_eig[2] * 100, 1)
    
    # Create ordination dataframe
    Ordination <- phyloseq::plot_ordination(ps_dist, ord_rclr, type = "samples", justDF = TRUE)
    Ordination$Time <- factor(Ordination$Time, levels = c("D-0", "D8", "D28", "D42", "D84"))
    levels(Ordination$Group)
    # Conditional formatting
    if (nlevels(Ordination$Group) == 3) {
        palette <- c("#bad9a2", "#80B39F", "#687967")
        Ordination$Group <- factor(Ordination$Group,
                                   levels = c("Rhizosphere_Carex_OSPW", "Sediments_Carex_OSPW",
                                              "Sediments_No_plant_OSPW"),
                                   labels = c("Carex rhizosphere", "Carex sediments", 
                                              "Unplanted sediments"))
    } else {
        palette <- c("#9DEE74", "#687967")
        Ordination$Group <- factor(Ordination$Group,
                                   levels = c("Carex_OSPW", "No_plant_OSPW"),
                                   labels = c("Carex OSPW", "Unplanted OSPW"))
    }
    
    # Create plots
    plot_stype <- ggplot(Ordination, aes(x = Axis.1, y = Axis.2, color = Group)) +
        geom_point(size = 4, aes(stroke = 3)) +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.box = "vertical",
              legend.spacing = unit(0.01, "cm"),
              text = element_text(size = 25)) +
        xlab(paste0("PCo1 (", PC1, "%)")) + 
        ylab(paste0("PCo2 (", PC2, "%)")) + 
        labs(title = "" ) + #Add space above graph
        scale_color_manual(values = palette, name = "") +
        guides(color = guide_legend(nrow = 1))
    
    plot_time <- ggplot(Ordination, aes(x = Axis.1, y = Axis.2, color = Time, shape = Group)) +
        geom_point(size = 4, aes(stroke = 3)) +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.box = "vertical",
              legend.spacing = unit(0.01, "cm"),
              text = element_text(size = 25)) +
        xlab(paste0("PCo1 (", PC1, "%)")) + 
        ylab(paste0("PCo2 (", PC2, "%)")) + 
        labs(title = "" ) + #Add space above graph
        scale_shape_manual(values = c(1, 2, 3, 4), name = "") +
        scale_color_manual(values = c("#0D0887FF", "#7E03A8FF", "#CC4678FF", 
                                      "#F89441FF", "#F0F921FF"), name = "") +
        guides(shape = guide_legend(nrow = 1))
    
    return(list(plot_stype = plot_stype, plot_time = plot_time))
}

source(here::here("src", "libraries.R"))

# Process all datasets
Bac_solid <- Beta_div_Ordination("data/16s_sed/", "ps_16S_Sed.RData")
Bac_water <- Beta_div_Ordination("data/16s_wat/", "ps_16S_Wat.RData")
Fun_solid <- Beta_div_Ordination("data/its_sed/", "ps_fungi.RData")
Euk_water <- Beta_div_Ordination("data/18s_wat/", "ps_18S_Wat.RData")

# Arrange all the figures in one plot
Beta_stype <- ggarrange(plot(Bac_solid$plot_stype),
                        plot(Bac_water$plot_stype),
                        plot(Fun_solid$plot_stype),
                        plot(Euk_water$plot_stype),labels="AUTO")


ggsave(Beta_stype,device='svg',height = 10, width = 18 ,
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/Beta_stype2.svg")



###____####
#### Figure 3 - Relative Abundance ####

source(here::here("src", "libraries.R"))

Relative_Abundance <- function(data_path, file_name,label,dataset) {
    
    # Load phyloseq object
    load(here::here(data_path, file_name))
    
    # Data processing pipeline
    exclude_groups <-  c("Roots_Carex_OSPW",
                         "Sediments_Carex_OSPW",
                         "Sediments_No_plant_Artificial_OSPW",
                         "No_plant_Artificial_OSPW")
    
    # Groups that are not of interest for the plot are removed 
    ps_sub <- prune_samples(!(sample_data(ps)$Group %in% exclude_groups), ps)
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
    
    # Select taxonomic level to group by
    tax_lvl <- "Class"
    ps_glom <- tax_glom(ps_sub, taxrank = tax_lvl)
    
    # Transform to long format for plotting purposes
    ps_long <- psmelt(ps_glom)
    
    # Group by different columns and compute the mean and relative mean (mean across replicates)
    ps_mean <-  ps_long %>%
        group_by(Group, Time,Compartment,Water_type, Sample_type, !!sym(tax_lvl)) %>%
        summarise(Mean = mean(Abundance, na.rm = TRUE))%>%
        mutate(Relative_Mean = Mean / sum(Mean))
    
    if (file_name=="ps_18S_Wat.RData"){
        # Select the top  most abundant taxa for each sample type (group) x Time combination 
        # For Eukaryotes select the top 5 as there is a lot of variation
        Top <- ps_mean %>%
            group_by(Group, Time, !!sym(tax_lvl)) %>%
            arrange(Group, Time, desc(Relative_Mean)) %>%
            group_by(Group, Time) %>%
            slice_head(n = 5) %>%
            pull(!!sym(tax_lvl)) %>%
            unique()
        # For the other dataset select the top10 
    } else {
        Top <- ps_mean %>%
        group_by(Group, Time, !!sym(tax_lvl)) %>%
        arrange(Group, Time, desc(Relative_Mean)) %>%
        group_by(Group, Time) %>%
        slice_head(n = 10) %>%
        pull(!!sym(tax_lvl)) %>%
        unique()
    }
    
    # Classify all the other taxa that remain that are not in the top10 nor NA as "Other"
    ps_mean_top <-  ps_mean %>%
        group_by(Group, Time,Compartment,Water_type, Sample_type, !!sym(tax_lvl)) %>%
        mutate(!!sym(tax_lvl) := ifelse(!!sym(tax_lvl) %in% Top | !!sym(tax_lvl) == "NA", !!sym(tax_lvl),"Other"))

     
    # Rename categories for plotting purposes
    
    ps_mean_top$Sample_type <- factor(ps_mean_top$Sample_type, 
                                      levels = c( "Carex","No_plant"), 
                                      labels = c("Carex","Unplanted"))
    
    ps_mean_top$Water_type <- factor(ps_mean_top$Water_type, 
                                     levels = c( "OSPW","Artificial_OSPW"), 
                                     labels = c("OSPW","LPW"))
    
    if (dataset == "Sediments") {
        # For the sediments datasets
        ps_mean_top$Group <- factor(ps_mean_top$Group,
                                    levels=c("Roots_Carex_OSPW","Rhizosphere_Carex_OSPW","Sediments_Carex_OSPW","Sediments_No_plant_OSPW","Sediments_No_plant_Artificial_OSPW"),
                                    labels=c("Carex Ro.","Carex Rh.","Carex Sed.","Unplanted Sed.","LPW Sed."))
        
    } else {
        # For the water datasets
        ps_mean_top$Group <- factor(ps_mean_top$Group,
                                    levels=c("Carex_OSPW","No_plant_OSPW","No_plant_Artificial_OSPW"),
                                    labels=c("Carex OSPW","Unplanted OSPW","LPW"))
        
    }
    
    
    # Put Other category and NA at the end of the legend
    ps_mean_top[[tax_lvl]] <- fct_relevel(ps_mean_top[[tax_lvl]], "Other", after = Inf)
    ps_mean_top[[tax_lvl]] <- fct_relevel(ps_mean_top[[tax_lvl]], "NA", after = Inf)
    
    levels(ps_mean_top[[tax_lvl]])
    
    # Bac Water
    if (file_name == "ps_16S_Wat.RData") {
        Class=c("#61a5c2","#003049","#d62828","#f77f00","#fcbf49",
                "#b8f2e6","#00a5cf","#7ae582","#aacc00","#dddf00",
                "#d38b5d","#e26d5c","#ffe1a8","#b56576",
                "#9d4edd","#9f86c0","#5a189a",
                "beige","grey40")
        
   # Class = c("midnightblue", "#776dd5","#1aaaf3","#4b2a12","#cc77dc",
   #           "#fcbf49", "#941e12","#6bde45","#FFF700","#781bd3",
   #           "#e4d999","turquoise4", "forestgreen", "#57d6e0","#7d7d99",
   #              "#81F2a3", "#8b328a","beige","grey40")
   # 
    
    } else if (file_name == "ps_16S_Sed.RData") {
        Class=c("#003049","#d62828","#f77f00","#fcbf49","#eae2b7",
                "#00a5cf","#7ae582","#2c5530","#739e82","#d38b5d",
                "#99621e","#b56576","#9d4edd","#5a189a",
                "beige","grey40")
        
        #99621e, #d38b5d, #f3ffb6, #739e82, #2c5530
       #Class = c("#776dd5","#1aaaf3","#4b2a12","#cc77dc","#cedf42",
       #          "#941e12","#6bde45","#1e2d21" ,"#438737","#e4d999",
       #          "#E28BA5", "#57d6e0","#7d7d99","#8b328a",
       #                     "beige","grey40")
    } else if (file_name == "ps_fungi.RData") {
        Class = c("#001219ff","#0a9396ff","#94d2bdff","#e9d8a6ff",
                  "#ee9b00ff","#ca6702ff","#a4161a","#613dc1",
                  "#8edf34","#80c423","#509724","#1dd3b0",
                  "beige","grey40")
        
    } else if (file_name == "ps_18S_Wat.RData"){
        Class = c("#eff7cf","#bad9b5","#aba361","#732c2c","#420c14",
                  "#444b6e","#9ab87a","#f8f991","#f9c74f","#720977",
                  "#277da1","#43aa8b",
                  "beige","grey40")
        
    } else print("Dataset not recognized")
    
    Barplot <- ggplot(ps_mean_top, aes(x = Group, y = Relative_Mean, fill = .data[[tax_lvl]])) +
        geom_bar(stat = "identity") +
        theme_bw() +
        facet_grid(~ Time , scales = "fixed")+
        theme(title = element_text(size = 22),
              #ggh4x.facet.nestline = element_line(color = "black"), 
              strip.text = element_text(face="bold",size=15),
              strip.background = element_rect(color="black",fill="white"),
              axis.title = element_text(size = 16),
              axis.text.x = element_text(size = 14, angle = 40,hjust=1),
              axis.text.y = element_text(size = 15, angle = 90, hjust = 0.5),
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 15, face = "italic"),
              legend.background = element_rect(color = "black"),
        )+
        guides(fill=guide_legend(ncol=1))+
        labs(fill = label, x="")+
        scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.01))) +
        scale_fill_manual(values = Class) +  # Use the dynamically selected palette
        ylab("Relative abundance")
    
    Barplot
    
    return(list(Barplot = Barplot))
    
}



Bac_solid <- Relative_Abundance("data/16s_sed/", "ps_16S_Sed.RData","Bacteria","Sediments")
Bac_water <- Relative_Abundance("data/16s_wat/", "ps_16S_Wat.RData","Bacteria","Water")
Fun_solid <- Relative_Abundance("data/its_sed/", "ps_fungi.RData","Fungi","Sediments")
Euk_water <- Relative_Abundance("data/18s_wat/", "ps_18S_Wat.RData","Eukaryotes","Water")


# Arrange all the figures in one plot
Relab_plot <- ggarrange(plot(Bac_solid$Barplot),
                        plot(Bac_water$Barplot),
                        plot(Fun_solid$Barplot),
                        plot(Euk_water$Barplot),labels="AUTO")

ggsave(Relab_plot,device='svg',height = 15, width = 15 ,
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/Relab_plot2.svg")




###____####

#### Figure 4  - Taxa correlation with final NA concentration####

Cor_taxa <- function(data_path, file_name,taxa,taxo_lvl) {
    
    load(here::here(data_path, file_name))
    
    metadata <- read.csv(file = paste0(data_path,"/NA_predict.csv"), 
                         dec = ".", header = T, row.names = 1, sep = ";", 
                         comment.char = "") 
    
    exclude_groups <-  c("Roots_Carex_OSPW",
                         #"Rhizosphere_Carex_OSPW",
                         "Sediments_Carex_OSPW",
                         "Sediments_No_plant_Artificial_OSPW",
                         "No_plant_Artificial_OSPW")
   
    ps_sub <- prune_samples(!(sample_data(ps)$Group %in% exclude_groups), ps)
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
    
    ps_rclr <- microbiome::transform(ps_sub, "rclr") 
    ps <- ps_rclr
    sample_data(ps) <- metadata
    
    # Get logical vector of taxa matching your criteria
    taxa_to_keep <- taxa_names(ps)[tax_table(ps)[, taxo_lvl] == taxa]
    
    # Subset phyloseq object
    ps_tax <- prune_taxa(taxa_to_keep, ps)
    ps_tax <- tax_glom(ps_tax, taxrank = taxo_lvl, NArm = F)
    
    plot_comp<- psmelt(ps_tax)
    
    plot_comp$Time <- factor(plot_comp$Time, levels = c("D-0","D8","D28","D42","D84"))
    
    plot_ps_mean <- plot_comp%>%
        group_by(Group, Time) %>%
        summarise(Mean_NA=mean(NA_concentration,na.rm = TRUE))
    
    plot_ps_mean <- merge(plot_comp,plot_ps_mean,by=c("Group","Time"))
    
    
    if(any(grepl("Rhizosphere", plot_ps_mean$Group))){
        label_1 <- "Carex rhizosphere"
        label_2 <- "Unplanted sediments"
    }else{
        label_1 <- "Carex OSPW"
        label_2 <- "Unplanted OSPW"
    }
    
    Cor_plot <- ggplot(plot_ps_mean, aes(x = Abundance, y = Final_NA_conc)) +
        theme_bw()+
        theme(axis.title=element_text(size=20),
              axis.text.x = element_text(size=18,angle=45,hjust=1.2,vjust=1.2),
              axis.text.y=element_text(size=18),
              legend.title = element_text(size = 22),
              legend.text = element_text(size = 20),
              legend.background = element_rect(color = "black"),
              legend.position = "bottom",
              strip.text = element_text(face="bold",size=20),
              strip.background = element_rect(color="black",fill="white"))+
        facet_grid(.~Time,scales="free")+
        geom_point(aes(fill=Group),size=2.2,shape=21)+
        stat_cor(label.y = 90,method="spearman",cor.coef.name="rho",label.sep="\n",size=7) +
        stat_regline_equation(label.y = -10,size=6.7)+
        geom_smooth(method = "lm", se = TRUE,color="black",linetype=5)+ 
        scale_fill_manual(labels= c(label_1,label_2),
                          values=c("#bad9a2", "#687967"),
                          name="Sample type :")+
        ylim(c(-20,110))+
        labs(title="\n",
             x=bquote(italic(.(taxa))~"- RCLR-transformed abundance"), 
             y = "[NAFC] (mg/L) at D84")
    
    return(list(Cor_plot = Cor_plot))
    
}

##### » Bac sed ####

#Clostridia <- Cor_taxa("data/16s_sed/", "ps_16S_Sed.RData","Clostridia","Class")
# Signif at D28

#Clostridiales <- Cor_taxa("data/16s_sed/", "ps_16S_Sed.RData","Clostridiales","Order")
# Non signif

Hungate <- Cor_taxa("data/16s_sed/", "ps_16S_Sed.RData","Hungateiclostridiaceae","Family")
# Non signif, almost at D8

#Christensenellaceae <- Cor_taxa("data/16s_sed/", "ps_16S_Sed.RData","Christensenellaceae","Family")
# Signif at D28

#Clostridiaceae <- Cor_taxa("data/16s_sed/", "ps_16S_Sed.RData","Clostridiaceae","Family")
# Non-signif

Dechlo <- Cor_taxa("data/16s_sed/", "ps_16S_Sed.RData","Dechloromonas","Genus")
# signif D8 - D28 - D42

#Anaero <- Cor_taxa("data/16s_sed/", "ps_16S_Sed.RData","Anaerovorax","Genus")
# Non-signif


##### » Fungi sed ####
Morti <- Cor_taxa("data/its_sed/", "ps_fungi.RData","Mortierellaceae","Family")
# signif D28 - D42

Didymellaceae <- Cor_taxa("data/its_sed/", "ps_fungi.RData","Didymellaceae","Family")
# signif D42

#Calophoma <- Cor_taxa("data/its_sed/", "ps_fungi.RData","Calophoma","Genus")
# signif D42

#Oidiodendron <- Cor_taxa("data/its_sed/", "ps_fungi.RData","Oidiodendron","Genus")
# signif D28

#Leptodontidium <- Cor_taxa("data/its_sed/", "ps_fungi.RData","Leptodontidium","Genus")
# non-signif

#Pezoloma <- Cor_taxa("data/its_sed/", "ps_fungi.RData","Pezoloma","Genus")
#Phialocephala <- Cor_taxa("data/its_sed/", "ps_fungi.RData","Phialocephala","Genus")
#Oculimacula <- Cor_taxa("data/its_sed/", "ps_fungi.RData","Oculimacula","Genus")
#Clathrosphaerina <- Cor_taxa("data/its_sed/", "ps_fungi.RData","Clathrosphaerina","Genus")

#Oculimacula <- Cor_taxa("data/its_sed/", "ps_fungi.RData","Oculimacula","Genus")


##### » Bac wat ####

Gemmataceae <- Cor_taxa("data/16s_wat/", "ps_16S_Wat.RData","Gemmataceae","Family")
#D-0 and D8

Isosphaeraceae <- Cor_taxa("data/16s_wat/", "ps_16S_Wat.RData","Isosphaeraceae","Family")
# D8

# Arrange all the figures in one plot
Bac_sed_cor <- ggarrange(plot(Hungate$Cor_plot),plot(Dechlo$Cor_plot),
                         common.legend = T,
                         labels=c("A","B"),
                         font.label = list(size=20),
                         legend = "bottom")

Bac_wat_cor <- ggarrange(plot(Gemmataceae$Cor_plot),plot(Isosphaeraceae$Cor_plot),
                         common.legend = T,
                         labels=c("C","D"),
                         font.label = list(size=20),
                         legend = "bottom")

Fun_sed_cor <- ggarrange(plot(Morti$Cor_plot),plot(Didymellaceae$Cor_plot),
                         common.legend = T,
                         labels=c("E","F"),
                         font.label = list(size=20),
                         legend = "bottom")

Figure_4 <- ggarrange(Bac_sed_cor,
                      Bac_wat_cor,
                      Fun_sed_cor,
                      nrow=3)


ggsave(Figure_4,device='svg',height = 20, width = 20 ,
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/Cor_taxa_v2.svg")



###____####
#### Figure 5 - Community Change (D8 vs D-0) ####

source(here::here("src", "libraries.R"))
data_path <- here("data/16s_sed/")
load(here(data_path,"ps_16S_Sed.RData"))

data_path <- here("data/16s_wat/")
load(here(data_path,"ps_16S_Wat.RData"))

data_path <- here("data/its_sed/")
load(here(data_path,"ps_fungi.RData"))


Community_change <- function(data_path, file_name,tax_lvl,Stype) {

load(here::here(data_path, file_name))
    
ps_glom <- tax_glom(ps,taxrank=tax_lvl,NArm = F)

Sdate1 <- "D-0"# Analogous to control (Will be on the left side of the volcano)
Sdate2 <- "D8"# Analogous to treatment (Will be on the right side of the volcano)

ps_comp <- prune_samples(sample_data(ps_glom)$Group == Stype & sample_data(ps_glom)$Time%in%c(Sdate1,Sdate2),ps_glom)
ps_comp <- prune_taxa(taxa_sums(ps_comp)>0, ps_comp)

#Keep only those taxa that have at least some non-zero abundance in at least 75% of the samples.
ps_comp <- prune_taxa(rowSums(otu_table(ps_comp) == 0) 
                      < ncol(otu_table(ps_comp)) * 0.75, ps_comp)

# ps_comp_long <- psmelt(ps_comp)
# ps_comp_long1 <-select(ps_comp_long ,Abundance,Time,Group,Mesocosm,tax_lvl)


diagdds = phyloseq_to_deseq2(ps_comp, ~ Time)
diagdds$Time <- relevel(diagdds$Time, ref = Sdate1)
# Sdate1 as the reference 

diagdds = estimateSizeFactors(diagdds, type = "poscounts")

diagdds = DESeq(diagdds,
                test = "Wald",
                fitType="local")
plotDispEsts(diagdds) # local has a better fit


resultsNames(diagdds)
res <- results(diagdds, contrast=c("Time", Sdate2, Sdate1), alpha=0.01)
# This means that 
# Right side of the volcano plot (positive log2FoldChange): Higher abundance in Group2
# Left side of the volcano plot (negative log2FoldChange): Higher abundance in Group1

plotMA(res)
res2 = res[order(res$padj, na.last=NA), ] # ordered and removes the ASVs which have padj = NAs 
plotMA(res2)


# Volcano plots
res2 <- as.data.frame(res2)
tax <- data.frame(tax_table(ps_glom))

# Add taxonomy to each taxa scores
res2_taxa <- merge(res2,tax,by='row.names')

res2_taxa$neg_log10_pvalue <- -log10(res2_taxa$padj)

# Following function allows to find the p_value cutoff 
# So that the chances of finding 
find_power_of_ten <- function(x) {
    y <- 0
    while (dim(x)[1] / 10^y >= 1) {
        y <- y + 1
    }
    return(y)
}
neg_log10_pval <- find_power_of_ten(res2_taxa)

# Filter signif_taxa based on desired cutoffs 
filtered_signif_taxa <- res2_taxa[(abs(res2_taxa$log2FoldChange) > 1) & 
                                      (res2_taxa$neg_log10_pvalue >= neg_log10_pval), ]


# Get labels
labels <- res2_taxa[[tax_lvl]]
labels_noNA <- labels[labels!="NA"]


if (file_name=="ps_16S_Sed.RData"){ 
    if (Stype=="Rhizosphere_Carex_OSPW"){
    Fig_title <- "Carex rhizosphere - Bacteria"
    
        } else if(Stype=="Sediments_Carex_OSPW"){
            Fig_title <- "Carex sediments - Bacteria"
            
        } else if(Stype=="Sediments_No_plant_OSPW"){
        Fig_title <- "Unplanted sediments - Bacteria"}
    
} else if(file_name=="ps_fungi.RData"){ 
    if (Stype=="Rhizosphere_Carex_OSPW"){
        Fig_title <- "Carex rhizosphere - Fungi"
        
    } else if(Stype=="Sediments_Carex_OSPW"){
        Fig_title <- "Carex sediments - Fungi"
        
    } else if(Stype=="Sediments_No_plant_OSPW"){
        Fig_title <- "Unplanted sediments - Fungi"}
        
} else if(file_name=="ps_16S_Wat.RData"){ 
    if (Stype=="Carex_OSPW"){
        Fig_title <- "Carex OSPW - Bacteria"
        
    } else if(Stype=="No_plant_OSPW"){
        Fig_title <- "Unplanted OSPW - Bacteria"}
}
    
# Volcano plot with dynamic labeling / text and p_value cutoff
volcano_plot<- EnhancedVolcano(res2_taxa, 
                               lab = labels, 
                               selectLab = labels_noNA,
                               x = 'log2FoldChange', 
                               y = "padj",
                               pCutoff =  10^(-neg_log10_pval),
                               title = Fig_title,
                               titleLabSize = 18,
                               subtitle = bquote(.(Sdate2) ~italic("vs") ~ .(Sdate1)   ~ .(paste("-",tax_lvl,"level"))),
                               subtitleLabSize = 15,
                               pointSize = 3,
                               labSize = 5.5,
                               labFace = "italic",
                               colAlpha = 0.5,
                               legendPosition = 'bottom',
                               legendLabSize = 14,
                               legendIconSize = 4.0,
                               drawConnectors = TRUE,
                               widthConnectors = 0.6,
                               max.overlaps = 10,
                               border="full")+
    theme(legend.box.spacing = unit(0, "pt"),
          legend.text = element_text(size=18),
          legend.margin = margin(10, 0, 0, 0),
          plot.caption = element_text(hjust=1),
          plot.title = element_text(face="plain",margin=margin(10,0,4,0)),
          plot.subtitle = element_text(margin=margin(0,0,3,0)))


}

### Bac in solid samples  ###
Bac_rhizo <- Community_change("data/16s_sed/", "ps_16S_Sed.RData",
                              tax_lvl = "Family",Stype = "Rhizosphere_Carex_OSPW")

Bac_NoPlant <- Community_change("data/16s_sed/", "ps_16S_Sed.RData",
                                tax_lvl = "Family",Stype = "Sediments_No_plant_OSPW")


### Bac in water column  ###
Bac_Carex_OSPW <- Community_change("data/16s_wat/", "ps_16S_Wat.RData",
                              tax_lvl = "Family",Stype = "Carex_OSPW")

Bac_NoPlant_OSPW <- Community_change("data/16s_wat/", "ps_16S_Wat.RData",
                              tax_lvl = "Family",Stype = "No_plant_OSPW")

### Fungi in solid samples ###
Fun_rhizo <- Community_change("data/its_sed/", "ps_fungi.RData",
                              tax_lvl = "Family",Stype = "Rhizosphere_Carex_OSPW")

Fun_NoPlant <- Community_change("data/its_sed/", "ps_fungi.RData",
                              tax_lvl = "Family",Stype = "Sediments_No_plant_OSPW")





Primary_response <- ggarrange(Bac_rhizo,Bac_NoPlant,
                              Bac_Carex_OSPW,Bac_NoPlant_OSPW,
                              Fun_rhizo,Fun_NoPlant,
                              labels = "AUTO",
                              ncol=2,nrow=3)

ggsave(Primary_response,device='svg',height = 19, width = 15,
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/AAA_community_change.svg")



Primary_response2 <- ggarrange(Bac_rhizo,Bac_NoPlant,
                               Bac_Carex_OSPW,Bac_NoPlant_OSPW,
                               labels = "AUTO",
                               ncol=2,nrow=2)

ggsave(Primary_response2,device='svg',height = 15, width = 15,
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/AAB_community_change.svg")


####____####
#### SUPPLEMENTARY FIGURES #####
####____####



#### Figure S1 to S8 - Alpha div ####
Alpha.div.Stype <- function(data_path, file_name, rarefaction=FALSE) {
    
    load(here::here(data_path, file_name))
    
    exclude_groups <-  c("Sediments_No_plant_Artificial_OSPW",
                         "No_plant_Artificial_OSPW")
    
    ps_sub <- prune_samples(!(sample_data(ps)$Group %in% exclude_groups), ps)
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
    
    if (rarefaction){
        set.seed(12)
        ps_sub <- rarefy_even_depth(ps_sub,
                                    sample.size = min(sample_sums(ps)),
                                    rngseed = FALSE, replace = FALSE, 
                                    trimOTUs = TRUE, verbose = TRUE)
        } else {
            ps_sub <- ps_sub
        }
    alpha.div<- plot_richness(ps_sub, x="Group", measures=c("Observed","Simpson", "Shannon"))
        
    plot <- select(alpha.div$data,c("samples","Time","Group","variable","value"))
    plot$Time <- factor(plot$Time, levels = c("D-0","D8","D28","D42","D84"))
        
    if(any(grepl("Rhizosphere", plot$Group))){
        plot$Group <- factor(plot$Group,levels=c("Roots_Carex_OSPW","Rhizosphere_Carex_OSPW","Sediments_Carex_OSPW","Sediments_No_plant_OSPW","Sediments_No_plant_Artificial_OSPW"),
                             labels=c("Carex roots","Carex rhizosphere","Carex sediments","Unplanted sediments","LPW sediments"))
        palette <- c("#9DEE74","#bad9a2", "#80B39F", "#687967", "#7CCAEC")
    }else{plot$Group <- factor(plot$Group, levels = c("Carex_OSPW","No_plant_OSPW","No_plant_Artificial_OSPW"), labels = c("Carex OSPW","Unplanted OSPW","Unplanted LPW"))
          palette <- c("#9DEE74","#687967", "#7CCAEC")
        }

Alpha_div <- ggplot(plot,aes(x=variable,y=value,color=Group))+
    geom_boxplot(outlier.shape = NA,size = 1.3,alpha=0.5,show.legend = T)+
    theme_bw()+
    facet_wrap(~variable,scales = "free")+
    theme(#title = element_text(size=20),
    axis.text.x =element_blank(),
    axis.text.y =element_text(size=12),
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size=15),
    legend.text.align = 0,
    axis.ticks.x = element_blank(),
    legend.text=element_text(size=15),
    legend.title=element_text(size=20),
    legend.position = "bottom",
    legend.box="vertical",
    legend.spacing = unit(0.01,"cm"),
    strip.text = element_text(size=15),
    strip.background = element_rect(fill='white'))+
    geom_pwc(method="dunn_test",hide.ns=T,
           p.adjust.method="bonferroni",
           label="p.adj.format",
           label.size = 4,
           position="identity",
           step.increase = 0.1,
           linetype = 1,
           size = 1)+
    scale_color_manual(values=palette,name=element_blank())+
    ylab("Alpha diversity score")
    guides(color=guide_legend(order = 1))
            
                
    Alpha_stat <- Alpha_div$data %>%
        group_by(Group,variable) %>% 
        summarize(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE))%>% 
        arrange(by=variable)
    
    table <- as.data.frame(Alpha_stat)
    table[c("mean", "sd")] <- lapply(table[c("mean", "sd")], round, digits = 3)    
    table <- ggtexttable(table,theme = ttheme(base_style = "light",
                                               colnames.style = colnames_style(color = "black", fill = "white"),
                                               tbody.style = tbody_style(color = "black", fill = "white"),
                                               rownames.style = rownames_style(fill="white")))
                                   
Alpha_div_plot <- ggarrange(Alpha_div,table,
                            widths=c(1.5,1),
                            ncol=2)

}
Sed_Bac_nonRar <- Alpha.div.Stype("data/16s_sed/", "ps_16S_Sed.RData",rarefaction=F)
Sed_Bac_Rar <- Alpha.div.Stype("data/16s_sed/", "ps_16S_Sed.RData",rarefaction=T)

ggsave(Sed_Bac_nonRar,device='svg',height = 8, width = 12 ,
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/Alpha_Bac_Sed.svg")

Sed_Fun_nonRar <- Alpha.div.Stype("data/its_sed/", "ps_fungi.RData",rarefaction=F)
Sed_Fun_Rar <- Alpha.div.Stype("data/its_sed/", "ps_fungi.RData",rarefaction=T)

ggsave(Sed_Fun_nonRar,device='svg',height = 8, width = 12,
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/Alpha_Fun_Sed.svg")

Wat_Bac_nonRar <- Alpha.div.Stype("data/16s_wat/", "ps_16S_Wat.RData",rarefaction=F)
Wat_Bac_Rar <- Alpha.div.Stype("data/16s_wat/", "ps_16S_Wat.RData",rarefaction=T)

ggsave(Wat_Bac_nonRar,device='svg',height = 8, width = 12,
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/Alpha_Bac_Wat.svg")


Wat_Euk_nonRar <- Alpha.div.Stype("data/18s_wat/", "ps_18S_Wat.RData",rarefaction=F)
Wat_Euk_Rar <- Alpha.div.Stype("data/18s_wat/", "ps_18S_Wat.RData",rarefaction=T)

ggsave(Wat_Euk_nonRar,device='svg',height = 8, width = 12,
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/Alpha_Euk_Wat.svg")


Alpha.div.Time <- function(data_path, file_name, rarefaction=FALSE) {
    
    load(here::here(data_path, file_name))
    
    exclude_groups <-  c("Sediments_No_plant_Artificial_OSPW",
                         "No_plant_Artificial_OSPW")
    
    ps_sub <- prune_samples(!(sample_data(ps)$Group %in% exclude_groups), ps)
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
    
    if (rarefaction){
        set.seed(12)
        ps_sub <- rarefy_even_depth(ps_sub,
                                    sample.size = min(sample_sums(ps)),
                                    rngseed = FALSE, replace = FALSE, 
                                    trimOTUs = TRUE, verbose = TRUE)
    } else {
        ps_sub <- ps_sub
    }
    alpha.div<- plot_richness(ps_sub, x="Group", measures=c("Observed","Simpson", "Shannon"))
    
    plot <- select(alpha.div$data,c("samples","Time","Group","variable","value"))
    plot$Time <- factor(plot$Time, levels = c("D-0","D8","D28","D42","D84"))
    
    if(any(grepl("Rhizosphere", plot$Group))){
        plot$Group <- factor(plot$Group,levels=c("Roots_Carex_OSPW","Rhizosphere_Carex_OSPW","Sediments_Carex_OSPW","Sediments_No_plant_OSPW","Sediments_No_plant_Artificial_OSPW"),
                             labels=c("Carex roots","Carex rhizosphere","Carex sediments","Unplanted sediments","LPW sediments"))
        }else{plot$Group <- factor(plot$Group, levels = c("Carex_OSPW","No_plant_OSPW","No_plant_Artificial_OSPW"), labels = c("Carex OSPW","Unplanted OSPW","Unplanted LPW"))
        }
    
    palette <- c("#0D0887FF","#7E03A8FF", "#CC4678FF", "#F89441FF", "#F0F921FF")
    
    Alpha_div <- ggplot(plot,aes(x=Time,y=value,color=Time))+
        #geom_boxplot(aes(x=variable,y=value,color=Group,fill=Time),position=position_dodge(width=0.75),size=0.5,outlier.shape = NA)+
        #geom_point(aes(x=variable,y=value,color=Group,shape=Time,fill=Time,stroke = 0.7),alpha=0.3,position=position_dodge(width=0.75,preserve = "total"),size=2.5)+
        geom_boxplot(outlier.shape = NA,size = 1,alpha=0.5,show.legend = T)+
        theme_bw()+
        facet_grid(variable ~ Group,scales = "free_y")+
        theme(#title = element_text(size=20),
            axis.text.x =element_blank(),
            axis.text.y =element_text(size=12),
            axis.title.x = element_blank(), 
            axis.title.y = element_text(size=15),
            legend.text.align = 0,
            axis.ticks.x = element_blank(),
            legend.text=element_text(size=15),
            legend.title=element_text(size=20),
            legend.position = "bottom",
            legend.box="vertical",
            legend.spacing = unit(0.01,"cm"),
            strip.text = element_text(size=15),
            strip.background = element_rect(fill='white'))+
        geom_pwc(method="dunn_test",hide.ns=T,
                 p.adjust.method="bonferroni",
                 label="p.adj.format",
                 label.size = 4,
                 #tip.length = 0, 
                 position="identity",
                 step.increase = 0.1,
                 linetype = 1,
                 size = 1)+
        scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
        scale_color_manual(values=palette,name=element_blank())+
        ylab("Alpha diversity score")+
        guides(color=guide_legend(order = 1))
    Alpha_div
    
    
    Alpha_stat <- Alpha_div$data %>%
        group_by(Time,variable) %>% 
        summarize(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE))%>% 
        arrange(by=variable)
    
    table <- as.data.frame(Alpha_stat)
    table[c("mean", "sd")] <- lapply(table[c("mean", "sd")], round, digits = 3)    
    table <- ggtexttable(table,theme = ttheme(base_style = "light",
                                              colnames.style = colnames_style(color = "black", fill = "white"),
                                              tbody.style = tbody_style(color = "black", fill = "white"),
                                              rownames.style = rownames_style(fill="white")))
    
    Alpha_div_plot <- ggarrange(Alpha_div,table,
                                widths=c(2.5,1),
                                ncol=2)
    
}


Sed_Bac_nonRar <- Alpha.div.Time("data/16s_sed/", "ps_16S_Sed.RData",rarefaction=F)
Sed_Bac_Rar <- Alpha.div.Time("data/16s_sed/", "ps_16S_Sed.RData",rarefaction=T)

ggsave(Sed_Bac_nonRar,device='svg',height = 8, width = 12 ,
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/Alpha_Bac_Sed_Time.svg")

Sed_Fun_nonRar <- Alpha.div.Time("data/its_sed/", "ps_fungi.RData",rarefaction=F)
Sed_Fun_Rar <- Alpha.div.Time("data/its_sed/", "ps_fungi.RData",rarefaction=T)

ggsave(Sed_Fun_nonRar,device='svg',height = 8, width = 12,
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/Alpha_Fun_Sed_Time.svg")

Wat_Bac_nonRar <- Alpha.div.Time("data/16s_wat/", "ps_16S_Wat.RData",rarefaction=F)
Wat_Bac_Rar <- Alpha.div.Time("data/16s_wat/", "ps_16S_Wat.RData",rarefaction=T)

ggsave(Wat_Bac_nonRar,device='svg',height = 8, width = 12,
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/Alpha_Bac_Wat_Time.svg")


Wat_Euk_nonRar <- Alpha.div.Time("data/18s_wat/", "ps_18S_Wat.RData",rarefaction=F)
Wat_Euk_Rar <- Alpha.div.Time("data/18s_wat/", "ps_18S_Wat.RData",rarefaction=T)

ggsave(Wat_Euk_nonRar,device='svg',height = 8, width = 12,
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/Alpha_Euk_Wat_Time.svg")


#### Figure S9 and S10 - Beta div with roots ####

Beta_div_Ordination <- function(data_path, file_name) {
    
    # Load phyloseq object
    load(here::here(data_path, file_name))
    
    # Data processing pipeline
    exclude_groups <-  c("Sediments_No_plant_Artificial_OSPW",
                         "No_plant_Artificial_OSPW")
    # Groups that are not of interest for the plot are removed 
    ps_sub <- prune_samples(!(sample_data(ps)$Group %in% exclude_groups), ps)
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
    # 
    
    # Ordination analysis
    ps_dist <- microbiome::transform(ps_sub, "rclr")
    ord_rclr <- phyloseq::ordinate(ps_dist, "PCoA", distance = "euclidean")
    
    # Calculate variance explained
    PC1 <- round(ord_rclr$values$Relative_eig[1] * 100, 1)
    PC2 <- round(ord_rclr$values$Relative_eig[2] * 100, 1)
    
    # Create ordination dataframe
    Ordination <- phyloseq::plot_ordination(ps_dist, ord_rclr, type = "samples", justDF = TRUE)
    Ordination$Time <- factor(Ordination$Time, levels = c("D-0", "D8", "D28", "D42", "D84"))

    
    palette <- c("#9DEE74","#bad9a2", "#80B39F", "#687967")
    Ordination$Group <- factor(Ordination$Group,
                                   levels = c("Roots_Carex_OSPW","Rhizosphere_Carex_OSPW", "Sediments_Carex_OSPW",
                                              "Sediments_No_plant_OSPW"),
                                   labels = c("Carex roots","Carex rhizosphere", "Carex sediments", 
                                              "Unplanted sediments"))
    
    # Create plots
    plot_stype <- ggplot(Ordination, aes(x = Axis.1, y = Axis.2, color = Group)) +
        geom_point(size = 4, aes(stroke = 3)) +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.box = "vertical",
              legend.spacing = unit(0.01, "cm"),
              text = element_text(size = 25)) +
        xlab(paste0("PCo1 (", PC1, "%)")) + 
        ylab(paste0("PCo2 (", PC2, "%)")) + 
        labs(title = "" ) + #Add space above graph
        scale_color_manual(values = palette, name = "") +
        guides(color = guide_legend(nrow = 1))
    
    plot_time <- ggplot(Ordination, aes(x = Axis.1, y = Axis.2, color = Time, shape = Group)) +
        geom_point(size = 4, aes(stroke = 3)) +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.box = "vertical",
              legend.spacing = unit(0.01, "cm"),
              text = element_text(size = 25)) +
        xlab(paste0("PCo1 (", PC1, "%)")) + 
        ylab(paste0("PCo2 (", PC2, "%)")) + 
        labs(title = "" ) + #Add space above graph
        scale_shape_manual(values = c(1, 2, 3, 4), name = "") +
        scale_color_manual(values = c("#0D0887FF", "#7E03A8FF", "#CC4678FF", 
                                      "#F89441FF", "#F0F921FF"), name = "") +
        guides(shape = guide_legend(nrow = 1))
    
    return(list(plot_stype = plot_stype, plot_time = plot_time))
    
    
}

source(here::here("src", "libraries.R"))

# Process all datasets
Bac_solid <- Beta_div_Ordination("data/16s_sed/", "ps_16S_Sed.RData")
Fun_solid <- Beta_div_Ordination("data/its_sed/", "ps_fungi.RData")

Bac_solid$plot_stype

ggsave(Bac_solid$plot_stype,device='svg',height = 10, width = 18 ,
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/Beta_Bac_roots.svg")


Fun_solid$plot_stype


ggsave(Fun_solid$plot_stype,device='svg',height = 10, width = 18 ,
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/Beta_Fun_roots.svg")



#### Figure S11 to S13 - Diff ab per time  ####
#### 
#### 

#Group1 <- "Sediments_No_plant_OSPW" # Analogous to control (Will be on the left side of the volcano)
#Group2 <- "Rhizosphere_Carex_OSPW" # Analogous to treatment (Will be on the right side of the volcano)


Change.over.time <- function(data_path, file_name,Group1,Group2) {

load(here::here(data_path, file_name))


metadata <- sample_data(ps)
metadata$Time <- factor(metadata$Time, levels = c("D-0", "D8", "D28", "D42", "D84"))
plot_list <- list()
i <- 
for (i in levels(metadata$Time)){

ps_comp <- prune_samples(sample_data(ps)$Time == i & 
                                 sample_data(ps)$Group%in%c(Group1,Group2),ps)
    
#ps_comp<-  subset_samples(ps,Time==i & 
#                              (Group==Group1 | Group ==Group2))

ps_comp <- prune_taxa(taxa_sums(ps_comp)>0, ps_comp)
ps_comp <- prune_taxa(rowSums(otu_table(ps_comp) == 0) 
                      < ncol(otu_table(ps_comp)) * 0.75, ps_comp)
#Keeps only those taxa that have at least some non-zero abundance in at least 75% of the samples.

#ps_comp_long <- psmelt(ps_comp)
#ps_comp_long1 <-select(ps_comp_long ,Abundance,Time,Group,tax_lvl)


diagdds = phyloseq_to_deseq2(ps_comp, ~ Group)
levels(diagdds$Group)  # This will show you the levels and their order
#levels(diagdds$Time)
diagdds$Group <- relevel(diagdds$Group, ref = Group1)


diagdds = estimateSizeFactors(diagdds, type = "poscounts")

diagdds = DESeq(diagdds,
                test = "Wald",
                fitType="local")
plotDispEsts(diagdds) # local has a better fit


resultsNames(diagdds)
res <- results(diagdds, contrast=c("Group", Group2, Group1), alpha=0.01)
# This means that 
# Right side of the volcano plot (positive log2FoldChange): Higher abundance in Group2
# Left side of the volcano plot (negative log2FoldChange): Higher abundance in Group1

plotMA(res)

res2 = res[order(res$padj, na.last=NA), ] # ordered and removes the ASVs which have padj = NAs 
plotMA(res2)


#dds1 <- DESeq(phyloseq_to_deseq2(ps_comp, design = ~ Group + Time))
#dds2 <- DESeq(phyloseq_to_deseq2(ps_comp, design = ~ Group * Time))
#results_lrt <- results(dds2, test = "LRT", reduced = ~ Group + Time)

# Volcano plots
res2 <- as.data.frame(res2)
tax <- data.frame(tax_table(ps))
res2_taxa <- merge(res2,tax,by='row.names')

res2_taxa$neg_log10_pvalue <- -log10(res2_taxa$padj)

# This function allows to find the p_value cutoff 
# So that the chances of finding 
find_power_of_ten <- function(x) {
    y <- 0
    while (dim(x)[1] / 10^y >= 1) {
        y <- y + 1
    }
    return(y)
}
neg_log10_pval <- find_power_of_ten(res2_taxa)

# Filter signif_taxa based on cutoffs to see the output
filtered_signif_taxa <- res2_taxa[(abs(res2_taxa$log2FoldChange) > 1) & (res2_taxa$neg_log10_pvalue >= neg_log10_pval), ]


if (Group1=="Sediments_No_plant_OSPW"){
    Label1 <- "Unplanted sediments"
} else if(Group1=="No_plant_OSPW"){
    Label1 <- "Unplanted OSPW"}
    
if(Group2=="Rhizosphere_Carex_OSPW"){
    Label2 <- "Carex rhizosphere"
} else if (Group2=="Carex_OSPW"){
    Label2 <- "Carex OSPW"}

# Volcano plot with dynamic labeling / text and p_value cutoff
volcano_plot<- EnhancedVolcano(res2_taxa, 
                               lab=res2_taxa$Row.names,
                               selectLab=res2_taxa$Row.names,
                               x = 'log2FoldChange', 
                               y = "padj",
                               pCutoff =  10^(-neg_log10_pval),
                               title = i,
                               titleLabSize = 18,
                               subtitle = bquote(.(Label2) ~italic("vs") ~ .(Label1)),
                               subtitleLabSize = 15,
                               pointSize = 3,
                               labSize = 5.5,
                               labFace = "italic",
                               colAlpha = 0.5,
                               legendPosition = 'bottom',
                               legendLabSize = 14,
                               legendIconSize = 4.0,
                               drawConnectors = TRUE,
                               widthConnectors = 0.6,
                               max.overlaps = 10,
                               border="full")+
    theme(legend.box.spacing = unit(0, "pt"),
          legend.text = element_text(size=18),
          legend.margin = margin(10, 0, 0, 0),
          plot.caption = element_text(hjust=1),
          plot.title = element_text(face="plain",margin=margin(10,0,4,0)),
          plot.subtitle = element_text(margin=margin(0,0,3,0)))

plot_list[[i]] <- volcano_plot
}
return(plot_list)
}

BacSolidCompTime <- Change.over.time(data_path="data/16s_sed/",
                                     file_name="ps_16S_Sed.RData",
                                     Group1="Sediments_No_plant_OSPW",
                                     Group2="Rhizosphere_Carex_OSPW")

Bac_solid_Change_Time <- ggarrange(BacSolidCompTime$`D-0`,BacSolidCompTime$D8,
                                   BacSolidCompTime$D28,BacSolidCompTime$D42,
                                   BacSolidCompTime$D84,ncol=2,nrow=3)
     

ggsave(Bac_solid_Change_Time,device='svg', height = 15, width = 20, 
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/Change_time_bac_solid.svg")




BacWaterCompTime <- Change.over.time(data_path="data/16s_wat/",
                                     file_name="ps_16S_Wat.RData",
                                     Group1="No_plant_OSPW",
                                     Group2="Carex_OSPW")

Bac_Water_Change_Time <- ggarrange(BacWaterCompTime$`D-0`,BacWaterCompTime$D8,
                                   BacWaterCompTime$D28,BacWaterCompTime$D42,
                                   BacWaterCompTime$D84,ncol=2,nrow=3)


ggsave(Bac_Water_Change_Time,device='svg', height = 15, width = 20, 
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/Change_time_bac_Water.svg")





FunSolidCompTime <- Change.over.time(data_path="data/its_sed/",
                                     file_name="ps_fungi.RData",
                                     Group1="Sediments_No_plant_OSPW",
                                     Group2="Rhizosphere_Carex_OSPW")

Fun_solid_Change_Time <- ggarrange(FunSolidCompTime$`D-0`,FunSolidCompTime$D8,
                                   FunSolidCompTime$D28,FunSolidCompTime$D42,
                                   FunSolidCompTime$D84,ncol=2,nrow=3)


ggsave(Fun_solid_Change_Time,device='svg', height = 15, width = 20, 
       filename="/Users/Simon/OneDrive - INRS/Documents/INRS/Research_projects/Projet_GROW/Meso1_metabarcoding/Manuscripts/Figures/Change_time_Fun_solid.svg")







diff_ab_ASVs <- paste0("ASV",c("104",
                               "103","82","80",
                               "4","122","47","139","146",
                               "48","90","21","189"))
diff_ab_tax <- subset(tax,row.names(tax)%in%diff_ab_ASVs)




#### Extra Figure - Root RelAb ####

Root_Relab_Abundance <- function(data_path, file_name,label,dataset) {
    
    # Load phyloseq object
    load(here::here(data_path, file_name))
    
    # Data processing pipeline
    Root_group <-  "Roots_Carex_OSPW"
    
    # Groups that are not of interest for the plot are removed 
    ps_sub <- prune_samples((sample_data(ps)$Group %in% Root_group), ps)
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
    
    # Select taxonomic level to group by
    tax_lvl <- "Class"
    ps_glom <- tax_glom(ps_sub, taxrank = tax_lvl)
    
    # Transform to long format for plotting purposes
    ps_long <- psmelt(ps_glom)
    
    # Group by different columns and compute the mean and relative mean (mean across replicates)
    ps_mean <-  ps_long %>%
        group_by(Group, Time,Compartment,Water_type, Sample_type, !!sym(tax_lvl)) %>%
        summarise(Mean = mean(Abundance, na.rm = TRUE))%>%
        mutate(Relative_Mean = Mean / sum(Mean))
    
        Top <- ps_mean %>%
            group_by(Group, Time, !!sym(tax_lvl)) %>%
            arrange(Group, Time, desc(Relative_Mean)) %>%
            group_by(Group, Time) %>%
            slice_head(n = 10) %>%
            pull(!!sym(tax_lvl)) %>%
            unique()
    
    # Classify all the other taxa that remain that are not in the top10 nor NA as "Other"
    ps_mean_top <-  ps_mean %>%
        group_by(Group, Time,Compartment,Water_type, Sample_type, !!sym(tax_lvl)) %>%
        mutate(!!sym(tax_lvl) := ifelse(!!sym(tax_lvl) %in% Top | !!sym(tax_lvl) == "NA", !!sym(tax_lvl),"Other"))
    
    
    # Rename categories for plotting purposes
    
    ps_mean_top$Sample_type <- factor(ps_mean_top$Sample_type, 
                                      levels = c( "Carex","No_plant"), 
                                      labels = c("Carex","Unplanted"))
    
    ps_mean_top$Water_type <- factor(ps_mean_top$Water_type, 
                                     levels = c( "OSPW","Artificial_OSPW"), 
                                     labels = c("OSPW","LPW"))
    
    if (dataset == "Sediments") {
        # For the sediments datasets
        ps_mean_top$Group <- factor(ps_mean_top$Group,
                                    levels=c("Roots_Carex_OSPW","Rhizosphere_Carex_OSPW","Sediments_Carex_OSPW","Sediments_No_plant_OSPW","Sediments_No_plant_Artificial_OSPW"),
                                    labels=c("Carex Ro.","Carex Rh.","Carex Sed.","Unplanted Sed.","LPW Sed."))
        
    } else {
        # For the water datasets
        ps_mean_top$Group <- factor(ps_mean_top$Group,
                                    levels=c("Carex_OSPW","No_plant_OSPW","No_plant_Artificial_OSPW"),
                                    labels=c("Carex OSPW","Unplanted OSPW","LPW"))
        
    }
    
    
    # Put Other category and NA at the end of the legend
    ps_mean_top[[tax_lvl]] <- fct_relevel(ps_mean_top[[tax_lvl]], "Other", after = Inf)
    ps_mean_top[[tax_lvl]] <- fct_relevel(ps_mean_top[[tax_lvl]], "NA", after = Inf)
    
    levels(ps_mean_top[[tax_lvl]])
    
    # Bac Water
    if (file_name == "ps_16S_Sed.RData") {
        Class=c("#003049","#d62828","#f77f00","#fcbf49","#eae2b7",
                "#00a5cf","#7ae582","#aacc00","#2c5530","#BED981",
                "#739e82","#d38b5d","#99621e","#F9DCEB",
                "#b56576","#9f86c0","#5a189a",
                "beige","grey40")
        

 
    } else if (file_name == "ps_fungi.RData") {
        Class = c("#001219ff","#0a9396ff","#94d2bdff","#e9d8a6ff",
                  "#ee9b00ff","#ca6702ff","#a4161a","#520005","#613dc1",
                  "#8A6A86","#B7ABED",
                  "#8edf34","#80c423","#509724","#1dd3b0",
                  "beige","grey40")
        
        
    } else print("Dataset not recognized")
    
    Barplot <- ggplot(ps_mean_top, aes(x = Group, y = Relative_Mean, fill = .data[[tax_lvl]])) +
        geom_bar(stat = "identity") +
        theme_bw() +
        facet_grid(~ Time , scales = "fixed")+
        theme(title = element_text(size = 22),
              #ggh4x.facet.nestline = element_line(color = "black"), 
              strip.text = element_text(face="bold",size=15),
              strip.background = element_rect(color="black",fill="white"),
              axis.title = element_text(size = 16),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 15, angle = 90, hjust = 0.5),
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 15, face = "italic"),
              legend.background = element_rect(color = "black"),
        )+
        guides(fill=guide_legend(ncol=1))+
        labs(fill = label, x="")+
        scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.01))) +
        scale_fill_manual(values = Class) +  # Use the dynamically selected palette
        ylab("Relative abundance")
    
    Barplot
    
    return(list(Barplot = Barplot))
    
}

Root_Bac <- Root_Relab_Abundance("data/16s_sed/", "ps_16S_Sed.RData","Bacteria","Sediments")
Root_Fun <- Root_Relab_Abundance("data/its_sed/", "ps_fungi.RData","Fungi","Sediments")
