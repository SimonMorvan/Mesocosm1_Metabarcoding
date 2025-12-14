# Installs the libraries needed

#Installs pacman if needed
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

# Load  pacman
library(pacman)

# Forces the CRAN officiel CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# --------
# CRAN packages
# --------
cran_pkgs <- c("agricolae", "colorspace", "cowplot", "dunn.test",
               "glue", "ggdist", "ggforce", "ggh4x",
               "ggpubr", "ggtext", "forcats", "here", "labdsv",
               "MASS","psych", "tidyverse", "roperators", "vegan")


# Install and load CRAN packages
pacman::p_load(char = cran_pkgs, install = TRUE)


# --------
# Bioconductor packages
# --------
bioc_pkgs <- c("DESeq2", "BiocGenerics","phyloseq","EnhancedVolcano",
               "MicEco","pairwiseAdonis")  # BiocGenerics est souvent requis par DESeq2

# Install and load  Bioconductor packages
for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
    library(pkg, character.only = TRUE)
}

# --------
# GitHub packages
# --------

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
library(devtools)

# MicEco
if (!requireNamespace("MicEco", quietly = TRUE)) {
    devtools::install_github("Russel88/MicEco")
}
library(MicEco)

# pairwiseAdonis
if (!requireNamespace("pairwiseAdonis", quietly = TRUE)) {
    devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
}
library(pairwiseAdonis)


