
rm(list = ls()) # Clear all
# SETUP ####

# the following command sets the working directory to the folder where this
# script is located (similar to RStudio menu "Session - Set Working Directory -
# To Source File Location"):
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# LOAD PACKAGES ####
library(openxlsx) 
library(terra)  
library(gridExtra)
library(ncdf4)
library(raster)
library(lattice)
library(ggplot2)
library(ENMTools) # niche breadth
library(fields)
library(maps)
library(sf)
library(ENMeval)
library(dplyr)
library(ggpubr)
library(ape)
library(phangorn)
library(SYNCSA)
library(phyr)
library(openxlsx)
library(lme4)
library(lmerTest)
library(phytools)
library(pglmm)
library(rr2)
library(TreeTools)


getwd()
my_path<- "" # The directory where outputs will be saved
getwd()
species_list <- c("Ursus.spelaeus", "Ursus.arctos", "Panthera.spelaea",  "Panthera.leo", "Panthera.pardus", 
                   "Crocuta.crocuta","Canis.lupus", "Meles.meles",
                  "Vulpes.vulpes", "Vulpes.lagopus", "Felis.silvestris",
                   "Homo_neand","Homo.sapiens")

#### MEAN NICHE BREADTH ####
# Before running this code, it is necessary to have run the previous 1_SpeciesDistributionModels.R code.
# The code below computes the mean Levin's B2 for each species during the MIS3.
# You can compute the Levin's B1 by replacing "B2" with "B1" below:

niche_beadth_results <- data.frame(species = species_list, niche_breadth = NA)

for (i in 1:length(species_list)) {
  species <- species_list[i]
  path_species <- file.path(my_path, species, "Layers", "TIF")
  tif_files <- list.files(path_species, pattern = "f\\.tif$", full.names = TRUE)
  
  raster_list <- lapply(tif_files, raster)
  raster_stack <- stack(raster_list)
  raster_stack_1 <- calc(raster_stack, mean)
  # Plot mean
  plot(raster_stack_1, main = species)
  species_niche_breadth <- calc.B2(raster_stack_1)
  species_result <- data.frame(species = species, niche_breadth = species_niche_breadth)
  niche_beadth_results[nrow(niche_beadth_results) + 1,] <- c(species, species_niche_breadth)
  }
niche_beadth_results <- na.omit(niche_beadth_results)
niche_beadth_results$breadth <- as.numeric(niche_beadth_results$niche_breadth)

ggplot(data=niche_beadth_results, aes(x=species, y= breadth ))+
  geom_bar(stat = "identity")


# SPECIFIC NICHE & SPATIAL OVERLAP THROUGHOUT THE MIS3 ####
# The code below computes the niche and spatial overlap for each species pair for every 1,000
# years between 55 and 27 kyr BP:

# Specify the directory with the TIFF files:
path_TIF <- "" 

# Here you can use the same directory:
my_path <- ""

# Function to calculate Schoener's D (niche overlap):
calculate_overlap_D <- function(raster_stack_1, raster_stack_2) {
  return(nicheOverlap(raster_stack_1, raster_stack_2, stat='D'))
}
# Function to compute Warren's I (niche overlap):
calculate_overlap_I <- function(raster_stack_1, raster_stack_2) {
  return(nicheOverlap(raster_stack_1, raster_stack_2, stat='I'))
}
# Function to compute the spatial overlap (geographic overlap):
geographic_overlap <- function(raster_stack_1, raster_stack_2) {
  # Obtener los valores de favorabilidad igual o mayor a 0.5 para cada especie
  raster1_0.5 <- raster_stack_1 >= 0.5
  raster2_0.5 <- raster_stack_2 >= 0.5
  raster1_0.5 <- as(raster1_0.5, "SpatRaster")
  raster2_0.5 <- as(raster2_0.5, "SpatRaster")
  plot(raster1_0.5)
  plot(raster2_0.5)
  e <- raster1_0.5 * raster2_0.5
  plot(e)
  overlap_cells <- sum(values(e, na.rm = TRUE)>0.5)
  return(overlap_cells)

}



# Get lists of TIFF files for each species
ref_species <- species_list[1]  # Let's use the first as reference
path_ref_species <- file.path(my_path, ref_species, "Layers", "TIF")
tif_files_ref <- list.files(path_ref_species, pattern = "f\\.tif$", full.names = TRUE)

# Create matrixes to save the overlap results:
O_D_matrix <- array(0, dim = c(length(species_list), length(species_list), length(tif_files_ref)))
dimnames(O_D_matrix) <- list(species_list, species_list, paste0("Age_", 1:length(tif_files_ref)))

O_I_matrix <- array(0, dim = c(length(species_list), length(species_list), length(tif_files_ref)))
dimnames(O_I_matrix) <- list(species_list, species_list, paste0("Age_", 1:length(tif_files_ref)))

O_G_matrix <- array(0, dim = c(length(species_list), length(species_list), length(tif_files_ref)))
dimnames(O_G_matrix) <- list(species_list, species_list, paste0("Age_", 1:length(tif_files_ref)))


# Iterate for each species and age
for (i in 1:length(species_list)) {
  species <- species_list[i]
  path_species <- file.path(my_path, species, "Layers", "TIF")
  tif_files <- list.files(path_species, pattern = "f\\.tif$", full.names = TRUE)
  
  if (length(tif_files) > 0) {
    for (tif_index in 1:length(tif_files)) {
      tif_file <- tif_files[tif_index]
      raster_stack_1 <- calc(raster(tif_file), mean)
      age <- as.numeric(gsub(".*?/(\\d+)f\\.tif$", "\\1", tif_file))
      
      for (j in 1:length(species_list)) {
        other_species <- species_list[j]
        if (other_species != species) {
          path_other_species <- file.path(my_path, other_species, "Layers", "TIF")
          tif_files_other <- list.files(path_other_species, pattern = paste0("^", age, "f\\.tif$"), full.names = TRUE)
          
          if (length(tif_files_other) > 0) {
            raster_list_other <- lapply(tif_files_other, raster)
            raster_stack_other <- stack(raster_list_other)
            raster_stack_2 <- calc(raster_stack_other, mean)
            
            # Calcular la superposición de nicho
            O_D <- calculate_overlap_D(raster_stack_1, raster_stack_2)
            O_I <- calculate_overlap_I(raster_stack_1, raster_stack_2)
            O_G <- geographic_overlap(raster_stack_1, raster_stack_2)
 
            
            # Almacenar el resultado en las matrices correspondientes
            O_D_matrix[i, j, tif_index] <- O_D
            O_I_matrix[i, j, tif_index] <- O_I
            O_G_matrix[i, j, tif_index] <- O_G

            
          }
        }
      }
    }
  }
}



# Create a list with the age names
age_names <- paste0("", 27:58)

# Create an empty dataframe to save the combined outputs
combined_df_D <- data.frame()
combined_df_I <- data.frame()
combined_df_G <- data.frame()

# Create a dataframe for each overlap variable
for (age_name in age_names) {
  age_df <- as.data.frame(O_D_matrix[,, age_name])
  age_df$Age <- age_name
  
  combined_df_D <- rbind(combined_df_D, age_df)
}

for (age_name in age_names) {
  age_df <- as.data.frame(O_I_matrix[,, age_name])
  age_df$Age <- age_name
  
  combined_df_I <- rbind(combined_df_I, age_df)
}

for (age_name in age_names) {
  age_df <- as.data.frame(O_G_matrix[,, age_name])
  age_df$Age <- age_name
  
  combined_df_G <- rbind(combined_df_G, age_df)
}


# Save:
write.csv(combined_df_D, "Outcomes/combined_df_D.csv", row.names = TRUE)
write.csv(combined_df_I, "Outcomes/combined_df_I.csv", row.names = TRUE)
write.csv(combined_df_G, "Outcomes/combined_df_G.csv", row.names = TRUE)


# Summary
Overalp_D<- read.csv("Outcomes/combined_df_D.csv", sep =";")
head(Overalp_D)
summary_stats <- Overlap_D %>%
  group_by(Age, Species) %>%
  summarise(
    mean_Panthera_spelaea = mean(Panthera.spelaea, na.rm = TRUE),
    sd_Panthera_spelaea = sd(Panthera.spelaea, na.rm = TRUE)
  )

########## Analyses ####

# BoxPlot Overlap####
getwd()
data<-read.xlsx("Outcomes/Geographic_Overlap.xlsx", rowNames=FALSE, 
                colNames=TRUE, sheet="Hoja1")
head(data)

data<- subset(data, Species != "Homo_neand")
data<- subset(data, Species != "Homo_sapiens")

data$Species <- factor(data$Species, levels = c( "Ursus.arctos","Ursus.spelaeus", "Panthera.leo",
    "Panthera.spelaea", "Panthera.pardus" ,  "Crocuta.crocuta",  "Vulpes.vulpes", "Vulpes.lagopus",
    "Canis.lupus", "Felis.silvestris",   "Meles.meles"))



boxplot<- ggplot(data, aes(x= Species, y=shared, fill= Homo))+
  geom_boxplot()+
  scale_fill_manual(values = c("3", "2"))+
  theme_classic()+
  ylab("Potential geographic range shared with humans (%)")


boxplot + stat_compare_means(aes(group = Homo), label="p.signif", method = "wilcox.test")


# Phylogenetic mixed Model ####

# Directory where "Complete_phylogeny.nex" (Obtained from Phylacine https://datadryad.org/stash/dataset/doi:10.5061/dryad.bp26v20) was located
getwd()
file <- "Complete_phylogeny.nex"

# Species list
species_list <- c("Ursus_spelaeus", "Ursus_arctos", "Panthera_spelaea",
                  "Panthera_leo", "Panthera_pardus", "Crocuta_crocuta",
                  "Canis_lupus", "Meles_meles", "Vulpes_vulpes",
                  "Vulpes_lagopus", "Felis_silvestris", "Homo_neanderthalensis",
                  "Homo_sapiens")

# Load phylogenetic trees
phylo_tree <- read.nexus(file = file)
str(phylo_tree)
# Get labels for tree nodes: 
node_labels <- phylo_tree$tip.label
# Print nodes
print(node_labels)

#phylo_tree_1 <- phylo_tree[[1]]
#e<- consensus



# Filter the tree to include only the species of the list:

filtered_tree <- keep.tip(phylo_tree, tip = species_list)
# Build a consensus tree
consesus_tree<- phytools::consensus.edges(filtered_tree)
# Plot consensus tree
plot(consesus_tree, edge.color = ifelse(consesus_tree$edge.length > 0.1, "red", "black"))

# PGLM:
getwd()
data_NicheBreadth<-read.xlsx("Outcomes/PLMM.xlsx", rowNames=FALSE, 
                             colNames=TRUE, sheet="Outcomes_Breadth_Co_occ")
head(data_NicheBreadth)
nrow(data_NicheBreadth)
data_NicheBreadth_ag <- subset(data_NicheBreadth, Pgt_Count > 0) # Aggregations
data_NicheBreadth_seg <- subset(data_NicheBreadth, P_lt_Count > 0) # Segregations


head(data_NicheBreadth)

## Covariance matrix with phylogenetic structure:
phylo_cov <- vcv.phylo(consesus_tree)
phylo_cov


###### TOTAL NUMBER OF CO-OCCURRENCES ####
Total_Obs_cooccur_NicheB_B1 <- pglmm(Total_Obs_cooccur ~ NicheB1 + (1 | Species), 
                                     data = data_NicheBreadth, 
                                     cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                                     REML = TRUE, 
                                     verbose = FALSE)

summary(Total_Obs_cooccur_NicheB_B1)
round(rr2::R2(Total_Obs_cooccur_NicheB_B1[[1]]), 2)


Total_Obs_cooccur_NicheB_B2 <- pglmm(Total_Obs_cooccur ~ NicheB2 + (1 | Species), 
                                     data = data_NicheBreadth, 
                                     cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                                     REML = TRUE, 
                                     verbose = FALSE)

summary(Total_Obs_cooccur_NicheB_B2)
round(rr2::R2(Total_Obs_cooccur_NicheB_B2[[1]]), 2)

Total_Obs_cooccur_Favorable_area <- pglmm(Total_Obs_cooccur ~ Favorable_area + (1 | Species), 
                                          data = data_NicheBreadth, 
                                          cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                                          REML = TRUE, 
                                          verbose = FALSE)

summary(Total_Obs_cooccur_Favorable_area)
round(rr2::R2(Total_Obs_cooccur_Favorable_area), 2)

Total_Obs_cooccur_Overlap_G <- pglmm(Total_Obs_cooccur ~ Overlap_G + (1 | Species), 
                                     data = data_NicheBreadth, 
                                     cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                                     REML = TRUE, 
                                     verbose = FALSE)

summary(Total_Obs_cooccur_Overlap_G)
round(rr2::R2(Total_Obs_cooccur_Overlap_G)[[1]], 2)

Total_Obs_cooccur_Overlap_D <- pglmm(Total_Obs_cooccur ~ Overlap_D + (1 | Species), 
                                     data = data_NicheBreadth, 
                                     cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                                     REML = TRUE, 
                                     verbose = FALSE)

summary(Total_Obs_cooccur_Overlap_D)
round(rr2::R2(Total_Obs_cooccur_Overlap_D)[[1]],2)

Total_Obs_cooccur_Overlap_I <- pglmm(Total_Obs_cooccur ~ Overlap_I + (1 | Species), 
                                     data = data_NicheBreadth, 
                                     cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                                     REML = TRUE, 
                                     verbose = FALSE)

summary(Total_Obs_cooccur_Overlap_I)
round(rr2::R2(Total_Obs_cooccur_Overlap_I)[[1]], 2)



###### CO-OCCURRENCE PROBABILITY ####

Avg_Prob_Coocurrence_NicheB_B1 <- pglmm(Avg_Prob_Coocurrence ~ NicheB1 + (1 | Species), 
                                        data = data_NicheBreadth, 
                                        cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                                        REML = TRUE, 
                                        verbose = FALSE)

summary(Avg_Prob_Coocurrence_NicheB_B1)
round(rr2::R2(Avg_Prob_Coocurrence_NicheB_B1)[[1]],2)

Avg_Prob_Coocurrence_NicheB_B2 <- pglmm(Avg_Prob_Coocurrence ~ NicheB2 + (1 | Species), 
                                        data = data_NicheBreadth, 
                                        cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                                        REML = TRUE, 
                                        verbose = FALSE)

summary(Avg_Prob_Coocurrence_NicheB_B2)
round(rr2::R2(Avg_Prob_Coocurrence_NicheB_B2)[[1]],2)

Avg_Prob_Coocurrence_Favorable_area <- pglmm(Avg_Prob_Coocurrence ~ Favorable_area + (1 | Species), 
                                             data = data_NicheBreadth, 
                                             cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                                             REML = TRUE, 
                                             verbose = FALSE)

summary(Avg_Prob_Coocurrence_Favorable_area)
round(rr2::R2(Avg_Prob_Coocurrence_Favorable_area)[[1]],2)

Avg_Prob_Coocurrence_Overlap_G <- pglmm(Avg_Prob_Coocurrence ~ Overlap_G + (1 | Species), 
                                        data = data_NicheBreadth, 
                                        cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                                        REML = TRUE, 
                                        verbose = FALSE)

summary(Avg_Prob_Coocurrence_Overlap_G)
round(rr2::R2(Avg_Prob_Coocurrence_Overlap_G)[[1]],2)

Avg_Prob_Coocurrence_Overlap_D <- pglmm(Avg_Prob_Coocurrence ~ Overlap_D + (1 | Species), 
                                        data = data_NicheBreadth, 
                                        cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                                        REML = TRUE, 
                                        verbose = FALSE)

summary(Avg_Prob_Coocurrence_Overlap_D)
round(rr2::R2(Avg_Prob_Coocurrence_Overlap_D)[[1]],2)

Avg_Prob_Coocurrence_Overlap_I <- pglmm(Avg_Prob_Coocurrence ~ Overlap_I + (1 | Species), 
                                        data = data_NicheBreadth, 
                                        cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                                        REML = TRUE, 
                                        verbose = FALSE)

summary(Avg_Prob_Coocurrence_Overlap_I)
round(rr2::R2(Avg_Prob_Coocurrence_Overlap_I)[[1]],2)



###### AGGREGATIONS ####
Pgt_Count_NicheB_B1 <- pglmm(Pgt_Count ~ NicheB1 + (1 | Species), 
                             data = data_NicheBreadth_ag, 
                             cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                             REML = TRUE, 
                             verbose = FALSE)

summary(Pgt_Count_NicheB_B1)
round(rr2::R2(Pgt_Count_NicheB_B1)[[1]],2)

Pgt_Count_NicheB_B2 <- pglmm(Pgt_Count ~ NicheB2 + (1 | Species), 
                             data = data_NicheBreadth_ag, 
                             cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                             REML = TRUE, 
                             verbose = FALSE)

summary(Pgt_Count_NicheB_B2)
round(rr2::R2(Pgt_Count_NicheB_B2)[[1]],2)

Pgt_Count_Favorable_area <- pglmm(Pgt_Count ~ Favorable_area + (1 | Species), 
                                  data = data_NicheBreadth_ag, 
                                  cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                                  REML = TRUE, 
                                  verbose = FALSE)

summary(Pgt_Count_Favorable_area)
round(rr2::R2(Pgt_Count_Favorable_area)[[1]],2)

Pgt_Count_Overlap_G <- pglmm(Pgt_Count ~ Overlap_G + (1 | Species), 
                             data = data_NicheBreadth_ag, 
                             cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                             REML = TRUE, 
                             verbose = FALSE)

summary(Pgt_Count_Overlap_G)
round(rr2::R2(Pgt_Count_Overlap_G)[[1]],2)

Pgt_Count_Overlap_D <- pglmm(Pgt_Count ~ Overlap_D + (1 | Species), 
                             data = data_NicheBreadth_ag, 
                             cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                             REML = TRUE, 
                             verbose = FALSE)

summary(Pgt_Count_Overlap_D)
round(rr2::R2(Pgt_Count_Overlap_D)[[1]],2)

Pgt_Count_Overlap_I <- pglmm(Pgt_Count ~ Overlap_I + (1 | Species), 
                             data = data_NicheBreadth_ag, 
                             cov_ranef = list(Species = as(as(as(phylo_cov, "dMatrix"), "generalMatrix"), "TsparseMatrix")), 
                             REML = TRUE, 
                             verbose = FALSE)

summary(Pgt_Count_Overlap_I)
round(rr2::R2(Pgt_Count_Overlap_I)[[1]],2)


#### SEGREGATIONS ####

P_lt_Count_NicheB_B <- pglmm(P_lt_Count ~ NicheB1+ (1 | Species), 
                             data = data_NicheBreadth_seg, 
                             cov_ranef = list(Species = phylo_cov), 
                             REML = TRUE, 
                             verbose = FALSE)

summary(P_lt_Count_NicheB_B)
round(rr2::R2(P_lt_Count_NicheB_B)[[1]],2)

P_lt_Count_NicheB2 <- pglmm(P_lt_Count ~ NicheB2+ (1 | Species), 
                            data = data_NicheBreadth_seg, 
                            cov_ranef = list(Species = phylo_cov), 
                            REML = TRUE, 
                            verbose = FALSE)

summary(P_lt_Count_NicheB2)
round(rr2::R2(P_lt_Count_NicheB2)[[1]],2)

P_lt_Count_Favorable_area <- pglmm(P_lt_Count ~ Favorable_area+ (1 | Species), 
                                   data = data_NicheBreadth_seg, 
                                   cov_ranef = list(Species = phylo_cov), 
                                   REML = TRUE, 
                                   verbose = FALSE)

summary(P_lt_Count_Favorable_area)
round(rr2::R2(P_lt_Count_Favorable_area)[[1]],2)

P_lt_Count_Overlap_G <- pglmm(P_lt_Count ~ Overlap_G+ (1 | Species), 
                              data = data_NicheBreadth_seg, 
                              cov_ranef = list(Species = phylo_cov), 
                              REML = TRUE, 
                              verbose = FALSE)

summary(P_lt_Count_Overlap_G)
round(rr2::R2(P_lt_Count_Overlap_G)[[1]],2)


P_lt_Count_Overlap_I <- pglmm(P_lt_Count ~ Overlap_I+ (1 | Species), 
                              data = data_NicheBreadth_seg, 
                              cov_ranef = list(Species = phylo_cov), 
                              REML = TRUE, 
                              verbose = FALSE)

summary(P_lt_Count_Overlap_I)
round(rr2::R2(P_lt_Count_Overlap_I)[[1]],2)

P_lt_Count_Overlap_D <- pglmm(P_lt_Count ~ Overlap_D+ (1 | Species), 
                              data = data_NicheBreadth_seg, 
                              cov_ranef = list(Species = phylo_cov), 
                              REML = TRUE, 
                              verbose = FALSE)

summary(P_lt_Count_Overlap_D)
round(rr2::R2(P_lt_Count_Overlap_D)[[1]],2)



