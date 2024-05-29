
rm(list = ls()) # Clear all
# SETUP ####

# the following command sets the working directory to the folder where this
# script is located (similar to RStudio menu "Session - Set Working Directory -
# To Source File Location"):
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
library(cooccur)
library(igraph)
library(openxlsx) # Excel
library(dplyr)
library(viridis)
library(tidyr)
library(readxl)
library(ggplot2)
library(iNEXT)



getwd()


############ RAREFACTION ANALYSIS ####

# PHASE I
Matrix_Phase_I <- read.xlsx("Data.xlsx", rowNames=FALSE,
                     colNames=TRUE, sheet="50_45")
# Re-configurate: species as rows, assemblages as columns:
Matrix_Phase_I<- t(Matrix_Phase_I[,-c(1:4)])
head(Matrix_Phase_I)

SpeciesDataset = list(Phase_I= Matrix_Phase_I)
out.raw <- iNEXT(SpeciesDataset, datatype="incidence_raw", endpoint=300, nboot=500)
ggiNEXT(out.raw, type = 1 ,
        se = TRUE,
        facet.var = "None",
        color.var = "Order.q",
        grey = FALSE)+
  theme_classic() #PLOT
print(out.raw)

# PHASE II
Matrix_Phase_II <- read.xlsx("Data.xlsx", rowNames=FALSE,
                            colNames=TRUE, sheet="45_40")
# Re-configurate: species as rows, assemblages as columns:
Matrix_Phase_II<- t(Matrix_Phase_II[,-c(1:4)])
head(Matrix_Phase_II)

SpeciesDataset = list(Phase_II= Matrix_Phase_II)
out.raw <- iNEXT(SpeciesDataset, datatype="incidence_raw", endpoint=300, nboot=500)
ggiNEXT(out.raw, type = 1 ,
        se = TRUE,
        facet.var = "None",
        color.var = "Order.q",
        grey = FALSE)+
  theme_classic() #PLOT
print(out.raw)

# PHASE III
Matrix_Phase_III <- read.xlsx("Data.xlsx", rowNames=FALSE,
                             colNames=TRUE, sheet="40_35")
# Re-configurate: species as rows, assemblages as columns:
Matrix_Phase_III<- t(Matrix_Phase_III[,-c(1:4)])
head(Matrix_Phase_III)

SpeciesDataset = list(Phase_III= Matrix_Phase_III)
out.raw <- iNEXT(SpeciesDataset, datatype="incidence_raw", endpoint=300, nboot=500)
ggiNEXT(out.raw, type = 1 ,
        se = TRUE,
        facet.var = "None",
        color.var = "Order.q",
        grey = FALSE)+
  theme_classic() #PLOT
print(out.raw)

# PHASE IV
Matrix_Phase_IV <- read.xlsx("Data.xlsx", rowNames=FALSE,
                              colNames=TRUE, sheet="35_30")
# Re-configurate: species as rows, assemblages as columns:
Matrix_Phase_IV<- t(Matrix_Phase_IV[,-c(1:4)])
head(Matrix_Phase_IV)

SpeciesDataset = list(Phase_IV= Matrix_Phase_IV)
out.raw <- iNEXT(SpeciesDataset, datatype="incidence_raw", endpoint=300, nboot=500)
ggiNEXT(out.raw, type = 1 ,
        se = TRUE,
        facet.var = "None",
        color.var = "Order.q",
        grey = FALSE)+
  theme_classic() #PLOT
print(out.raw)


############ CO-OCCURRENCES ANALYSIS ####

# Phase I (50-45 kyr BP) ####
species_matrix <- read.xlsx("Data.xlsx", rowNames=FALSE,
                            colNames=TRUE, sheet = "50_45")
species_matrix <- species_matrix[, -(1:4)]
head(species_matrix)
# Change matrix configuration (species as rows, assemblages as columns): 
species_matrix <- t(species_matrix)
# Co-occurrence:
cooc_matrix <- cooccur(as.matrix(species_matrix), type = "spp_site", thresh = FALSE, eff_matrix = TRUE, spp_names = TRUE, prob = "comb", true_rand_classifier = FALSE)
plot(cooc_matrix)
cooc_results <- prob.table(cooc_matrix)
# Prepare and save outcomes in a new dataframe
summary_df_1 <- cooc_results %>%
  group_by(sp1_name) %>%
  dplyr::summarise(
    Total_Obs_cooccur = sum(obs_cooccur),
    Avg_Prob_Coocurrence = mean(prob_cooccur),
    P_lt_Count = sum(p_lt < 0.05),
    Pgt_Count = sum(p_gt < 0.05)
  )
summary_df_1
summary_df_2 <- cooc_results %>%
  group_by(sp2_name) %>%
  dplyr::summarise(
    Total_Obs_cooccur = sum(obs_cooccur),
    Avg_Prob_Coocurrence = mean(prob_cooccur),
    P_lt_Count = sum(p_lt < 0.05),
    Pgt_Count = sum(p_gt < 0.05)
  )
print(summary_df_2)

# Combine in a dataframe:
combined_summary <- bind_rows(
  mutate(summary_df_1, Species = sp1_name),
  mutate(summary_df_2, Species = sp2_name)
) %>%
  group_by(Species) %>%
  dplyr::summarise(
    Total_Obs_cooccur = sum(Total_Obs_cooccur),
    Avg_Prob_Coocurrence = mean(Avg_Prob_Coocurrence),
    P_lt_Count = sum(P_lt_Count),
    Pgt_Count = sum(Pgt_Count)
  )
# Save results:
write.csv(combined_summary, "Outcomes/50_45_co_occ.csv")
write.csv(cooc_results, "Outcomes/50_45_co_occ_details.csv")


# Phase II (45-40 kyr BP) ####
species_matrix <- read.xlsx("Data.xlsx", rowNames=FALSE,
                            colNames=TRUE, sheet = "45_40")
species_matrix <- species_matrix[, -(1:4)]
head(species_matrix)
# Change matrix configuration (species as rows, assemblages as columns): 
species_matrix <- t(species_matrix)
# Co-occurrence:
cooc_matrix <- cooccur(as.matrix(species_matrix), type = "spp_site", thresh = FALSE, eff_matrix = TRUE, spp_names = TRUE, prob = "comb", true_rand_classifier = FALSE)
plot(cooc_matrix)
cooc_results <- prob.table(cooc_matrix)
# Prepare and save outcomes in a new dataframe
summary_df_1 <- cooc_results %>%
  group_by(sp1_name) %>%
  dplyr::summarise(
    Total_Obs_cooccur = sum(obs_cooccur),
    Avg_Prob_Coocurrence = mean(prob_cooccur),
    P_lt_Count = sum(p_lt < 0.05),
    Pgt_Count = sum(p_gt < 0.05)
  )
summary_df_1
summary_df_2 <- cooc_results %>%
  group_by(sp2_name) %>%
  dplyr::summarise(
    Total_Obs_cooccur = sum(obs_cooccur),
    Avg_Prob_Coocurrence = mean(prob_cooccur),
    P_lt_Count = sum(p_lt < 0.05),
    Pgt_Count = sum(p_gt < 0.05)
  )
print(summary_df_2)

# Combine in a dataframe:
combined_summary <- bind_rows(
  mutate(summary_df_1, Species = sp1_name),
  mutate(summary_df_2, Species = sp2_name)
) %>%
  group_by(Species) %>%
  dplyr::summarise(
    Total_Obs_cooccur = sum(Total_Obs_cooccur),
    Avg_Prob_Coocurrence = mean(Avg_Prob_Coocurrence),
    P_lt_Count = sum(P_lt_Count),
    Pgt_Count = sum(Pgt_Count)
  )
# Save results:
write.csv(combined_summary, "Outcomes/45_40_co_occ.csv")
write.csv(cooc_results, "Outcomes/45_40_co_occ_details.csv")

# Phase III (40-35 kyr BP) ####
species_matrix <- read.xlsx("Data.xlsx", rowNames=FALSE,
                            colNames=TRUE, sheet = "40_35")
species_matrix <- species_matrix[, -(1:4)]
head(species_matrix)
# Change matrix configuration (species as rows, assemblages as columns): 
species_matrix <- t(species_matrix)
# Co-occurrence:
cooc_matrix <- cooccur(as.matrix(species_matrix), type = "spp_site", thresh = FALSE, eff_matrix = TRUE, spp_names = TRUE, prob = "comb", true_rand_classifier = FALSE)
plot(cooc_matrix)
cooc_results <- prob.table(cooc_matrix)
# Prepare and save outcomes in a new dataframe
summary_df_1 <- cooc_results %>%
  group_by(sp1_name) %>%
  dplyr::summarise(
    Total_Obs_cooccur = sum(obs_cooccur),
    Avg_Prob_Coocurrence = mean(prob_cooccur),
    P_lt_Count = sum(p_lt < 0.05),
    Pgt_Count = sum(p_gt < 0.05)
  )
summary_df_1
summary_df_2 <- cooc_results %>%
  group_by(sp2_name) %>%
  dplyr::summarise(
    Total_Obs_cooccur = sum(obs_cooccur),
    Avg_Prob_Coocurrence = mean(prob_cooccur),
    P_lt_Count = sum(p_lt < 0.05),
    Pgt_Count = sum(p_gt < 0.05)
  )
print(summary_df_2)

# Combine in a dataframe:
combined_summary <- bind_rows(
  mutate(summary_df_1, Species = sp1_name),
  mutate(summary_df_2, Species = sp2_name)
) %>%
  group_by(Species) %>%
  dplyr::summarise(
    Total_Obs_cooccur = sum(Total_Obs_cooccur),
    Avg_Prob_Coocurrence = mean(Avg_Prob_Coocurrence),
    P_lt_Count = sum(P_lt_Count),
    Pgt_Count = sum(Pgt_Count)
  )
# Save results:
write.csv(combined_summary, "Outcomes/40_35_co_occ.csv")
write.csv(cooc_results, "Outcomes/40_35_co_occ_details.csv")


# Phase IV (35-30 kyr BP) ####
species_matrix <- read.xlsx("Data.xlsx", rowNames=FALSE,
                            colNames=TRUE, sheet = "35_30")
species_matrix <- species_matrix[, -(1:4)]
head(species_matrix)
# Change matrix configuration (species as rows, assemblages as columns): 
species_matrix <- t(species_matrix)
# Co-occurrence:
cooc_matrix <- cooccur(as.matrix(species_matrix), type = "spp_site", thresh = FALSE, eff_matrix = TRUE, spp_names = TRUE, prob = "comb", true_rand_classifier = FALSE)
plot(cooc_matrix)
cooc_results <- prob.table(cooc_matrix)
# Prepare and save outcomes in a new dataframe
summary_df_1 <- cooc_results %>%
  group_by(sp1_name) %>%
  dplyr::summarise(
    Total_Obs_cooccur = sum(obs_cooccur),
    Avg_Prob_Coocurrence = mean(prob_cooccur),
    P_lt_Count = sum(p_lt < 0.05),
    Pgt_Count = sum(p_gt < 0.05)
  )
summary_df_1
summary_df_2 <- cooc_results %>%
  group_by(sp2_name) %>%
  dplyr::summarise(
    Total_Obs_cooccur = sum(obs_cooccur),
    Avg_Prob_Coocurrence = mean(prob_cooccur),
    P_lt_Count = sum(p_lt < 0.05),
    Pgt_Count = sum(p_gt < 0.05)
  )
print(summary_df_2)

# Combine in a dataframe:
combined_summary <- bind_rows(
  mutate(summary_df_1, Species = sp1_name),
  mutate(summary_df_2, Species = sp2_name)
) %>%
  group_by(Species) %>%
  dplyr::summarise(
    Total_Obs_cooccur = sum(Total_Obs_cooccur),
    Avg_Prob_Coocurrence = mean(Avg_Prob_Coocurrence),
    P_lt_Count = sum(P_lt_Count),
    Pgt_Count = sum(Pgt_Count)
  )
# Save results:
write.csv(combined_summary, "Outcomes/35_30_co_occ.csv")
write.csv(cooc_results, "Outcomes/35_30_co_occ_details.csv")





## PLOT CO-OCCURRENCE ####


combined_summary<- read.csv("Outcomes/50_45_co_occ.csv", row.names = 1) # Specify the phase you want to plot:
cooc_results <-read.csv("Outcomes/50_45_co_occ_details.csv", row.names = 1) 
# Phase I: combined_summary<- read.csv("Outcomes/50_45_co_occ.csv", row.names = 1)
#          cooc_results <-read.csv("Outcomes/50_45_co_occ_details.csv", row.names = 1)
# Phase II: combined_summary<- read.csv("Outcomes/45_40_co_occ.csv", row.names = 1)
#          cooc_results <-read.csv("Outcomes/45_40_co_occ_details.csv", row.names = 1)
# Phase III: combined_summary<- read.csv("Outcomes/40_35_co_occ.csv", row.names = 1)
#          cooc_results <-read.csv("Outcomes/40_35_co_occ_details.csv", row.names = 1)
# Phase VI: combined_summary<- read.csv("Outcomes/35_30_co_occ.csv", row.names = 1)
#          cooc_results <-read.csv("Outcomes/35_30_co_occ_details.csv", row.names = 1)
 
  # Create graph from the .csv files previously saved:
  g <- graph.data.frame(cooc_results[, c("sp1_name", "sp2_name")], directed = FALSE)
  # Get the species order in V(g)$name:
  species_order <- V(g)$name
  # Re-order combined_summary according to the order of species in V(g)$name:
  combined_summary <- combined_summary[match(species_order, combined_summary$Species), ]
  #Check that the order match and is correct: 
  head(combined_summary)
  # Get node sizes from combined_summary:
  node_sizes <- combined_summary$Avg_Prob_Coocurrence
  # Color for plot and legend:
  breaks <- seq(0, 0.3, by = 0.1)
  # Node colors according to the co-occurrence value:
  node_colors <- viridis(length(breaks) - 1)
  # Links according to the statistical significance (p-value) of the aggregations and segregations:
  V(g)$color <- node_colors[findInterval(node_sizes, breaks)]
  E(g)$color <- ifelse(round(cooc_results$p_lt, 2) < 0.05 , "red", 
                       ifelse(round(cooc_results$p_gt, 2) < 0.05, "blue", "grey"))
  # Associate the node sizes:
  V(g)$size <- node_sizes * 100
  # Plot:
  plot(g, layout = layout.circle, vertex.label = NA, vertex.size = V(g)$size, vertex.color = V(g)$color)
  # Labels
  text(x = layout.circle(g)[,1], y = layout.circle(g)[,2] - 0.1, labels = V(g)$name, pos = 1, cex = 0.7)
  # Legend
  legend_values <- breaks[-length(breaks)]  # Excluir el último valor para evitar problemas de superposición de colores
  legend("bottomright", legend = legend_values, title = "Total Obs Cooccur", pch = 16,
         col = viridis(length(breaks) - 1), cex = 0.8)
  
  
  
  
  
##### SENSITIVITY TEST ####
  # The following code repeats the previous co-occurrence analyses 1000 times, each time with a subsample of 70
  # assemblages per Phase:
  
  # Load data
  species_matrix <- read.xlsx("Data.xlsx", rowNames=FALSE, colNames=TRUE, sheet = "50_45")
  species_matrix <- species_matrix[, -(1:4)]
  species_matrix <- t(species_matrix)
  
  combined_summary_all <- data.frame()
  cooc_results_all <- data.frame()
  
  # Resampling and co-occurrence
  for (experiment in 1:1000) {
    # Submuestreo aleatorio de 70 sitios
    index_assemblage <- sample(1:ncol(species_matrix), 70)
    species_matrix_sub <- species_matrix[, index_assemblage]
    
    # Co-occurrence analysis
    cooc_matrix <- cooccur(as.matrix(species_matrix_sub), type = "spp_site", thresh = FALSE, eff_matrix = TRUE, spp_names = TRUE, prob = "comb", true_rand_classifier = FALSE)
    cooc_results <- prob.table(cooc_matrix)
    
    plot(cooc_matrix)

    # Summary
    summary_df_1 <- cooc_results %>%
      group_by(sp1_name) %>%
      dplyr::summarise(
        Total_Obs_cooccur = sum(obs_cooccur),
        Avg_Prob_Coocurrence = mean(prob_cooccur),
        P_lt_Count = sum(p_lt < 0.05),
        Pgt_Count = sum(p_gt < 0.05)
      )
    
    summary_df_2 <- cooc_results %>%
      group_by(sp2_name) %>%
      dplyr::summarise(
        Total_Obs_cooccur = sum(obs_cooccur),
        Avg_Prob_Coocurrence = mean(prob_cooccur),
        P_lt_Count = sum(p_lt < 0.05),
        Pgt_Count = sum(p_gt < 0.05)
      )
    
    combined_summary <- bind_rows(
      mutate(summary_df_1, Species = sp1_name),
      mutate(summary_df_2, Species = sp2_name)
    ) %>%
      group_by(Species) %>%
      dplyr::summarise(
        Total_Obs_cooccur = sum(Total_Obs_cooccur),
        Avg_Prob_Coocurrence = mean(Avg_Prob_Coocurrence),
        P_lt_Count = sum(P_lt_Count),
        Pgt_Count = sum(Pgt_Count)
      )
    
    plot(cooc_matrix)
    
    # Save results
    combined_summary_all <- rbind(combined_summary_all, combined_summary)
    cooc_results$experiment <- experiment
    cooc_results_all <- rbind(cooc_results_all, cooc_results)
    
  }
  
  
 # Save

  write.csv(cooc_results_all, "Outcomes/Sensitivity_Cooccurrence/cooc_results_all_50_45_.csv") # Remember you have to specify the phase (50_45, 45_40, 40_35 or 35_30)
  
  ### plot
  getwd()
  
  
 Df<- read.csv("Outcomes/Sensitivity_Cooccurrence/cooc_results_all_50_45_.csv")
  head(Df)


  # Resumen de los resultados
  summary_df_1 <- Df %>%
    group_by(sp1_name, experiment) %>%
    dplyr::summarise(
      Total_Obs_cooccur = sum(obs_cooccur),
      SD_Total_cooc= sd(obs_cooccur) ,
      
      Avg_Prob_Coocurrence = mean(prob_cooccur),
      SD_Prob= sd(prob_cooccur),
      P_lt_Count = sum(p_lt < 0.05),
      Pgt_Count = sum(p_gt < 0.05)
      
    )
  
  summary_df_1
  summary_df_2 <- Df %>%
    group_by(sp2_name, experiment) %>%
    dplyr::summarise(
      Total_Obs_cooccur = sum(obs_cooccur),
      SD_Total_cooc= sd(obs_cooccur) ,
      
      Avg_Prob_Coocurrence = mean(prob_cooccur),
      SD_Prob= sd(prob_cooccur),
      P_lt_Count = sum(p_lt < 0.05),
      Pgt_Count = sum(p_gt < 0.05)
      
    )
  summary_df_2
  combined_summary <- bind_rows(
    mutate(summary_df_1, Species = sp1_name),
    mutate(summary_df_2, Species = sp2_name)
  ) %>%
    group_by(Species, experiment) %>%
    dplyr::summarise(
      Total_Obs_cooccur = sum(Total_Obs_cooccur),
      Avg_Prob_Coocurrence = mean(Avg_Prob_Coocurrence),
      P_lt_Count = sum(P_lt_Count),
      Pgt_Count = sum(Pgt_Count)
    )
  write.csv( combined_summary, "Outcomes/Sensitivity_Cooccurrence/Sensitivity_combined_summary_50_45.csv")
  
  # Species order:
  species_list <- c("Ursus.spelaeus", "Ursus.arctos", "Panthera.spelaea", "Panthera.leo", "Panthera.pardus", 
                    "Crocuta.crocuta", "Canis.lupus", "Meles.meles",
                    "Vulpes.vulpes", "Vulpes.lagopus", "Felis.silvestris",
                    "Homo_neand", "Homo.sapiens")
  
  # Species as factor:
  combined_summary$Species <- factor(combined_summary$Species, levels = species_list)
  
  
  # Plot mean value with SD
Aggregations_50_45<-  ggplot(combined_summary, aes(x = Species, y = Pgt_Count)) +
  stat_summary(fun.data = "mean_sdl", color="blue",
               fun.args = list(
                 mult = 1.96
               ))+
    labs(title = "Aggregations (50-45 kya)",
         x = "Especie", y = "Aggregations") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_y_continuous(limits = c(-1, 10), breaks = seq(0, 10, by = 1))
  
Aggregations_50_45


Seggregations_50_45<-  ggplot(combined_summary, aes(x = Species, y = P_lt_Count)) +
  stat_summary(fun.data = "mean_sdl", color= "red",
               fun.args = list(
                 mult = 1.96
               ))+
  labs(title = "Segregations (50-45 kya)",
       x = "Especie", y = "Segregations") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_y_continuous(limits = c(-1, 10), breaks = seq(0, 10, by = 1))

Seggregations_50_45

grid.arrange(Aggregations_50_45, Seggregations_50_45)

# Repeat the sensitivity test for each Phase.



######### BODY MASS RATIO #####

Df<- read.csv("Outcomes/Sensitivity_Cooccurrence/cooc_results_all_50_45_.csv") # Select the Phase of interest. 

head(Df)

# Body masses according to the Phylacine Database:
body_mass <- data.frame(
  Species = c("Canis.lupus", "Crocuta.crocuta", "Felis.silvestris", "Homo.sapiens", 
              "Homo_neand", "Meles.meles", "Panthera.leo", "Panthera.pardus", 
              "Panthera.spelaea", "Ursus.arctos", "Ursus.spelaeus", "Vulpes.lagopus", 
              "Vulpes.vulpes"),
  BodyMass = c(32183.3, 62999.9, 5500, 59500, 76000, 13000, 161499.1, 54999.7, 
               380189.4, 180520.4, 390416.7, 4867.7, 5318.2)
)

# Merge:
Df <- merge(Df, body_mass, by.x = "sp1_name", by.y = "Species", all.x = TRUE)
Df <- merge(Df, body_mass, by.x = "sp2_name", by.y = "Species", all.x = TRUE, suffixes = c("_sp1", "_sp2"))

# Rename columns
names(Df)[names(Df) == "BodyMass_sp1"] <- "BM1"
names(Df)[names(Df) == "BodyMass_sp2"] <- "BM2"

# Compute the Body mass ration
Df$BM_ratio <- ifelse(Df$BM1 < Df$BM2, Df$BM1 / Df$BM2, Df$BM2 / Df$BM1)
# Identify aggregations and segregations:
Df$Segregations <- ifelse(Df$p_lt < 0.05, 1, 0)
Df$Aggregations <- ifelse(Df$p_gt < 0.05, 1, 0)

head(Df)
# Create three groups
Df$BM_group <- cut(Df$BM_ratio, breaks = c(0, 0.33, 0.66, Inf), labels = c("0.33", "0.66", "0.66+"))
# New variable for the BM ratio groups:
Df <- Df %>%
  mutate(species_group = ifelse(sp2_name == "Homo_neand" | sp1_name == "Homo_neand", "Neand",
                                ifelse(sp2_name == "Homo.sapiens" | sp1_name == "Homo.sapiens", "Sapiens", "Other")))
# save outputs
write.csv(Df, "Outcomes/Sensitivity_Cooccurrence/df_50_45_BM_RATIO.csv")

# Repeat for all phases



# PLOT

df_45_40<- read.csv("Outcomes/Sensitivity_Cooccurrence/df_45_40_BM_RATIO.csv")
head(df_45_40)  

head(Df)

Df$species_group
grouped_data <- Df %>%
  dplyr::group_by(BM_group, species_group) %>%
  dplyr::summarise(total_prob_cooccur = (prob_cooccur))
head(grouped_data)
# Crear el gráfico de barras con una barra para cada species_group
plot_50_45<- ggplot(data = grouped_data, aes(x = BM_group, y = total_prob_cooccur, fill = species_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  #geom_line(aes(group = species_group, color= species_group), position = "dodge") +
  scale_fill_manual(values = c("cyan3","black",   "grey")) +
  labs(x = "BM Ratio Group", y = "Total Probability Co-occurrence") +
  theme_minimal()

plot_50_45 # Repeat with all phases







