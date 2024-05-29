rm(list = ls()) # Clear all

# SETUP ####

# the following command sets the working directory to the folder where this
# script is located (similar to RStudio menu "Session - Set Working Directory -
# To Source File Location"):
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

getwd()
library(openxlsx)
library(ggplot2)
library(raster)
library(gridExtra)
library(rworldmap)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(viridis)
library(ggpubr)
library(lattice)



##### Load Data for Boxplots
Neand_df<-read.xlsx("Outcomes/Sensitivity_SDM/Output_Sensitivity_neand.xlsx", rowNames=FALSE, 
                         colNames=TRUE, sheet="Sheet 1")
Sap_df<-read.xlsx("Outcomes/Sensitivity_SDM/Output_Sensitivity_sapiens.xlsx", rowNames=FALSE, 
                    colNames=TRUE, sheet="Sheet 1")
Uarctos_df<-read.xlsx("Outcomes/Sensitivity_SDM/Output_Sensitivity_Ursus.arctos.xlsx", rowNames=FALSE, 
                  colNames=TRUE, sheet="Sheet 1")
Uspelaeus_df<-read.xlsx("Outcomes/Sensitivity_SDM/Output_Sensitivity_Ursus.spelaeus.xlsx", rowNames=FALSE, 
                      colNames=TRUE, sheet="Sheet 1")
Pleo_df<-read.xlsx("Outcomes/Sensitivity_SDM/Output_Sensitivity_Panthera.leo.xlsx", rowNames=FALSE, 
                      colNames=TRUE, sheet="Sheet 1")
Pspelaea_df<-read.xlsx("Outcomes/Sensitivity_SDM/Output_Sensitivity_Panthera.spelaea.xlsx", rowNames=FALSE, 
                   colNames=TRUE, sheet="Sheet 1")
Ppardus_df<-read.xlsx("Outcomes/Sensitivity_SDM/Output_Sensitivity_Panthera.pardus.xlsx", rowNames=FALSE, 
                       colNames=TRUE, sheet="Sheet 1")
Ccrocuta_df<-read.xlsx("Outcomes/Sensitivity_SDM/Output_Sensitivity_Crocuta.crocuta.xlsx", rowNames=FALSE, 
                      colNames=TRUE, sheet="Sheet 1")
Calupus_df<-read.xlsx("Outcomes/Sensitivity_SDM/Output_Sensitivity_Canis.lupus.xlsx", rowNames=FALSE, 
                       colNames=TRUE, sheet="Sheet 1")
Mmeles_df<-read.xlsx("Outcomes/Sensitivity_SDM/Output_Sensitivity_Meles.meles.xlsx", rowNames=FALSE, 
                       colNames=TRUE, sheet="Sheet 1")
Fsilvestris_df<-read.xlsx("Outcomes/Sensitivity_SDM/Output_Sensitivity_Felis.silvestris.xlsx", rowNames=FALSE, 
                       colNames=TRUE, sheet="Sheet 1")
Vvulpes_df<-read.xlsx("Outcomes/Sensitivity_SDM/Output_Sensitivity_Vulpes.vulpes.xlsx", rowNames=FALSE, 
                       colNames=TRUE, sheet="Sheet 1")
Vlagopus_df<-read.xlsx("Outcomes/Sensitivity_SDM/Output_Sensitivity_Vulpes.lagopus.xlsx", rowNames=FALSE, 
                      colNames=TRUE, sheet="Sheet 1")


# Load Data to compute z-scores. First, we have to create the co-occurrence matrix:
Presence_data<-read.xlsx("Data.xlsx", rowNames=FALSE, 
                         colNames=TRUE, sheet="F_df")
head(Presence_data)


##### BOXPLOT WITH AREA FAVOURABILITY AND Z-SCORE WITH CO-OCCURRENCE: ####
# Divide and round the age
Presence_data$Age <- round(Presence_data$Age / 1000)

# Identify the column species names
species_columns <- grep("^\\D+$", colnames(Presence_data))

# Numeric conversion
Presence_data[, species_columns] <- sapply(Presence_data[, species_columns], as.numeric)

#Sum column species
presence_sum <- colSums(Presence_data[, species_columns], na.rm = TRUE)

# Only select species with, at least, 30 occurrences
selected_species <- names(presence_sum[presence_sum >= 30])

# Replace empty values with 0
Presence_data[is.na(Presence_data)] <- 0
Presence_data <- Presence_data[, c("Site/Level", "x", "y", "Age", selected_species)]

# Eliminate empty rows:
species_matrix<- Presence_data[rowSums(Presence_data[, -c(1:8)]) != 0, ]
head(species_matrix)


# COMPUTE Z-SCORES AND PLOT
Sap_df <- Sap_df[Sap_df$Time_k_BP >= 30 & Sap_df$Time_k_BP <= 50, ]

gr_sapiens <-ggplot(data = Sap_df, aes(x = Time_k_BP, y = Area_km / 1000000, fill=Species, group = interaction(Time_k_BP, Species))) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 1, outlier.size = 1, outlier.color = "grey") +
  theme_classic() +
  xlab("Age (kya BP)") +
  ylab(expression(paste("Area "("10" ^{6}, paste("km"^{-2})))))+
  labs(fill = "Species") +
  scale_x_reverse()+
  scale_fill_manual(values = c("2")) +
  scale_y_continuous(limits = c(1, 3), breaks = seq(0.5, 3, by = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")

gr_sapiens


# Z-SCORE
select_Homo_sapiens <- species_matrix$Homo_sapiens == 1
species_matrix_with_Homo_sapiens <- species_matrix[select_Homo_sapiens, ]

# Sum the species' co-occurrences/Age
sum_by_age <- tapply(rowSums(species_matrix_with_Homo_sapiens[, -c(1:8)]), species_matrix_with_Homo_sapiens$Age, sum)

# Create a dataframe with output
data <- data.frame(Age = as.numeric(names(sum_by_age)),
                   cooccurrence = as.numeric(sum_by_age))
data <- data[data$Age >= 30 & data$Age <= 50, ]
data
# Bin size
bin_size <- 5
# Group by bin
data$Age_bin <- cut(data$Age, breaks = seq(50, 30, by = -bin_size))
data$Age_bin <- gsub("\\[|\\(|\\)|\\]", "", data$Age_bin)
data$Age_bin <- gsub(",", "-", data$Age_bin)
# Mean and SD
mean_value <- mean(data$cooccurrence)
standard_deviation <- sd(data$cooccurrence)
# Compute z-score
data$Z_score <- (data$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- aggregate(cooccurrence ~ Age_bin, data = data, FUN = mean, na.rm = TRUE)
mean_value_cooccurrence$Z_score <- (mean_value_cooccurrence$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- transform(mean_value_cooccurrence,
                                  Age_plot = sapply(strsplit(as.character(Age_bin), "-"), function(x) mean(as.numeric(x))))
# Plot
Sapiens_z<-ggplot(mean_value_cooccurrence, aes(x = Age_plot, y = Z_score, group = 1)) +
  geom_line(color = "2", linetype = "dashed") +
  geom_point(color = "2", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_reverse(limits = c(50, 30)) +
  theme_minimal() +
  labs(x = "Age", y = "Z-score of co-occurrence")
Sapiens_z



sapiens_box<- plot_grid(gr_sapiens,Sapiens_z, ncol = 1, align = "v", rel_heights = c(2,1))
sapiens_box


Neand_df <- Neand_df[Neand_df$Time_k_BP >= 30 & Neand_df$Time_k_BP <= 50, ]
gr_neand <-ggplot(data = Neand_df, aes(x = Time_k_BP, y = Area_km / 1000000, fill=Species, group = interaction(Time_k_BP, Species))) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 1, outlier.size = 1, outlier.color = "grey") +
  theme_classic() +
  xlab("Age (kya BP)") +
  ylab(expression(paste("Area "("10" ^{6}, paste("km"^{-2})))))+
  labs(fill = "Species") +
  scale_x_reverse()+
  scale_fill_manual(values = c("3")) +
  scale_y_continuous(limits = c(1, 3), breaks = seq(0.5, 3, by = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")

gr_neand
select_Homo_neand <- species_matrix$Homo_neand == 1
species_matrix_with_Homo_neand <- species_matrix[select_Homo_neand, ]

# Sum the species' co-occurrences/Age
sum_by_age <- tapply(rowSums(species_matrix_with_Homo_neand[, -c(1:8)]), species_matrix_with_Homo_neand$Age, sum)

# Create a dataframe with output
data <- data.frame(Age = as.numeric(names(sum_by_age)),
                   cooccurrence = as.numeric(sum_by_age))
data <- data[data$Age >= 30 & data$Age <= 50, ]
data
# Bin size
bin_size <- 5
# Group by bin
data$Age_bin <- cut(data$Age, breaks = seq(50, 30, by = -bin_size))
data$Age_bin <- gsub("\\[|\\(|\\)|\\]", "", data$Age_bin)
data$Age_bin <- gsub(",", "-", data$Age_bin)
# Mean and SD
mean_value <- mean(data$cooccurrence)
standard_deviation <- sd(data$cooccurrence)
# Compute z-score
data$Z_score <- (data$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- aggregate(cooccurrence ~ Age_bin, data = data, FUN = mean, na.rm = TRUE)
mean_value_cooccurrence$Z_score <- (mean_value_cooccurrence$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- transform(mean_value_cooccurrence,
                                  Age_plot = sapply(strsplit(as.character(Age_bin), "-"), function(x) mean(as.numeric(x))))
# Plot
neand_z<-ggplot(mean_value_cooccurrence, aes(x = Age_plot, y = Z_score, group = 1)) +
  geom_line(color = "3", linetype = "dashed") +
  geom_point(color = "3", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_reverse(limits = c(50, 30)) +
  labs(x = "Age", y = "Z-score of co-occurrence")
neand_z



neand_box<- plot_grid(gr_neand,neand_z, ncol = 1, align = "v", rel_heights = c(2,1))
neand_box


Ccrocuta_df <- Ccrocuta_df[Ccrocuta_df$Time_k_BP >= 30 & Ccrocuta_df$Time_k_BP <= 50, ]
gr_Ccrocuta <-ggplot(data = Ccrocuta_df, aes(x = Time_k_BP, y = Area_km / 1000000, fill=Species, group = interaction(Time_k_BP, Species))) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 1, outlier.size = 1, outlier.color = "grey") +
  theme_classic() +
  xlab("Age (kya BP)") +
  ylab(expression(paste("Area "("10" ^{6}, paste("km"^{-2})))))+
  labs(fill = "Species") +
  scale_x_reverse()+
  scale_fill_manual(values = c("4")) +
  scale_y_continuous(limits = c(1, 3), breaks = seq(0.5, 3, by = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")

gr_Ccrocuta
select_Crocuta_crocuta <- species_matrix$Crocuta.crocuta == 1
species_matrix_with_Crocuta_crocuta <- species_matrix[select_Crocuta_crocuta, ]

# Sum the species' co-occurrences/Age
sum_by_age <- tapply(rowSums(species_matrix_with_Crocuta_crocuta[, -c(1:8)]), species_matrix_with_Crocuta_crocuta$Age, sum)

# Create a dataframe with output
data <- data.frame(Age = as.numeric(names(sum_by_age)),
                   cooccurrence = as.numeric(sum_by_age))
data <- data[data$Age >= 30 & data$Age <= 50, ]
data
# Bin size
bin_size <- 5
# Group by bin
data$Age_bin <- cut(data$Age, breaks = seq(50, 30, by = -bin_size))
data$Age_bin <- gsub("\\[|\\(|\\)|\\]", "", data$Age_bin)
data$Age_bin <- gsub(",", "-", data$Age_bin)
# Mean and SD
mean_value <- mean(data$cooccurrence)
standard_deviation <- sd(data$cooccurrence)
# Compute z-score
data$Z_score <- (data$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- aggregate(cooccurrence ~ Age_bin, data = data, FUN = mean, na.rm = TRUE)
mean_value_cooccurrence$Z_score <- (mean_value_cooccurrence$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- transform(mean_value_cooccurrence,
                                  Age_plot = sapply(strsplit(as.character(Age_bin), "-"), function(x) mean(as.numeric(x))))
# Plot
crocuta_z<-ggplot(mean_value_cooccurrence, aes(x = Age_plot, y = Z_score, group = 1)) +
  geom_line(color = "4", linetype = "dashed") +
  geom_point(color = "4", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_reverse(limits = c(50, 30)) +
  labs(x = "Age", y = "Z-score of co-occurrence")
crocuta_z

crocuta_box<- plot_grid(gr_Ccrocuta,
          crocuta_z,
          ncol = 1, align = "v", rel_heights = c(2,1))
crocuta_box



Pleo_df <- Pleo_df[Pleo_df$Time_k_BP >= 30 & Pleo_df$Time_k_BP <= 50, ]
gr_Pleo <-ggplot(data = Pleo_df, aes(x = Time_k_BP, y = Area_km / 1000000, fill=Species, group = interaction(Time_k_BP, Species))) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 1, outlier.size = 1, outlier.color = "grey") +
  theme_classic() +
  xlab("Age (kya BP)") +
  ylab(expression(paste("Area "("10" ^{6}, paste("km"^{-2})))))+
  labs(fill = "Species") +
  scale_x_reverse()+
  scale_fill_manual(values = c("5")) +
  scale_y_continuous(limits = c(1, 3.5), breaks = seq(0.5, 3, by = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")

gr_Pleo
select_Panthera_leo <- species_matrix$Panthera.leo == 1
species_matrix_with_Panthera_leo <- species_matrix[select_Panthera_leo, ]

# Sum the species' co-occurrences/Age
sum_by_age <- tapply(rowSums(species_matrix_with_Panthera_leo[, -c(1:8)]), species_matrix_with_Panthera_leo$Age, sum)

# Create a dataframe with output
data <- data.frame(Age = as.numeric(names(sum_by_age)),
                   cooccurrence = as.numeric(sum_by_age))
data <- data[data$Age >= 30 & data$Age <= 50, ]
data
# Bin size
bin_size <- 5
# Group by bin
data$Age_bin <- cut(data$Age, breaks = seq(50, 30, by = -bin_size))
data$Age_bin <- gsub("\\[|\\(|\\)|\\]", "", data$Age_bin)
data$Age_bin <- gsub(",", "-", data$Age_bin)
# Mean and SD
mean_value <- mean(data$cooccurrence)
standard_deviation <- sd(data$cooccurrence)
# Compute z-score
data$Z_score <- (data$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- aggregate(cooccurrence ~ Age_bin, data = data, FUN = mean, na.rm = TRUE)
mean_value_cooccurrence$Z_score <- (mean_value_cooccurrence$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- transform(mean_value_cooccurrence,
                                  Age_plot = sapply(strsplit(as.character(Age_bin), "-"), function(x) mean(as.numeric(x))))
# Plot
leo_z<-ggplot(mean_value_cooccurrence, aes(x = Age_plot, y = Z_score, group = 1)) +
  geom_line(color = "5", linetype = "dashed") +
  geom_point(color = "5", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_reverse(limits = c(50, 30)) +
  labs(x = "Age", y = "Z-score of co-occurrence")
leo_z


pleo_box<- plot_grid(gr_Pleo,leo_z,ncol = 1, align = "v", rel_heights = c(2,1))
pleo_box



Pspelaea_df <- Pspelaea_df[Pspelaea_df$Time_k_BP >= 30 & Pspelaea_df$Time_k_BP <= 50, ]
gr_Pspelaea <-ggplot(data = Pspelaea_df, aes(x = Time_k_BP, y = Area_km / 1000000, fill=Species, group = interaction(Time_k_BP, Species))) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 1, outlier.size = 1, outlier.color = "grey") +
  theme_classic() +
  xlab("Age (kya BP)") +
  ylab(expression(paste("Area "("10" ^{6}, paste("km"^{-2})))))+
  labs(fill = "Species") +
  scale_x_reverse()+
  scale_fill_manual(values = c("6")) +
  scale_y_continuous(limits = c(1, 2.5), breaks = seq(0.5, 3, by = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")

gr_Pspelaea
select_Panthera_spelaea <- species_matrix$Panthera.spelaea == 1
species_matrix_with_Panthera_spelaea <- species_matrix[select_Panthera_spelaea, ]

# Sum the species' co-occurrences/Age
sum_by_age <- tapply(rowSums(species_matrix_with_Panthera_spelaea[, -c(1:8)]), species_matrix_with_Panthera_spelaea$Age, sum)

# Create a dataframe with output
data <- data.frame(Age = as.numeric(names(sum_by_age)),
                   cooccurrence = as.numeric(sum_by_age))
data <- data[data$Age >= 30 & data$Age <= 50, ]
data
# Bin size
bin_size <- 5
# Group by bin
data$Age_bin <- cut(data$Age, breaks = seq(50, 30, by = -bin_size))
data$Age_bin <- gsub("\\[|\\(|\\)|\\]", "", data$Age_bin)
data$Age_bin <- gsub(",", "-", data$Age_bin)
# Mean and SD
mean_value <- mean(data$cooccurrence)
standard_deviation <- sd(data$cooccurrence)
# Compute z-score
data$Z_score <- (data$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- aggregate(cooccurrence ~ Age_bin, data = data, FUN = mean, na.rm = TRUE)
mean_value_cooccurrence$Z_score <- (mean_value_cooccurrence$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- transform(mean_value_cooccurrence,
                                  Age_plot = sapply(strsplit(as.character(Age_bin), "-"), function(x) mean(as.numeric(x))))
# Plot
spelaea_z<-ggplot(mean_value_cooccurrence, aes(x = Age_plot, y = Z_score, group = 1)) +
  geom_line(color = "6", linetype = "dashed") +
  geom_point(color = "6", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_reverse(limits = c(50, 30)) +
  labs(x = "Age", y = "Z-score of co-occurrence")
spelaea_z

pspelaea_box<- plot_grid(gr_Pspelaea,spelaea_z, ncol = 1, align = "v", rel_heights = c(2,1))
pspelaea_box

Ppardus_df <- Ppardus_df[Ppardus_df$Time_k_BP >= 30 & Ppardus_df$Time_k_BP <= 50, ]
gr_Ppardus <-ggplot(data = Ppardus_df, aes(x = Time_k_BP, y = Area_km / 1000000, fill=Species, group = interaction(Time_k_BP, Species))) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 1, outlier.size = 1, outlier.color = "grey") +
  theme_classic() +
  xlab("Age (kya BP)") +
  ylab(expression(paste("Area "("10" ^{6}, paste("km"^{-2})))))+
  labs(fill = "Species") +
  scale_x_reverse()+
  scale_fill_manual(values = c("7")) +
  scale_y_continuous(limits = c(0.5, 2.5), breaks = seq(0.5, 3, by = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")

gr_Ppardus
select_Panthera_pardus <- species_matrix$Panthera.pardus == 1
species_matrix_with_Panthera_pardus <- species_matrix[select_Panthera_pardus, ]

# Sum the species' co-occurrences/Age
sum_by_age <- tapply(rowSums(species_matrix_with_Panthera_pardus[, -c(1:8)]), species_matrix_with_Panthera_pardus$Age, sum)

# Create a dataframe with outputs
data <- data.frame(Age = as.numeric(names(sum_by_age)),
                   cooccurrence = as.numeric(sum_by_age))
data <- data[data$Age >= 30 & data$Age <= 50, ]
data
# Bin size
bin_size <- 5
# Group by bin
data$Age_bin <- cut(data$Age, breaks = seq(50, 30, by = -bin_size))
data$Age_bin <- gsub("\\[|\\(|\\)|\\]", "", data$Age_bin)
data$Age_bin <- gsub(",", "-", data$Age_bin)
# Mean and SD
mean_value <- mean(data$cooccurrence)
standard_deviation <- sd(data$cooccurrence)
# Compute z-score
data$Z_score <- (data$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- aggregate(cooccurrence ~ Age_bin, data = data, FUN = mean, na.rm = TRUE)
mean_value_cooccurrence$Z_score <- (mean_value_cooccurrence$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- transform(mean_value_cooccurrence,
                                  Age_plot = sapply(strsplit(as.character(Age_bin), "-"), function(x) mean(as.numeric(x))))
# Plot
pardus_z<-ggplot(mean_value_cooccurrence, aes(x = Age_plot, y = Z_score, group = 1)) +
  geom_line(color = "7", linetype = "dashed") +
  geom_point(color = "7", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_reverse(limits = c(50, 30)) +
  labs(x = "Age", y = "Z-score of co-occurrence")
pardus_z


ppardus_box<- plot_grid(gr_Ppardus,pardus_z, ncol = 1, align = "v", rel_heights = c(2,1))   
ppardus_box


Vvulpes_df <- Vvulpes_df[Vvulpes_df$Time_k_BP >= 30 & Vvulpes_df$Time_k_BP <= 50, ]
gr_Vvulpes <-ggplot(data = Vvulpes_df, aes(x = Time_k_BP, y = Area_km / 1000000, fill=Species, group = interaction(Time_k_BP, Species))) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 1, outlier.size = 1, outlier.color = "grey") +
  theme_classic() +
  xlab("Age (kya BP)") +
  ylab(expression(paste("Area "("10" ^{6}, paste("km"^{-2})))))+
  labs(fill = "Species") +
  ggtitle(expression(italic("Vulpes vulpes")))+
  scale_x_reverse()+
  scale_fill_manual(values = c("8")) +
  scale_y_continuous(limits = c(1, 3), breaks = seq(0.5, 3, by = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")

gr_Vvulpes
select_Vulpes_vulpes <- species_matrix$Vulpes.vulpes == 1
species_matrix_with_Vulpes_vulpes <- species_matrix[select_Vulpes_vulpes, ]

# Sum the species' co-occurrences/Age
sum_by_age <- tapply(rowSums(species_matrix_with_Vulpes_vulpes[, -c(1:8)]), species_matrix_with_Vulpes_vulpes$Age, sum)

# Create a dataframe with output
data <- data.frame(Age = as.numeric(names(sum_by_age)),
                   cooccurrence = as.numeric(sum_by_age))
data <- data[data$Age >= 30 & data$Age <= 50, ]
data
# Bin size
bin_size <- 5
# Group by bin
data$Age_bin <- cut(data$Age, breaks = seq(50, 30, by = -bin_size))
data$Age_bin <- gsub("\\[|\\(|\\)|\\]", "", data$Age_bin)
data$Age_bin <- gsub(",", "-", data$Age_bin)
# Mean and SD
mean_value <- mean(data$cooccurrence)
standard_deviation <- sd(data$cooccurrence)
# Compute z-score
data$Z_score <- (data$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- aggregate(cooccurrence ~ Age_bin, data = data, FUN = mean, na.rm = TRUE)
mean_value_cooccurrence$Z_score <- (mean_value_cooccurrence$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- transform(mean_value_cooccurrence,
                                  Age_plot = sapply(strsplit(as.character(Age_bin), "-"), function(x) mean(as.numeric(x))))
# Plot
vulpes_z<-ggplot(mean_value_cooccurrence, aes(x = Age_plot, y = Z_score, group = 1)) +
  geom_line(color = "8", linetype = "dashed") +
  geom_point(color = "8", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_reverse(limits = c(50, 30)) +
  labs(x = "Age", y = "Z-score of co-occurrence")
vulpes_z

vvulpes_box<- plot_grid(gr_Vvulpes,vulpes_z,ncol = 1, align = "v", rel_heights = c(2,1))
vvulpes_box

Vlagopus_df <- Vlagopus_df[Vlagopus_df$Time_k_BP >= 30 & Vlagopus_df$Time_k_BP <= 50, ]
gr_Vlagopus <-ggplot(data = Vlagopus_df, aes(x = Time_k_BP, y = Area_km / 1000000, fill=Species, group = interaction(Time_k_BP, Species))) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 1, outlier.size = 1, outlier.color = "grey") +
  theme_classic() +
  xlab("Age (kya BP)") +
  ylab(expression(paste("Area "("10" ^{6}, paste("km"^{-2})))))+
  labs(fill = "Species") +
  scale_x_reverse()+
  scale_fill_manual(values = c("10")) +
  scale_y_continuous(limits = c(0.5, 3.5), breaks = seq(0.5, 3.5, by = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")

gr_Vlagopus
select_Vulpes_lagopus <- species_matrix$Vulpes.lagopus == 1
species_matrix_with_Vulpes_lagopus <- species_matrix[select_Vulpes_lagopus, ]

# Sum the species' co-occurrences/Age
sum_by_age <- tapply(rowSums(species_matrix_with_Vulpes_lagopus[, -c(1:8)]), species_matrix_with_Vulpes_lagopus$Age, sum)

# Create a dataframe with output
data <- data.frame(Age = as.numeric(names(sum_by_age)),
                   cooccurrence = as.numeric(sum_by_age))
data <- data[data$Age >= 30 & data$Age <= 50, ]
data
# Bin size
bin_size <- 5
# Group by bin
data$Age_bin <- cut(data$Age, breaks = seq(50, 30, by = -bin_size))
data$Age_bin <- gsub("\\[|\\(|\\)|\\]", "", data$Age_bin)
data$Age_bin <- gsub(",", "-", data$Age_bin)
# Mean and SD
mean_value <- mean(data$cooccurrence)
standard_deviation <- sd(data$cooccurrence)
# Compute z-score
data$Z_score <- (data$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- aggregate(cooccurrence ~ Age_bin, data = data, FUN = mean, na.rm = TRUE)
mean_value_cooccurrence$Z_score <- (mean_value_cooccurrence$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- transform(mean_value_cooccurrence,
                                  Age_plot = sapply(strsplit(as.character(Age_bin), "-"), function(x) mean(as.numeric(x))))
# Plot
lagopus_z<-ggplot(mean_value_cooccurrence, aes(x = Age_plot, y = Z_score, group = 1)) +
  geom_line(color = "10", linetype = "dashed") +
  geom_point(color = "10", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_reverse(limits = c(50, 30)) +
  labs(x = "Age", y = "Z-score of co-occurrence")
lagopus_z

vlagopus_box<-plot_grid(gr_Vlagopus,lagopus_z, ncol = 1, align = "v", rel_heights = c(2,1,1))
vlagopus_box


Calupus_df <- Calupus_df[Calupus_df$Time_k_BP >= 30 & Calupus_df$Time_k_BP <= 50, ]
gr_Calupus <-ggplot(data = Calupus_df, aes(x = Time_k_BP, y = Area_km / 1000000, fill=Species, group = interaction(Time_k_BP, Species))) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 1, outlier.size = 1, outlier.color = "grey") +
  theme_classic() +
  xlab("Age (kya BP)") +
  ylab(expression(paste("Area "("10" ^{6}, paste("km"^{-2})))))+
  labs(fill = "Species") +
  scale_x_reverse()+
  scale_fill_manual(values = c("11")) +
  scale_y_continuous(limits = c(1, 3.2), breaks = seq(0.5, 3, by = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")

gr_Calupus

select_Canis_lupus <- species_matrix$Canis.lupus == 1
species_matrix_with_Canis_lupus <- species_matrix[select_Canis_lupus, ]

# Sum the species' co-occurrences/Age
sum_by_age <- tapply(rowSums(species_matrix_with_Canis_lupus[, -c(1:8)]), species_matrix_with_Canis_lupus$Age, sum)

# Create a dataframe with output
data <- data.frame(Age = as.numeric(names(sum_by_age)),
                   cooccurrence = as.numeric(sum_by_age))
data <- data[data$Age >= 30 & data$Age <= 50, ]
data
# Bin size
bin_size <- 5
# Group by bin
data$Age_bin <- cut(data$Age, breaks = seq(50, 30, by = -bin_size))
data$Age_bin <- gsub("\\[|\\(|\\)|\\]", "", data$Age_bin)
data$Age_bin <- gsub(",", "-", data$Age_bin)
# Mean and SD
mean_value <- mean(data$cooccurrence)
standard_deviation <- sd(data$cooccurrence)
# Compute z-score
data$Z_score <- (data$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- aggregate(cooccurrence ~ Age_bin, data = data, FUN = mean, na.rm = TRUE)
mean_value_cooccurrence$Z_score <- (mean_value_cooccurrence$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- transform(mean_value_cooccurrence,
                                  Age_plot = sapply(strsplit(as.character(Age_bin), "-"), function(x) mean(as.numeric(x))))
# Plot
lupus_z<-ggplot(mean_value_cooccurrence, aes(x = Age_plot, y = Z_score, group = 1)) +
  geom_line(color = "11", linetype = "dashed") +
  geom_point(color = "11", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_reverse(limits = c(50, 30)) +
  labs(x = "Age", y = "Z-score of co-occurrence")
lupus_z

calupus_box<- plot_grid(gr_Calupus,lupus_z, ncol = 1, align = "v", rel_heights = c(2,1))
calupus_box

Uarctos_df <- Uarctos_df[Uarctos_df$Time_k_BP >= 30 & Uarctos_df$Time_k_BP <= 50, ]
gr_Uarctos <-ggplot(data = Uarctos_df, aes(x = Time_k_BP, y = Area_km / 1000000, fill=Species, group = interaction(Time_k_BP, Species))) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 1, outlier.size = 1, outlier.color = "grey") +
  theme_classic() +
  xlab("Age (kya BP)") +
  ylab(expression(paste("Area "("10" ^{6}, paste("km"^{-2})))))+
  labs(fill = "Species") +
  scale_x_reverse()+
  scale_fill_manual(values = c("12")) +
  scale_y_continuous(limits = c(1, 3.5), breaks = seq(0.5, 3.5, by = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")

gr_Uarctos
select_Ursus_arctos <- species_matrix$Ursus.arctos == 1
species_matrix_with_Ursus_arctos <- species_matrix[select_Ursus_arctos, ]

# Sum the species' co-occurrences/Age
sum_by_age <- tapply(rowSums(species_matrix_with_Ursus_arctos[, -c(1:8)]), species_matrix_with_Ursus_arctos$Age, sum)

# Create a dataframe with output
data <- data.frame(Age = as.numeric(names(sum_by_age)),
                   cooccurrence = as.numeric(sum_by_age))
data <- data[data$Age >= 30 & data$Age <= 50, ]
data
# Bin size
bin_size <- 5
# Group by bin
data$Age_bin <- cut(data$Age, breaks = seq(50, 30, by = -bin_size))
data$Age_bin <- gsub("\\[|\\(|\\)|\\]", "", data$Age_bin)
data$Age_bin <- gsub(",", "-", data$Age_bin)
# Mean and SD
mean_value <- mean(data$cooccurrence)
standard_deviation <- sd(data$cooccurrence)
# Compute z-score
data$Z_score <- (data$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- aggregate(cooccurrence ~ Age_bin, data = data, FUN = mean, na.rm = TRUE)
mean_value_cooccurrence$Z_score <- (mean_value_cooccurrence$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- transform(mean_value_cooccurrence,
                                  Age_plot = sapply(strsplit(as.character(Age_bin), "-"), function(x) mean(as.numeric(x))))
# Plot
arctos_z<-ggplot(mean_value_cooccurrence, aes(x = Age_plot, y = Z_score, group = 1)) +
  geom_line(color = "12", linetype = "dashed") +
  geom_point(color = "12", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_reverse(limits = c(50, 30)) +
  labs(x = "Age", y = "Z-score of co-occurrence")
arctos_z

uarctos_box<-plot_grid(gr_Uarctos,arctos_z, ncol = 1, align = "v", rel_heights = c(2,1))
uarctos_box

Uspelaeus_df <- Uspelaeus_df[Uspelaeus_df$Time_k_BP >= 30 & Uspelaeus_df$Time_k_BP <= 50, ]
gr_Uspelaeus <-ggplot(data = Uspelaeus_df, aes(x = Time_k_BP, y = Area_km / 1000000, fill=Species, group = interaction(Time_k_BP, Species))) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 1, outlier.size = 1, outlier.color = "grey") +
  theme_classic() +
  xlab("Age (kya BP)") +
  ylab(expression(paste("Area "("10" ^{6}, paste("km"^{-2})))))+
  labs(fill = "Species") +
  scale_x_reverse()+
  scale_fill_manual(values = c("13")) +
  scale_y_continuous(limits = c(1, 2.5), breaks = seq(0.5, 3.5, by = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")

gr_Uspelaeus
select_Ursus_spelaeus <- species_matrix$Ursus.spelaeus == 1
species_matrix_with_Ursus_spelaeus <- species_matrix[select_Ursus_spelaeus, ]

# Sum the species' co-occurrences/Age
sum_by_age <- tapply(rowSums(species_matrix_with_Ursus_spelaeus[, -c(1:8)]), species_matrix_with_Ursus_spelaeus$Age, sum)
Ursus_counts_por_edad <- tapply(species_matrix_with_Ursus_spelaeus$Ursus.spelaeus, species_matrix_with_Ursus_spelaeus$Age, sum)

# Create a dataframe with output
data <- data.frame(Age = as.numeric(names(sum_by_age)),
                   cooccurrence = as.numeric(sum_by_age))

data <- data[data$Age >= 30 & data$Age <= 50, ]

# Bin size
bin_size <- 5
# Group by bin
data$Age_bin <- cut(data$Age, breaks = seq(50, 30, by = -bin_size))
data$Age_bin <- gsub("\\[|\\(|\\)|\\]", "", data$Age_bin)
data$Age_bin <- gsub(",", "-", data$Age_bin)
# Mean and SD
mean_value <- mean(data$cooccurrence)
standard_deviation <- sd(data$cooccurrence)
# Compute z-score
data$Z_score <- (data$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- aggregate(cooccurrence ~ Age_bin, data = data, FUN = mean, na.rm = TRUE)
mean_value_cooccurrence$Z_score <- (mean_value_cooccurrence$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- transform(mean_value_cooccurrence,
                                  Age_plot = sapply(strsplit(as.character(Age_bin), "-"), function(x) mean(as.numeric(x))))
# Plot
spelaeus_z<-ggplot(mean_value_cooccurrence, aes(x = Age_plot, y = Z_score, group = 1)) +
  geom_line(color = "13", linetype = "dashed") +
  geom_point(color = "13", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_reverse(limits = c(50, 30)) +
  labs(x = "Age", y = "Z-score of co-occurrence")
spelaeus_z

uspelaeus_box<- plot_grid(gr_Uspelaeus,spelaeus_z,ncol = 1, align = "v", rel_heights = c(2,1))
uspelaeus_box


Fsilvestris_df <- Fsilvestris_df[Fsilvestris_df$Time_k_BP >= 30 & Fsilvestris_df$Time_k_BP <= 50, ]
gr_Fsilvestris <-ggplot(data = Fsilvestris_df, aes(x = Time_k_BP, y = Area_km / 1000000, fill=Species, group = interaction(Time_k_BP, Species))) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 1, outlier.size = 1, outlier.color = "grey") +
  theme_classic() +
  xlab("Age (kya BP)") +
  ylab(expression(paste("Area "("10" ^{6}, paste("km"^{-2})))))+
  labs(fill = "Species") +
  scale_x_reverse()+
  scale_fill_manual(values = c("14")) +
  scale_y_continuous(limits = c(1, 2.5), breaks = seq(0.5, 3.5, by = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")

gr_Fsilvestris
select_Felis_silvestris <- species_matrix$Felis.silvestris == 1
species_matrix_with_Felis_silvestris <- species_matrix[select_Felis_silvestris, ]

# Sum the species' co-occurrences/Age
sum_by_age <- tapply(rowSums(species_matrix_with_Felis_silvestris[, -c(1:8)]), species_matrix_with_Felis_silvestris$Age, sum)

# Create a dataframe with output
data <- data.frame(Age = as.numeric(names(sum_by_age)),
                   cooccurrence = as.numeric(sum_by_age))
data <- data[data$Age >= 30 & data$Age <= 50, ]
data
# Bin size
bin_size <- 5
# Group by bin
data$Age_bin <- cut(data$Age, breaks = seq(50, 30, by = -bin_size))
data$Age_bin <- gsub("\\[|\\(|\\)|\\]", "", data$Age_bin)
data$Age_bin <- gsub(",", "-", data$Age_bin)
# Calcular la mean_value y la desviación estándar de las co-ocurrencias de Felis silvestris
mean_value <- mean(data$cooccurrence)
standard_deviation <- sd(data$cooccurrence)
# Compute z-score
data$Z_score <- (data$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- aggregate(cooccurrence ~ Age_bin, data = data, FUN = mean, na.rm = TRUE)
mean_value_cooccurrence$Z_score <- (mean_value_cooccurrence$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- transform(mean_value_cooccurrence,
                                  Age_plot = sapply(strsplit(as.character(Age_bin), "-"), function(x) mean(as.numeric(x))))
# Plot
silvestris_z<-ggplot(mean_value_cooccurrence, aes(x = Age_plot, y = Z_score, group = 1)) +
  geom_line(color = "14", linetype = "dashed") +
  geom_point(color = "14", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_reverse(limits = c(50, 30)) +
  labs(x = "Age", y = "Z-score of co-occurrence")
silvestris_z

fsilvestris_box<- plot_grid(gr_Fsilvestris,silvestris_z, ncol = 1, align = "v", rel_heights = c(2,1))
fsilvestris_box

Mmeles_df <- Mmeles_df[Mmeles_df$Time_k_BP >= 30 & Mmeles_df$Time_k_BP <= 50, ]
gr_Mmeles <-ggplot(data = Mmeles_df, aes(x = Time_k_BP, y = Area_km / 1000000, fill=Species, group = interaction(Time_k_BP, Species))) +
  geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 1, outlier.size = 1, outlier.color = "grey") +
  theme_classic() +
  xlab("Age (kya BP)") +
  ylab(expression(paste("Area "("10" ^{6}, paste("km"^{-2})))))+
  labs(fill = "Species") +
  scale_x_reverse()+
  scale_fill_manual(values = c("15")) +
  scale_y_continuous(limits = c(1, 2.5), breaks = seq(0.5, 3.5, by = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")

gr_Mmeles
select_Meles_meles <- species_matrix$Meles.meles == 1
species_matrix_with_Meles_meles <- species_matrix[select_Meles_meles, ]

# Sum the species' co-occurrences/Age
sum_by_age <- tapply(rowSums(species_matrix_with_Meles_meles[, -c(1:8)]), species_matrix_with_Meles_meles$Age, sum)

# Create a dataframe with output
data <- data.frame(Age = as.numeric(names(sum_by_age)),
                   cooccurrence = as.numeric(sum_by_age))
data <- data[data$Age >= 30 & data$Age <= 50, ]
data
# Bin size
bin_size <- 5
# Group by bin
data$Age_bin <- cut(data$Age, breaks = seq(50, 30, by = -bin_size))
data$Age_bin <- gsub("\\[|\\(|\\)|\\]", "", data$Age_bin)
data$Age_bin <- gsub(",", "-", data$Age_bin)
# Calcular la mean_value y la desviación estándar de las co-ocurrencias de Meles meles
mean_value <- mean(data$cooccurrence)
standard_deviation <- sd(data$cooccurrence)
# Compute z-score
data$Z_score <- (data$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- aggregate(cooccurrence ~ Age_bin, data = data, FUN = mean, na.rm = TRUE)
mean_value_cooccurrence$Z_score <- (mean_value_cooccurrence$cooccurrence - mean_value) / standard_deviation
mean_value_cooccurrence <- transform(mean_value_cooccurrence,
                                  Age_plot = sapply(strsplit(as.character(Age_bin), "-"), function(x) mean(as.numeric(x))))
# Plot
meles_z<-ggplot(mean_value_cooccurrence, aes(x = Age_plot, y = Z_score, group = 1)) +
  geom_line(color = "15", linetype = "dashed") +
  geom_point(color = "15", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_reverse(limits = c(50, 30)) +
  labs(x = "Age", y = "Z-score of co-occurrence")
meles_z


mmeles_box<- plot_grid(gr_Mmeles,meles_z,ncol = 1, align = "v", rel_heights = c(2,1,1))
mmeles_box






##### HOVMOLLER PLOTS WITH AREA FAVOURABILITY: ####

# GEOGRAPHIC EXTENT FOR HOVMOLLER DIAGRAMS
mywindow <- c(-15, 50, 33, 60)  # xmin, xmax, ymin, ymax
ext <- extent(mywindow)

# COUNTRIES
mywindow <- c(-15, 60, 33, 75)  # xmin, xmax, ymin, ymax

# Extract the subset of the world map data
world_map_subset <- map_data("world") %>%
  filter(long >= mywindow[1] & long <= mywindow[2] & lat >= mywindow[3] & lat <= mywindow[4] &
           region != "Morocco" & region != "Algeria" & region != "Tunisia")


getwd()



######## SAPIENS ####

mi_paleta <- colorRampPalette(brewer.pal(11, "RdYlBu"))
clrs <- rev(mi_paleta(100))
species <- "Homo.sapiens"

my_path<- " " # Your directory
my_path<- paste0(my_path, "/", species, "/Layers/TIF")
layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)
Sapiens_layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)
# Create a raster stack from the list of raster files
Sapiens_rasters <- stack(Sapiens_layers)
z_values <- c("55","54","53","52","51","50","49","48",
  "47","46","45","44","43","42","41","40","39","38","37",
              "36","35","34","33","32", "31", "30", "29", "28", "27")
Sapiens_rasters <- setZ(Sapiens_rasters, z_values)

# Check the information of the raster stack
print(Sapiens_rasters)

# Extent of RasterStack
Sapiens_rasters_c <- crop(Sapiens_rasters, ext)

# Use the hovmoller function with the raster stack
dirLayer <- init(Sapiens_rasters_c, v='x')
z <- zonal(Sapiens_rasters_c, dirLayer, FUN='mean', digits=2)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)
limits <- seq(0, 0.8, by = 0.1)

# Plot
Longitude_sapiens<- contourplot(z ~ x * y, data = dat,
            xlab = 'Chronology (kya)', ylab = 'Longitude (º)',
            panel = panel.levelplot.raster,
            colorkey = TRUE, col.regions = viridis,
            scales = list(x = list(at = c(1, 6, 11, 16, 21, 26)))
            )

dirLayer <- init(Sapiens_rasters_c, v='y')
z <- zonal(Sapiens_rasters_c, dirLayer, FUN='mean', digits=2)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Latitude_sapiens<- contourplot(z ~ x*y, data=dat,
                              xlab='Chronology (kya)', ylab='Latitude (º)',
                              panel=panel.levelplot.raster,
                              colorkey = TRUE, col.regions = viridis,
                              scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))


# Difference between Phase I (50-45 kyr BP) and Phase III (40-35 kyr BP)

path_TIF <- "" # write the directory where the TIF files are located (e.g. path_TIF <- "F:/Data$Code/Species/Homo.sapiens/Layers/TIF/")

Sapiens_50_45 <- c(
  paste0(path_TIF, "50f.tif"),
  paste0(path_TIF, "49f.tif"),
  paste0(path_TIF, "48f.tif"),
  paste0(path_TIF, "47f.tif"),
  paste0(path_TIF, "46f.tif"),
  paste0(path_TIF, "45f.tif")
)
Sapiens_50_45 <- mean(stack(Sapiens_50_45))

Sapiens_40_35 <- c(
  paste0(path_TIF, "40f.tif"),
  paste0(path_TIF, "39f.tif"),
  paste0(path_TIF, "38f.tif"),
  paste0(path_TIF, "37f.tif"),
  paste0(path_TIF, "36f.tif"),
  paste0(path_TIF, "35f.tif")
)
Sapiens_40_35 <- mean(stack(Sapiens_40_35))
# Compute difference:
Dif_sapiens<-  Sapiens_40_35 - Sapiens_50_45
plot(Dif_sapiens)
# Load raster files
raster_stack_50 <- stack(Sapiens_50_45)
raster_stack_40 <- stack(Sapiens_40_35)
raster_stack_dif <- stack(Dif_sapiens)

# Convert raster
raster_df_s_50 <- as.data.frame(raster_stack_50, xy = TRUE)
raster_df_s_40 <- as.data.frame(raster_stack_40 , xy = TRUE)
raster_df_s_dif <- as.data.frame(raster_stack_dif, xy = TRUE)



# Plot
sapiens_b_MIS3<-ggplot(raster_df_s_50, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", linewidth=0.1)

sapiens_b_MIS3
sapiens_e_MIS3<-ggplot(raster_df_s_40, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")
sapiens_dif<-ggplot(raster_df_s_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), limits = c( -0.4, 0.5), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size=0.1)

plot(sapiens_dif)

sapiens_dif_m <- sapiens_dif + theme(plot.margin = unit(c(1,1,1,0.1), "cm"))

# Plot:
ggarrange( sapiens_box, Latitude_sapiens,Longitude_sapiens, sapiens_dif_m, ncol = 4)



#### NEANDERTHALS ####

species <- "Homo_neand"
my_path<- " "
my_path<- paste0(my_path, "/", species, "/Layers/TIF")
layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)
Neand_layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)
Neand_layers

# Create a raster stack from the list of raster files
Neand_rasters <- stack(Neand_layers)


Neand_rasters <- setZ(Neand_rasters, z_values)

Neand_rasters_c <- crop(Neand_rasters, ext)

# Use the hovmoller function with the raster stack

dirLayer <- init(Neand_rasters_c, v='x')
z <- zonal(Neand_rasters_c, dirLayer, FUN='mean', digits=2)
head(dat)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Longitude_neand<- contourplot(z ~ x*y, data=dat,
                              xlab='Chronology (kya)', ylab='Longitude (º)',
                              panel=panel.levelplot.raster,
                              colorkey = TRUE,
                              col.regions = viridis,
                              scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

dirLayer <- init(Neand_rasters_c, v='y')
z <- zonal(Neand_rasters_c, dirLayer, FUN='mean', digits=2)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Latitude_neand<- contourplot(z ~ x*y, data=dat,
                             xlab='Chronology (kya)', ylab='Latitude (º)',
                             panel=panel.levelplot.raster,
                             colorkey = TRUE,
                             col.regions = viridis,
                             scales = list(x = list(at = c(1, 6, 11, 16, 21, 26)))
                             )

# Difference between Phase I (50-45 kr BP) and Phase III (40-35 kyr BP)
path_TIF <- " " # write the directory where the TIF files are located (e.g. path_TIF <- "F:/Data$Code/Species/Homo.sapiens/Layers/TIF/")

Neand_50_45 <- c(
  paste0(path_TIF, "50f.tif"),
  paste0(path_TIF, "49f.tif"),
  paste0(path_TIF, "48f.tif"),
  paste0(path_TIF, "47f.tif"),
  paste0(path_TIF, "46f.tif"),
  paste0(path_TIF, "45f.tif")
)
Neand_50_45 <- mean(stack(Neand_50_45))

Neand_40_35 <- c(
  paste0(path_TIF, "40f.tif"),
  paste0(path_TIF, "39f.tif"),
  paste0(path_TIF, "38f.tif"),
  paste0(path_TIF, "37f.tif"),
  paste0(path_TIF, "36f.tif"),
  paste0(path_TIF, "35f.tif")
)
Neand_40_35 <- mean(stack(Neand_40_35))
Dif_neand<-  Neand_40_35 - Neand_50_45

# Load raster files
raster_stack_50 <- stack(Neand_50_45)
raster_stack_40 <- stack(Neand_40_35)
raster_stack_dif <- stack(Dif_neand)

# Convert data
raster_df_n_50 <- as.data.frame(raster_stack_50, xy = TRUE)
rraster_df_n_40 <- as.data.frame(raster_stack_40 , xy = TRUE)
raster_df_n_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

# Plot
neand_b_MIS3<-ggplot(raster_df_n_50, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")

neand_b_MIS3
neand_e_MIS3<-ggplot(rraster_df_n_40, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")
neand_dif<-ggplot(raster_df_n_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), limits = c(-0.4, 0.5), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", linewidth = 0.1)
neand_dif


neand_dif_m <- neand_dif + theme(plot.margin = unit(c(1,1,1,0.1), "cm"))

ggarrange(neand_box, Latitude_neand,Longitude_neand, neand_dif_m,ncol = 4)


#### URSUS SPELAEUS ####

species <- "Ursus.spelaeus"
my_path<- " " # Your directory
my_path<- paste0(my_path, "/", species, "/Layers/TIF")
layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)
Uspelaeus_layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)

# Create a raster stack from the list of raster files
Uspelaeus_rasters <- stack(Uspelaeus_layers)

Uspelaeus_rasters <- setZ(Uspelaeus_rasters, z_values)

Uspelaeus_rasters_c <- crop(Uspelaeus_rasters, ext)

# Use the hovmoller function with the raster stack

dirLayer <- init(Uspelaeus_rasters_c, v='x')
z <- zonal(Uspelaeus_rasters_c, dirLayer, FUN='mean', digits=2)
head(dat)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Longitude_Uspelaeus<- contourplot(z ~ x*y, data=dat,
                                  xlab='Chronology (kya)', ylab='Longitude (º)',
                                  panel=panel.levelplot.raster,
                                  colorkey = TRUE,
                                  col.regions = viridis,
                                  scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

dirLayer <- init(Uspelaeus_rasters_c, v='y')
z <- zonal(Uspelaeus_rasters_c, dirLayer, FUN='mean', digits=2)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Latitude_Uspelaeus<- contourplot(z ~ x*y, data=dat,
                                 xlab='Chronology (kya)', ylab='Latitude (º)',
                                 panel=panel.levelplot.raster,
                                 colorkey = TRUE,
                                 col.regions = viridis,
                                 scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

# Difference between Phase I (50-45 kr BP) and Phase III (40-35 kyr BP)
path_TIF <- " " # write the directory where the TIF files are located (e.g. path_TIF <- "F:/Data$Code/Species/Homo.sapiens/Layers/TIF/")


Uspelaeus_50_45 <- c(
  paste0(path_TIF, "50f.tif"),
  paste0(path_TIF, "49f.tif"),
  paste0(path_TIF, "48f.tif"),
  paste0(path_TIF, "47f.tif"),
  paste0(path_TIF, "46f.tif"),
  paste0(path_TIF, "45f.tif")
)
Uspelaeus_50_45 <- mean(stack(Uspelaeus_layers))

Uspelaeus_40_35 <- c(
  paste0(path_TIF, "40f.tif"),
  paste0(path_TIF, "39f.tif"),
  paste0(path_TIF, "38f.tif"),
  paste0(path_TIF, "37f.tif"),
  paste0(path_TIF, "36f.tif"),
  paste0(path_TIF, "35f.tif")
)
Uspelaeus_40_35 <- mean(stack(Uspelaeus_40_35))
Dif_Uspelaeus<-  Uspelaeus_40_35 - Uspelaeus_50_45

# Load raster files
raster_stack_50 <- stack(Uspelaeus_layers)
raster_stack_40 <- stack(Uspelaeus_40_35)
raster_stack_dif <- stack(Dif_Uspelaeus)

# Convert data
raster_df_n_50 <- as.data.frame(raster_stack_50, xy = TRUE)
rraster_df_n_40 <- as.data.frame(raster_stack_40 , xy = TRUE)
raster_df_n_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

# Plot
Uspelaeus_b_MIS3<-ggplot(raster_df_n_50, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")

Uspelaeus_b_MIS3
Uspelaeus_e_MIS3<-ggplot(rraster_df_n_40, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")
Uspelaeus_e_MIS3
Uspelaeus_dif<-ggplot(raster_df_n_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)),limits = c(-0.4, 0.5), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", linewidth = 0.1)

plot(Uspelaeus_dif)

Uspelaeus_dif_m <- Uspelaeus_dif + theme(plot.margin = unit(c(1,1,1,0.1), "cm"))

ggarrange( Latitude_Uspelaeus,Longitude_Uspelaeus, Uspelaeus_dif_m,ncol = 3)


#### URSUS ARCTOS ####

species <- "Ursus.arctos"
my_path<- " " # Your directory
my_path<- paste0(my_path, "/", species, "/Layers/TIF")
layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)
Uarctos_layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)

# Create a raster stack from the list of raster files
Uarctos_rasters <- stack(Uarctos_layers)

Uarctos_rasters <- setZ(Uarctos_rasters, z_values)

Uarctos_rasters_c <- crop(Uarctos_rasters, ext)

# Use the hovmoller function with the raster stack

dirLayer <- init(Uarctos_rasters_c, v='x')
z <- zonal(Uarctos_rasters_c, dirLayer, FUN='mean', digits=2)
head(dat)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Longitude_Uarctos<- contourplot(z ~ x*y, data=dat,
                                xlab='Chronology (kya)', ylab='Longitude (º)',
                                panel=panel.levelplot.raster,
                                 colorkey = TRUE,
                                col.regions = viridis,
                                scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

dirLayer <- init(Uarctos_rasters_c, v='y')
z <- zonal(Uarctos_rasters_c, dirLayer, FUN='mean', digits=2)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Latitude_Uarctos<- contourplot(z ~ x*y, data=dat,
                               xlab='Chronology (kya)', ylab='Latitude (º)',
                               panel=panel.levelplot.raster,
                                colorkey = TRUE,
                               col.regions = viridis,
                               scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

# Difference between Phase I (50-45 kr BP) and Phase III (40-35 kyr BP)
path_TIF <- " " # write the directory where the TIF files are located (e.g. path_TIF <- "F:/Data$Code/Species/Homo.sapiens/Layers/TIF/")

Uarctos_50_45 <- c(
  paste0(path_TIF, "50f.tif"),
  paste0(path_TIF, "49f.tif"),
  paste0(path_TIF, "48f.tif"),
  paste0(path_TIF, "47f.tif"),
  paste0(path_TIF, "46f.tif"),
  paste0(path_TIF, "45f.tif")
)
Uarctos_50_45 <- mean(stack(Uarctos_50_45))

Uarctos_40_35 <- c(
  paste0(path_TIF, "40f.tif"),
  paste0(path_TIF, "39f.tif"),
  paste0(path_TIF, "38f.tif"),
  paste0(path_TIF, "37f.tif"),
  paste0(path_TIF, "36f.tif"),
  paste0(path_TIF, "35f.tif")
)
Uarctos_40_35 <- mean(stack(Uarctos_40_35))
Dif_Uarctos<-  Uarctos_40_35 - Uarctos_50_45

# Load raster files
raster_stack_50 <- stack(Uarctos_50_45)
raster_stack_40 <- stack(Uarctos_40_35)
raster_stack_dif <- stack(Dif_Uarctos)

# Convert data
raster_df_n_50 <- as.data.frame(raster_stack_50, xy = TRUE)
rraster_df_n_40 <- as.data.frame(raster_stack_40 , xy = TRUE)
raster_df_n_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

# Plot
Uarctos_b_MIS3<-ggplot(raster_df_n_50, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")

Uarctos_b_MIS3
Uarctos_e_MIS3<-ggplot(rraster_df_n_40, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")
Uarctos_dif<-ggplot(raster_df_n_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), limits = c(-0.4, 0.5), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", linewidth = 0.1)

plot(Uarctos_dif)

Uarctos_dif_m <- Uarctos_dif + theme(plot.margin = unit(c(1,1,1,0.1), "cm"))

ggarrange( Latitude_Uarctos,Longitude_Uarctos, Uarctos_dif_m,ncol = 3)


#### PANTHERA LEO ####

species <- "Panthera.leo"
my_path<- ""
my_path<- paste0(my_path, "/", species, "/Layers/TIF")
layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)
Pleo_layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)

# Create a raster stack from the list of raster files
Pleo_rasters <- stack(Pleo_layers)

Pleo_rasters <- setZ(Pleo_rasters, z_values)

Pleo_rasters_c <- crop(Pleo_rasters, ext)

# Use the hovmoller function with the raster stack

dirLayer <- init(Pleo_rasters_c, v='x')
z <- zonal(Pleo_rasters_c, dirLayer, FUN='mean', digits=2)
head(dat)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Longitude_Pleo<- contourplot(z ~ x*y, data=dat,
                             xlab='Chronology (kya)', ylab='Longitude (º)',
                             panel=panel.levelplot.raster,
                              colorkey = TRUE,
                             col.regions = viridis,
                             scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

dirLayer <- init(Pleo_rasters_c, v='y')
z <- zonal(Pleo_rasters_c, dirLayer, FUN='mean', digits=2)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Latitude_Pleo<- contourplot(z ~ x*y, data=dat,
                            xlab='Chronology (kya)', ylab='Latitude (º)',
                            panel=panel.levelplot.raster,
                             colorkey = TRUE,
                            col.regions = viridis,
                            scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

# Difference between Phase I (50-45 kr BP) and Phase III (40-35 kyr BP)
path_TIF <- "" # write the directory where the TIF files are located (e.g. path_TIF <- "F:/Data$Code/Species/Homo.sapiens/Layers/TIF/")

Pleo_50_45 <- c(
  paste0(path_TIF, "50f.tif"),
  paste0(path_TIF, "49f.tif"),
  paste0(path_TIF, "48f.tif"),
  paste0(path_TIF, "47f.tif"),
  paste0(path_TIF, "46f.tif"),
  paste0(path_TIF, "45f.tif")
)
Pleo_50_45 <- mean(stack(Pleo_50_45))

Pleo_40_35 <- c(
  paste0(path_TIF, "40f.tif"),
  paste0(path_TIF, "39f.tif"),
  paste0(path_TIF, "38f.tif"),
  paste0(path_TIF, "37f.tif"),
  paste0(path_TIF, "36f.tif"),
  paste0(path_TIF, "35f.tif")
)
Pleo_40_35 <- mean(stack(Pleo_40_35))
Dif_Pleo<-  Pleo_40_35 - Pleo_50_45
plot(Dif_Pleo)
# Load raster files
raster_stack_50 <- stack(Pleo_50_45)
raster_stack_40 <- stack(Pleo_40_35)
raster_stack_dif <- stack(Dif_Pleo)

# Convert data
raster_df_n_50 <- as.data.frame(raster_stack_50, xy = TRUE)
rraster_df_n_40 <- as.data.frame(raster_stack_40 , xy = TRUE)
raster_df_n_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

# Plot
Pleo_b_MIS3<-ggplot(raster_df_n_50, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")

Pleo_b_MIS3
Pleo_e_MIS3<-ggplot(rraster_df_n_40, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")
Pleo_dif<-ggplot(raster_df_n_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)),limits = c(-0.4, 0.5), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", linewidth = 0.1)

plot(Pleo_dif)

Pleo_dif_m <- Pleo_dif + theme(plot.margin = unit(c(1,1,1,0.1), "cm"))

ggarrange( Latitude_Pleo,Longitude_Pleo, Pleo_dif_m,ncol = 3)


#### PANTHERA SPELAEA ####

species <- "Panthera.spelaea"
my_path<- ""
my_path<- paste0(my_path, "/", species, "/Layers/TIF")
layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)
Pspelaea_layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)

# Create a raster stack from the list of raster files
Pspelaea_rasters <- stack(Pspelaea_layers)

Pspelaea_rasters <- setZ(Pspelaea_rasters, z_values)

Pspelaea_rasters_c <- crop(Pspelaea_rasters, ext)

# Use the hovmoller function with the raster stack

dirLayer <- init(Pspelaea_rasters_c, v='x')
z <- zonal(Pspelaea_rasters_c, dirLayer, FUN='mean', digits=2)
head(dat)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Longitude_Pspelaea<- contourplot(z ~ x*y, data=dat,
                                 xlab='Chronology (kya)', ylab='Longitude (º)',
                                 panel=panel.levelplot.raster,
                                  colorkey = TRUE,
                                 col.regions = viridis,
                                 scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

dirLayer <- init(Pspelaea_rasters_c, v='y')
z <- zonal(Pspelaea_rasters_c, dirLayer, FUN='mean', digits=2)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Latitude_Pspelaea<- contourplot(z ~ x*y, data=dat,
                                xlab='Chronology (kya)', ylab='Latitude (º)',
                                panel=panel.levelplot.raster,
                                 colorkey = TRUE,
                                col.regions = viridis,
                                scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

# Difference between Phase I (50-45 kr BP) and Phase III (40-35 kyr BP)
path_TIF <- "" # write the directory where the TIF files are located (e.g. path_TIF <- "F:/Data$Code/Species/Homo.sapiens/Layers/TIF/")

Pspelaea_50_45 <- c(
  paste0(path_TIF, "50f.tif"),
  paste0(path_TIF, "49f.tif"),
  paste0(path_TIF, "48f.tif"),
  paste0(path_TIF, "47f.tif"),
  paste0(path_TIF, "46f.tif"),
  paste0(path_TIF, "45f.tif")
)
Pspelaea_50_45 <- mean(stack(Pspelaea_50_45))

Pspelaea_40_35 <- c(
  paste0(path_TIF, "40f.tif"),
  paste0(path_TIF, "39f.tif"),
  paste0(path_TIF, "38f.tif"),
  paste0(path_TIF, "37f.tif"),
  paste0(path_TIF, "36f.tif"),
  paste0(path_TIF, "35f.tif")
)
Pspelaea_40_35 <- mean(stack(Pspelaea_40_35))
Dif_Pspelaea<-  Pspelaea_40_35 - Pspelaea_50_45

# Load raster files
raster_stack_50 <- stack(Pspelaea_50_45)
raster_stack_40 <- stack(Pspelaea_40_35)
raster_stack_dif <- stack(Dif_Pspelaea)

# Convert data
raster_df_n_50 <- as.data.frame(raster_stack_50, xy = TRUE)
rraster_df_n_40 <- as.data.frame(raster_stack_40 , xy = TRUE)
raster_df_n_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

# Plot
Pspelaea_b_MIS3<-ggplot(raster_df_n_50, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")

Pspelaea_b_MIS3
Pspelaea_e_MIS3<-ggplot(rraster_df_n_40, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")

Pspelaea_dif<-ggplot(raster_df_n_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)),limits = c(-0.4, 0.5), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", linewidth = 0.1)

plot(Pspelaea_dif)

Pspelaea_dif_m <- Pspelaea_dif + theme(plot.margin = unit(c(1,1,1,0.1), "cm"))

ggarrange( Latitude_Pspelaea,Longitude_Pspelaea, Pspelaea_dif_m,ncol = 3)


#### PANTHERA PARDUS ####

species <- "Panthera.pardus"
my_path<- ""
my_path<- paste0(my_path, "/", species, "/Layers/TIF")
layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)
Ppardus_layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)

# Create a raster stack from the list of raster files
Ppardus_rasters <- stack(Ppardus_layers)

Ppardus_rasters <- setZ(Ppardus_rasters, z_values)

Ppardus_rasters_c <- crop(Ppardus_rasters, ext)

# Use the hovmoller function with the raster stack

dirLayer <- init(Ppardus_rasters_c, v='x')
z <- zonal(Ppardus_rasters_c, dirLayer, FUN='mean', digits=2)
head(dat)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)
limits <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
Longitude_Ppardus<- contourplot(z ~ x*y, data=dat,
                                xlab='Chronology (kya)', ylab='Longitude (º)',
                                panel=panel.levelplot.raster,
                                 colorkey = TRUE,
                                col.regions = viridis,
                                scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))),
                                at = limits)


dirLayer <- init(Ppardus_rasters_c, v='y')
z <- zonal(Ppardus_rasters_c, dirLayer, FUN='mean', digits=2)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Latitude_Ppardus<- contourplot(z ~ x*y, data=dat,
                               xlab='Chronology (kya)', ylab='Latitude (º)',
                               panel=panel.levelplot.raster,
                                colorkey = TRUE,
                               col.regions = viridis,
                               scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

# Difference between Phase I (50-45 kr BP) and Phase III (40-35 kyr BP)
path_TIF <- "" # write the directory where the TIF files are located (e.g. path_TIF <- "F:/Data$Code/Species/Homo.sapiens/Layers/TIF/")

Ppardus_50_45 <- c(
  paste0(path_TIF, "50f.tif"),
  paste0(path_TIF, "49f.tif"),
  paste0(path_TIF, "48f.tif"),
  paste0(path_TIF, "47f.tif"),
  paste0(path_TIF, "46f.tif"),
  paste0(path_TIF, "45f.tif")
)
Ppardus_50_45 <- mean(stack(Ppardus_50_45))

Ppardus_40_35 <- c(
  paste0(path_TIF, "40f.tif"),
  paste0(path_TIF, "39f.tif"),
  paste0(path_TIF, "38f.tif"),
  paste0(path_TIF, "37f.tif"),
  paste0(path_TIF, "36f.tif"),
  paste0(path_TIF, "35f.tif")
)
Ppardus_40_35 <- mean(stack(Ppardus_40_35))
Dif_Ppardus<-  Ppardus_40_35 - Ppardus_50_45

# Load raster files
raster_stack_50 <- stack(Ppardus_50_45)
raster_stack_40 <- stack(Ppardus_40_35)
raster_stack_dif <- stack(Dif_Ppardus)

# Convert data
raster_df_n_50 <- as.data.frame(raster_stack_50, xy = TRUE)
rraster_df_n_40 <- as.data.frame(raster_stack_40 , xy = TRUE)
raster_df_n_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

# Plot
Ppardus_b_MIS3<-ggplot(raster_df_n_50, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")

Ppardus_b_MIS3
Ppardus_e_MIS3<-ggplot(rraster_df_n_40, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")
Ppardus_dif<-ggplot(raster_df_n_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", linewidth = 0.1)

plot(Ppardus_dif)

Ppardus_dif_m <- Ppardus_dif + theme(plot.margin = unit(c(1,1,1,0.1), "cm"))

ggarrange( Latitude_Ppardus,Longitude_Ppardus, Ppardus_dif_m,ncol = 3)


#### CROCUTA CROCUTA ####

species <- "Crocuta.crocuta"
my_path<- ""
my_path<- paste0(my_path, "/", species, "/Layers/TIF")
layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)
Ccrocuta_layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)


# Create a raster stack from the list of raster files
Ccrocuta_rasters <- stack(Ccrocuta_layers)

Ccrocuta_rasters <- setZ(Ccrocuta_rasters, z_values)

Ccrocuta_rasters_c <- crop(Ccrocuta_rasters, ext)

# Use the hovmoller function with the raster stack

dirLayer <- init(Ccrocuta_rasters_c, v='x')
z <- zonal(Ccrocuta_rasters_c, dirLayer, FUN='mean', digits=2)
head(dat)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Longitude_Ccrocuta<- contourplot(z ~ x*y, data=dat,
                                 xlab='Chronology (kya)', ylab='Longitude (º)',
                                 panel=panel.levelplot.raster,
                                  colorkey = TRUE,
                                 col.regions = viridis,
                                 scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

dirLayer <- init(Ccrocuta_rasters_c, v='y')
z <- zonal(Ccrocuta_rasters_c, dirLayer, FUN='mean', digits=2)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Latitude_Ccrocuta<- contourplot(z ~ x*y, data=dat,
                                xlab='Chronology (kya)', ylab='Latitude (º)',
                                panel=panel.levelplot.raster,
                                 colorkey = TRUE,
                                col.regions = viridis,
                                scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

# Difference between Phase I (50-45 kr BP) and Phase III (40-35 kyr BP)
path_TIF <- "" # write the directory where the TIF files are located (e.g. path_TIF <- "F:/Data$Code/Species/Homo.sapiens/Layers/TIF/")

Ccrocuta_50_45 <- c(
  paste0(path_TIF, "50f.tif"),
  paste0(path_TIF, "49f.tif"),
  paste0(path_TIF, "48f.tif"),
  paste0(path_TIF, "47f.tif"),
  paste0(path_TIF, "46f.tif"),
  paste0(path_TIF, "45f.tif")
)
Ccrocuta_50_45 <- mean(stack(Ccrocuta_50_45))

Ccrocuta_40_35 <- c(
  paste0(path_TIF, "40f.tif"),
  paste0(path_TIF, "39f.tif"),
  paste0(path_TIF, "38f.tif"),
  paste0(path_TIF, "37f.tif"),
  paste0(path_TIF, "36f.tif"),
  paste0(path_TIF, "35f.tif")
)
Ccrocuta_40_35 <- mean(stack(Ccrocuta_40_35))
Dif_Ccrocuta<-  Ccrocuta_40_35 - Ccrocuta_50_45

# Load raster files
raster_stack_50 <- stack(Ccrocuta_50_45)
raster_stack_40 <- stack(Ccrocuta_40_35)
raster_stack_dif <- stack(Dif_Ccrocuta)

# Convert data
raster_df_n_50 <- as.data.frame(raster_stack_50, xy = TRUE)
rraster_df_n_40 <- as.data.frame(raster_stack_40 , xy = TRUE)
raster_df_n_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

# Plot
Ccrocuta_b_MIS3<-ggplot(raster_df_n_50, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")

Ccrocuta_b_MIS3
Ccrocuta_e_MIS3<-ggplot(rraster_df_n_40, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")
Ccrocuta_dif<-ggplot(raster_df_n_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)),limits = c(-0.4, 0.5), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", linewidth = 0.1)

plot(Ccrocuta_dif)

Ccrocuta_dif_m <- Ccrocuta_dif + theme(plot.margin = unit(c(1,1,1,0.1), "cm"))

ggarrange( Latitude_Ccrocuta,Longitude_Ccrocuta, Ccrocuta_dif_m,ncol = 3)


#### CANIS LUPUS ####

species <- "Canis.lupus"
my_path<- ""
my_path<- paste0(my_path, "/", species, "/Layers/TIF")
layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)
Calupus_layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)


# Create a raster stack from the list of raster files
Calupus_rasters <- stack(Calupus_layers)


Calupus_rasters <- setZ(Calupus_rasters, z_values)

Calupus_rasters_c <- crop(Calupus_rasters, ext)

# Use the hovmoller function with the raster stack

dirLayer <- init(Calupus_rasters_c, v='x')
z <- zonal(Calupus_rasters_c, dirLayer, FUN='mean', digits=2)
head(dat)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Longitude_Calupus<- contourplot(z ~ x*y, data=dat,
                                xlab='Chronology (kya)', ylab='Longitude (º)',
                                panel=panel.levelplot.raster,
                                 colorkey = TRUE,
                                col.regions = viridis,
                                scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

dirLayer <- init(Calupus_rasters_c, v='y')
z <- zonal(Calupus_rasters_c, dirLayer, FUN='mean', digits=2)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Latitude_Calupus<- contourplot(z ~ x*y, data=dat,
                               xlab='Chronology (kya)', ylab='Latitude (º)',
                               panel=panel.levelplot.raster,
                                colorkey = TRUE,
                               col.regions = viridis,
                               scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

# Difference between Phase I (50-45 kr BP) and Phase III (40-35 kyr BP)
path_TIF <- "" # write the directory where the TIF files are located (e.g. path_TIF <- "F:/Data$Code/Species/Homo.sapiens/Layers/TIF/")

Calupus_50_45 <- c(
  paste0(path_TIF, "50f.tif"),
  paste0(path_TIF, "49f.tif"),
  paste0(path_TIF, "48f.tif"),
  paste0(path_TIF, "47f.tif"),
  paste0(path_TIF, "46f.tif"),
  paste0(path_TIF, "45f.tif")
)
Calupus_50_45 <- mean(stack(Calupus_50_45))

Calupus_40_35 <- c(
  paste0(path_TIF, "40f.tif"),
  paste0(path_TIF, "39f.tif"),
  paste0(path_TIF, "38f.tif"),
  paste0(path_TIF, "37f.tif"),
  paste0(path_TIF, "36f.tif"),
  paste0(path_TIF, "35f.tif")
)
Calupus_40_35 <- mean(stack(Calupus_40_35))
Dif_Calupus<-  Calupus_40_35 - Calupus_50_45

# Load raster files
raster_stack_50 <- stack(Calupus_50_45)
raster_stack_40 <- stack(Calupus_40_35)
raster_stack_dif <- stack(Dif_Calupus)

# Convert data
raster_df_n_50 <- as.data.frame(raster_stack_50, xy = TRUE)
rraster_df_n_40 <- as.data.frame(raster_stack_40 , xy = TRUE)
raster_df_n_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

# Plot
Calupus_b_MIS3<-ggplot(raster_df_n_50, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")

Calupus_b_MIS3
Calupus_e_MIS3<-ggplot(rraster_df_n_40, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")
Calupus_dif<-ggplot(raster_df_n_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)),limits = c(-0.4, 0.5), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", linewidth = 0.1)

plot(Calupus_dif)

Calupus_dif_m <- Calupus_dif + theme(plot.margin = unit(c(1,1,1,0.1), "cm"))

ggarrange( Latitude_Calupus,Longitude_Calupus, Calupus_dif_m,ncol = 3)


#### VULPES VULPES ####

species <- "Vulpes.vulpes"
my_path<- ""
my_path<- paste0(my_path, "/", species, "/Layers/TIF")
layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)
Vvulpes_layers <-list.files(path=my_path, pattern='f.tif$', full.names=TRUE)


# Create a raster stack from the list of raster files
Vvulpes_rasters <- stack(Vvulpes_layers)

Vvulpes_rasters <- setZ(Vvulpes_rasters, z_values)

Vvulpes_rasters_c <- crop(Vvulpes_rasters, ext)

# Use the hovmoller function with the raster stack

dirLayer <- init(Vvulpes_rasters_c, v='x')
z <- zonal(Vvulpes_rasters_c, dirLayer, FUN='mean', digits=2)
head(dat)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Longitude_Vvulpes<- contourplot(z ~ x*y, data=dat,
                                xlab='Chronology (kya)', ylab='Longitude (º)',
                                panel=panel.levelplot.raster,
                                 colorkey = TRUE,
                                col.regions = viridis,
                                scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

dirLayer <- init(Vvulpes_rasters_c, v='y')
z <- zonal(Vvulpes_rasters_c, dirLayer, FUN='mean', digits=2)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Latitude_Vvulpes<- contourplot(z ~ x*y, data=dat,
                               xlab='Chronology (kya)', ylab='Latitude (º)',
                               panel=panel.levelplot.raster,
                                colorkey = TRUE,
                               col.regions = viridis,
                               scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

# Difference between Phase I (50-45 kr BP) and Phase III (40-35 kyr BP)
path_TIF <- "" # write the directory where the TIF files are located (e.g. path_TIF <- "F:/Data$Code/Species/Homo.sapiens/Layers/TIF/")

Vvulpes_50_45 <- c(
  paste0(path_TIF, "50f.tif"),
  paste0(path_TIF, "49f.tif"),
  paste0(path_TIF, "48f.tif"),
  paste0(path_TIF, "47f.tif"),
  paste0(path_TIF, "46f.tif"),
  paste0(path_TIF, "45f.tif")
)
Vvulpes_50_45 <- mean(stack(Vvulpes_50_45))

Vvulpes_40_35<- c(
  paste0(path_TIF, "40f.tif"),
  paste0(path_TIF, "39f.tif"),
  paste0(path_TIF, "38f.tif"),
  paste0(path_TIF, "37f.tif"),
  paste0(path_TIF, "36f.tif"),
  paste0(path_TIF, "35f.tif")
)
Vvulpes_40_35<- mean(stack(Vvulpes_40_35))
Dif_Vvulpes<-  Vvulpes_40_35- Vvulpes_50_45

# Load raster files
raster_stack_50 <- stack(Vvulpes_50_45)
raster_stack_40 <- stack(Vvulpes_40_35)
raster_stack_dif <- stack(Dif_Vvulpes)

# Convert data
raster_df_n_50 <- as.data.frame(raster_stack_50, xy = TRUE)
rraster_df_n_40 <- as.data.frame(raster_stack_40 , xy = TRUE)
raster_df_n_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

# Plot
Vvulpes_b_MIS3<-ggplot(raster_df_n_50, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")

Vvulpes_b_MIS3
Vvulpes_e_MIS3<-ggplot(rraster_df_n_40, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")
Vvulpes_dif<-ggplot(raster_df_n_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), limits = c(-0.4, 0.5),na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", linewidth = 0.1)

plot(Vvulpes_dif)

Vvulpes_dif_m <- Vvulpes_dif + theme(plot.margin = unit(c(1,1,1,0.1), "cm"))

ggarrange( Latitude_Vvulpes,Longitude_Vvulpes, Vvulpes_dif_m,ncol = 3)


#### VULPES LAGOPUS ####

species <- "Vulpes.lagopus"
my_path<- " "
my_path<- paste0(my_path, "/", species, "/Layers/TIF")
layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)
Vlagopus_layers <-list.files(path=my_path, pattern='f.tif$', full.names=TRUE)

# Create a raster stack from the list of raster files
Vlagopus_rasters <- stack(Vlagopus_layers)

Vlagopus_rasters <- setZ(Vlagopus_rasters, z_values)

Vlagopus_rasters_c <- crop(Vlagopus_rasters, ext)

# Use the hovmoller function with the raster stack

dirLayer <- init(Vlagopus_rasters_c, v='x')
z <- zonal(Vlagopus_rasters_c, dirLayer, FUN='mean', digits=2)
head(dat)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Longitude_Vlagopus<- contourplot(z ~ x*y, data=dat,
                                 xlab='Chronology (kya)', ylab='Longitude (º)',
                                 panel=panel.levelplot.raster,
                                  colorkey = TRUE,
                                 col.regions = viridis,
                                 scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

dirLayer <- init(Vlagopus_rasters_c, v='y')
z <- zonal(Vlagopus_rasters_c, dirLayer, FUN='mean', digits=2)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Latitude_Vlagopus<- contourplot(z ~ x*y, data=dat,
                                xlab='Chronology (kya)', ylab='Latitude (º)',
                                panel=panel.levelplot.raster,
                                 colorkey = TRUE,
                                col.regions = viridis,
                                scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

# Difference between Phase I (50-45 kr BP) and Phase III (40-35 kyr BP)
path_TIF <- "" # write the directory where the TIF files are located (e.g. path_TIF <- "F:/Data$Code/Species/Homo.sapiens/Layers/TIF/")

Vlagopus_50_45 <- c(
  paste0(path_TIF, "50f.tif"),
  paste0(path_TIF, "49f.tif"),
  paste0(path_TIF, "48f.tif"),
  paste0(path_TIF, "47f.tif"),
  paste0(path_TIF, "46f.tif"),
  paste0(path_TIF, "45f.tif")
)
Vlagopus_50_45 <- mean(stack(Vlagopus_50_45))

Vlagopus_40_35 <- c(
  paste0(path_TIF, "40f.tif"),
  paste0(path_TIF, "39f.tif"),
  paste0(path_TIF, "38f.tif"),
  paste0(path_TIF, "37f.tif"),
  paste0(path_TIF, "36f.tif"),
  paste0(path_TIF, "35f.tif")
)
Vlagopus_40_35 <- mean(stack(Vlagopus_40_35))
Dif_Vlagopus<-  Vlagopus_40_35 - Vlagopus_50_45

# Load raster files
raster_stack_50 <- stack(Vlagopus_50_45)
raster_stack_40 <- stack(Vlagopus_40_35)
raster_stack_dif <- stack(Dif_Vlagopus)

# Convert data
raster_df_n_50 <- as.data.frame(raster_stack_50, xy = TRUE)
rraster_df_n_40 <- as.data.frame(raster_stack_40 , xy = TRUE)
raster_df_n_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

# Plot
Vlagopus_b_MIS3<-ggplot(raster_df_n_50, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")

Vlagopus_b_MIS3
Vlagopus_e_MIS3<-ggplot(rraster_df_n_40, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")
Vlagopus_dif<-ggplot(raster_df_n_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)),limits = c(-0.4, 0.5), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", linewidth = 0.1)

plot(Vlagopus_dif)

Vlagopus_dif_m <- Vlagopus_dif + theme(plot.margin = unit(c(1,1,1,0.1), "cm"))

ggarrange( Latitude_Vlagopus,Longitude_Vlagopus, Vlagopus_dif_m,ncol = 3)

#### FELIS SILVESTRIS ####

species <- "Felis.silvestris"
my_path<- ""
my_path<- paste0(my_path, "/", species, "/Layers/TIF")
layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)
Fsilvestris_layers <-list.files(path=my_path, pattern='f.tif$', full.names=TRUE)


# Create a raster stack from the list of raster files
Fsilvestris_rasters <- stack(Fsilvestris_layers)


Fsilvestris_rasters <- setZ(Fsilvestris_rasters, z_values)

Fsilvestris_rasters_c <- crop(Fsilvestris_rasters, ext)

# Use the hovmoller function with the raster stack

dirLayer <- init(Fsilvestris_rasters_c, v='x')
z <- zonal(Fsilvestris_rasters_c, dirLayer, FUN='mean', digits=2)
head(dat)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Longitude_Fsilvestris<- contourplot(z ~ x*y, data=dat,
                                    xlab='Chronology (kya)', ylab='Longitude (º)',
                                    panel=panel.levelplot.raster,
                                     colorkey = TRUE,
                                    col.regions = viridis,
                                    scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

dirLayer <- init(Fsilvestris_rasters_c, v='y')
z <- zonal(Fsilvestris_rasters_c, dirLayer, FUN='mean', digits=2)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Latitude_Fsilvestris<- contourplot(z ~ x*y, data=dat,
                                   xlab='Chronology (kya)', ylab='Latitude (º)',
                                   panel=panel.levelplot.raster,
                                    colorkey = TRUE,
                                   col.regions = viridis,
                                   scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

# Difference between Phase I (50-45 kr BP) and Phase III (40-35 kyr BP)
path_TIF <- "" # write the directory where the TIF files are located (e.g. path_TIF <- "F:/Data$Code/Species/Homo.sapiens/Layers/TIF/")

Fsilvestris_50_45 <- c(
  paste0(path_TIF, "50f.tif"),
  paste0(path_TIF, "49f.tif"),
  paste0(path_TIF, "48f.tif"),
  paste0(path_TIF, "47f.tif"),
  paste0(path_TIF, "46f.tif"),
  paste0(path_TIF, "45f.tif")
)
Fsilvestris_50_45 <- mean(stack(Fsilvestris_50_45))

Fsilvestris_40_35 <- c(
  paste0(path_TIF, "40f.tif"),
  paste0(path_TIF, "39f.tif"),
  paste0(path_TIF, "38f.tif"),
  paste0(path_TIF, "37f.tif"),
  paste0(path_TIF, "36f.tif"),
  paste0(path_TIF, "35f.tif")
)
Fsilvestris_40_35 <- mean(stack(Fsilvestris_40_35))
Dif_Fsilvestris<-  Fsilvestris_40_35 - Fsilvestris_50_45

# Load raster files
raster_stack_50 <- stack(Fsilvestris_50_45)
raster_stack_40 <- stack(Fsilvestris_40_35)
raster_stack_dif <- stack(Dif_Fsilvestris)

# Convert data
raster_df_n_50 <- as.data.frame(raster_stack_50, xy = TRUE)
rraster_df_n_40 <- as.data.frame(raster_stack_40 , xy = TRUE)
raster_df_n_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

# Plot
Fsilvestris_b_MIS3<-ggplot(raster_df_n_50, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")

Fsilvestris_b_MIS3
Fsilvestris_e_MIS3<-ggplot(rraster_df_n_40, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")
Fsilvestris_dif<-ggplot(raster_df_n_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)),limits = c(-0.6, 0.6), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", linewidth = 0.1)

plot(Fsilvestris_dif)

Fsilvestris_dif_m <- Fsilvestris_dif + theme(plot.margin = unit(c(1,1,1,0.1), "cm"))

ggarrange( Latitude_Fsilvestris,Longitude_Fsilvestris, Fsilvestris_dif_m,ncol = 3)

#### MELES MELES ####

species <- "Meles.meles"
my_path<- ""
my_path<- paste0(my_path, "/", species, "/Layers/TIF")
layers <- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)
Mmeles_layers<- list.files(path=my_path, pattern='f.tif$', full.names=TRUE)


# Create a raster stack from the list of raster files
Mmeles_rasters <- stack(Mmeles_layers)

Mmeles_rasters <- setZ(Mmeles_rasters, z_values)

Mmeles_rasters_c <- crop(Mmeles_rasters, ext)

# Use the hovmoller function with the raster stack

dirLayer <- init(Mmeles_rasters_c, v='x')
z <- zonal(Mmeles_rasters_c, dirLayer, FUN='mean', digits=2)
head(dat)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Longitude_Mmeles<- contourplot(z ~ x*y, data=dat,
                               xlab='Chronology (kya)', ylab='Longitude (º)',
                               panel=panel.levelplot.raster,
                                colorkey = TRUE,
                               col.regions = viridis,
                               scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

dirLayer <- init(Mmeles_rasters_c, v='y')
z <- zonal(Mmeles_rasters_c, dirLayer, FUN='mean', digits=2)
dat <- expand.grid(y=z[,1], x=z_values)
dat$z <- as.vector(z[,-1], mode='numeric')
dat$x<- rev(dat$x)

Latitude_Mmeles<- contourplot(z ~ x*y, data=dat,
                              xlab='Chronology (kya)', ylab='Latitude (º)',
                              panel=panel.levelplot.raster,
                               colorkey = TRUE,
                              col.regions = viridis,
                              scales = list(x = list(at = c(1, 6, 11, 16, 21, 26))))

# Difference between Phase I (50-45 kr BP) and Phase III (40-35 kyr BP)
path_TIF <- "" # write the directory where the TIF files are located (e.g. path_TIF <- "F:/Data$Code/Species/Homo.sapiens/Layers/TIF/")

Mmeles_50_45 <- c(
  paste0(path_TIF, "50f.tif"),
  paste0(path_TIF, "49f.tif"),
  paste0(path_TIF, "48f.tif"),
  paste0(path_TIF, "47f.tif"),
  paste0(path_TIF, "46f.tif"),
  paste0(path_TIF, "45f.tif")
)
Mmeles_50_45 <- mean(stack(Mmeles_50_45))

Mmeles_40_35 <- c(
  paste0(path_TIF, "40f.tif"),
  paste0(path_TIF, "39f.tif"),
  paste0(path_TIF, "38f.tif"),
  paste0(path_TIF, "37f.tif"),
  paste0(path_TIF, "36f.tif"),
  paste0(path_TIF, "35f.tif")
)
Mmeles_40_35 <- mean(stack(Mmeles_40_35))
Dif_Mmeles<-  Mmeles_40_35 - Mmeles_50_45

# Load raster files
raster_stack_50 <- stack(Mmeles_50_45)
raster_stack_40 <- stack(Mmeles_40_35)
raster_stack_dif <- stack(Dif_Mmeles)

# Convert data
raster_df_n_50 <- as.data.frame(raster_stack_50, xy = TRUE)
rraster_df_n_40 <- as.data.frame(raster_stack_40 , xy = TRUE)
raster_df_n_dif <- as.data.frame(raster_stack_dif, xy = TRUE)

# Plot
Mmeles_b_MIS3<-ggplot(raster_df_n_50, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")

Mmeles_b_MIS3
Mmeles_e_MIS3<-ggplot(rraster_df_n_40, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)), na.value = "white") +
  theme_classic() +
  ylab("")+
  xlab("")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black")
Mmeles_dif<-ggplot(raster_df_n_dif, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colours = rev(mi_paleta(100)),limits = c(-0.4, 0.5), na.value = "white") +
  theme_classic() +
  ylab("Latitude (º)")+
  xlab("Longitude (º)")+
  coord_cartesian(xlim = c(-10, 45), ylim = c(36, 60)) +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", linewidth = 0.1)

plot(Mmeles_dif)

Mmeles_dif_m <- Mmeles_dif + theme(plot.margin = unit(c(1,1,1,0.1), "cm"))

ggarrange( Latitude_Mmeles,Longitude_Mmeles, Mmeles_dif_m,ncol = 3)





##### PLOTS ####


plot_arrange <- grid.arrange(sapiens_box, Latitude_sapiens, Longitude_sapiens, sapiens_dif_m,
                             neand_box, Latitude_neand, Longitude_neand, neand_dif_m,
                             pleo_box, Latitude_Pleo, Longitude_Pleo, Pleo_dif_m,
                             pspelaea_box, Latitude_Pspelaea, Longitude_Pspelaea, Pspelaea_dif_m,
                             ppardus_box, Latitude_Ppardus, Longitude_Ppardus, Ppardus_dif_m,
                             uspelaeus_box, Latitude_Uspelaeus, Longitude_Uspelaeus, Uspelaeus_dif_m,
                             uarctos_box, Latitude_Uarctos, Longitude_Uarctos, Uarctos_dif_m,
                             crocuta_box, Latitude_Ccrocuta, Longitude_Ccrocuta, Ccrocuta_dif_m,
                             calupus_box, Latitude_Calupus, Longitude_Calupus, Calupus_dif_m,
                             vvulpes_box, Latitude_Vvulpes, Longitude_Vvulpes, Vvulpes_dif_m,
                             vlagopus_box, Latitude_Vlagopus, Longitude_Vlagopus, Vlagopus_dif_m,
                             fsilvestris_box, Latitude_Fsilvestris, Longitude_Fsilvestris, Fsilvestris_dif_m,
                             mmeles_box, Latitude_Mmeles, Longitude_Mmeles, Mmeles_dif_m,
                             ncol = 4, nrow = 13)




hominin <- grid.arrange(sapiens_box, Latitude_sapiens, Longitude_sapiens, sapiens_dif_m,
                        neand_box, Latitude_neand, Longitude_neand, neand_dif_m, ncol=4, nrow=2)
library(egg)

plot_grid(sapiens_box, Latitude_sapiens, Longitude_sapiens, sapiens_dif_m, 
                   ncol = 4, heights = c(2, 0.2, 0.2, 0.2))+ theme(plot.margin = unit(c(0, 0,-20 , 0), "lines"))

neand + theme(plot.margin = unit(c(0, 0, -0.5, 0), "lines"))

  sapiens <- grid.arrange(sapiens_box, Latitude_sapiens, Longitude_sapiens, sapiens_dif_m,ncol=4, nrow=1)
  neand <- grid.arrange(neand_box, Latitude_neand, Longitude_neand, neand_dif_m,ncol=4, nrow=1)
  pleo <- grid.arrange(pleo_box, Latitude_Pleo, Longitude_Pleo, Pleo_dif_m, ncol=4, nrow=1)
  pspelaea <- grid.arrange(pspelaea_box, Latitude_Pspelaea, Longitude_Pspelaea, Pspelaea_dif_m,ncol=4, nrow=1)
  ppardus <- grid.arrange(ppardus_box, Latitude_Ppardus, Longitude_Ppardus, Ppardus_dif_m,ncol=4, nrow=1)
  fsilvestris <- grid.arrange(fsilvestris_box, Latitude_Fsilvestris, Longitude_Fsilvestris, Fsilvestris_dif_m,ncol=4, nrow=1)
  uarctos <- grid.arrange(uarctos_box, Latitude_Uarctos, Longitude_Uarctos, Uarctos_dif_m,ncol=4, nrow=1)
  uspelaeus <- grid.arrange(uspelaeus_box, Latitude_Uspelaeus, Longitude_Uspelaeus, Uspelaeus_dif_m,ncol=4, nrow=1)
  crocuta <- grid.arrange(crocuta_box, Latitude_Ccrocuta, Longitude_Ccrocuta, Ccrocuta_dif_m,ncol=4, nrow=1)
  calupus <- grid.arrange(calupus_box, Latitude_Calupus, Longitude_Calupus, Calupus_dif_m, ncol=4, nrow=1)
  vvulpes<- grid.arrange(vvulpes_box, Latitude_Vvulpes, Longitude_Vvulpes, Vvulpes_dif_m,ncol=4, nrow=1)
  vlagopus <- grid.arrange(vlagopus_box, Latitude_Vlagopus, Longitude_Vlagopus, Vlagopus_dif_m, ncol=4, nrow=1)
  mmeles<- grid.arrange(mmeles_box, Latitude_Mmeles, Longitude_Mmeles, Mmeles_dif_m, ncol=4, nrow=1)



