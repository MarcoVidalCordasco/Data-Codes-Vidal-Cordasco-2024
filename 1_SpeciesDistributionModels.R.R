################## The script below was written in  Rstudio, so it is  recommended
############ to be run within RStudio.
rm(list = ls()) # Clear all
# SETUP ####
# the following command sets the working directory to the folder where this
# script is located (similar to RStudio menu "Session - Set Working Directory -
# To Source File Location"):
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Check the directory :
getwd()

# LOAD PACKAGES ####
library(openxlsx) 
library(terra)  
library(geodata) 
library(sdmpredictors)  
library(fuzzySim)  
library(pastclim) 
library(gridExtra)
library(raster)
library(lattice)
library(gam)  
library(maxnet)  
library(randomForest)  
library(gbm)  
library(embarcadero)  
library(beanplot)
library(ENMTools) 
library(tidyverse)
library(magick)
library(plyr)
library(rworldmap)
library(modEvA)  
library(blockCV)  
library(landscapemetrics)
library(ncdf4)
library(fields)
library(maps)
library(sf)
library(car)
library(pROC)
library(geosphere)
library(caret)
library(fitMaxnet)
library(ggplot2)
library(tidyr)
library(RColorBrewer)


# Specify where you want to save the outputs
my_path<- " "



#_________________________________________________________________________
################# OPEN DATA AND SELECT SPECIES
#

Presence_data<-read.xlsx("Data.xlsx", rowNames=FALSE, 
          colNames=TRUE, sheet="F_df")
head(Presence_data)

# for each species, just write exactly the same name below:

species <- "Panthera.spelaea" # specify the species
db <- subset(Presence_data,Panthera.spelaea==1) #specify the species

head(db)
nrow(db)

# Create a new data frame

Presence_data<- data.frame(db$`Site/Level`)
Presence_data$x <- db$x
Presence_data$y<- db$y
Presence_data$Age <- db$Age
Presence_data$SD <- db$SD
Presence_data$presence <- 1
Presence_data

# Create a new folder for each species. Here will be saved all outcomes obtained:
pres_folder <- paste0(species)
pres_folder
dir.create(pres_folder)  # '../' goes up one level from the current working directory, so this creates the 'species' folder 


# map the species occurrence data:
# Define the extent for plotting
mywindow <- c(-15, 60, 33, 75)  # xmin, xmax, ymin, ymax

# Extract the subset of the world map data
world_map_subset <- map_data("world") %>%
  filter(long >= mywindow[1] & long <= mywindow[2] & lat >= mywindow[3] & lat <= mywindow[4] &
           region != "Morocco" & region != "Algeria" & region != "Tunisia")

# Plot the subset of the world map
ggplot() +
  geom_polygon(data = world_map_subset, aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "black") +
  geom_point(data = Presence_data, aes(x = x, y = y), color = "red")+
  theme_classic()




#######___________________________________________________________________________

# CHECK ALL PREDICTOR VARIABLES ####

get_vars_for_dataset(dataset="Krapp2021", details = TRUE)

# if this is the first time you run the code, run the following line chunk:
 set_data_path(" ") # introduce your own directory within the " ". There the .nc files will be stored

 # You have to run the next line chunk only the first time:
 
 download_dataset(c("bio01","bio04", "bio05", "bio06","bio07", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                   "bio17", "bio18", "bio19", "npp", "lai", "altitude", "rugosity"),data="Krapp2021")

# Take into account that these NetCDF files are directly dowloaded with a 0.5 spatial resolution. Each file was downscaled with
# bi-linear interpolation to 0.25 resolution with the Get_Env_F.R file. The bilinear interpolated maps are available in the same repository where 
# you found this code.

# Select all variables:
bio_variables = c("bio01","bio04", "bio05", "bio06", "bio07", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                  "bio17", "bio18", "bio19", "npp", "lai", "altitude", "rugosity")



Presence_data$time_bp<- Presence_data$Age * -1 # in pastclim, "Age" should be negative, so multiply it by * -1 
head(Presence_data)

# GET VARIABLES IN ALL PRESENCE DATA
# THIS STEP IS ONLY NECESSARY THE FIRST TIME YOU RUN THE CODE

# First, change name of x and y columns because some packages require "longitude" and "latitude" to be specifically stated:

colnames(Presence_data)[2] = "longitude"
colnames(Presence_data)[3] = "latitude"

# get the values of each variable for each location with the presence of the selected species:

Variables<- location_slice( x = Presence_data, bio_variables = c("bio01","bio04", "bio05", "bio06", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                                 "bio17", "bio18", "bio19",  "npp",  "lai", "altitude", "rugosity"), dataset = "Krapp2021", nn_interpol = TRUE, buffer=TRUE)

head(Presence_data)
head(Variables)


# Complete the information of "Variables" with that of "Presence_data":
Variables <- Variables %>%
  mutate(
    db..Site.Level.= coalesce(db..Site.Level., Presence_data$db..Site.Level),
    longitude = coalesce(longitude, Presence_data$longitude),
    latitude = coalesce(latitude, Presence_data$latitude),
    Age = coalesce(Age, Presence_data$Age),
    SD = coalesce(SD, Presence_data$SD),
    presence = coalesce(presence, Presence_data$presence),
    time_bp = coalesce(time_bp, Presence_data$time_bp)
  )
Variables$time_bp_slice <- round(Variables$time_bp, -3)

# Check it:
head(Variables)


# save outputs in your species' folder

write.xlsx(Variables, paste0( species, "/Variables_Krapp2021222_", species, ".xlsx"), colNames = TRUE, rowNames = FALSE, detectDates = FALSE)
getwd()

# I you already have run the code and you have this file ("Variables"), just load it:

Presence_data <-read.xlsx(paste0( species, "/Variables_Krapp2021_", species, ".xlsx"), rowNames=FALSE, 
                          colNames=TRUE, sheet="Sheet 1")
head(Presence_data)
nrow(Presence_data)

Presence_data<- Presence_data[complete.cases(Presence_data), ]


#_________________________________________________________________
# REPEAT THE SAME PROCEDURE, BUT FILTERING PRESENCE POINTS:
#
#___________________________________________________________
# Sample bias in species occurrence data has long been a recognized issue in SDM. 
# However, environmental filtering of observation data can improve model predictions 
# by reducing redundancy in environmental (e.g. climatic) hyper-space (Varela et al. 2014). 
# The following command removes duplicates and thin the points: only one occurrence per pixel and time-step


 # _________
 #
 # BUFFER
 #
 #__________
 # Before obtaining climate variables from points, we set an area enclosing all the fossil 
 # localities where the species was present, and then we create a buffer
 # around the polygon with a radius equal to 10% of the maximum distance between actual 
 # fossil occurrences. This is the calibration area:

countries <- readRDS("gadm36_adm0_r5_pk.rds")

 pres_points <- vect(Presence_data, geom = c("longitude", "latitude"), keepgeom = TRUE, crs = "epsg:4326")  
 buff_dist <- max(distance(pres_points) * 0.1 , na.rm = TRUE)  # can take time if there are many points!
 buff_dist
 pres_buff <- aggregate(buffer(pres_points, width = buff_dist))
 plot(pres_buff)
 plot(pres_points, col = "blue", add = TRUE)
 plot(countries, border = "grey", add = TRUE)
 
 
# Get the spatial resolution of one climatic/environmental layer for thinning:

 climate_layer <- region_slice(
   time_bp = -40000,
   bio_variables = bio_variables[[1]],
   dataset = "Krapp2021"
 )
 
# plot(climate_layer)

 plot(climate_layer[[1]])
 layers_cut <- crop(climate_layer[[1]], pres_buff, mask = TRUE)
 plot(layers_cut)
 # Create empty dataset where outcomes will be saved
  Outcomes_v_thinned <- data.frame()

  # The following command thins presence point to avoid over-representation. There will be only 1 presence/pixel/age bin (1k yr)
  for (step in unique(Presence_data$time_bp_slice)){
  dat <- subset(Presence_data, time_bp_slice== step) # select only points with same chronology
  if (nrow(dat)> 0) 
  {
  dat <- gridRecords(rst = layers_cut[[1]], pres.coords = dat[ , c("longitude", "latitude")], plot = TRUE) 
  dat<- subset(dat, presence==1)
  colnames(dat)[2] = "longitude"
  colnames(dat)[3] = "latitude"
  dat$time_bp <- step
  dat<- dat[ , !(names(dat) %in% "bio01")] # the raster with the variable bios01 was used for thinning, now we remove this column from dat
  Variables_thinned<- location_slice( x = dat, bio_variables = c("bio01","bio04", "bio05", "bio06", "bio07", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                                                       "bio17", "bio18", "bio19", "npp", "lai", "altitude", "rugosity"), dataset = "Krapp2021", nn_interpol = TRUE, buffer=TRUE)
 head(Variables_thinned)
 head(dat)
   Variables_thinned <- Variables_thinned %>%
    mutate(
      presence= coalesce(presence, dat$presence),
      longitude = coalesce(longitude, dat$longitude),
      latitude = coalesce(latitude, dat$latitude),
      cells = coalesce(cells, dat$cells),
      time_bp = coalesce(time_bp, dat$time_bp),
      time_bp_slice = coalesce(time_bp_slice, dat$time_bp_slice)
   
    )
  
  
  Outcomes_v_thinned <- rbind(Outcomes_v_thinned, Variables_thinned)
  }
}


#PLOT RAW PRESENCE POINTS AND THINNED PRESENCES
plot(countries, ext = mywindow)
plot(countries, xlim = mywindow[1:2], ylim = mywindow[3:4])
points(Presence_data[ , c("longitude", "latitude")], col = "blue")
points(Outcomes_v_thinned[ , c("longitude", "latitude")], col = "red")
nrow(Presence_data)
nrow(Outcomes_v_thinned)


# save outputs in your species' folder
write.xlsx(Outcomes_v_thinned, paste0( species, "/Outcomes_v_thinned_Krapp2021_", species, ".xlsx"))


# I you already have run the code previously and you have this file ("Outcomes_v_thinned"), just load it:

Presence_data_thinned <-read.xlsx(paste0( species, "/Outcomes_v_thinned_Krapp2021_", species, ".xlsx"), rowNames=FALSE, 
                          colNames=TRUE, sheet="Sheet 1")

head(Presence_data_thinned)


# PSEUDO-ABSENCE POINTS ####
# CYCLE OVER 50 * n observations TIME STEPS AND SAMPLE THE ENVIRONMENTAL CONDITIONS
# THIS STEP IS ONLY NECESSARY THE FIRST TIME YOU RUN THE CODE

# First, we create an empty list:
climate_list<-list()

# Extract variables from pseudoabsences within the calibration area:
for (step in unique(Presence_data_thinned$time_bp)){
  # extract climate and environmental variables
  step_climate <- region_slice(
    time_bp = step,
    bio_variables = bio_variables,
    dataset = "Krapp2021")
  # subset calibration area:
  crs(step_climate) <- crs(pres_buff)
  step_climate <- terra::mask(step_climate, pres_buff)# select calibration area (press_buff) 
  # Select number of observations per time step / chronology:
  n <- subset(Presence_data_thinned, time_bp== step) # select only points with same chronology
  # Select pesudoabsences:
  step_values <- terra::spatSample(step_climate, 50 * nrow(n), na.rm=TRUE, cells=TRUE, xy=TRUE) 
  step_values$time_bp<-step 
  climate_list[[as.character(step)]]<-step_values
}

# combine data into a single matrix and include "presence", "x" (longitude) and "y" (latitude) variables
baseline_climate <- do.call(rbind, climate_list)
baseline_climate$presence<- 0
baseline_climate$longitude<- baseline_climate$x
baseline_climate$latitude<- baseline_climate$y
nrow(baseline_climate)
# save output
write.xlsx(baseline_climate, paste0( species, "/Baseline_allvars", species, ".xlsx"))

# merge datasets
head(baseline_climate)
head(Presence_data_thinned)

db<- rbind.fill(Presence_data_thinned, baseline_climate)

colnames(db)

# Now, lets eliminate the empty or unnecesary columns for the next steps:
db <- db[, !(colnames(db) %in% c("cells", "time_bp_slice", "cell", "x", "y"))]
colnames(db)
# Save:
write.xlsx(db, paste0( species, "/ALL", species, ".xlsx"))

# I you have already run the previous code and you have this file ("ALL+species"), just load it:

db<- read.xlsx(paste0( species, "/ALL", species, ".xlsx"), rowNames = FALSE, colNames = TRUE, sheet="Sheet 1")
head(db)
nrow(db)
colnames(db)

# PLOT PSEUDOABSENCE AND PRESENCE DATA: ###
 plot(countries, ext = mywindow)
 # Plot the 'countries' data within the specified extent
 plot(countries, xlim = mywindow[1:2], ylim = mywindow[3:4])
 points(subset(db, db$presence == 1, select = c("longitude", "latitude")), pch = 20, cex = 0.7, col="red")
 points(subset(db, db$presence == 0, select = c("longitude", "latitude")), pch = 20, cex = 0.3, col="blue")



# SPECIES DISTRIBUTION MODELS ####
 
 # First, load all variables for each 1 k years from  between 55 and 27 k yrs BP. 
 # (These files are obtained from the file "Get_Env_F.R")
 
 path_env_f <- "Bilinear_interpolated_variables" # Check that the "Bilinear_interpolated_variables" folder is in the same script directory.
 
 Climate_55K_BP <- readRDS(paste0( path_env_f, "/Climate_55.rds"))
 Climate_54K_BP <- readRDS(paste0( path_env_f, "/Climate_54.rds"))
 Climate_53K_BP <- readRDS(paste0( path_env_f, "/Climate_53.rds"))
 Climate_52K_BP <- readRDS(paste0( path_env_f, "/Climate_52.rds"))
 Climate_51K_BP <- readRDS(paste0( path_env_f, "/Climate_51.rds"))
 Climate_50K_BP <- readRDS(paste0( path_env_f, "/Climate_50.rds"))
 Climate_49K_BP <- readRDS(paste0( path_env_f, "/Climate_49.rds"))
 Climate_48K_BP <- readRDS(paste0( path_env_f, "/Climate_48.rds"))
 Climate_47K_BP <- readRDS(paste0( path_env_f, "/Climate_47.rds"))
 Climate_46K_BP <- readRDS(paste0( path_env_f, "/Climate_46.rds"))
 Climate_45K_BP <- readRDS(paste0( path_env_f, "/Climate_45.rds"))
 Climate_44K_BP <- readRDS(paste0( path_env_f, "/Climate_44.rds"))
 Climate_43K_BP <- readRDS(paste0( path_env_f, "/Climate_43.rds"))
 Climate_42K_BP <- readRDS(paste0( path_env_f, "/Climate_42.rds"))
 Climate_41K_BP <- readRDS(paste0( path_env_f, "/Climate_41.rds"))
 Climate_40K_BP <- readRDS(paste0( path_env_f, "/Climate_40.rds"))
 Climate_39K_BP <- readRDS(paste0( path_env_f, "/Climate_39.rds"))
 Climate_38K_BP <- readRDS(paste0( path_env_f, "/Climate_38.rds"))
 Climate_37K_BP <- readRDS(paste0( path_env_f, "/Climate_37.rds"))
 Climate_36K_BP <- readRDS(paste0( path_env_f, "/Climate_36.rds"))
 Climate_35K_BP <- readRDS(paste0( path_env_f, "/Climate_35.rds"))
 Climate_34K_BP <- readRDS(paste0( path_env_f, "/Climate_34.rds"))
 Climate_33K_BP <- readRDS(paste0( path_env_f, "/Climate_33.rds"))
 Climate_32K_BP <- readRDS(paste0( path_env_f, "/Climate_32.rds"))
 Climate_31K_BP <- readRDS(paste0( path_env_f, "/Climate_31.rds"))
 Climate_30K_BP <- readRDS(paste0( path_env_f, "/Climate_30.rds"))
 Climate_29K_BP <- readRDS(paste0( path_env_f, "/Climate_29.rds"))
 Climate_28K_BP <- readRDS(paste0( path_env_f, "/Climate_28.rds"))
 Climate_27K_BP <- readRDS(paste0( path_env_f, "/Climate_27.rds"))
 

 # COMPUTE SPECIES DISTRIBUTION MODELS ###
 
 ######################### Select the predictive variables
 
 # Define the modelling columns:

 names(db)
 dat <- db
 spc_col <- "presence"  # species presence/absence column is named "presence" in the dataset
 var_cols <-  names(db[,5:25]) 
 var_cols
 

    
  # Function to compute the Root Mean Square Error (RMSE):
  
  calculate_rmse <- function(model) {
    predicted_values <- predict(model)
    observed_values <- dat$presence
    rmse <- sqrt(mean((observed_values - predicted_values)^2))
    return(rmse)
  }
  
  # Function to compute the Variance Inflation Factor (VIF):
  
  calculate_vif <- function(model) {
    vif_values <- tryCatch({
      vif(model)
    }, error = function(e) {
      return(6) 
    })
    return(vif_values)
  }
  
  all_predictors = c("bio01","bio04", "bio05", "bio06", "bio07", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                    "bio17", "bio18", "bio19", "npp", "lai", "altitude", "rugosity")
  
  # Check highly correlated variables (r^2 > 0.8)
  head(dat)[,5:25]
  correlation_matrix <- cor(dat[,5:25])

  # Function to select variables with correlation coefficients < 0.8:
  
  select_variables <- function(matrix) {
    selected_variables <- character(0)
    remaining_variables <- colnames(matrix)
    
    while (length(remaining_variables) > 0) {
      reference_variable <- remaining_variables[1]
      selected_variables <- c(selected_variables, reference_variable)
      
      # Filter variable with correlation coefficients lower than 0.8
      filtered_variables <- remaining_variables[-1][abs(matrix[reference_variable, remaining_variables[-1]]) < 0.8]
      
      remaining_variables <- filtered_variables
    }
    
    return(selected_variables)
  }
  
  # Get the list of selected variables: 
  selected_vars <- select_variables(correlation_matrix)
  cat("Selected variables:", selected_vars, "\n")
  
  
  
  
  
  
  
  ## GENERALIZED LINEAR MODEL (GLM) ####

 form_glm <- as.formula(paste(spc_col, "~", paste(selected_vars, collapse = "+")))
 form_glm
 mod_glm <- glm(formula = form_glm, family = binomial, data = dat)
 summary(mod_glm)
 
 # Initialize the criteria used to select the best combination of predictive variables:
 best_AUC_value <- -Inf
 best_rmse_value <- Inf
 best_model <- NULL
 best_predictors <- NULL
 counter<-0
 
 dev.off()
 
 for (i in 2:length(selected_vars)) {
   combinations <- combn(selected_vars, i)
   
   for (j in 1:ncol(combinations)) {
     selected_predictors <- combinations[, j] # for each combination, we run a GLM:
     
     form_glm <- as.formula(paste(spc_col, "~", paste(selected_predictors, collapse = "+")))
     mod_glm <- glm(formula = form_glm, family = binomial, data = dat)
     summary(mod_glm)
     head(dat)
     # and get the p-value, AUC, VIG and RME:
     p_values <-summary(mod_glm)$coefficients[, "Pr(>|z|)"]
     AUC_value<- AUC(mod_glm)$AUC
     vif_values <- calculate_vif(mod_glm)
     rmse_value<-  calculate_rmse(mod_glm)
     
         # In any VIF > 5, the combination is discarded (Vif_exclusion == 1):
         if (any(vif_values > 5) ) {
           Vif_exclusion <- 1
         } else {
           Vif_exclusion <- 0
         }
         
         # Check whether the obtained model is the best model so far:
         if (Vif_exclusion == 0 && AUC_value > best_AUC_value || (AUC_value == best_AUC_value && rmse_value < best_rmse_value)) 
           {
             best_AUC_value <- AUC_value
             print(best_AUC_value)
             best_rmse_value <- rmse_value
             print(best_rmse_value)
             best_model_glm <- mod_glm
             best_predictors_glm <- selected_predictors
             print(best_predictors_glm)
           }
         }
       }
     
 
 
 # Outcomes:
 cat("Final best AUC: ", best_AUC_value, "\n")
 cat("Final best RMSE: ", best_rmse_value, "\n")
 cat("Final best predictors: ", best_predictors_glm, "\n")
 print(best_predictors_glm)

 # Variables used for each species to run the definitive model:
 form_glm <- as.formula(paste(spc_col, "~", paste(best_predictors_glm, collapse = "+")))
 form_glm
 mod_glm <- glm(formula = form_glm, family = binomial, data = dat)
 AUC(mod_glm)$AUC
 
 
 ## GENERALIZED ADDITIVE MODEL (GAM) ####

 form_gam <- as.formula(paste(spc_col, "~", paste0("s(",  selected_vars, ")", collapse = "+"))) 
 form_gam
 
 mod_gam <- gam(formula = form_gam, family = binomial, data = dat)
 summary(mod_gam)

 # Initialize the criteria used to select the best combination of predictive variables:
 best_AUC_value <- -Inf
 best_rmse_value <- Inf
 best_model <- NULL
 best_predictors <- NULL
 counter<-0
 dev.off()
 
 for (i in 2:length(selected_vars)) {
   combinations <- combn(selected_vars, i)# for each combination, we run a GAM:
   
   for (j in 1:ncol(combinations)) {
     selected_predictors <- combinations[, j]
     form_gam <- as.formula(paste(spc_col, "~", paste0("s(",  selected_predictors, ")", collapse = "+")))  
     mod_gam <- gam(formula = form_gam, family = binomial, data = dat)
     summary(mod_gam)
     head(dat)
     # and get the p-value, AUC, VIG and RME:
     p_values <-summary(mod_gam)$coefficients[, "Pr(>F)"]
     AUC_value<- AUC(mod_gam)$AUC
     vif_values <- calculate_vif(mod_gam)
     rmse_value<-  calculate_rmse(mod_gam)
     
     # In any VIFs > 5 = Vif_exclusion == 1
     if (any(vif_values > 5) ) {
       Vif_exclusion <- 1
     } else {
       Vif_exclusion <- 0
     }
     
     # Check whether this is the best model obtained so far
     if (Vif_exclusion == 0 && AUC_value > best_AUC_value || (AUC_value == best_AUC_value && rmse_value < best_rmse_value)) {
       
         best_AUC_value <- AUC_value
         print(best_AUC_value)
         best_rmse_value <- rmse_value
         print(best_rmse_value)
         best_model_gam <- mod_gam
         best_predictors_gam <- selected_predictors
         print(best_predictors_gam)
       }
     }
   }
 
 

 # Outcomes:
 cat("Final best AUC: ", best_AUC_value, "\n")
 cat("Final best RMSE: ", best_rmse_value, "\n")
 cat("Final best predictors: ", best_predictors_gam, "\n")
 
 print(best_predictors_gam)
 # Check multicollinearity:
 corSelect(data = dat, var.cols =  best_predictors_gam, cor.thres = 0.8, select = "VIF")
 
 selected_vars_gam <- best_predictors_gam
 
 # Run the definitive model
 form_gam <- as.formula(paste(spc_col, "~", paste0("s(",  selected_vars_gam, ")", collapse = "+")))  # GAM with smoothing splines ('s')
 form_gam
 mod_gam <- gam(formula = form_gam, family = binomial, data = dat)
 AUC(mod_gam)$AUC
 
 
 ## MAXIMUM ENTROPY (MAXENT) ####
 

 mod_mxt <- maxnet(p = dat[ , spc_col], data = dat[ , selected_vars], f = maxnet.formula(dat[ , spc_col], dat[ , selected_vars])) 
 summary (mod_mxt)

# VIF cannot be computed from MAXENT with the previous loop, so lets test first the VIF manually:
 
 corSelect(data = dat, var.cols =   selected_vars, cor.thres = 0.8, select = "VIF")
 # Now introduce into option1 the variables with VIF < 5:
 option1<- c(" ") 
 # Now check VIF < 5 in the selected variables:
 corSelect(data = dat, var.cols =   option1, cor.thres = 0.8, select = "VIF")
 selected_vars<- option1
 selected_vars


# From here, we can apply the same loop previously used:
best_AUC_value <- -Inf
best_rmse_value <- Inf
best_model <- NULL
best_predictors <- NULL
counter<-0

for (i in 2:length(selected_vars)) {
  combinations <- combn(selected_vars, i)
  cat("Combinación ", counter, ": ", combinations, "\n")
  
  for (j in 1:ncol(combinations)) {
    selected_predictors <- combinations[, j]
    mod_mxt <- maxnet(p = dat[ , spc_col], data = dat[ , selected_predictors], f = maxnet.formula(dat[ , spc_col], dat[ , selected_predictors])) 
    summary(mod_mxt)
    # Get predictions from the model
    predicted_scores <- predict(mod_mxt, newdata = dat[ , selected_predictors])
    
    # Compute AUC and RMSE
    roc_obj <- roc(dat[ , spc_col], predicted_scores)
    AUC_value <- auc(roc_obj)
    plot(AUC_value, main=paste0( species, " :", AUC_value), col="red")
    rmse_value<-  sqrt(mean((predicted_scores - dat[, spc_col])^2))
    
    # Check whether this is the best model so far
    if (AUC_value > best_AUC_value || (AUC_value == best_AUC_value && rmse_value < best_rmse_value)) 
      {
        best_AUC_value <- AUC_value
        print(best_AUC_value)
        best_rmse_value <- rmse_value
        print(best_rmse_value)
        best_model_mxt<- mod_mxt
        best_predictors_mxt <- selected_predictors
        print(best_predictors_mxt)
      }
    }
  }


# Resultados del mejor modelo encontrado al final del bucle
cat("Final best AUC: ", best_AUC_value, "\n")
cat("Final best RMSE: ", best_rmse_value, "\n")
cat("Final best predictors: ", best_predictors_mxt, "\n")
print(best_predictors_mxt)
selected_vars
# Check multicolinearity
corSelect(data = dat, var.cols =  best_predictors_mxt, cor.thres = 0.8, select = "VIF")

selected_vars_maxt <-best_predictors_mxt

# Select the definitive model:

mod_mxt <- maxnet(p = dat[ , spc_col], data = dat[ , selected_vars_maxt], f = maxnet.formula(dat[ , spc_col], dat[ , selected_vars_maxt])) 
summary (mod_mxt)

# Get predictions
predicted_scores <- predict(mod_mxt, newdata = dat[ , selected_vars])

# Compute AUC
roc_obj <- roc(dat[ , spc_col], predicted_scores)
auc <- auc(roc_obj)
auc


# If you have changed selected_vars because of the VIF (in lines 600-606), now select again the original selected_vars:

selected_vars <- select_variables(correlation_matrix) # = line 443

# Print
cat("Selected variables:", selected_vars, "\n")


 ## BAYESIAN ADDITIVE REGRESSION TREES (BART) ####
 
 # BART selects the minimal subset of relevant variables. BART runs the model with different variable combinations
 # and select those with the lowest RMSE.
  set.seed(123)
  n_iter <- 50 # if it takes too long, reduce the number.
  varselect_bart <- variable.step(x.data = dat[ , selected_vars], y.data = dat[ , spc_col], n.trees = 10, iter= n_iter)
  varselect_bart

  # Now we check the VIF of the selected variables:
  corSelect(data = dat, var.cols =  varselect_bart, cor.thres = 0.8, select = "VIF")
  
 # If VIF > 5 and the RMSE is reduced meaningfully,and use "varlselect_bart" in the mod_bart_final (below),
 # If the RMSE remains identical or similar, or VIF > 5, we keep "selected_vars"
  
 mod_bart <- dbarts::bart(x.train = dat[ , varselect_bart], y.train = dat[ , "presence"], keeptrees = TRUE)
 summary( mod_bart)
 # If you want to use this BART model in future R sessions, you need to explicitly
 # ask for the full information to be included when you next save the model object (see "Saving" section in ?bart):
 invisible(mod_bart$fit$state)
 
 # For the predictive variables in BART we have used the RMSE and VIF, but unlike with the previous alorithms, AUC values
 # have not been  used. Therefore, we repeat the procedure to check the AUC value obtained from all combination possible
 # of predictive variables:
 

 # Loop:
 best_AUC_value <- -Inf
 best_model <- NULL
 best_predictors <- NULL
 counter<-0
 
 for (i in 2:length(varselect_bart)) {
   combinations <- combn(varselect_bart, i)
   cat("Combinación ", counter, ": ", combinations, "\n")
   
   for (j in 1:ncol(combinations)) {
     selected_predictors <- combinations[, j]
     mod_bart <- bart(x.train = dat[ , selected_predictors], y.train = dat[ , spc_col], keeptrees = TRUE)
     summary(mod_bart)
   # Get AUC
     match_auc <- regmatches(capture.output(summary(mod_bart)), regexpr("AUC = [0-9.]+", capture.output(summary(mod_bart))))
     AUC_value <- as.numeric(gsub("AUC = ", "", match_auc))
     
     # Check whether this the best AUC so far
     if (AUC_value > best_AUC_value ) {
       best_AUC_value <- AUC_value
       print(best_AUC_value)
       best_model_bart<- mod_bart
       best_predictors_bart <- selected_predictors
       print(best_predictors_bart)
     }
   }
 }
 
 summary(best_model_bart)
 varselect_bart<- c(" ") # Include the predictive variables
 mod_bart <- dbarts::bart(x.train = dat[ , varselect_bart], y.train = dat[ , "presence"], keeptrees = TRUE)
 
 
 # SAVE ALL INDIVIDUAL SPECIES DISTRIBUTION MODELS ####
 
 models <- list(mod_glm ,mod_gam, mod_mxt, mod_bart)
 names(models) <- c("glm", "gam", "mxt", "bart")
 pred_folder <- paste0( species, "/predictions")
 dir.create(pred_folder, recursive = TRUE)
 saveRDS(models, file = paste0(pred_folder, "/models.rds"))
 
 # COMPUTE AND COMPARE THE AUC, BOYCE INDEX & VARIABLE IMPORTANCE ####
 
 # IF YOU HAVE PREVIOUSLY RUN THE CODE AND SAVED THE GAM, MXT AND BART MODELS, JUST LOAD IT:
 
 # Load models:
 setwd(" ") # write here the directory
 getwd()

 species <-"Panthera.spelaea" # Write here the neame of the species
 models <- readRDS(paste0(species, "/predictions/models.rds"))
 summary(models)

# Ara Under the Curve (AUC) values

 round(AUC(models$glm)$AUC, 2)
 round(AUC(models$gam)$AUC, 2)
 round(AUC(models$bart)$AUC, 2)
 # for the MAXENT AUC, we need to compute it:
 species
 db<- read.xlsx(paste0( species, "/ALL", species, ".xlsx"), rowNames = FALSE, colNames = TRUE, sheet="Sheet 1")
 head(db)
 nrow(db)
 dat<-db
 predicted_scores <- predict(models$mxt, newdata = dat)
 roc_obj <- roc(dat[ , "presence"], predicted_scores)
 auc(roc_obj)
 round(auc(roc_obj), 2)
 
 
# Boyce Index values and Variable Importance:

dat<- db
head(dat)
# Boyce Index for GLM
predicted <- predict(models$glm, newdata = dat, type = "response")
glm_fav <- Fav(pred = predicted, sample.preval = prevalence(model = models$glm))
Boyce(obs = dat$presence, pred = glm_fav)$Boyce # If you get the warning message saying that some bins have
# less than 30 values, then the Boyce Index value is msileading. To fix it, include "bin.width = 0.2" in the chunk and run it again

# Variable Importance for GLM
caret::varImp(models$glm) # for details, see ?caret::varImp



# Boyce Index for GAM
predicted <- predict(models$gam, newdata = dat, type = "response")
gam_fav <- Fav(pred = predicted, sample.preval = prevalence(model = models$gam))
Boyce(obs = dat$presence, pred = gam_fav)$Boyce

# Variable Importance for GAM
caret::varImp(models$gam)# for details, see ?caret::varImp

# Boyce Index for MAXENT
predicted <- predict(models$mxt, newdata = dat, clamp = FALSE, type = "cloglog")
mxt_fav <- models$mxt
Boyce(obs = dat$presence, pred = predicted)$Boyce

# Variable Importance for MAXENT
varImportance(mod_mxt, occSWD =dat[dat$presence == "1", ], bkgSWD=dat[dat$presence == "0", ])
?varImportance
# Boyce Index for BART
summary(models$bart)

bart_selected_vars <- c(" ") # Include here the variables used by the model (available in "Predictor list")


model_bart<- bart(y.train = dat[ , "presence"], x.train = dat[ , bart_selected_vars], keeptrees = TRUE, verbose = FALSE)
predicted<- colMeans(stats::predict(model_bart, dat))
bart_fav <- Fav(pred = predicted, sample.preval = prevalence(model = models$bart))
Boyce(obs = dat$presence, pred = bart_fav)$Boyce

# Variable Importance for BART
varimp(models$bart) # for details, see ?varimp

# PLOT VARIABLE IMPORTANCE:
setwd(" ") # write directory where all outcomes are available
getwd()

VI_df<-read.xlsx("F:/SPD_NP_v2_final/SubmitPaper/Outcomes/Varimp.xlsx", rowNames=FALSE, 
                         colNames=TRUE, sheet="Var_Imp_Ensemble")

head(VI_df)


# Define the variables order:
order_vars <- c("Bio01", "Bio04", "Bio05", "Bio06", "Bio07", 
                "Bio08", "Bio09", "Bio10", "Bio11", "Bio12", 
                "Bio13", "Bio14", "Bio15", "Bio16", "Bio17", 
                "Bio18", "Bio19", "npp", "lai", "altitude", "rugosity")

# Filter species to get your target species
VI_df_filtered <- VI_df %>% 
  filter(Species %in% c("Homo neanderthalensis", "Homo sapiens"))

# Format data
VI_df_filtered_long <- gather(VI_df_filtered, variable, value, -Species)
# Variables order
VI_df_filtered_long$variable <- factor(VI_df_filtered_long$variable, levels = order_vars)
# Plot
plot_Homo<- ggplot(VI_df_filtered_long, aes(x = variable, y = value, color = Species, group=Species, fill = Species)) +
 # geom_line() +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
  #geom_point()+
  labs(title = "Homo",
       x = "Variable",
       y = "Importance",
       color = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_Homo

# PANTHERA
VI_df
VI_df_filtered <- VI_df %>% 
  filter(Species %in% c("Panthera leo", "Panthera spelaea", "Panthera pardus", "Felis silvestris"))

VI_df_filtered_long <- gather(VI_df_filtered, variable, value, -Species)

VI_df_filtered_long$variable <- factor(VI_df_filtered_long$variable, levels = order_vars)

plot_Felidae<- ggplot(VI_df_filtered_long, aes(x = variable, y = value, color = Species, group=Species, fill = Species)) +
  # geom_line() +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
  #geom_point()+
  labs(title = "Felidae",
       x = "Variable",
       y = "Importance",
       color = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_Felidae

# URSIDAE
VI_df
VI_df_filtered <- VI_df %>% 
  filter(Species %in% c("Ursus arctos", "Ursus spelaeus"))

VI_df_filtered_long <- gather(VI_df_filtered, variable, value, -Species)

VI_df_filtered_long$variable <- factor(VI_df_filtered_long$variable, levels = order_vars)

plot_Ursidae<- ggplot(VI_df_filtered_long, aes(x = variable, y = value, color = Species, group=Species, fill = Species)) +
  # geom_line() +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
  #geom_point()+
  labs(title = "Ursidae",
       x = "Variable",
       y = "Importance",
       color = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_Ursidae

# CANIDAE
VI_df
VI_df_filtered <- VI_df %>% 
  filter(Species %in% c("Canis lupus", "Vulpes vulpes", "Vulpes lagopus"))

VI_df_filtered_long <- gather(VI_df_filtered, variable, value, -Species)

VI_df_filtered_long$variable <- factor(VI_df_filtered_long$variable, levels = order_vars)

plot_Canidae<- ggplot(VI_df_filtered_long, aes(x = variable, y = value, color = Species, group=Species, fill = Species)) +
  # geom_line() +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
  #geom_point()+
  labs(title = "Canidae",
       x = "Variable",
       y = "Importance",
       color = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_Canidae

# HYAENA AND MELES
VI_df
VI_df_filtered <- VI_df %>% 
  filter(Species %in% c("Crocuta crocuta", "Meles meles"))

VI_df_filtered_long <- gather(VI_df_filtered, variable, value, -Species)

VI_df_filtered_long$variable <- factor(VI_df_filtered_long$variable, levels = order_vars)

plot_Others<- ggplot(VI_df_filtered_long, aes(x = variable, y = value, color = Species, group=Species, fill = Species)) +
  # geom_line() +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
  #geom_point()+
  labs(title = "Others",
       x = "Variable",
       y = "Importance",
       color = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_Others

grid.arrange(plot_Homo, plot_Ursidae, plot_Felidae, plot_Canidae, plot_Others, ncol=2)


 # CROSS-BLOCK VALIDATION ####

  

  
  # DIVIDE STUDY AREA INTO SPATIAL BLOCKS ###
  
  # first, convert 'dat' to a spatial points object with its CRS (mind that you
  # need to ascertain the correct CRS to which your coordinates refer:
  dat <- db
  nrow(dat)
  names(dat)
  head(dat)
  dat_points <- vect(dat, geom = c("longitude", "latitude"), crs = "epsg:4326")
  
  # the size (or range) of the spatial blocks is sometimes defined as the range of
  # spatial autocorrelation in the variables (which you can compute with function
  # 'spatialAutoRange'), but this is now discouraged (see Wadoux et al. 2021,
  # https://doi.org/10.1016/j.ecolmodel.2021.109692) here we'll use a size we
  # consider reasonable for our species and study area: 500-km blocks. However, for each species,
  # you should check that all class blocks (identified with colors or numbers) have occurrences (points):
  layers_cut <- crop(climate_layer[[1]], mywindow, mask = TRUE)
  plot(layers_cut)
  blocks <- spatialBlock(speciesData = as(dat_points, "Spatial"), 
                         rasterLayer = raster(layers_cut[[1]]), theRange =  500000, # 500000 / 1000000
                         k = 5, seed = 2023)  # you can also use an optional additional argument, species = spc_col; see ?spatialBlock
  
  # add the spatial block ID to the dataset:
  dat$foldID <- blocks$foldID
  
  # map blocks and presences:
  dev.off()
  plot(vect(blocks$blocks), "folds", col = hcl.colors(5, "geyser"), main = species)
  plot(countries, lwd = 1, add = TRUE)
  plot(subset(dat_points, dat_points$presence == 1), col = "blue", cex = 0.8, add = TRUE)
  
  # GET PREDICTIONS FOR CROSS-VALIDATION ###
  
  # (i.e., build models again, leaving out each fold in turn)
  
  folds <- sort(unique(dat$foldID))
  names(dat)
  spc_col <- "presence"  
  
  # Specify the specific predictive variables for the MAXENT and BART models:
  summary(models$bart)
  bart_selected_vars<- c(" ") # Introduce the predictive variables
  models$mxt$varmin
  mxt_selected_vars<- c(" ") # Panthera spelaea
  

# BLOCK VALIDATION:
  for (f in folds) {
    cat("modelling outside fold", f, "...\n")  # inform of progress
    dat_train <- subset(dat, foldID != f)
    
    mod_glm_fold <- glm(formula = models$glm$formula, family = binomial, data = dat_train)
    dat[ , paste0("glm_fold", f, "_p")] <- predict(mod_glm_fold, newdata = dat, type = "response")
    
    mod_gam_fold <- gam(formula = models$gam$formula, family = binomial, data = dat_train)
    dat[ , paste0("gam_fold", f, "_p")] <- predict(mod_gam_fold, dat, type = "response")
 

    mod_mxt_fold <- maxnet(p = dat_train[ , spc_col], data = dat_train[ , mxt_selected_vars], f = maxnet.formula(dat_train[ , spc_col], dat_train[ , mxt_selected_vars]))  # you can add to 'maxnet.formula' e.g. classes="lq", to use only linear ('l') and quadratic ('q') features
    dat[ , paste0("mxt_fold", f, "_p")] <- predict(mod_mxt_fold, dat, clamp = FALSE, type = "cloglog")
    

    mod_bart_fold <- bart(y.train = dat_train[ , spc_col], x.train = dat_train[ , bart_selected_vars], keeptrees = TRUE, verbose = FALSE)
    dat[ , paste0("bart_fold", f, "_p")] <- colMeans(stats::predict(mod_bart_fold, dat))
  
    
    gc()  # cleanup memory before next loop iteration
  }  
  
  # see the new predictions added to the data frame:
  head(dat)
  

  # EVALUATE EACH MODEL ON ITS VALIDATION FOLD ###
  
  fold_cols <- grep("_fold", names(dat))
  names(dat)[fold_cols]
  
  # choose e.g. 3 metrics with complementary information (discrimination, classification, calibration):
  metrics <- c("AUC")
  
  # create an empty table to receive the cross-validation results:
  crossval <- as.data.frame(matrix(nrow = length(folds), ncol = length(metrics) * length(names(models))))
  colnames(crossval) <- c(outer(1, metrics, FUN = paste, sep = "_"))
  crossval  # for now it's only filled with NAs
  

  for (m in names(models))  for (f in folds) {
    fold_name <- paste0("fold", f)
    fold_col <- names(dat)[grep(paste0(m, "_fold", f), names(dat))]
    fold_dat <- subset(dat, foldID == f)
    crossval[f, paste(m, "AUC", sep = "_")] <- AUC(obs = fold_dat[ , spc_col], pred = fold_dat[ , fold_col], simplif = TRUE, plot = TRUE, main = paste(m, "AUC"))
  }
  
  
  # get more information with boxplots of the cross-validation metrics:

  # remove empty columns
  empty_columns <- sapply(crossval, function(x) all(is.na(x) | x == ""))
  dev.off()
  e<- crossval[, !empty_columns]
  e
  boxplot(e, col = grey(10:1/10), each = length(names(models)), las = 2, main = species)
  abline(h = 0.7, col = "darkred", lty = 1, lwd = 2)  

  # save outputs
  write.csv(crossval, paste0("write_here_your_directory", species, "/performance.csv"))



 # MAKE PREDICTIONS AND GET THE RESULTS FOR EACH TIME STEP OF THE MIS3 ####
  
  # the next command creates a folder for output files:
  # IMPORTANT! Check your directory (it should be outside the species folder):
  getwd()
  # If that is not the correct directory, change it here: 
  setwd(" ") # e.g. setwd("F:/Data&Code")
  dir.create(paste0(species, "/Layers") )  # '../' goes up one level from the current working directory, so this creates the 'outputs' folder just outside the 'scripts' folder
  # Now two additional folders are created, one for the .TIF files and other one for the .JPEG files
  dir.create(paste0(species, "/Layers/TIF") )
  dir.create(paste0(species, "/Layers/JPEG") ) 
  
  # Specify the directory for the SDMs' projections:
  my_path<-"F:/SPD_NP_v2_final/SubmitPaper/Data&Code" # e.g. my_path<-"F:/Data&Code"
  path_TIF <- paste0(my_path, "/", species, "/Layers/TIF")
  path_JPEG <- paste0(my_path, "/", species, "/Layers/JPEG") 
 
 
 # Now, for visual representation purposes, we will include the coast line
 getwd()
 coast_file <- read.csv("coastline__120m.csv")
 
 # Convert coordinates from WKT to geometric objects
 coast_geom <- st_as_sfc(coast_file$the_geom, crs = 4326)
 
 # Create object sf
 coast_sf <- st_sf(coast_geom, data = coast_file)

 # create a list with all raster files:
 f <- list.files(path=path_env_f, pattern='rds$', full.names=TRUE) 
 # Set colors
 mi_paleta <- colorRampPalette(brewer.pal(11, "RdYlBu"))
 clrs <- rev(mi_paleta(100))

 dev.off()
 
 
 
 # Before running the next chunk, check you have the "performance_alldata.csv" file with the AUC values
 # of each model to compute the weighted mean value in the ensemble projection. This file can be created
 # with the performance.csv file saved in line 1061, or with the AUC values obtained
 # from the models in lines 766-778.
 # IMPTORTANT! If any model (GLM, GAM, MAXENT or BART) has a mean AUC or Boyce Index value < than 0.7,
 # go to line 1159 before running this loop:



 for (i in 1:length(f)) {
   raster <- readRDS(f[[i]])
  
   glm_p <- predict(stack(raster), mod_glm, type = "response")
   glm_fav <- Fav(pred = glm_p, sample.preval = prevalence(model = mod_glm))
   plot(glm_fav, col = clrs, range = c(0, 1), main = "GLM")
   glm_fav <- as(glm_fav, "SpatRaster")
   
   gam_p <- predict(stack(raster), mod_gam, type = "response")
   gam_fav <- Fav(pred = gam_p, sample.preval = prevalence(model = mod_gam))
   plot(gam_fav, col = clrs, range = c(0, 1), main = "GAM")
   gam_fav <- as(gam_fav, "SpatRaster")
   
   mxt_p <- rast(predict(stack(raster), mod_mxt, type = "cloglog"))  # raster format conversions needed because 'predict' for MX is not working correctly with 'terra' raster format
   plot(mxt_p, col = clrs, range = c(0, 1), main = "MXT")
   plot(coastsCoarse, add=T, col="black")
   mxt_fav <- as(mxt_p, "SpatRaster")
   
   bart_p <- rast(predict(model_bart, x.layers = stack(raster)))
   bart_fav <- Fav(pred = bart_p, sample.preval = prevalence(model = mod_bart))
   plot(bart_fav,  range = c(0, 1),col=clrs, main = paste0("", (26 + i), "bart k yrs BP"))
   bart_fav <- as(bart_fav, "SpatRaster")
 
   # Ensemble during the loop
   # read csv with cross-validation outcomes:
   
   ens_df<- read.csv(paste0(my_path, "/", species, "/performance_alldata.csv"), sep = ";")

   AUC_GLM <- mean(ens_df$glm_AUC)
   AUC_GAM <- mean(ens_df$gam_AUC)
   AUC_MXT <- mean(ens_df$mxt_AUC)
   AUC_BART <- mean(ens_df$bart_AUC)
   

   # Ensemble favourability
   glm_fav <- as(glm_fav, "SpatRaster")
   gam_fav <- as(gam_fav, "SpatRaster")
   
   # Ensure that below, "b" and "weighted.mean()" only include the models with AUC and BI >=0.7
    b <- brick(x = c(glm_fav, gam_fav, mxt_p, bart_fav))
    ensemble_f <- weighted.mean(b, w= c(AUC_GLM, AUC_GAM, AUC_MXT, AUC_BART))

   plot(ensemble_f, col = clrs, range = c(0, 1), main = "ENSEMBLE f")
  
   # Save rasters as tif
   writeRaster(ensemble_f,paste0(path_TIF, "/", (26+ i), "f.tif"),options=c('TFW=YES'))
   
   # Save rasters as jpeg to make a video
   jpeg(paste0(path_JPEG, "/", (40 - i), "p.jpeg"),  units="in", width=9, height=7, res=300)

   dev.off()
   
  
 }
 


############# COMPUTE NICHE BREADTH ########
path_TIF<- NULL
species<- "Panthera.spelaea" # Specify the species name
path_TIF <- paste0(my_path, "/", species, "/Layers/TIF")
path_TIF


# Get the .tif files
f <- list.files(path_TIF,
                         pattern = "\\.tif$",
                         recursive = FALSE,
                         all.files = FALSE,
                         full.names = TRUE)
# Select the favourability files
f <- f[grep("f\\.tif$", f)]


# First, we create a new data frame to save the outputs obtained from the niche breadth computations
Time_k_BP <- seq(from = 27, to = 55, by = 1)
Output <- data.frame(Time_k_BP)
Output$NichBrd_1_f <- NA
Output$NichBrd_2_f<- NA
Output$Favourable_area<- NA
Output$Favourable_area_abs<- NA
Output$Area_km<- NA
Output$Species <- species



head(Output)

# Get the numbers
file_numbers <- as.integer(gsub("[^0-9]", "", basename(f)))

# Order the numbers
order_index <- order(file_numbers, decreasing = FALSE)

# Order the f list accordingly 
f_sorted <- f[order_index]

# Get results and store data
results <- list()
for (i in 1:length(f_sorted)) {
  archivo_raster <- f[i]
  
  # Carga el raster desde el archivo
  raster <- raster::raster(archivo_raster)
  ensemble_f <- raster

  ensemble_f<- as(ensemble_f, "SpatRaster")
  NicheBreadth_f<- raster.breadth(ensemble_f,verbose=FALSE) # Measures the spatial heterogeneity of the distribution of suitability scores: Levin's B1 and B2
  Output$NichBrd_1_f[[i]]<- NicheBreadth_f$B1
  Output$NichBrd_2_f[[i]]<- NicheBreadth_f$B2
  raster <- ensemble_f >= 0.5 # Pixels with favourability values =>0.5 are considered as favourables.
  plot(raster)
  Output$Favourable_area_abs[[i]]<- sum(raster[] == 1, na.rm = TRUE)
  absence <- sum(raster[] == 0, na.rm = TRUE)
  Output$Favourable_area[[i]]<-sum(raster[] == 1, na.rm = TRUE) / sum(raster[] == 0, na.rm = TRUE)
  # area in km
  r <- raster(raster)
  polys <- rasterToPolygons(r, fun=function(x){x==1}, dissolve=TRUE)
  Output$Area_km[[i]] <- sum(areaPolygon(polys)) / 1E6  # en km2
  
}

Output
getwd()
# save data
my_path
write.xlsx(Output, file = paste0(my_path, "/", species, "/Output_NB_f", species, ".xlsx") )









# SENSITIVITY TEST 1: AGE UNCERTAINTY ####
# This sensitivity tests does exactly the same procedures as above to run the SDMs
# but, in this case, environmental and climatic variables are not obtained from the
# median age of each archaeo-paleontological assemblage, but from 100 random ages
# obtained from a normal distribution around the mean age:
 
Presence_data<-read.xlsx("Data.xlsx", rowNames=FALSE, 
                         colNames=TRUE, sheet="F_df")

head(Presence_data)

# for each species just write exactly the same name below where you see "specify the species"
species <- "Panthera.spelaea" # specify the species
db <- subset(Presence_data,Panthera.spelaea==1) #specify the species

head(db)
nrow(db)
# Create as new data frame

Presence_data<- data.frame(db$`Site/Level`)
Presence_data$x <- db$x
Presence_data$y<- db$y
Presence_data$Age <- db$Age
Presence_data$SD <- db$SD
Presence_data$presence <- 1
Presence_data




# Create a new folder for each species. Here will be saved all outcomes obtained.
pres_folder <- paste0(species, "/Sensitivity") 
pres_folder
dir.create(pres_folder)  # '../' goes up one level from the current working directory, so this creates the 'species' folder 


# Select all variables
bio_variables = c("bio01","bio04", "bio05", "bio06", "bio07", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                  "bio17", "bio18", "bio19", "npp", "lai", "altitude", "rugosity")


#Specify the specific predictive variables
selected_vars_glm <- c( " " ) 
selected_vars_gam <- c(" ") 
selected_vars_maxt <- c(" ") 
varselect_bart <- c(" ") 

 
 
 # New Dataframe

 Output_AUC <- data.frame("auc_bart" = c(NA),
                          "AUC_maxent" = c(NA),
                          "AUC_GAM"= c(NA),
                          "AUC_GLM"= c(NA))
 # First, we create a new data frame to save the outputs obtained from the niche breadth computations
 Time_k_BP <- NA
 Output <- data.frame(Time_k_BP)
 Output$Experiment <- NA
 Output$NichBrd_1_f <- NA
 Output$NichBrd_2_f<- NA
 Output$Favourable_area<- NA
 Output$Favourable_area_abs<- NA
 Output$Area_km<- NA
 Output$Species <- species
 
 head(Output)
 ?rnorm
 
 # Create a function to randomly select climate variables according to the 95% CI
 # of the archaeo-paleontological levels:
 
 Sensitivity_Age <- function(x) 
 {

   for(a in 1:100) {                                          
     n <- 1
   # Resample chronology following a normal distribution and taking into account the mean and standard deviation of each chronology
     func1 <- function(x, n) rnorm(n, mean = x[1], sd = x[2])
     df <- apply(Presence_data[, 4:5], 1, FUN = func1, n = n)
     Presence_data$df <- df
     head(Presence_data)
     Presence_data$time_bp <- (Presence_data$df * -1) # in pastclim Age  should be negative, so multiply it by * -1 
     head(Presence_data)
    
    
     # Get varibles
     colnames(Presence_data)[2] = "longitude"
     colnames(Presence_data)[3] = "latitude"
     
     # get the values of each variable for each location with the presence of the selected species:
     
     Variables<- location_slice( x = Presence_data, bio_variables = c("bio01","bio04", "bio05", "bio06", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                                                                      "bio17", "bio18", "bio19",  "npp",  "lai", "altitude", "rugosity"), dataset = "Krapp2021", nn_interpol = TRUE, buffer=TRUE)
     
     head(Presence_data)
     head(Variables)
   
     # Replace NA values in Variables with Presence_data values
     Variables <- Variables %>%
       mutate(
         db..Site.Level.= coalesce(db..Site.Level., Presence_data$db..Site.Level),
         longitude = coalesce(longitude, Presence_data$longitude),
         latitude = coalesce(latitude, Presence_data$latitude),
         Age = coalesce(Age, Presence_data$Age),
         SD = coalesce(SD, Presence_data$SD),
         presence = coalesce(presence, Presence_data$presence),
         time_bp = coalesce(time_bp, Presence_data$time_bp)
       )
     
     Variables$time_bp_slice <- round(Variables$time_bp, -3)
     head(Variables)
     
     
     Presence_data <-Variables
     Presence_data <- subset(Presence_data, select = -df)
     head(Presence_data)
     nrow(Presence_data)
     
     Presence_data<- Presence_data[complete.cases(Presence_data), ]
     
     
     # _________
     #
     # BUFFER
     #
     #__________
     
     countries <- readRDS("gadm36_adm0_r5_pk.rds")
     
     
     pres_points <- vect(Presence_data, geom = c("longitude", "latitude"), keepgeom = TRUE, crs = "epsg:4326")  
     buff_dist <- max(distance(pres_points) * 0.1 , na.rm = TRUE)  # can take time if there are many points!
     buff_dist
     pres_buff <- aggregate(buffer(pres_points, width = buff_dist))
     plot(pres_buff)
     plot(pres_points, col = "blue", add = TRUE)
     plot(countries, border = "grey", add = TRUE)
     
     plot(climate_layer[[1]])
     layers_cut <- crop(climate_layer[[1]], pres_buff, mask = TRUE)
     
     Outcomes_v_thinned <- data.frame()
     
     # The following command thins presence point to avoid over-representation of points. There will be only 1 presence/pixel/age
     for (step in unique(Presence_data$time_bp_slice)){
       dat <- subset(Presence_data, time_bp_slice== step) # select only points with same chronology
       if (nrow(dat)> 0) 
       {
         dat <- gridRecords(rst = layers_cut[[1]], pres.coords = dat[ , c("longitude", "latitude")], plot = TRUE) # thin, for details see ?gridRecords
         dat<- subset(dat, presence==1)
         nrow(dat)
         colnames(dat)[2] = "longitude"
         colnames(dat)[3] = "latitude"
         dat$time_bp <- step
         dat<- dat[ , !(names(dat) %in% "bio01")] # the raster with the variable bios01 was used for thinning, now we remove this column from dat
         Variables_thinned<- location_slice( x = dat, bio_variables = c("bio01","bio04", "bio05", "bio06", "bio07", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16",
                                                                        "bio17", "bio18", "bio19", "npp", "lai", "altitude", "rugosity"), dataset = "Krapp2021", nn_interpol = TRUE, buffer=TRUE)
         head(Variables_thinned)
         head(dat)
         Variables_thinned <- Variables_thinned %>%
           mutate(
             presence= coalesce(presence, dat$presence),
             longitude = coalesce(longitude, dat$longitude),
             latitude = coalesce(latitude, dat$latitude),
             cells = coalesce(cells, dat$cells),
             time_bp = coalesce(time_bp, dat$time_bp),
             time_bp_slice = coalesce(time_bp_slice, dat$time_bp_slice)
             
           )
         
         
         Outcomes_v_thinned <- rbind(Outcomes_v_thinned, Variables_thinned)
       }
     }
     
     
     #PLOT RAW PRESENCE POINTS AND FILTERED PRESENCES
     plot(countries, ext = mywindow)
     plot(countries, xlim = mywindow[1:2], ylim = mywindow[3:4])
     points(Presence_data[ , c("longitude", "latitude")], col = "blue")
     points(Outcomes_v_thinned[ , c("longitude", "latitude")], col = "red")
     nrow(Presence_data)
     nrow(Outcomes_v_thinned)
     
     Outcomes_v_thinned
     
     Presence_data_thinned <-Outcomes_v_thinned
     
     head(Presence_data_thinned)
     
     
# BASELINE ENVIRONMENTAL CONDITIONS
     
     # First, we create an empty list:
     climate_list<-list()
     
     for (step in unique(Presence_data_thinned$time_bp)){
       # extract climate for the world
       step_climate <- region_slice(
         time_bp = step,
         bio_variables = bio_variables,
         dataset = "Krapp2021"
       )
       
       # subset calibration area:
       crs(step_climate) <- crs(pres_buff)
       step_climate <- terra::mask(step_climate, pres_buff)# select calibration area (press_buff) 
       # Select number of observations per time step / chronology:
       n <- subset(Presence_data_thinned, time_bp== step) # select only points with same chronology
       nrow(n) 
       # Select pesudoabsences:
       step_values <- terra::spatSample(step_climate, 50 * nrow(n), na.rm=TRUE, cells=TRUE, xy=TRUE) # Change 50. Amend it
       step_values$time_bp<-step # Time step
       climate_list[[as.character(step)]]<-step_values
     }
     
     # combine all into a single matrix and include "presence", "x" (longitude) and "y" (latitude) variables
     baseline_climate <- do.call(rbind, climate_list)
     baseline_climate$presence<- 0
     baseline_climate$longitude<- baseline_climate$x
     baseline_climate$latitude<- baseline_climate$y
     nrow(baseline_climate)

     
     # merge datasets
     
     head(baseline_climate)
     head(Presence_data_thinned)
     
     db<- rbind.fill(Presence_data_thinned, baseline_climate)
     colnames(db)
     
     # Now, lets eliminate the empty or unnecesary columns for the next steps:
     db <- db[, !(colnames(db) %in% c("cells", "time_bp_slice", "cell", "x", "y"))]
     colnames(db)

     # SPECIES DISTRIBUTION MODELS // SENSITIVITY TEST #####
     
     
     # define the modelling columns:
     
     names(db)
     dat <- db
     spc_col <- "presence"  # species presence/absence column is named "presence" in the dataset
     var_cols <-  names(db[,5:25]) 
     var_cols
     
     ## GENERALIZED LINEAR MODEL (GLM) // SENSITIVITY ####
     
     form_glm <- as.formula(paste(spc_col, "~", paste(selected_vars_glm, collapse = "+")))
     form_glm
     mod_glm <- glm(formula = form_glm, family = binomial, data = dat)
     summary(mod_glm)
     AUC_GLM<- AUC(mod_glm)$AUC
     
     
     ## GENERALIZED ADDITIVE MODEL (GAM) // SENSITIVITY ####

     form_gam <- as.formula(paste(spc_col, "~", paste0("s(",  selected_vars_gam, ")", collapse = "+")))  # GAM with smoothing splines ('s')
     form_gam
     mod_gam <- gam(formula = form_gam, family = binomial, data = dat)
     summary(mod_gam)
     AUC_GAM<- AUC(mod_gam)$AUC
     
     
     
     
     ## MAXIMUM ENTROPY (MAXENT) // SENSITIVITY ####

     
     mod_mxt <- maxnet(p = dat[ , spc_col], data = dat[ , selected_vars_maxt], f = maxnet.formula(dat[ , spc_col], dat[ , selected_vars_maxt])) 
     summary (mod_mxt)
     # Compute AUC
     # Get predictions from the model
     predicted_scores <- predict(mod_mxt, newdata = dat[ , selected_vars_maxt])
     roc_obj <- roc(dat[ , spc_col], predicted_scores)
     AUC_maxent <- auc(roc_obj)
     
     
     ## BAYESIAN ADDITIVE REGRESSION TREES (BART) // SENSITIVITY ####
    
     mod_bart <- bart(x.train = dat[ , varselect_bart], y.train = dat[ , spc_col], keeptrees = TRUE)
     
     # Extraer el valor de AUC desde el resumen
     match_auc <- regmatches(capture.output(summary(mod_bart)), regexpr("AUC = [0-9.]+", capture.output(summary(mod_bart))))
     # Extraer el número de la cadena coincidente
     auc_bart <- as.numeric(gsub("AUC = ", "", match_auc))
     
     ## Save AUC values ####
     Output_AUC[nrow(Output_AUC) + 1,] <- c(AUC_GLM, AUC_GAM, AUC_maxent, auc_bart)
     
     
     
     # Save models
     models <- list(mod_glm ,mod_gam, mod_mxt, mod_bart)
     names(models) <- c("glm", "gam", "mxt", "bart")
     pred_folder <- paste0( species, "/Sensitivity")
     dir.create(pred_folder, recursive = TRUE)
     saveRDS(models, file = paste0(pred_folder, "/models.rds"))
     # Load models:
     models <- readRDS(paste0(pred_folder, "/models.rds"))
     summary(models)
     
     
     # MAKE PREDICTIONS AND GET THE RESULTS FOR EACH TIME STEP OF THE MIS3 // SENSITIVITY ####
     # the next command creates a folder for output files:

     dir.create(paste0(pred_folder, "/Layers") )  # '../' goes up one level from the current working directory, so this creates the 'outputs' folder just outside the 'scripts' folder
     # Now two additional folders are created, one for the .TIF files and other one for the .JPEG files
     dir.create(paste0(pred_folder, "/Layers/TIF") )
     dir.create(paste0(pred_folder, "/Layers/JPEG") ) 
     
     path_TIF <- paste0(my_path, "/", pred_folder, "/Layers/TIF")
     path_JPEG <- paste0(my_path, "/", pred_folder, "/Layers/JPEG") 
     
     # Now, for visual representation purposes, we will include the coast line
     coast_file <- read.csv("coastline__120m.csv")
     # Convert coordinates from WKT to geometric objects
     coast_geom <- st_as_sfc(coast_file$the_geom, crs = 4326)
     # Create object sf
     coast_sf <- st_sf(coast_geom, data = coast_file)
     
     f <- list.files(path=path_env_f, pattern='rds$', full.names=TRUE) # create a list with all raster files to make predictions
     # clrs <- hcl.colors(n = 100, palette = "Lajolla")
     library(RColorBrewer)
     mi_paleta <- colorRampPalette(brewer.pal(11, "RdYlBu"))
     clrs <- rev(mi_paleta(100))
     

     
    ##### MAKE PREDICTIONS // SENSITIVITY####
     
     for (i in 1:length(f)) {
       raster <- readRDS(f[[i]])
       
       glm_p <- predict(stack(raster), mod_glm, type = "response")
       glm_fav <- Fav(pred = glm_p, sample.preval = prevalence(model = mod_glm))
       glm_fav <- as(glm_fav, "SpatRaster")
       
       gam_p <- predict(stack(raster), mod_gam, type = "response")
       gam_fav <- Fav(pred = gam_p, sample.preval = prevalence(model = mod_gam))
       gam_fav <- as(gam_fav, "SpatRaster")
       
       mxt_p <- rast(predict(stack(raster), mod_mxt, type = "cloglog"))  # raster format conversions needed because 'predict' for MX is not working correctly with 'terra' raster format
       mxt_fav <- as(mxt_p, "SpatRaster")
       bart_p <- rast(predict(mod_bart, x.layers = stack(raster)))
       bart_fav <- Fav(pred = bart_p, sample.preval = prevalence(model = mod_bart))
       bart_fav <- as(bart_fav, "SpatRaster")
       
       # Ensemble during the loop
       # read csv with cross-validation outcomes:
       
       AUC_GLM <- AUC_GLM
       AUC_GAM <- AUC_GAM
       AUC_MXT <- AUC_maxent
       AUC_BART <- auc_bart
       
       # Ensemble "response"
       glm_p <- as(glm_p, "SpatRaster")
       gam_p <- as(gam_p, "SpatRaster")
       
       crs(bart_p) <- crs(mxt_p)
       b <- brick(x = c(glm_p, gam_p, mxt_p, bart_p))
 
       
       ensemble_r <- weighted.mean(b, w= c(AUC_GLM, AUC_GAM, AUC_MXT, AUC_BART))

       plot(ensemble_r, col = clrs, range = c(0, 1), main = paste0("", (26 + i), "Ensemble.r k yrs BP"))
       plot(coastsCoarse, add=T, col="black")
       
       
       # Ensemble "favourability"
       glm_fav <- as(glm_fav, "SpatRaster")
       gam_fav <- as(gam_fav, "SpatRaster")
       
       crs(bart_fav) <- crs(mxt_p)
       
       # Ensure that you exclude models with AUC or BI values lower than 0.7:
       b <- brick(x = c(glm_fav, gam_fav, mxt_p, bart_fav))
       ensemble_f <- weighted.mean(b, w= c(AUC_GLM, AUC_GAM, AUC_MXT, AUC_BART))
       
       plot(ensemble_f, col = clrs, range = c(0, 1),  main = paste0("", (26 + i), "Ensemble.f k yrs BP"))
       plot(coastsCoarse, add=T, col="black")

       ## End ensemble
       # Save rasters as tif
       writeRaster(ensemble_f,paste0(path_TIF, "/", (26+ i), "f.tif"),options=c('TFW=YES'), overwrite=TRUE)

       dev.off()
       
       
     }
     
     
     ############# COMPUTE NICHE BREADTH // SENSITIVITY ########

     path_TIF
     
     # Get the .tif files
     f <- list.files(path_TIF,
                     pattern = "\\.tif$",
                     recursive = FALSE,
                     all.files = FALSE,
                     full.names = TRUE)
     # Select the favourability files
     f <- f[grep("f\\.tif$", f)]

     # Get the numbers
     file_numbers <- as.integer(gsub("[^0-9]", "", basename(f)))
     
     # Order the numbers
     order_index <- order(file_numbers, decreasing = FALSE)
     
     # Order the f list accordingly 
     f_sorted <- f[order_index]
     
     # Get and Store data
     resultados <- list()
     for (i in 1:length(f_sorted)) {
       archivo_raster <- f_sorted[i]
       
       # Carga el raster desde el archivo
       raster <- raster::raster(archivo_raster)
       plot(raster)
       
       ensemble_f <- raster
       # ensemble_f <- weighted.mean(b, w= c(AUC_GAM, AUC_MXT, AUC_BART))
       # plot(ensemble_f, range = c(0, 1), main = "ENSEMBLE f")
       ensemble_f<- as(ensemble_f, "SpatRaster")
       NicheBreadth_f<- raster.breadth(ensemble_f,verbose=FALSE) # Measures the spatial heterogeneity of the distribution of suitability scores
     
       
       Time <- as.integer(gsub("[^0-9]", "", basename(archivo_raster)))
       bloque <- 29
       bloque_actual <- ((i - 1) %/% bloque) + 1
       Experiment <- bloque_actual
       NichBrd_1_f <- NicheBreadth_f$B1
       NichBrd_2_f<- NicheBreadth_f$B2
       raster <- ensemble_f >= 0.5 # Pixels with favourability values =>0.5 = presence.
       plot(raster)
       absence <- sum(raster[] == 0, na.rm = TRUE)
       Favourable_area<- sum(raster[] == 1, na.rm = TRUE) / sum(raster[] == 0, na.rm = TRUE)
       Favourable_area_abs<- sum(raster[] == 1, na.rm = TRUE)
       Species <- species
       # area in km
       r <- raster(raster)
       polys <- rasterToPolygons(r, fun=function(x){x==1}, dissolve=TRUE)
       Area_km <- sum(areaPolygon(polys)) / 1E6  # en km2
       
  
       #Save outputs
       Output[nrow(Output) + 1,] <- c(Time, Experiment, NichBrd_1_f, NichBrd_2_f, 
                                      Favourable_area, Favourable_area_abs, Area_km, Species)

     }
     

     
   }
   print(Output)
   # save data
   write.xlsx(Output, file = paste0(pred_folder, "/Output_Sensitivity.xlsx") )
   write.xlsx(Output_AUC, file = paste0(pred_folder, "/Output_AUC_Sensitivity.xlsx") )
   
 }
 
 

 # Run the sensitivity test. It will need time. 
 Sensitivity_Age()
 


 
 
