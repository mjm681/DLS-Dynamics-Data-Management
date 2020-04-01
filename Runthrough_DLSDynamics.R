##### Creating Environment #####

setwd("~/M6 Project/Final Data/")

library(stringr)
library(ggplot2)
library(tidyverse)
library(stats)
library(factoextra)
library(ggpubr)
library(naniar)
library(stats)
library(factoextra)
library(dplyr)
library(purrr)
library(devtools)



rm(list = ls()) # Removes everything from the environemnt
# check they are removed by doing
ls()

##### Load Data #####

data_csv <- read.csv("Data.csv", stringsAsFactors = FALSE, na.strings = c("", "--", "NULL", "NA"))

load_data <- function(data) {
  # Set correct headers
  peak_column_names <- c("Meas Num","Item","Series","Time Stamp","Date","Time (s)","Measurement ID","Image","Intensity (Cnt/s)",
                         "Normalized Intensity (Cnt/s)","Number Acqs","Number Unmarked Acqs","Number Marked Acqs","% Acqs Unmarked",
                         "Acq Time (s)","DLS Temp (C)","Attenuation Level (%)","Laser Power (%)","Set Temp (C)","Ramp Rate (C/min)",
                         "Well","Row","Col","Configuration","Radius (nm)","Amplitude","Diffusion Coefficient (cm^2/s)",
                         "Diameter (nm)","Polydispersity (nm)","%PD","PD Index","Mw-R (kDa)","Baseline","SOS","Viscosity from Sample Rh (cP)",
                         "Lambda (1/s)","Sigma (1/s^2)","D10 (nm)","D50 (nm)","D90 (nm)","Span (D90 - D10)/D50","RMS Error",
                         "Number Unfiltered Peaks (I)","Number Filtered Peaks (I)","Pre-correction Internal Standard Rh (I) (nm)",
                         "Internal Standard Viscosity (I) (cP)","Mw-S (Da)","1/Mw-S (1/Da)","Particle Concentration (1/mL)",
                         "Particle Concentration Calculation","Particle Material","Particle Shape","Particle Core Real RI",
                         "Particle Core Imaginary RI","Sample","Concentration (mg/mL)","Mw-R Model","dn/dc (mL/g)",
                         "Conformation Model for SLS","A2 (mol mL/g^2)","Viscosity Calculation: Use Internal Standard",
                         "Internal Standard Rh (nm)","Internal Standard Rh Range Minimum (nm)","Internal Standard Rh Range Maximum (nm)",
                         "Solvent","Rfr Idx @ 589nm & 20C","Viscosity (cP)","Temp Model","Peak Num","Peak Radius (nm) (I)",
                         "Peak Diameter (nm) (I)","Peak Diffusion Coefficient (cm^2/s) (I)","Peak Mw-R (kDa) (I)","Peak %Pd (I)",
                         "Peak Pd (nm) (I)","Peak Pd Index (I)","Peak %Intensity (I)","Peak %Mass (I)","Peak %Number (I)",
                         "Membrane_Type", "Protein_Type", "Protein_Conc", "Lipid_Type", "Lipid_Conc", "Polymer_Type",
                         "Polymer_Percentage", "Buffer_Type", "Buffer_Conc", "pH", "Salt_Conc", "Polymer_Mr", "Latex_Nanospheres",
                         "Experiment_ID"
  )
  colnames(data) <- peak_column_names # Give dataframe correct column names
  # Replace -- with NA
  data <- replace_with_na_all(data, condition = ~.x == "--") 
  data <- replace_with_na_all(data, condition = ~.x == "NULL")
  data <- replace_with_na_all(data, condition = ~.x == "")
  # Make Item a combination of the conditions
  data$Conditions <- paste(data$Membrane_Type, data$Protein_Type,
                     data$Protein_Conc, data$Lipid_Type,
                     data$Lipid_Conc, data$Polymer_Type,
                     data$Polymer_Percentage, data$Buffer_Type,
                     data$Buffer_Conc, data$pH,
                     data$Salt_Conc, data$Polymer_Mr,
                     data$Latex_Nanospheres,
                     sep = " - ")
  # Change columns into correct data type
  data$Item <- as.factor(data$Item)
  data$`Peak Num` <- as.factor(data$`Peak Num`)
  data$`Peak Diameter (nm) (I)` <- as.numeric(data$`Peak Diameter (nm) (I)`)
  data$Membrane_Type <- as.factor(data$Membrane_Type)
  data$Protein_Type <- as.factor(data$Protein_Type)
  data$Lipid_Type <- as.factor(data$Lipid_Type)
  data$Lipid_Conc <- as.factor(data$Lipid_Conc)
  data$Polymer_Type <- as.factor(data$Polymer_Type)
  data$Buffer_Type <- as.factor(data$Buffer_Type)
  data$Latex_Nanospheres <- as.factor(data$Latex_Nanospheres)
  data$Conditions <- as.factor(data$Conditions)
  data$Experiment_ID <- as.factor(data$Experiment_ID)
  data$Protein_Conc <- as.numeric(data$Protein_Conc)
  data$Polymer_Percentage <- as.numeric(data$Polymer_Percentage)
  data$Buffer_Conc <- as.numeric(data$Buffer_Conc)
  data$pH <- as.numeric(data$pH)
  data$Salt_Conc <- as.numeric(data$Salt_Conc)
  data$Polymer_Mr <- as.numeric(data$Polymer_Mr)
  data$Image <- as.factor(data$Image)
  # Log transform the diameter and PD
  data$`Log Diameter` <- log10(data$`Peak Diameter (nm) (I)`)
  data$`Log Polydispersity` <- log10(data$`Peak Pd (nm) (I)`)
  return(data)
}

data <- load_data(data_csv)

Items_Conditions <- function(data) {
  items_conditions <- data[,c(93, 94, 2)]
  return(items_conditions)
}

items_conditions <- Items_Conditions(data = data)

sort_data <- function(data, sorting){
  # Sort order based on item
  attach(data)
  data <- data[order(sorting),]
  detach(data)
  return(data)
}

sorting <- data$Conditions

sorted_data <- sort_data(data = data, sorting = sorting)

##### Plot logD #####

log_diameter_boxplot <- function(data, xvar = Conditions, variable, target = T, original = F, indent = 1) {
  xvar <- enquo(xvar)
  variable <- enquo(variable)
  plot <- ggplot(data, aes(x = !!xvar, y = `Log Diameter`, col = !!variable)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Log Diameter for Single Peak") +
    labs(colour = "Variable of Interest")
  if ( target == TRUE ) {
    plot <- plot + geom_hline(yintercept = 1, col = "black") +
      annotate("text", label = "10nm", x = indent, y = 1, size = 3, colour = "black")
  }
  if ( original == TRUE ) {
    plot <- plot + geom_hline(yintercept = 2.3, col = "blue") +
      annotate("text", label = "200nm", x = indent, y = 2.3, size = 3, colour = "blue")
  }
  return(plot)
}

log_diameter_boxplot(sorted_data, xvar = pH, variable = Buffer_Type, target = T, original = T, indent = 2)

##### Split logD #####

split_peaks <- function(data, peak){
  split_peaks <- data %>% filter(`Peak Num` == peak)
  return(split_peaks)
}

split_peak_variable <- function(data, peak, variable) {
  variable <- enquo(variable) # Allows calling of variable as variable
  split_data <- split_peaks(data, peak)
  split_variable <- select(split_data, !!variable)
  return(unlist(split_variable)) # Return the variable as a vector
}

split_peak_xvar <- function(data, peak, xvar) {
  xvar <- enquo(xvar) # Allows calling of xvar as variable
  split_data <- split_peaks(data, peak)
  split_xvar <- select(split_data, !!xvar)
  return(unlist(split_xvar)) # Return the variable as a vector
}

split_peak_1 <- split_peaks(sorted_data, 1)
split_var <- split_peak_variable(sorted_data, 1, Buffer_Type)
split_xvar <- split_peak_xvar(sorted_data, 1, Conditions)

##### Plot Split #####

log_diameter_peak_boxplot <- function(data, peak, xvar, variable, target, original, indent) {
  variable <- enquo(variable) # Allows calling of variable as variable
  xvar <- enquo(xvar) # Allows calling of xvar as variable
  split_peaks <- split_peaks(data, peak)
  split_variable <- split_peak_variable(data, peak, !!variable)
  split_xvar <- split_peak_xvar(data, peak, !!xvar)
  log_diameter_boxplot(split_peaks, split_xvar, split_variable, target, original, indent)
}

log_diameter_peak_boxplot(data, 1, Experiment_ID, Buffer_Type, target = T, original = T, indent = 1)

##### Plot Logs #####

plot_logD_v_logPD <- function(data, variable) {
  variable <- enquo(variable)
  ggplot(data, aes(x = `Log Polydispersity`, y = `Log Diameter`, col = !!variable)) +
    geom_point() +
    labs(xlab = "Log Polydispersity", ylab = "Log Diameter", title = "Log Diameter vs Log Polydispersity")
}

plot_logD_v_logPD(data, Lipid_Conc)

plot_logD_v_logPD_peak <- function(data, peak, variable) {
  variable <- enquo(variable) # Allows calling of variable as variable
  split_peaks <- split_peaks(data, peak)
  split_variable <- split_peak_variable(data, peak, !!variable)
  plot_logD_v_logPD(split_peaks, split_variable)
}

plot_logD_v_logPD_peak(sorted_data, 2, Polymer_Percentage)

##### plot logPD #####

log_PD_boxplot <- function(data, xvar, variable, indent) {
  xvar <- enquo(xvar)
  variable <- enquo(variable)
  plot <- ggplot(data, aes(x = !!xvar, y = `Log Polydispersity`, col = !!variable)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Log Polydispersity for Single Peak") +
    labs(colour = "Variable of Interest")
  return(plot)
}

log_PD_boxplot(sorted_data, Polymer_Percentage, Buffer_Type, indent = 1)

##### Quartiles #####

quartile_extraction <- function (data, grouping) {
  fct_explicit_na(data$Buffer_Type, na_level = "(Missing)")
  p <- c(0, 0.25, 0.5, 0.75, 1)
  p_names <- map_chr(p, ~paste0(.x*100, "%"))
  p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>%
    set_names(nm = p_names)
  grouping <- enquo(grouping)
  both <- data %>%
    group_by(!!grouping) %>%
    summarize_at(vars(`Log Diameter`, `Log Polydispersity`), funs(!!!p_funs)) %>%
    select(!!grouping, contains("Log Diameter"), contains("Log Polydispersity"))
}
