

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
  data <- naniar::replace_with_na_all(data, condition = ~.x == "--")
  data <- naniar::replace_with_na_all(data, condition = ~.x == "NULL")
  data <- naniar::replace_with_na_all(data, condition = ~.x == "")
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
  data$`Log % Intensity` <- log10(data$`Peak %Intensity (I)`)
  return(data)
}

Items_Conditions <- function(data) {
  items_conditions <- data[,c(93, 94, 2)]
  return(items_conditions)
}

sort_data <- function(data, sorting){
  # Sort order based on item
  attach(data)
  data <- data[order(sorting),]
  detach(data)
  return(data)
}


log_diameter_boxplot <- function(data, xvar = Conditions, variable, target = T, original = F, indent = 1) {
  xvar <- dplyr::enquo(xvar)
  variable <- dplyr::enquo(variable)
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = !!xvar, y = `Log Diameter`, col = !!variable)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(title = "Log Diameter") +
    ggplot2::labs(colour = "Variable of Interest")
  if ( target == TRUE ) {
    plot <- plot + ggplot2::geom_hline(yintercept = 1, col = "black") +
      ggplot2::annotate("text", label = "10nm", x = indent, y = 1, size = 3, colour = "black")
  }
  if ( original == TRUE ) {
    plot <- plot + ggplot2::geom_hline(yintercept = 2.3, col = "blue") +
      ggplot2::annotate("text", label = "200nm", x = indent, y = 2.3, size = 3, colour = "blue")
  }
  return(plot)
}


split_peaks <- function(data, peak){
  if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr")
  library(dplyr)
  split_peaks <- data  %>% filter(`Peak Num` == peak)
  return(split_peaks)
}

split_peak_variable <- function(data, peak, variable) {
  variable <- dplyr::enquo(variable) # Allows calling of variable as variable
  split_data <- split_peaks(data, peak)
  split_variable <- select(split_data, !!variable)
  return(unlist(split_variable)) # Return the variable as a vector
}

split_peak_xvar <- function(data, peak, xvar) {
  xvar <- dplyr::enquo(xvar) # Allows calling of xvar as variable
  split_data <- split_peaks(data, peak)
  split_xvar <- select(split_data, !!xvar)
  return(unlist(split_xvar)) # Return the variable as a vector
}


log_diameter_peak_boxplot <- function(data, peak, xvar, variable, target, original, indent) {
  variable <- dplyr::enquo(variable) # Allows calling of variable as variable
  xvar <- dplyr::enquo(xvar) # Allows calling of xvar as variable
  split_peaks <- split_peaks(data, peak)
  split_variable <- split_peak_variable(data, peak, !!variable)
  split_xvar <- split_peak_xvar(data, peak, !!xvar)
  log_diameter_boxplot(split_peaks, split_xvar, split_variable, target, original, indent)
}


plot_logD_v_logPD <- function(data, variable) {
  variable <- dplyr::enquo(variable)
  ggplot2::ggplot(data, ggplot2::aes(x = `Log Polydispersity`, y = `Log Diameter`, col = !!variable)) +
    ggplot2::geom_point() +
    ggplot2::labs(xlab = "Log Polydispersity", ylab = "Log Diameter", title = "Log Diameter vs Log Polydispersity")
}


plot_logD_v_logPD_peak <- function(data, peak, variable) {
  variable <- dplyr::enquo(variable) # Allows calling of variable as variable
  split_peaks <- split_peaks(data, peak)
  split_variable <- split_peak_variable(data, peak, !!variable)
  plot_logD_v_logPD(split_peaks, split_variable)
}


log_PD_boxplot <- function(data, xvar, variable, indent) {
  xvar <- dplyr::enquo(xvar)
  variable <- dplyr::enquo(variable)
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = !!xvar, y = `Log Polydispersity`, col = !!variable)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(title = "Log Polydispersity for Single Peak") +
    ggplot2::labs(colour = "Variable of Interest")
  return(plot)
}

plot_logD_v_logintensity <- function(data, variable) {
  variable <- dplyr::enquo(variable)
  ggplot2::ggplot(data, ggplot2::aes(x = `Log Diameter`, y = `Log % Intensity`, col = !!variable)) +
    ggplot2::geom_point() +
    ggplot2::labs(xlab = "% Intensity", ylab = "Log Diameter", title = "Log Diameter vs Log % Intensity")
}


plot_logD_v_logintensity_peak <- function(data, peak, variable) {
  variable <- dplyr::enquo(variable) # Allows calling of variable as variable
  split_peaks <- split_peaks(data, peak)
  split_variable <- split_peak_variable(data, peak, !!variable)
  plot_logD_v_logintensity(split_peaks, split_variable)
}
