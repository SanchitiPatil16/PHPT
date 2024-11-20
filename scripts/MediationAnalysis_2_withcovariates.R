# Mediation Analysis 
# With covariates
# Part2

rm(list=ls())


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(openxlsx)
library(mediation)
library(gt)
#LBRAU-0039,GACK-0041,MGETZ-0040

# Set working directory
setwd("/Users/amnahsiddiqa/Documents/Colabs/Hamid/coulmbiaData/Sep08")

data<-read.xlsx("/Users/amnahsiddiqa/Documents/Colabs/Hamid/coulmbiaData/Sep08/Columbia  data for Amnah 08 29 24.xlsx", check.names = FALSE)
metadata<-read.xlsx("/Users/amnahsiddiqa/Documents/Colabs/Hamid/coulmbiaData/NEJM/November2024/Patient Age and gender.xlsx", check.names = FALSE)
data<-left_join(metadata, data)
rm(metadata)
colnames(data)
med=data[,c("Patient.ID",
            "Age" ,
            "Gender" , 
             "Bifidobacterium.longum" ,
             "CD4+CD3+TNF+.T.cell.%",
             "CD4+CD3+TNF+.T.cell.MFI" ,
             "PBMC.Tnf.mRNA.(relative.expression)" ,
             "CD4+CD3+-Th17.%",
             "CD4+CD3+-Th17.MFI" ,
             "PBMC.Il17a.mRNA.(relative.expression)" ,
             "Tibia.Tb.vBMD","Tibia.Ct.vBMD","Tibia.BV/TV",
             "UD.Radius.Tb.vBMD", "UD.Radius.Ct.vBMD"  ,"UD.Radius.BV/TV"   ,
             "1/3.Distal.Radius.Tb.vBMD" ,"1/3.Distal.Radius.Ct.vBMD" ,"1/3.Distal.Radius.BV/TV" )]
med<-na.omit(med)
dim(med)
colnames(med)
data=med
rm(med)

set.seed(1190)

colnames(data)<-c("PID",
                  "BL", "Age" ,
                  "Gender" , 
                  "perc_tnfPTcell", "tnfPTcell_MFI", "tnf_mrna","perc_th17", "MFI_th17", "Il17a_mrna",
                  "TibiaTbvBMD", "TibiaCtvBMD", "TibiaBV_TV",
                  "UDRadiusTbvBMD", "UD.Radius.Ct.vBMD", "UDRadiusBV_TV",
                  "X1_3DistalRadiusTbvBMD", "X1_3.DistalRadiusCt_vBMD", "X1_3DistalRadiusBV_TV"
                  )


# Define mediators, outcomes, and exposure
mediators <- c("perc_th17", "MFI_th17", "Il17a_mrna", "perc_tnfPTcell", "tnfPTcell_MFI", "tnf_mrna")
exposure <- "BL"
outcomes_tibia <- c("TibiaTbvBMD", "TibiaCtvBMD", "TibiaBV_TV")
outcomes_radius <- c("UDRadiusTbvBMD", "UD.Radius.Ct.vBMD", "UDRadiusBV_TV")
outcomes_distal_radius <- c("X1_3DistalRadiusTbvBMD", "X1_3.DistalRadiusCt_vBMD", "X1_3DistalRadiusBV_TV")

# Function to perform mediation analysis for a given set of outcomes
perform_mediation_analysis <- function(outcomes) {
  results <- data.frame(
    Mediator = character(),
    `Beta a` = numeric(),
    `p-value a` = numeric(),
    `Beta b` = numeric(),
    `p-value b` = numeric(),
    `ACME Beta` = numeric(),
    `ACME p-value` = numeric(),
    `ADE Beta` = numeric(),
    `ADE p-value` = numeric(),
    Outcome = character()
  )
  
  for (mediator in mediators) {
    for (outcome in outcomes) {
      set.seed(1190)
      mediator_model <- lm(as.formula(paste(mediator, "~", exposure, "+Age+Gender")), data = data)
      set.seed(1190)
      outcome_model <- lm(as.formula(paste(outcome, "~", mediator, "+", exposure, "+Age+Gender")), data = data)
      set.seed(1190)
      mediation_result <- mediate(mediator_model, outcome_model, treat = exposure, mediator = mediator)
      
      # Extract beta coefficients and p-values
      beta_a <- summary(mediator_model)$coefficients[exposure, "Estimate"]
      p_value_a <- summary(mediator_model)$coefficients[exposure, "Pr(>|t|)"]
      beta_b <- summary(outcome_model)$coefficients[mediator, "Estimate"]
      p_value_b <- summary(outcome_model)$coefficients[mediator, "Pr(>|t|)"]
      beta_acme <- mediation_result$d0
      p_value_acme <- summary(mediation_result)$d0.p
      beta_ade <- mediation_result$z0
      p_value_ade <- summary(mediation_result)$z0.p
      
      # Append results to the data frame
      results <- rbind(results, data.frame(
        Mediator = mediator,
        `Beta a` = beta_a,
        `p-value a` = p_value_a,
        `Beta b` = beta_b,
        `p-value b` = p_value_b,
        `ACME Beta` = beta_acme,
        `ACME p-value` = p_value_acme,
        `ADE Beta` = beta_ade,
        `ADE p-value` = p_value_ade,
        Outcome = outcome
      ))
    }
  }
  
  return(results)
}

# Perform mediation analysis for each anatomical site
results_tibia <- perform_mediation_analysis(outcomes_tibia)
results_radius <- perform_mediation_analysis(outcomes_radius)
results_distal_radius <- perform_mediation_analysis(outcomes_distal_radius)


results<-rbind(results_distal_radius, results_radius, results_tibia)



head(results)
# Function to add stars based on p-value
add_stars <- function(p) {
  ifelse(p < 0.001, paste0(p, "***"),
         ifelse(p < 0.01, paste0(p, "**"),
                ifelse(p < 0.05, paste0(p, "*"), as.character(p))))
}

# Apply the function to the relevant p-value columns
results$p.value.a <- add_stars(results$p.value.a)
results$p.value.b <- add_stars(results$p.value.b)
results$ACME.p.value <- add_stars(results$ACME.p.value)
results$ADE.p.value<- add_stars(results$ADE.p.value)

# View the updated dataset
head(results)


desired_mediator_order <- c("tnf_mrna",
                            "perc_tnfPTcell" ,
                            "tnfPTcell_MFI",
                            "Il17a_mrna",
                            "perc_th17", 
                            "MFI_th17")      


mediator_mapping <- c(
  "tnf_mrna" = "PBMC Tnf mRNA",
  "perc_tnfPTcell" = "TNF+ T cells (% CD4+ T cells)",
  "tnfPTcell_MFI" = "TNF+ T cell (MFI)",
  "Il17a_mrna" = "PBMC Il17 mRNA",
  "perc_th17" = "Th17 cells (%CD4+ T cells)",
  "MFI_th17" = "Th17 cells (MFI)"
)

outcomes=unique(results$Outcome)
outcomes
# Generate formatted tables for each anatomical site and save as CSV
for (anoutcome in outcomes) {
  outcomes_sub <- results%>%filter(Outcome%in%anoutcome)
  print(outcomes_sub)
  # Convert Mediator to a factor with desired order
  outcomes_sub$Mediator <- factor(outcomes_sub$Mediator, levels = desired_mediator_order)
  
  # Reorder rows based on the ordered factor levels
  outcomes_sub <- outcomes_sub  %>%
    arrange(Mediator)
  
  outcomes_sub<-outcomes_sub %>%
    mutate(Mediator = mediator_mapping[Mediator])
  
  
  # Save the formatted table as CSV
  write.csv(outcomes_sub, paste0("/Users/amnahsiddiqa/Documents/Colabs/Hamid/coulmbiaData/NEJM/November2024/Mediation_results_covariates/", anoutcome, "_mediation_analysis.csv"), row.names = FALSE)
}

