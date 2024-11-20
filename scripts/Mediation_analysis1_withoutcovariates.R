# Mediation Analysis 
# Without covariates
# Part1

#clear env
rm(list=ls())
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(openxlsx)
library(mediation)
library(gt)

setwd("/Users/amnahsiddiqa/Documents/Colabs/Hamid/coulmbiaData/Sep08")

#read Metadata
data<-read.xlsx("/Users/amnahsiddiqa/Documents/Colabs/Hamid/coulmbiaData/Sep08/Columbia  data for Amnah 08 29 24.xlsx", check.names = FALSE)


# Mediation Analysis  -----------------------------------------------------
colnames(data)
med=data[,c("Patient.ID",
            "Bifidobacterium.longum" ,
            "Spine.BMD",
            "1/3.Radius.BMD",
            "CD4+CD3+TNF+.T.cell.%",
            "CD4+CD3+TNF+.T.cell.MFI" ,
            "PBMC.Tnf.mRNA.(relative.expression)" ,
            "CD4+CD3+-Th17.%",
            "CD4+CD3+-Th17.MFI" ,
            "PBMC.Il17a.mRNA.(relative.expression)"  )]

colnames(med)

colnames(med)<-c("Patient.ID",
                 "BL" ,
                 "SpBMD",
                 "rBMD",
                 "perc_tnfPTcell",
                 "tnfPTcell_MFI" ,
                 "tnf_mrna" ,
                 "perc_th17",
                 "MFI_th17" ,
                 "X1l17_mrna"  )


colnames(med)

data=med
rm(med)
data<-na.omit(data)

# Define mediators and outcomes
mediators <- c("perc_tnfPTcell", "tnfPTcell_MFI", "tnf_mrna", "perc_th17", "MFI_th17", "X1l17_mrna")
outcomes <- c("rBMD", "SpBMD")
# Initialize the results data frame with additional columns for beta values
results <- data.frame(
  Mediator = character(),
  `p-value a` = numeric(),
  `Beta a` = numeric(),
  `p-value b (1/3 Radius BMD)` = numeric(),
  `Beta b (1/3 Radius BMD)` = numeric(),
  `ACME p-value (1/3 Radius BMD)` = numeric(),
  `ACME Beta (1/3 Radius BMD)` = numeric(),
  `ADE p-value (1/3 Radius BMD)` = numeric(),
  `ADE Beta (1/3 Radius BMD)` =numeric(),
  `p-value b (Spine BMD)` = numeric(),
  `Beta b (Spine BMD)` = numeric(),
  `ACME p-value (Spine BMD)` = numeric(),
  `ACME Beta (Spine BMD)` = numeric(),
  `ADE p-value (Spine BMD)` = numeric(),
  `ADE Beta (Spine BMD)` =numeric()
  
)

# Perform mediation analysis for each mediator and outcome
for (mediator in mediators) {
  for (outcome in outcomes) {
    
    print(mediator)
    print(outcome)
    

    # Fit the mediator and outcome models
    mediator_model <- lm(as.formula(paste(mediator, "~ BL")), data = data)
    outcome_model <- lm(as.formula(paste(outcome, "~", mediator, "+ BL")), data = data)
    mediator_model
    outcome_model
    # Perform mediation analysis
    mediation_result <- mediate(mediator_model, outcome_model, treat = "BL", mediator = mediator)
    
    # Extract p-values and beta coefficients
    p_value_a <- summary(mediator_model)$coefficients["BL", "Pr(>|t|)"]
    beta_a <- summary(mediator_model)$coefficients["BL", "Estimate"]

    
    p_value_b <- summary(outcome_model)$coefficients[mediator, "Pr(>|t|)"]
    beta_b <- summary(outcome_model)$coefficients[mediator, "Estimate"]

 
    p_value_acme <- summary(mediation_result)$d0.p
    beta_acme <- mediation_result$d0
    
    p_value_ade <- summary(mediation_result)$z0.p
    beta_ade <- mediation_result$z0
    
    # Add to results data frame based on outcome
    if (outcome == "rBMD") {
      results <- rbind(results, data.frame(
        Mediator = mediator,
        `p-value a` = p_value_a,
        `Beta a` = beta_a,
        `p-value b (1/3 Radius BMD)` = p_value_b,
        `Beta b (1/3 Radius BMD)` = beta_b,
        `ACME p-value (1/3 Radius BMD)` = p_value_acme,
        `ACME Beta (1/3 Radius BMD)` = beta_acme,
        `ADE p-value (1/3 Radius BMD)` = p_value_ade,
        `ADE Beta (1/3 Radius BMD)` = beta_ade,
        `p-value b (Spine BMD)` = NA,
        `Beta b (Spine BMD)` = NA,
        `ACME p-value (Spine BMD)` = NA,
        `ACME Beta (Spine BMD)` = NA,
        `ADE p-value (Spine BMD)` = NA,
        `ADE Beta (Spine BMD)` = NA
      ))
    } else {
      results <- rbind(results, data.frame(
        Mediator = mediator,
        `p-value a` = p_value_a,
        `Beta a` = beta_a,
        `p-value b (1/3 Radius BMD)` = NA,
        `Beta b (1/3 Radius BMD)` = NA,
        `ACME p-value (1/3 Radius BMD)` = NA,
        `ACME Beta (1/3 Radius BMD)` = NA,
        `ADE p-value (1/3 Radius BMD)` =NA,
        `ADE Beta (1/3 Radius BMD)` = NA,
        `p-value b (Spine BMD)` = p_value_b,
        `Beta b (Spine BMD)` = beta_b,
        `ACME p-value (Spine BMD)` = p_value_acme,
        `ACME Beta (Spine BMD)` = beta_acme,
        `ADE p-value (Spine BMD)` = p_value_ade,
        `ADE Beta (Spine BMD)` = beta_ade
      ))
    }
  }
}
# Define the function to add stars based on p-value
# Define the function to add stars based on p-value
add_stars <- function(p) {
  ifelse(
    is.na(p), NA,
    ifelse(p < 0.001, paste0(formatC(p, format = "e", digits = 3), "***"),
           ifelse(p < 0.01, paste0(formatC(p, format = "e", digits = 3), "**"),
                  ifelse(p < 0.05, paste0(formatC(p, format = "e", digits = 3), "*"),
                         formatC(p, format = "e", digits = 3)))))
}



# Format the table for ppt ------------------------------------------------
# Define the desired order for the mediators
desired_order <- c(
  "perc_tnfPTcell", 
  "tnfPTcell_MFI", 
  "tnf_mrna", 
  "perc_th17", 
  "MFI_th17", 
  "X1l17_mrna"
)

# Combine rows by filling NAs with non-NA values for each Mediator
results_combined <- results %>%
  group_by(Mediator) %>%
  summarise(across(everything(), ~ first(na.omit(.), default = NA))) %>%
  ungroup()


# Reorder the mediators based on the desired order
results_combined <- results_combined %>%
  mutate(Mediator = factor(Mediator, levels = desired_order)) %>%
  arrange(Mediator)

# Apply `add_stars` to p-value columns and format other numeric columns
results_combined <- results_combined %>%
  mutate(across(contains("p.value"), add_stars)) %>%
  mutate(across(where(is.numeric), ~ as.numeric(format(., scientific = FALSE, trim = TRUE))))

mediator_mapping <- c(
  "tnf_mrna" = "PBMC Tnf mRNA",
  "perc_tnfPTcell" = "TNF+ T cells (% CD4+ T cells)",
  "tnfPTcell_MFI" = "TNF+ T cell (MFI)",
  "Il17a_mrna" = "PBMC Il17 mRNA",
  "perc_th17" = "Th17 cells (%CD4+ T cells)",
  "MFI_th17" = "Th17 cells (MFI)"
)
results_combined <-results_combined  %>%
  mutate(Mediator = mediator_mapping[Mediator])

# Correct mapping by converting Mediator to a character vector first
results_combined <- results_combined %>%
  mutate(Mediator = recode(as.character(Mediator), !!!mediator_mapping))


# Save the final table to Excel
write.table(results_combined, "/Users/amnahsiddiqa/Documents/Colabs/Hamid/coulmbiaData/NEJM/November2024/Mediation_results_withoutcovariates/mediation_analysis_1p3_spineradius.csv", sep=",", row.names=FALSE)

