library(brms)
library(tidyverse)

setwd('/Users/qifanwu/Documents/research/huron/code')

# proprocessing
df <- read_csv("../data/processed/WQ_BE_CLM_By_Month.csv")

colnames(df)[1] <- "ID"

df1 <- df[, c("ID", "Site_ID", "TP", "TSS", "eCol", "Cond", "DO", "Discharge")]
df1 <- df1 %>%
  remove_rownames() %>%
  column_to_rownames(var = 'ID')

variables <- c('tmean','ppt',
              'RoadDensity', 'LaneDensity', 
               # 'GravelDensity', 'AsphaltDensity', 'ConcreteDensity', # delete road pavement variables
               'StateRdDensity', 'CountyRdDensity', 'CityRdDensity', 'ParcelCtDensity',
               'ParcelAvgArea', 'AgPercent', 'CommercialPercent', 'GreenPercent',
               'IndustryPercent', 'ServicePercent', 'ResidentialPercent',
               'VaccantPercent', 'OpenDevPercent', 'LowDevPercent', 'MedDevPercent',
               'HighDevPercent', 'MedianFloors', 'MedianBuiltYr', 'RangeBuiltYr',
               'FootprintDensity', 'UnitDensity', 'SoilBPercent', 'SoilCPercent',
               'SoilDPercent', 'SoilDepth', 'Slope')

df2 <- df[, c("ID", variables)]
df2 <- df2 %>%
  remove_rownames() %>%
  column_to_rownames(var = 'ID')
df2 <- scale(df2)


d <- merge(df1, df2, by="row.names")

# variables selected from knockoff
ecol_v <- c("tmean","ppt",
            "LowDevPercent","SoilDPercent","AgPercent","MedDevPercent","ParcelCtDensity",
            "MedianFloors","StateRdDensity","ResidentialPercent","MedianBuiltYr",
            "ServicePercent","GreenPercent","ParcelAvgArea","UnitDensity","HighDevPercent","VaccantPercent",
            "SoilCPercent","SoilBPercent","CommercialPercent","CountyRdDensity","IndustryPercent","LaneDensity"
)

tp_v <- c("tmean","ppt",
          "FootprintDensity","SoilDPercent","VaccantPercent",
          "SoilDepth","GreenPercent","ParcelCtDensity"
)

tss_v <- c("tmean","ppt",
           "HighDevPercent","MedianBuiltYr","ParcelCtDensity","ServicePercent","ParcelAvgArea",
           "SoilDepth","VaccantPercent","GreenPercent","SoilCPercent","IndustryPercent","SoilDPercent",
           "RoadDensity","CityRdDensity","OpenDevPercent","FootprintDensity"
)

discharge_v <- c("tmean","ppt",
                 "MedianBuiltYr","ParcelAvgArea","Slope","GreenPercent",
                 "RoadDensity","SoilCPercent","OpenDevPercent","SoilBPercent","MedianFloors",
                 "UnitDensity","ParcelCtDensity","CountyRdDensity","IndustryPercent","VaccantPercent",
                 "CommercialPercent","CityRdDensity","FootprintDensity")
                 
# =========
# mixed linear model
# library('lme4')
# lmer(TSS ~ CountyRdDensity + ParcelAvgArea + CommercialPercent + ServicePercent + VaccantPercent + ppt + ppt*SoilCPercent + ppt*CityRdDensity
#      + (1 | Site_ID), data=d)
# =======

# ---------------------------
# run bayessian
# tp
bayesian_mixed_tp_gaussian = brm(
  TP ~ FootprintDensity + SoilDPercent + VaccantPercent + SoilDepth + GreenPercent + ParcelCtDensity + ppt + (1 | Site_ID),
  data = d, family = gaussian(), prior = set_prior("normal(-1,2000)")
)
# save the summary
sink("../output/bayesian/tp/bayesian_tp_gaussian.txt")
summary(bayesian_mixed_tp_gaussian)
sink()
# plot the result
png("../output/bayesian/tp/bayesian_mixed_tp_plot_%d.png", width = 800, height = 600)
plot(bayesian_mixed_tp_gaussian, ask = FALSE)
dev.off()

tp_bt_v <- c("FootprintDensity", "VaccantPercent", "GreenPercent", "ParcelCtDensity")
# tp_bt_v <- c("FootprintDensity")

# Loop through each variable to create, summarize, and plot each model
for (var in tp_bt_v) {
  # Create a dynamic model name
  model_name <- paste0("bayesian_tp_", var)
  formula <- as.formula(paste(
    "TP ~ FootprintDensity + SoilDPercent + VaccantPercent + SoilDepth + GreenPercent + ParcelCtDensity + ppt +",
    paste0(var, " * ppt + (1 | Site_ID)")
  ))
  
  model <- brm(
    formula = formula,
    data = d,
    family = gaussian(),
    prior = set_prior("normal(-1,2000)")
  )
  
  assign(model_name, model)  # Save the model in the environment
  
  # Save the summary to a file
  summary_file <- paste0("../output/bayesian/tp/", model_name, ".txt")
  summary_text <- capture.output(summary(model))
  writeLines(summary_text, summary_file)
  print(summary_text)
  
  plot_file <- paste0("../output/bayesian/tp/", model_name, "_plot_%d.png")
  png(plot_file, width = 800, height = 600)
  plot(get(model_name), ask = FALSE)
  dev.off()
}

# ---------------------------
# RUN bayesian model for TSS
bayesian_mixed_tss_gaussian = brm(
  TSS ~ HighDevPercent + MedianBuiltYr + ParcelCtDensity + ServicePercent + ParcelAvgArea + SoilDepth + VaccantPercent + GreenPercent + 
    SoilCPercent + IndustryPercent + SoilDPercent + RoadDensity + CityRdDensity + OpenDevPercent + FootprintDensity + ppt + (1 | Site_ID),
  data = d, family = gaussian(), prior = set_prior("normal(-1,2000)")
)

# save the summary
sink("../output/bayesian/tss/bayesian_tss_gaussian.txt")
summary(bayesian_mixed_tss_gaussian)
sink()
# plot the result
png("../output/bayesian/tss/bayesian_mixed_tss_plot_%d.png", width = 800, height = 600)
plot(bayesian_mixed_tss_gaussian, ask = FALSE)
dev.off()

tss_bt_v <- c("HighDevPercent", "ServicePercent", "VaccantPercent", "GreenPercent", "IndustryPercent", "RoadDensity", "CityRdDensity", "OpenDevPercent")

# Loop through each variable to create, summarize, and plot each model
for (var in tss_bt_v) {
  # Create a dynamic model name
  model_name <- paste0("bayesian_tss_", var)
  formula <- as.formula(paste(
    "TSS ~ HighDevPercent + MedianBuiltYr + ParcelCtDensity + ServicePercent + ParcelAvgArea + SoilDepth + VaccantPercent + GreenPercent + 
    SoilCPercent + IndustryPercent + SoilDPercent + RoadDensity + CityRdDensity + OpenDevPercent + FootprintDensity + ppt +",
    paste0(var, " * ppt + (1 | Site_ID)")
  ))
  
  model <- brm(
    formula = formula,
    data = d,
    family = gaussian(),
    prior = set_prior("normal(-1,2000)")
  )
  
  assign(model_name, model)  # Save the model in the environment
  
  # Save the summary to a file
  summary_file <- paste0("../output/bayesian/tss/", model_name, ".txt")
  summary_text <- capture.output(summary(model))
  writeLines(summary_text, summary_file)
  print(summary_text)
  
  plot_file <- paste0("../output/bayesian/tss/", model_name, "_plot_%d.png")
  png(plot_file, width = 800, height = 600)
  plot(get(model_name), ask = FALSE)
  dev.off()
}

# ---------------------------
# RUN bayesian model for eCol
bayesian_mixed_eCol_gaussian = brm(
  eCol ~ LowDevPercent + SoilDPercent + AgPercent + MedDevPercent + ParcelCtDensity + MedianFloors + StateRdDensity + ResidentialPercent + MedianBuiltYr + ServicePercent + GreenPercent + ParcelAvgArea + UnitDensity + HighDevPercent + VaccantPercent + SoilCPercent + SoilBPercent + 
    CommercialPercent + CountyRdDensity + IndustryPercent + LaneDensity + ppt + (1 | Site_ID),
  data = d, family = gaussian(), prior = set_prior("normal(-1,2000)")
)

# save the summary
sink("../output/bayesian/eCol/bayesian_eCol_gaussian.txt")
summary(bayesian_mixed_eCol_gaussian)
sink()
# plot the result
png("../output/bayesian/eCol/bayesian_mixed_eCol_plot_%d.png", width = 800, height = 600)
plot(bayesian_mixed_eCol_gaussian, ask = FALSE)
dev.off()

eCol_bt_v <- c("LowDevPercent", "AgPercent", "MedDevPercent", "StateRdDensity", "ResidentialPercent", "ServicePercent", "GreenPercent", "HighDevPercent", "VaccantPercent", "CommercialPercent", "CountyRdDensity", "IndustryPercent", "LaneDensity")

# Loop through each variable to create, summarize, and plot each model
for (var in eCol_bt_v) {
  # Create a dynamic model name
  model_name <- paste0("bayesian_eCol_", var)
  formula <- as.formula(paste(
    "eCol ~ LowDevPercent + SoilDPercent + AgPercent + MedDevPercent + ParcelCtDensity + MedianFloors + StateRdDensity + ResidentialPercent + MedianBuiltYr + ServicePercent + GreenPercent + ParcelAvgArea + UnitDensity + HighDevPercent + VaccantPercent + SoilCPercent + SoilBPercent + 
    CommercialPercent + CountyRdDensity + IndustryPercent + LaneDensity + ppt +",
    paste0(var, " * ppt + (1 | Site_ID)")
  ))
  
  model <- brm(
    formula = formula,
    data = d,
    family = gaussian(),
    prior = set_prior("normal(-1,2000)")
  )
  
  assign(model_name, model)  # Save the model in the environment
  
  # Save the summary to a file
  summary_file <- paste0("../output/bayesian/eCol/", model_name, ".txt")
  summary_text <- capture.output(summary(model))
  writeLines(summary_text, summary_file)
  
  plot_file <- paste0("../output/bayesian/eCol/", model_name, "_plot_%d.png")
  png(plot_file, width = 800, height = 600)
  plot(get(model_name), ask = FALSE)
  dev.off()
}

# ---------------------------
# RUN bayesian model for discharge
discharge_v <- c("MedianBuiltYr","ParcelAvgArea","Slope","GreenPercent",
                 "RoadDensity","SoilCPercent","OpenDevPercent","SoilBPercent","MedianFloors",
                 "UnitDensity","ParcelCtDensity","CountyRdDensity","IndustryPercent","VaccantPercent",
                 "CommercialPercent","CityRdDensity","FootprintDensity")
discharge_bt_v <- c("GreenPercent", "RoadDensity", "OpenDevPercent", "CountyRdDensity", "IndustryPercent", "VaccantPercent", "CommercialPercent", "CityRdDensity")


wq_lists <- list(Discharge=discharge_v)
wq_bt_list <- list(Discharge=discharge_bt_v)

for (wq_i in names(wq_lists)) {
  # with individual variables
  variables <- wq_lists[[wq_i]]
  formula_string <- paste(wq_i, "~", paste(variables, collapse = " + "), "+ ppt + (1 | Site_ID)")
  formula <- as.formula(formula_string)
  
  model_name <- paste0("bayesian_mixed_", tolower(wq_i))
  
  model <- brm(
    formula = formula,
    data = d,
    family = gaussian(),
    prior = set_prior("normal(-1,2000)")
  )
  
  summary_file <- paste0("../output/bayesian/", wq_i, "/", model_name, "_summary.txt")
  summary_text <- capture.output(summary(model))
  writeLines(summary_text, summary_file)
  
  plot_file <- paste0("../output/bayesian/", wq_i, "/", model_name, "_plot_%d.png")
  png(plot_file, width = 800, height = 600)
  plot(model, ask = FALSE)
  dev.off()
  
  # with one interaction at a time
  wq_bt_v <- wq_bt_list[[wq_i]]
  
  # Loop through each variable to create, summarize, and plot each model
  for (var in wq_bt_v) {
    # Create a dynamic model name
    model_name <- paste0("bayesian_", wq_i, "_", var)
    formula_string <- paste(wq_i, "~", paste(variables, collapse = " + "), "+ ppt +",  paste0(var, " * ppt + (1 | Site_ID)"))
    
    formula <- as.formula(formula_string)
    
    model <- brm(
      formula = formula,
      data = d,
      family = gaussian(),
      prior = set_prior("normal(-1,2000)")
    )
    
    assign(model_name, model)  # Save the model in the environment
    
    # Save the summary to a file
    summary_file <- paste0("../output/bayesian/", wq_i, "/", model_name, ".txt")
    summary_text <- capture.output(summary(model))
    writeLines(summary_text, summary_file)
    
    plot_file <- paste0("../output/bayesian/", wq_i, "/", model_name, "_plot_%d.png")
    png(plot_file, width = 800, height = 600)
    plot(get(model_name), ask = FALSE)
    dev.off()
  }
}


