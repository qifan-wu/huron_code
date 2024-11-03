# if doMC is not installed, please run
install.packages("doMC")

library(tidyverse)
library(dplyr)
library(knockoff)

#set working directory
setwd('/Users/qifanwu/Documents/research/huron/code')

# read dataset
df <- read.csv("../data/processed/WQ_BE_CLM_By_Month.csv")
basic_vars <- c('Site_ID', 'year', 'month')

water_quality_vars <- c('Cond', 'DO', 'Discharge', 'TP', 'TSS','eCol')
# select variables
variables <- c('RoadDensity', 'LaneDensity', 
               # 'GravelDensity', 'AsphaltDensity', 'ConcreteDensity', # delete road pavement variables
               'StateRdDensity', 'CountyRdDensity', 'CityRdDensity', 'ParcelCtDensity',
               'ParcelAvgArea', 'AgPercent', 'CommercialPercent', 'GreenPercent',
               'IndustryPercent', 'ServicePercent', 'ResidentialPercent',
               'VaccantPercent', 'OpenDevPercent', 'LowDevPercent', 'MedDevPercent',
               'HighDevPercent', 'MedianFloors', 'MedianBuiltYr', 'RangeBuiltYr',
               'FootprintDensity', 'UnitDensity', 'SoilBPercent', 'SoilCPercent',
               'SoilDPercent', 'SoilDepth', 'Slope')

df <- df[,c(basic_vars, water_quality_vars, variables)] 

# setup variables and response
pos_start <- which(names(df) == 'RoadDensity')
pos_cols <- seq.int(pos_start, ncol(df))

# normalize variables
for (i in (pos_start:max(pos_cols))){
  df[i] <- scale(df[i]) #(df[i] - mean(df[i])) / sd(df[i])
}

# ==== RUN START HERE FOR DIFFERENT WATER QUALITY INDICATORS === 
# variable matrix
X_df <- as.matrix(df[,pos_cols])

# response dataframe
dependent <- 'DO' # dependent variable: 'TP', 'TSS','eCol', 'Cond', 'DO', 'Discharge'
Y_df <- df[,c(dependent)]
Y_df <- as.data.frame(Y_df)

# making a column all positive by adding a constant
min_value <- min(Y_df[1], na.rm = TRUE)
if (min_value <= 0) {
  Y_df[1] <- Y_df[1] + abs(min_value) + 0.0001
}

# drop the rows where Y is na
combined_df <- cbind(Y_df, X_df)
clean_combined_df <- na.omit(combined_df)
Y_df <- clean_combined_df[, 1, drop = FALSE]
X_df <- clean_combined_df[, -1]  

# setup selection function
knockoff_ <- function (X, y, q) {
  # Log-transform the drug resistance measurements.
  y <- log(y)
  fdr <- q
  # Run the knockoff filter.
  result = knockoff.filter(X, y, fdr=fdr)
  # oudfut result
  list(result)
}

# a function for extracting variable names and lasso importance
knockoff_extraction <- function(result){
  # names and importances of selections
  selected_variable <- names(result[["Y_df"]][[1]][["selected"]])
  variable_importance <- result[["Y_df"]][[1]][["statistic"]]
  lasso_importance <- c()
  for (i in (1:length(selected_variable))){
    index <- result[["Y_df"]][[1]][["selected"]][[i]]
    lasso_importance <- c(lasso_importance, variable_importance[index]) 
  }
  # return variable names and lasso importances
  list(selected_variable, lasso_importance)
}


####### MANUALLY TEST FDR ###########
set.seed(123)
loop_time <- 500
fdr <- 0.2
res_list <- list() ##
variables_list <- c()
for (j in (1:loop_time)){
  print(paste("Running", j, "/", loop_time))
  results <- lapply(Y_df, function(y) knockoff_(X_df, y, fdr))
  if(length(results[["Y_df"]][[1]][["selected"]]) > 0){
    df_result <- knockoff_extraction(results)
    variables_list <- c(variables_list, df_result[[1]])
    
    res_list <- c(res_list, df_result) ##

  }
  else{
    variables_list <- variables_list
  }
}

selected_vars <- table(variables_list)
selected_vars_df <- as.data.frame(selected_vars)
selected_vars_df <- selected_vars_df[order(selected_vars_df$Freq, decreasing = TRUE), ]

write.csv(selected_vars_df, paste0("../output/knockoff/selected_vars_",dependent,".csv"), row.names = FALSE)

########## DO NOT RUN ##########
#====== Result, Do Not Run ======
# (results are the same given same seed)
# -------- random seed = 123
# fdr = 0.1, No result with selected variable in 200 loops

# fdr = 0.15, 8 results with selected variables in 200 loops
# one of the results: "GravelDensity""GreenPercent""VaccantPercent""LowDevPercent""MedianFloors""FootprintDensity""UnitDensity""SoilBPercent""SoilDPercent"    
# with lasso importance: 0.062277768 0.050077756 0.052674780 0.042914930 0.099367614 0.006998181 0.007847308 0.039388260 0.036984879
# the unique variable names are: "GravelDensity"    "AsphaltDensity"   "GreenPercent"     "VaccantPercent"   "MedianFloors"    
#                               "FootprintDensity" "UnitDensity"      "SoilBPercent"     "LowDevPercent"    "SoilDPercent"    

# fdr = 0.2, 14 results with selected variables in 200 loops
# one of the results: "GravelDensity""GreenPercent""VaccantPercent""LowDevPercent""MedianFloors""FootprintDensity" "UnitDensity""SoilBPercent""SoilDPercent"    
# with lasso importance: 0.07905548 0.05536728 0.06811788 0.04120136 0.09367562 0.01962384 0.02681137 0.03699014 0.03320929
# the unique variable names are: "GravelDensity"    "GreenPercent"     "VaccantPercent"   "LowDevPercent"    "MedianFloors" "FootprintDensity" "UnitDensity" "SoilBPercent""SoilDPercent"     "AsphaltDensity"  

# fdr = 0.25, 34 results with selected variables in 200 loops
# one of the results: "GreenPercent" "MedianFloors" "SoilBPercent" "SoilDPercent"
#  with lasso importance: 0.01231887 0.09068350 0.02873076 0.02625675
# the unique variable names are: "GravelDensity""GreenPercent""VaccantPercent""LowDevPercent""MedianFloors""FootprintDensity" "UnitDensity" "SoilDPercent"  "SoilBPercent" "AsphaltDensity" "ServicePercent"   

# -------- random seed = 456
# fdr = 0.2, 16 results with selected variables in 200 loops
# one of the results: "GravelDensity""GreenPercent""LowDevPercent""MedianFloors""FootprintDensity"
# with lasso importance: 0.10880786 0.05244399 0.06288588 0.09884661 0.05089595
# another result "GravelDensity""GreenPercent""VaccantPercent""LowDevPercent""MedianFloors""FootprintDensity" "UnitDensity""SoilBPercent""SoilDPercent" 
# with lasso importance: 0.06535814 0.02312273 0.05295628 0.03766965 0.09833257 0.01612993 0.01716026 0.04197344 0.03288625
# the unique variable names are: "GravelDensity""GreenPercent""VaccantPercent""MedianFloors""UnitDensity""SoilDPercent""LowDevPercent""FootprintDensity" "SoilBPercent"    

# -------- random seed = 789
# fdr = 0.2, 11 results with selected variables in 200 loops
# one of the results: "GravelDensity""GreenPercent""VaccantPercent" "LowDevPercent" "MedianFloors" "UnitDensity"    "SoilBPercent"   "SoilDPercent"  
# with lasso importance: 0.07469990 0.05483042 0.06197495 0.04273060 0.10180103 0.02472870 0.03686275 0.02002413
# another result "GravelDensity"  "GreenPercent"   "VaccantPercent" "LowDevPercent"  "MedianFloors"  
# with lasso importance: 0.06563165 0.04830392 0.05170203 0.04980743 0.08100649
# the unique variable names are: "GravelDensity""GreenPercent""VaccantPercent""MedianFloors""SoilBPercent""SoilDPercent""LowDevPercent""UnitDensity""ServicePercent""FootprintDensity""AsphaltDensity"  
# 
selected_vars <- table(variables_list)
# AsphaltDensity FootprintDensity    GravelDensity     GreenPercent    LowDevPercent     MedianFloors 
#       2                2               11               10                3               10 
# ServicePercent     SoilBPercent     SoilDPercent      UnitDensity   VaccantPercent 
#       1                5                5                6               11 
################################

#### DO NOT RUN BELOW #######
# test different fdr
# setup initial false discovery rate (FDR) 
initial_fdr <-  0.1
max_fdr <- 0.25
step <- 0.05

res_list <- c() ##

# get results based on different fdr
for (i in 1:((max_fdr-0.1)/step)){
  fdr <- initial_fdr
  # initialize variable list for saving selected variables in loops
  variables_list <- c()
  for (j in (1:500)){
    results <- lapply(Y_df, function(y) knockoff_(X_df, y, fdr))
    if(length(results[["Y_df"]][[1]][["selected"]]) > 0){
      res_list <- c(res_list, results) ##
      df_result <- knockoff_extraction(results)
      variables_list <- c(variables_list, df_result[[1]])
    }
    else{
      variables_list <- variables_list
    }
  }
  # calculate frequency for each selected variable collected in loops
  assign(paste0("variables_frequency_", i), table(variables_list))
  initial_fdr <- initial_fdr + step
}


