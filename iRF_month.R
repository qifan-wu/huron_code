# For using 'iFR' package, please update R and R Studio to latest version
# If using Mac, run in terminal: curl -OL http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
#                                sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
# And run in R Studio: devtools::install_github("karlkumbier/iRF2.0")


# returns in the 'interaction' of function 'iRF':
# @return int (string): interactions. names of variables. 
# The + and − indicate whether the interaction is characterized 
# by high or low levels of the indicated feature.
# @return prevalence (numeric): mean of interaction prevalence for class-1 leaf nodes
# @return precision (numeric): precision of leaf nodes for which an interaction is active
# @return cpe (numeric): difference in prevalence between class-0 and class-1 leaf nodes
# @return sta.cpe (numeric): the proportion of times (across n.bootsrap samples) 
# an interaction is more prevalent among class-1 leaf nodes

# more details: https://www.stat.berkeley.edu/~kkumbier/vignette.html
# https://github.com/sumbose/iRF/blob/master/R/stabilityScore.R 
# https://cran.r-project.org/web/packages/iRF/iRF.pdf 

library(iRF)
library(tidyverse)

setwd('/Users/qifanwu/Documents/research/huron/code')
# import event data
all <- read.csv("../data/processed/WQ_BE_CLM_By_Month.csv")

# Rename variables to no _ format
names(all)[names(all) == "ppt_rx1"] <- "pptRX1"

# input variable names of data
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
                 "CommercialPercent","CityRdDensity","FootprintDensity"
                )
variables <- list(eCol = ecol_v, TP = tp_v, TSS = tss_v, Discharge = discharge_v)

# choose water quality parameter
wq_i <- 'Discharge' # 'TP' 'TSS' 'eCol', 'Discharge'


# extract variables from event dataset
wq_vars <- variables[[wq_i]]
wq_X <- all[wq_vars]

# setup response of event data with scale
# wq_Y <- as.numeric(scale(all[,c(wq_i)]))
# setup response of event data without scale
wq_Y <- as.numeric(all[,c(wq_i)])

# drop the rows that has na values
combined_df <- cbind(wq_Y, wq_X)
clean_combined_df <- na.omit(combined_df)
wq_Y <- clean_combined_df[, 1, drop = FALSE]
wq_X <- clean_combined_df[, -1] 

# # scale
# num_column <- ncol(wq_X)
# for (j in (1:num_column)){
#     wq_X[j] <- scale(wq_X[j])
# }

# train set and test set
split_train_test <- function(data, split_ratio = 0.7) {
  set.seed(123)
  
  num_rows <- nrow(data)
  train_size <- floor(num_rows * split_ratio)
  
  train_indices <- sample(1:num_rows, size = train_size, replace = FALSE)
  test_indices <- setdiff(1:num_rows, train_indices)
  
  return(list(train_indices = train_indices, test_indices = test_indices))
}
  

indices <- split_train_test(wq_X, split_ratio = 0.7)

trainID <- indices$train_indices
testID <- indices$test_indices
X_matrix <- as.matrix(wq_X)

# set seed for parallel
library(doParallel)
set.seed(456)

# core_num <- 4
# cl <- makeCluster(core_num)  # Number of cores
# registerDoParallel(cl)
# set.seed(123, "L'Ecuyer-CMRG")

# run iRF
fit <- iRF(x=X_matrix[trainID,],
           y=wq_Y[trainID,],
           xtest=X_matrix[testID,],
           ytest=wq_Y[testID,],
           n.iter=5,
           ntree=1000,
           n.core=-1,
           interactions.return = TRUE,
           select.iter = TRUE,
           n.bootstrap=10)
  

intr <- as.data.frame(fit$interaction)

# write result to csv
write.csv(intr, paste0("irf_interaction_by_month_noscale_", wq_i, ".csv"), row.names = FALSE)


# write result to plot
# filtered_intr <- subset(intr, Freq > 0.5)
jpeg(file=paste0("irf_interaction_by_month_", wq_i, ".jpeg"), width = 800, height = 1200)
dotchart(as.matrix(intr['stability']),
         labels = as.matrix(intr['int']),
         xlab='stability score', xlim=c(0.4, 1),
         main=paste0('Prevalent Features/Interactions of ', wq_i))
dev.off()

#calculate model r-square
y_true <- wq_Y[trainID,]
y_pred <- fit$rf.list$predicted

residuals <- y_pred - y_true
ss_total <- sum((y_true - mean(y_pred))^2)
ss_res <- sum(residuals^2)
r_squared <- 1 - (ss_res / ss_total)
print(r_squared)
# with scale
# TP: -0.02687085 
# TSS: (seed123) -0.06887118 (seed456)-0.08132526
# eCol: 0.02192829, -0.05735514
# discharge: (seed123):0.6440919  (seed456) 0.6474106

# without scale (seed456)
# TP: 0.004427107
# TSS: -0.1047996
# eCol: -0.06439614
# discharge: 0.6448002


# ======= Random forest model ======
library(randomForest)

set.seed(123)
rf_fit <- randomForest(x=X_matrix[trainID,],
       y=wq_Y[trainID,],
       xtest=X_matrix[testID,],
       ytest=wq_Y[testID,],
       ntree=1000)
print(rf_fit)

# library(varImp)
# i_scores <- varImp(rf_fit, conditional=TRUE)
# i_scores <- i_scores %>% tibble::rownames_to_column("var") 
# i_scores$var<- i_scores$var %>% as.factor()
randomForest::varImpPlot(rf_fit, 
                         sort=FALSE, 
                         main="Variable Importance Plot")
# TP: Mean of squared residuals: 1.065577
      # % Var explained: -0.78
      # Test set MSE: 1.22
      # % Var explained: -3.71
# TSS: Mean of squared residuals: 0.6215183
      # % Var explained: -10.13
      # Test set MSE: 0.38
      # % Var explained: 12.17
# eCol: Mean of squared residuals: 1.132352
      # % Var explained: -7.79
      # Test set MSE: 0.81
      # % Var explained: 8.15
# Discharge: Mean of squared residuals: 0.4475618
#             % Var explained: 64.07
#             Test set MSE: 0.29
#             % Var explained: 44.81