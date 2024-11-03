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
all <- read.csv("../data/processed/WQ_BE_CLM_By_Day.csv")

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

variables <- list(eCol = ecol_v, TP = tp_v, TSS = tss_v)

# choose water quality parameter
wq_i <- 'TP' # 'TP' 'TSS' 'eCol'


# extract variables from event dataset
wq_vars <- variables[[wq_i]]
wq_X <- all[wq_vars]

# setup response of event data
wq_Y <- as.numeric(scale(all[,c(wq_i)]))


# drop the rows that has na values
combined_df <- cbind(wq_Y, wq_X)
clean_combined_df <- na.omit(combined_df)
wq_Y <- clean_combined_df[, 1, drop = FALSE]
wq_X <- clean_combined_df[, -1] 

num_column <- ncol(wq_X)
for (j in (1:num_column)){
    wq_X[j] <- scale(wq_X[j])
}

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
set.seed(123)

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
write.csv(intr, paste0("irf_interaction_by_day_", wq_i, ".csv"), row.names = FALSE)


# write result to plot
# filtered_intr <- subset(intr, Freq > 0.5)
jpeg(file=paste0("irf_interaction_by_day_", wq_i, ".jpeg"), width = 800, height = 1200)
dotchart(as.matrix(intr['stability']),
         labels = as.matrix(intr['int']),
         xlab='stability score', xlim=c(0.4, 1),
         main=paste0('Prevalent Features/Interactions of ', wq_i))
dev.off()
