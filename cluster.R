# cluster
library(tidyverse)
library(dplyr)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms

setwd('/Users/qifanwu/Documents/research/huron/code')

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


# library(ggsignif)