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
 
# import event data
ecol_event <- readr::read_csv("Desktop/IRF/ecol_complete.csv")
tp_event <- readr::read_csv("Desktop/IRF/tp_complete.csv")
tss_event <- readr::read_csv("Desktop/IRF/tss_complete.csv")
# import station data
ecol_station <- readr::read_csv("Desktop/IRF/ecol_nfc_wtshd.csv")
tp_station <- readr::read_csv("Desktop/IRF/tp_nfc_wtshd.csv")
tss_station <- readr::read_csv("Desktop/IRF/tss_nfc_wtshd.csv")

# input variable names of data
ecol_event_v <- c('ppt', 'tmean', 'med_high', 
                  'pstaterd', 'pcountyrd', 'rdtype_county', 'psoilA',
                  'pgreen', 'pres', 'pbar', 'pgrass', 
                  'm_ppt', 'AREA_MD_Low impervious', 'LPI_Tree', 'NP_High impervious',
                  'PARA_MD_Grass', 'PARA_MD_Tree', 'PD_Low impervious')
tp_event_v <- c('ppt', 'tmean', 'pbldgunit', 
                'sf_asphalt', 'sf_brick', 'psoilD', 'avg_depth',
                'pvac', 'PD_Tree')
tss_event_v <- c('ppt', 'tmean', 'psfconcrete', 
                 'sf_gravel', 'psoilC', 'psoilD', 'pgreen',
                 'pvac', 'NP_Tree', 'PARA_MD_Medium impervious', 'PD_Tree')

# extract variables from event dataset
ecol_event_X <- ecol_event[ecol_event_v]
tp_event_X <- tp_event[tp_event_v]
tss_event_X <- tss_event[tss_event_v]
event_list_X <- list(ecol_event_X, tp_event_X, tss_event_X)

# extract variables from station dataset
ecol_station_X <- ecol_station[ecol_event_v]
tp_station_X <- tp_station[tp_event_v]
tss_station_X <- tss_station[tss_event_v]
station_list_X <- list(ecol_station_X, tp_station_X, tss_station_X)

# normalize variables and convert variable dataframe to matrix
for (i in (1:3)) {
  num_column <- ncol(event_list_X[[i]])
  for (j in (1:num_column)){
    event_list_X[[i]][j] <- scale(event_list_X[[i]][j])
  }
  event_list_X[[i]] <- as.matrix(event_list_X[[i]])
}

for (i in (1:3)) {
  num_column <- ncol(station_list_X[[i]])
  for (j in (1:num_column)){
    station_list_X[[i]][j] <- scale(station_list_X[[i]][j])
  }
  station_list_X[[i]] <- as.matrix(station_list_X[[i]])
}

# setup response of event data
ecol_event_Y <- as.numeric(scale(ecol_event$results))
tp_event_Y <- as.numeric(scale(tp_event$results))
tss_event_Y <- as.numeric(scale(tss_event$results))
event_list_Y <- list(ecol_event_Y, tp_event_Y, tss_event_Y)
# setup response of station data 
ecol_station_Y <- as.numeric(scale(ecol_station$results))
tp_station_Y <- as.numeric(scale(tp_station$results))
tss_station_Y <- as.numeric(scale(tss_station$results))
station_list_Y <- list(ecol_station_Y, tp_station_Y, tss_station_Y)

# setup training set and testing set
datasplit <- function(response){
  # response should be event_list_X/station_list_X
  # 70% for training/30% for testing
  id_list <- list()
  for (i in (1:length(response))) {
    total_num <- length(response[[i]][,1])
    train_break <- round(0.7 * total_num)
    test_break <- train_break + 1
    train_id <- 1:train_break
    test_id <- test_break:total_num
    #deparse(substitute(event_list))
    id_list <- c(id_list, list(list(train_id, test_id)))  
  }
  return(id_list)
}

# train set and test set
event_id <- datasplit(event_list_X)
station_id <- datasplit(station_list_X)

# calculate stable interactions
stable_inter <- function(var, res, id){
  trainID <- id[[1]]
  testID <- id[[2]]
  fit <- iRF(x=var[trainID,], 
             y=res[trainID], 
             xtest=var[testID,], 
             ytest=res[testID], 
             n.iter=5, 
             n.core=4,
             interactions.return = TRUE,
             select.iter = TRUE,
             n.bootstrap=10)
  #fit[["interaction"]]
  return(fit)
}
 
# get stable interactions
#event_ineractoin <- stable_inter(event_list_X[[1]], event_list_Y[[1]], event_id[[1]])
for (i in (1:length(event_id))) {
  event_ineractoin <- stable_inter(event_list_X[[i]], event_list_Y[[i]], event_id[[i]])
  assign(paste0('event_', i), event_ineractoin)
}

for (i in (1:length(station_id))) {
  station_ineractoin <- stable_inter(station_list_X[[i]], station_list_Y[[i]], station_id[[i]])
  assign(paste0('station_', i), station_ineractoin)
}

# plot results
plot_result <- function(data, title = '', path) {
  #name <- deparse(substitute(data))
  # remove interaction with stability lower than 0.5
  data <- data$interaction
  data <- data[data$stability >= 0.5]
  data <- data[order(data$stability),]
  toplot <- data$stability
  names(toplot) <- data$int
  if(length(toplot) > 45){
    height <- 1200
    }else{height <- 800}
  jpeg(file=paste0(path, title, ".jpeg"), width = 800, height = height)
  dotchart(toplot, 
           xlab='stability score', xlim=c(0.4, 1),
           main=paste0('Prevalent Features/Interactions of ', title))
  dev.off()
}

# export all results as csv
save_csv <- function(data, name, path){
  data <- as.data.frame(data$interaction)
  write.csv(data, paste0(path, name, ".csv"), row.names = FALSE)
}

# save plot as image in jpg format
save_path <- '/Users/yangxiaohao/Desktop/IRF/'
result_list <- list(event_1, event_2, event_3)
name_list <- list('ecol_event', 'tp_event', 'tss_event')
for (i in (1:length(result_list))) {
  title_ = name_list[[i]]
  plot_result(result_list[[i]], title = title_, path = save_path)
}
# save csv
for (i in (1:length(result_list))) {
  name <- name_list[[i]]
  save_csv(result_list[[i]], name = name, path = save_path)
}

result_list <- list(station_1, station_2, station_3)
name_list <- list('ecol_station', 'tp_station', 'tss_station')
for (i in (1:length(result_list))) {
  title_ = name_list[[i]]
  plot_result(result_list[[i]], title = title_, path = save_path)
}
# save csv
for (i in (1:length(result_list))) {
  name <- name_list[[i]]
  save_csv(result_list[[i]], name = name, path = save_path)
}


