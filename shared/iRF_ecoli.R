
setwd("~/Desktop/summer project/Huron project/data analysis/manuscript")
library(iRF)

tp_event <- readr::read_csv("ecoli_log.csv")

## if change variables, should resave the ecoli_log files in the RF jupyter notebook, and then replace the files in this folder
variables <- c('m_ppt','m_tmean','med_high','medyr_built','pag','pcountyrd','pdev_low','ppt',
               'pres','pgreen','pservice','psoilA','psoilC','pstaterd','tmean')


tp_event_X <- tp_event[variables]

#convert variable dataframe to matrix

tp_event_X <- as.matrix(tp_event_X)


# log response of event data
tp_event_Y <- as.numeric(tp_event$results)


## train and test split

# drop Inf
drop_inf <- function(y, id){
  drop <- which(y == -Inf)
  id <- id[!(id %in% as.vector(drop))]
  return(id)
}

# setup training set and testing set
set.seed(123)
dt = sort(sample(nrow(tp_event), nrow(tp_event)*.9))
train<-tp_event[dt,]
test<-tp_event[-dt,]


# calculate stable interactions
stable_inter <- function(var, res, id){
  trainID <- id[[1]]
  testID <- id[[2]]
  fit <- iRF(x=var[trainID,],
             y=res[trainID],
             xtest=var[testID,],
             ytest=res[testID],
             n.iter=5,
             ntree=1000,
             n.core=-1,
             interactions.return = TRUE,
             select.iter = TRUE,
             n.bootstrap=10)
  #fit[["interaction"]]
  return(fit)
}


fit <- iRF(x=train[variables],
           y=train$results,
           xtest=test[variables],
           ytest=test$results,
           n.iter=5,
           ntree=1000,
           n.core=-1,
           interactions.return = TRUE,
           select.iter = TRUE,
           n.bootstrap=10)

write.csv(fit$interaction,"ecoli_irf.csv")