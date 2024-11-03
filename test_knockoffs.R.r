# if doMC is not installed, please run
install.packages("doMC")

# read dataset
tp<- read_csv("Desktop/knockoffs/tp_nfc_wtshd.csv")
# select variables
variables <- c('sum_hsunit','med_high','medyr_built','pftprint',
               'pbldgunit','struc_area','psfgravel','psfasphalt', 
               'psfconcrete', 'psfbrick','pstaterd','pcountyrd', 
               'pcityrd','road_tl','lanes_tl','road_dens',
               'lanes_dens','sf_gravel','sf_asphalt','sf_conc',
               'sf_brick','rdtype_state','rdtype_county','rdtype_city', 
               'ppccount','par_avgarea','psoilA','psoilB',
               'psoilC','psoilD','avg_depth','pag',
               'pcomm', 'pgreen','pind','pservice',
               'pres','pvac','pdev_open','pdev_low',
               'pdev_med','pdev_high','pbar','pfor',
               'pshrub','pgrass','pcrop')
# setup variables and response
pos_start <- which(names(tp) == 'ppt')
pos_cols <- seq.int(pos_start, ncol(tp))
# normalize variables
for (i in (pos_start:max(pos_cols))){
  tp[i] <- scale(tp[i]) #(tp[i] - mean(tp[i])) / sd(tp[i])
}
# variable matrix
X_tp <- as.matrix((tp[,pos_cols])[variables])
#X_tp <- as.matrix(tp[,pos_cols])
X_tp <- X_tp[,colSums(X_tp) != 0]
# response dataframe
Y_tp <- tp$results
Y_tp <- as.data.frame(Y_tp)

# setup selection function
knockoff_ <- function (X, y, q) {
  # Log-transform the drug resistance measurements.
  y <- log(y)
  fdr <- q
  # Run the knockoff filter.
  result = knockoff.filter(X, y, fdr=fdr)
  # output result
  list(result)
}

# a function for extracting variable names and lasso importance
knockoff_extraction <- function(result){
  # names and importances of selections
  selected_variable <- names(result[["Y_tp"]][[1]][["selected"]])
  variable_importance <- result[["Y_tp"]][[1]][["statistic"]]
  lasso_importance <- c()
  for (i in (1:length(selected_variable))){
    index <- result[["Y_tp"]][[1]][["selected"]][[i]]
    lasso_importance <- c(lasso_importance, variable_importance[index]) 
  }
  # return variable names and lasso importances
  list(selected_variable, lasso_importance)
}

########## DO NOT RUN ##########
# get results
# set false discovery rate
fdr = 0.2
results = lapply(Y_tp, function(y) knockoff_(X_tp, y, fdr))

variable_importance <- rrr[[2]][["statistic"]]
lasso_importance <- c()
for (i in (1:length(rrr[[2]][["selected"]]))){
  index <- rrr[[2]][["selected"]][[i]]
  lasso_importance <- c(lasso_importance, variable_importance[index]) 
}
################################

# test different fdr
# setup initial false discovery rate (FDR) 
initial_fdr <-  0.1
max_fdr <- 0.25
step <- 0.05

# get results based on different fdr
for (i in 1:((max_fdr-0.1)/step)){
  fdr <- initial_fdr
  # initialize variable list for saving selected variables in loops
  variables_list <- c()
  for (j in (1:500)){
    results <- lapply(Y_tp, function(y) knockoff_(X_tp, y, fdr))
    if(length(results[["Y_tp"]][[1]][["selected"]]) > 0){
      df_result <- knockoff_extraction(results)
      variables_list <- c(variables_list, df_result[[1]])
    }
    else{
      variables_list <- variables_list
    }
  }
  # calculate frequency for each selected variable collected in loops
  assign(paste0("variables_frenquency_", i), table(variables_list))
  initial_fdr <- initial_fdr + step
}


