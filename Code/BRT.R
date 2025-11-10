library(gbm)
library(dismo)

#1- Parametrization 
#Aim is to find the combination of parameters (lr, tc and nt) that achieves
#minimum predictive error (minimum error for predictions to independent samples).

order_data<-repl_P_PCA_B_T_HM %>% dplyr::select(Occurrence, decimallatitude, decimallongitude, S, PA_N, BEN, ISO, everything())
order_data_s1<-order_data %>% filter(S=="S1")
order_data_s2<-order_data %>% filter(S=="S2")

learningRate<-c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005)
bagFraction<-c(0.4,0.5,0.6, 0.7, 0.8) # In Elith et al., 2008, they say that best 0.5-0.75
treeComplexity<-c(1, 2, 3, 4, 6, 8, 10) #in Roura-Pascual said that best is the number of variables of a given model

df_BRTparamCombi <- expand_grid(
  lr = learningRate,
  bg = bagFraction, 
  tc = treeComplexity)

# Storage for results
results_S1 <- list()
results_S2 <- list()


ApplyingAllCombGBMStep <- function(all_comb, data, rep_id = NULL) {
  
  resultAllComb <- pbapply::pblapply(1:nrow(all_comb), function(i) {
    
    resultOneComb <- tryCatch({
      gbm.step(
        data = data, 
        gbm.x = 8:15, 
        gbm.y = 1,
        family = "bernoulli", 
        tree.complexity = all_comb[i, 3],
        learning.rate = all_comb[i, 1],
        bag.fraction = all_comb[i, 2],
        n.folds = 10,
        silent = TRUE,
        plot.main = FALSE
      )
    }, error = function(e) NULL)
    
    if(is.null(resultOneComb)) {
      return(data.frame(
        nTrees = NA,
        treeComplexity = all_comb[i, 3],
        learningRate = all_comb[i, 1],
        bagFraction = all_comb[i, 2],
        meanDev = NA,
        seDev = NA,
        PA_N = na.omit(data$PA_N)[1],
        S = na.omit(data$S)[1],
        rep = ifelse(is.null(rep_id), NA, rep_id),
        stringsAsFactors = FALSE
      ))
    }
    
    return(data.frame(
      nTrees = resultOneComb$n.trees,
      treeComplexity = all_comb[i, 3],
      learningRate = all_comb[i, 1],
      bagFraction = all_comb[i, 2],
      meanDev = resultOneComb$cv.statistics$deviance.mean,
      seDev = resultOneComb$cv.statistics$deviance.se,
      PA_N = na.omit(data$PA_N)[1],
      S = na.omit(data$S)[1],
      rep = ifelse(is.null(rep_id), NA, rep_id),
      stringsAsFactors = FALSE
    ))
  })
  
  return(dplyr::bind_rows(resultAllComb))
}


# ========================================
# FOR S1: 20 reps × 10 PA sets = 200 total
# ========================================
for(pa_set in 1:10) {
  
  # Extract data for this PA set
  data_current <- order_data_s1 %>% dplyr::filter(PA_N %in% c(NA, pa_set))  # Assuming you have a list
  
  # Run 20 reps for this PA set
  results_pa <- pblapply(1:20, function(rep) {
    set.seed(pa_set * 1000 + rep)
    result <- ApplyingAllCombGBMStep(df_BRTparamCombi, data_current, rep_id = rep)
    #result$PA_set <- pa_set
    return(result)
  })
  
  results_S1[[pa_set]] <- dplyr::bind_rows(results_pa)
}

# Combine all S1 results
results_S1_combined <- dplyr::bind_rows(results_S1)

##store the results
write.csv(results_S1_combined,"~XXXXXX/results_parametrization_S1") ##assign the path 


# ========================================
# FOR S2: 200 reps × 1 PA sets = 200 total
# ========================================
results_S2 <- pblapply(1:200, function(rep) {
  set.seed(rep)
  result <- ApplyingAllCombGBMStep(df_BRTparamCombi, data_current, rep_id = rep)
  return(result)
})

results_S2_combined <- dplyr::bind_rows(results_S2)

write.csv(results_S2_combined,"~XXXXXX/results_parametrization_S2") ##assign the path 



Sys.time()
# Test with just 5 reps
results_test <- pblapply(1:1, function(i) {
  set.seed(i)
  ApplyingAllCombGBMStep(df_BRTparamCombi, order_data_s2, rep_id = i)
})
Sys.time()
