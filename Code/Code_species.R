################################################################################
## Project: Intership- Relative influence invasion factor
## Author: Aitana Ralda Corfas
## Date: 05/11/2025
## Description: Pipline of project
################################################################################
library(dplyr)
library(tidyr)

source("Internship UNIL/Code/functions_internship.R")

################################################################################
# 1- Data loading and variable preparation ####
################################################################################
    ##A- Load data, cleaning and thinning ####

#setting name and path for the species 

sp_name<-"L_humile"
sp_path<-"LH"

#Loading occurrences of sp from GBIF in 2025 

sp<-read.csv(paste0("~/Internship UNIL/Data/species/", sp_path, "/", sp_name, ".csv"), sep= "\t")  ###modify for next

#rename some columns

sp<-sp %>% 
  rename(
    decimallatitude = decimalLatitude,
    decimallongitude = decimalLongitude,
    countrycode= countryCode)


##### Cleaning      ####

# We are going to clean the data using Olivia bates and Sebatien Ollier to 
# clean data removing Unvalid and duplicated, occurences, occurrences in lakes,
# sea etc. Information in cleaning.docx

res <- occ_cleaning(sp, details = FALSE, lon = "decimallongitude", lat = "decimallatitude", sea = FALSE, lak = FALSE, map = FALSE)
occ_cle <- sp[res,]
occ_cle<-occ_cle %>% drop_na(year) #remove data without a year
occ_cle$S[occ_cle$year <=1970]<-"S1" #assign the dataset 1 or 2 depending on year
occ_cle$S[occ_cle$year > 1970]<-"S2"
occ_cle$S<-as.factor(occ_cle$S)
table(occ_cle$S)

save(occ_cle, file=paste0("Internship UNIL/Data/species/", sp_path, "/intermediate_", sp_path, "/occ_cle", sp_path,".RData"))  #save intermediate results

### Thinning        ####

# Thinning of the data using code from Olivia bates and Sebatien Ollier. Using a 
# new function which applies thinning for data every 20 years. Which allows to 
# account for spatial autocorrelation of observation while keeping the temporal 
# dimension. More information in thnning.docx

occ_thin20<-thin_breaks(occ_cle, 20)
table(r$S) #new 197 1742 vs old w same method #195 1747 

#reduce number of columns

load("~/Internship UNIL/Data/species/occ.RData")
occ_cle_thin20 = subset(occ_thin20, select = c(colnames(occ), "S"))

save(occ_cle_thin20, file=paste0("Internship UNIL/Data/species/", sp_path, "/intermediate_", sp_path, "/occ_cle_thin", sp_path,".RData"))  #save intermediate results



##D- Trade ####
load("~/Internship UNIL/PREV 05-11/C-Trade/sumimport_LH.RData")
load("~/Internship UNIL/PREV 05-11/LH/i_data/P_PA_PCA_B.RData")


#select only 10 replicates of PA for S1 and 1 for S2
repl_P_PCA_Bentity<- P_PA_PCA_benity %>% filter(S== "S1" | S=="S2" & PA_N %in% c(NA, 1))
View(repl_P_PCA_Bentity)

#Re-assign some of the countries name 
repl_P_PCA_Bentity$ISO[repl_P_PCA_Bentity$ISO == "TUSC"]<-"ITA"
repl_P_PCA_Bentity$ISO[repl_P_PCA_Bentity$ISO == "NA"]<-NA

#done previously 
# > P_PA_PCA_benity$ISO[P_PA_PCA_benity$BEN== "BEN20190"]<-"IDN" #Malaku Islands as Indonesia
# > P_PA_PCA_benity$ISO[P_PA_PCA_benity$BEN== "BEN20594"]<-"IDN" #Riau Islands as Indonesia
# > P_PA_PCA_benity$ISO[P_PA_PCA_benity$BEN== "BEN20570"]<-"IDN" #Lesser Sunda Islands as Indonesia
# > P_PA_PCA_benity$ISO[P_PA_PCA_benity$BEN== "BEN20568"]<-"RUS" #Kuril islands as Russian Empire USSR
# > P_PA_PCA_benity$ISO[P_PA_PCA_benity$BEN== "BEN20560"]<-"GRC" #ionian silands as greece
# > P_PA_PCA_benity$ISO[P_PA_PCA_benity$BEN== "BEN20561"]<-"MEX" #ISlas Marias as mexico
# > P_PA_PCA_benity$ISO[P_PA_PCA_benity$BEN== "BEN20345"]<-"PNG" #Bismark archipelago as papua new guinea
# > P_PA_PCA_benity$ISO[P_PA_PCA_benity$BEN== "BEN20607"]<-"PNG" #Trobriand islands aS papua new guinea
# > P_PA_PCA_benity$ISO[P_PA_PCA_benity$BEN== "BEN20569"]<-"LA" #Lesser Antilles islands
# > P_PA_PCA_benity$ISO[P_PA_PCA_benity$BEN== "BEN20126"]<-"NA" #Europa Island as french western africa


#check the data 
sum(nrow(repl_P_PCA_Bentity[repl_P_PCA_Bentity$S=="S1" & repl_P_PCA_Bentity$Occurrence==0,])) #P= 196  PA=1960 Total: 2156 ---> Nb:Replicates: 10
sum(nrow(repl_P_PCA_Bentity[repl_P_PCA_Bentity$S=="S2"& repl_P_PCA_Bentity$Occurrence==1,])) #P=1713   PA=1731 Total: 3462----> Nb Replicates: 1

#apply trade information 
repl_P_PCA_B_T<-assign_trade(repl_P_PCA_Bentity, all_c_year_sumimport)


#E- Human modification ####
load("Internship UNIL/Data/global_rasters/hm_1970.RData")
load("Internship UNIL/Data/global_rasters/hm_2025.RData")

# Extract the PCA scores for each occurrence with terra (optimized approach):


repl_P_PCA_B_T_HM<-load_HM(repl_P_PCA_B_T)
View(repl_P_PCA_B_T_HM[is.na(repl_P_PCA_B_T_HM$HM),])

