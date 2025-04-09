
### LIBRARIES ###
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(gridExtra)
library(SoilR)
library(raster)


# DEFINING C INPUT ITERATIVELY RUNNING SPIN UP AND WARM UP AND COMPARING WITH MEAN SOC OF FIELD

### GLOBAL CONSTANTS ###
pE = 1.0 #Evaporation coefficient .Since we use ET-timeseries
soil.thick = 50 #Soil thickness (organic layer topsoil), in cm
bare = FALSE

### DATA ###
path<-"C:/Users/a_fer/Documents/SoilWatch/EcoGaia/field_campaign/"
#Load data
SOCmeasurements <- read.csv("C:/Users/a_fer/Documents/SoilWatch/EcoGaia/field_campaign/soc_winsorized_0-50.csv")
SOCmeasurements$Strata<-1 # we are not using by stratum now, so all the same

# separate calib and valid data
set.seed(123) 
calib_proportion <- 0.7
calib_indices <- sample(seq_len(nrow(SOCmeasurements)), size = floor(calib_proportion * nrow(SOCmeasurements)))
df_calibration <- SOCmeasurements[calib_indices, ] #calibration SOC data
df_validation <- SOCmeasurements[-calib_indices, ] #validation SOC data

df_NPP <-  read.csv2("C:/Users/a_fer/Documents/SoilWatch/EcoGaia/RothC/NPP_annualMean_SD.csv")
df_NPP$Strata<-1# we are not using by stratum now, so all the same
climate<-read.csv2("C:/Users/a_fer/Documents/SoilWatch/EcoGaia/RothC/climate_monthly_covar.csv")
climate$Strata<-1# we are not using by stratum now, so all the same
df_climate_spinup <- climate[which(climate$period=="jan2001-dec2019"),]  
df_climate_warmup <- climate[which(climate$period=="jan2020-dec2024"),]

df_clay<-SOCmeasurements[,c("Clay","Strata")]
df_clay$Strata<-1# we are not using by stratum now, so all the same
Clay<-mean(df_clay$Clay)



########## ROTH C MODEL FUNCTION###############
Roth_C<-function(Cinputs,years,DPMptf, RPMptf, BIOptf, HUMptf, FallIOM,Temp,Precip,Evp,soil_cover,soil.thick,SOC,clay,DR,bare1)
{
  #Temperature effects per month
  fT=fT.RothC(Temp[,2]) 
  #Moisture effects per month .
  fw1func<-function(P, E, S.Thick = 50, pClay = clay, pE = 1, bare) 
  { M = P - E * pE
  Acc.TSMD = NULL
  for (i in 2:length(M)) {
    B = ifelse(bare[i] == FALSE, 1, 1.8)
    Max.TSMD = -(20 + 1.3 * pClay - 0.01 * (pClay^2)) * (S.Thick/23) * (1/B)
    Acc.TSMD[1] = ifelse(M[1] > 0, 0, M[1])
    if (Acc.TSMD[i - 1] + M[i] < 0) {
      Acc.TSMD[i] = Acc.TSMD[i - 1] + M[i]
    }
    else (Acc.TSMD[i] = 0)
    if (Acc.TSMD[i] <= Max.TSMD) {
      Acc.TSMD[i] = Max.TSMD
    }
  }
  b = ifelse(Acc.TSMD > 0.444*Max.TSMD, 1,(0.2+0.8*((Max.TSMD-Acc.TSMD)/(Max.TSMD-0.444* Max.TSMD))))
  b<-clamp(b,lower=0.2)
  return(data.frame(Max.TSMD,Acc.TSMD,b,bare))	
  }
  fW_2<- fw1func(P=(Precip[,2]), E=(Evp[,2]), S.Thick = soil.thick, pClay = clay, pE = 1, bare=bare1)$b 
  #Vegetation Cover effects 
  fC<-soil_cover[1]
  # Set the factors frame for Model calculations
  xi.frame=data.frame(years,rep(fT*fW_2*fC,length.out=length(years)))
  # RUN THE MODEL FROM SOILR
  Model3_spin=RothCModel(t=years,C0=c(DPMptf, RPMptf, BIOptf, HUMptf, FallIOM),In=Cinputs,DR=DR,clay=clay,xi=xi.frame, pass=TRUE)
  Ct3_spin=getC(Model3_spin)
  return(Ct3_spin)
}



### USING ROTH C function #### 
CPOOL<-as.data.frame(NULL)


for (stratum in 1:1) {
  data_sample<-df_calibration$SOC_ESM_0.50[which(df_calibration$Strata==stratum)]
  SOC <- mean(data_sample)
  clay = as.numeric(df_clay$Clay[df_clay$Strata==stratum])#as.numeric(df_clay)
  clay = mean(clay)
  
  # IOM calculation using Falloon method
  FallIOM <- 0.049 * SOC^(1.139)
  
  C_in <- seq(1/10, 15, by = 1/10)
  C_loop <- NULL
  
  for (b in 1:length(C_in)) {
    # spin up
    n <- 500 # 500 years to find equilibrium
    years <- seq(1/12, n, by = 1/12) # 6000 rows (12 months * 500 years)
    soil_cover <- rep(1, length(years)) # fully covered soil
    
    spinup_data = df_climate_spinup[df_climate_spinup$Strata==stratum,]
    TempSpinup = df_climate_spinup[,c(1,5)]
    TempSpinup$mean_TempERA<-as.numeric(TempSpinup$mean_TempERA)
    PrecipSpinup = df_climate_spinup[,c(1,2)]
    PrecipSpinup$mean_Prec<-as.numeric(PrecipSpinup$mean_Prec)
    EvpSpinup = df_climate_spinup[,c(1,4)]
    EvpSpinup$mean_ET<-as.numeric(EvpSpinup$mean_ET)
    
    Ct1 <- Roth_C(Cinputs = C_in[b], years = years, DPMptf = 0, RPMptf = 0, BIOptf = 0, HUMptf = 0, 
                  FallIOM = FallIOM, Temp = TempSpinup, Precip = PrecipSpinup, Evp = EvpSpinup, soil_cover = soil_cover, soil.thick = soil.thick, 
                  SOC = SOC, clay = clay, bare1 = rep(F, 12), 
                  DR = 1.44)
    poolSize3_spin <- as.numeric(tail(Ct1, 1))
    Ct1_t <- sum(poolSize3_spin)
    Ct1_t <- as.data.frame(Ct1_t)
    colnames(Ct1_t) <- "SOCmod_spin"
    Ct1_t$Cinp_s<-C_in[b]
    Ct1_t$Cinput_loop <- b
    Ct1_t$stratum <- stratum
    Ct1_t$clay <- clay
    Ct1_t$DPMs <- poolSize3_spin[1]
    Ct1_t$RPMs <- poolSize3_spin[2]
    Ct1_t$BIOs <- poolSize3_spin[3]
    Ct1_t$HUMs <- poolSize3_spin[4]
    Ct1_t$IOMs <- poolSize3_spin[5]
    
    
    # warm up
    n <- 5 
    years <- seq(1/12, n, by = 1/12)
    
    #Warm-up phase data
    warmup_period = c(2019,2024)
    warmup_data = df_climate_warmup[df_climate_warmup$Strata==stratum,]
    Temp = df_climate_warmup[,c(1,5)]
    Temp$mean_TempERA<-as.numeric(Temp$mean_TempERA)
    Precip = df_climate_warmup[,c(1,2)]
    Precip$mean_Prec<-as.numeric(Precip$mean_Prec)
    Evp =  df_climate_warmup[,c(1,4)]
    Evp$mean_ET<-as.numeric(Evp$mean_ET)
    
    Ct2 <- Roth_C(Cinputs = Ct1_t$Cinp_s, years = years, DPMptf = Ct1_t$DPMs, RPMptf = Ct1_t$RPMs, BIOptf = Ct1_t$BIOs, 
                  HUMptf = Ct1_t$HUMs, FallIOM = FallIOM, Temp = Temp, Precip = Precip, Evp = Evp, 
                  soil_cover = soil_cover, soil.thick = soil.thick, SOC = SOC, clay = clay, bare1 = rep(F, 12), 
                  DR = 1.44)
    poolSize2 <- as.numeric(tail(Ct2, 1))
    Ct2_t <- sum(poolSize2)
    Ct1_t$SOCmod_warm<-Ct2_t
    Ct1_t$SOCsample <- SOC
    Ct1_t$DPMw <- poolSize2[1]
    Ct1_t$RPMw <- poolSize2[2]
    Ct1_t$BIOw <- poolSize2[3]
    Ct1_t$HUMw <- poolSize2[4]
    Ct1_t$IOMw <- poolSize2[5]
    C_loop <- rbind(C_loop, Ct1_t)
  }
  
  C_best <- C_loop$Cinput_loop[which.min(abs(SOC - C_loop$SOCmod_warm))]
  dfBEST<- C_loop[C_best,]
  
  CPOOL <- rbind(CPOOL, dfBEST)
}
# write.csv(CPOOL, "C:/Users/a_fer/Documents/SoilWatch/EcoGaia/RothC/RothC_CPOOL_meanclayfield_thick50_dr144_v2.csv")

# ### VALIDATION #### 
# CPOOL<-read.csv("C:/Users/a_fer/Documents/SoilWatch/EcoGaia/RothC/RothC_CPOOL_meanclayfield_thick50_dr144_v2.csv")
# 
# df_validation<-SOCmeasurements[-calib_indices, ] #validation SOC data#SOCmeasurements
# 
# validation <- data.frame()
# C_valid<- data.frame()
# stratum=1
# #for (stratum in 0:9) {
#  # if (stratum != 7) { # stratum 7 doesn't exist
#     ### Validation -- COMPARE PREDICTION OF SOC WITH VALID. SAMPLES RUNNED WITH MEAN BEST CINPUT AND CPOOLS BY STRATA
#     # validation_SOC is a vector of observed SOC values from your validation samples
#     # model_SOC is a vector of SOC values predicted by the RothC model for the same periods
#     data_sample <- df_validation[which(df_validation$Strata == stratum),]
# 
#     for(sample in 1:length(data_sample$X))
#     {
# 
#       # to run model for each validations sample use mean of Cinput and Cpools for the stratum
#       Cinput<- mean(CPOOL$Cinp_s[which(CPOOL$stratum == stratum)])#from SPIN UP
#       DPM<-mean(CPOOL$DPMs[which(CPOOL$stratum == stratum)])#from SPIN UP
#       RPM<-mean(CPOOL$RPMs[which(CPOOL$stratum == stratum)])#from SPIN UP
#       HUM<-mean(CPOOL$HUMs[which(CPOOL$stratum == stratum)])#from SPIN UP
#       BIO<-mean(CPOOL$BIOs[which(CPOOL$stratum == stratum)])#from SPIN UP
# 
#       # warm up
#       n <- 5
#       years <- seq(1/12, n, by = 1/12)
#       #Warm-up phase data
#       warmup_period = c(2019,2024)
#       warmup_data = df_climate_warmup[df_climate_warmup$Strata==stratum,]
#       Temp = df_climate_warmup[,c(1,5)]
#       Temp$mean_TempERA<-as.numeric(Temp$mean_TempERA)
#       Precip = df_climate_warmup[,c(1,2)]
#       Precip$mean_Prec<-as.numeric(Precip$mean_Prec)
#       Evp =  df_climate_warmup[,c(1,4)]
#       Evp$mean_ET<-as.numeric(Evp$mean_ET)
# 
#       SOC <- data_sample$SOC_ESM_0.50[sample]
#       soil_cover <- rep(1, length(years)) # fully covered soil
#       clay <-mean(as.numeric(df_clay$Clay[which(df_clay$Strata == stratum)]))
# 
#       #FallIOM <- 0.049 * SOC^(1.139) ## FallIOM from actual SOC valid. sample
#       ### use value from calibration as petrus did
# 
#       val <- Roth_C(Cinputs = Cinput, years = years, DPMptf = DPM, RPMptf = RPM, BIOptf = BIO,
#                     HUMptf = HUM, FallIOM = FallIOM, Temp = Temp, Precip = Precip, Evp = Evp,
#                     soil_cover = soil_cover, soil.thick = soil.thick, SOC = SOC, clay = clay,
#                     bare1 = rep(F, 12), DR = 1.44)
# 
#       poolSize <- as.numeric(tail(val, 1))
#       val_t <- as.data.frame(sum(poolSize))
#       colnames(val_t)<-"SOCvalidat"
#       val_t$DPMvalidat <- poolSize[1]
#       val_t$RPMvalidat <- poolSize[2]
#       val_t$BIOvalidat <- poolSize[3]
#       val_t$HUMvalidat <- poolSize[4]
#       val_t$IOMvalidat <- poolSize[5]
#       val_t$X<- data_sample$X[sample]
#       val_t$SOCsample<-data_sample$SOC_ESM_0.50[sample]
#       val_t$stratum<-data_sample$Strata[sample]
# 
# 
#       C_valid <- rbind(C_valid, val_t)
# 
#     }
#   #}
# #}
# C_valid_allstrat<-C_valid
# #write.csv(C_valid_allstrat, paste0(path,"sample_vs_modvalid-all_v3_500.csv"))
# 
# 
# #for (stratum in 0:9) {
# #  if (stratum != 7) {
# stratum=1
#    C_valid<-C_valid_allstrat[which(C_valid_allstrat$stratum == stratum),]
#     model_final_SOC <- C_valid$SOCvalidat
#     validation_SOC <-C_valid$SOCsample
# 
#       differences <- model_final_SOC - validation_SOC
# 
#       # Calculate validation metrics
#       MAE <- mean(abs(differences))
#       RMSE <- sqrt(mean(differences^2))
#       Bias <- mean(differences)
#       MSE <- mean(differences^2)
#       R2 <- 1 - (sum((validation_SOC - model_final_SOC)^2) / sum((validation_SOC - mean(validation_SOC))^2))
# 
#       validation_results <- data.frame(
#         Metric = c("Mean Absolute Error (MAE)-wALL",
#                    "Root Mean Squared Error (RMSE)-wALL",
#                    "Bias-wALL",
#                    "Mean Squared Error (MSE)-wALL",
#                    "Coefficient of Determination (R²)-wALL"),
#         Value = c(MAE, RMSE, Bias, MSE, R2)
#       )
# 
#       ### COMPARE TO MEAN VALUE
#       # Calculate mean of observed SOC values
#       mean_validation_SOC <- mean(validation_SOC)
# 
#       # Calculate mean of predicted SOC values
#       mean_model_final_SOC <- mean(model_final_SOC)
# 
#       # Calculate differences
#       difference <- mean_model_final_SOC - mean_validation_SOC
# 
#       # Calculate validation metrics
#       MAE <- abs(difference)
#       RMSE <- sqrt(difference^2)
#       Bias <- difference
#       MSE <- difference^2
#       #R2 <- 1 - ((mean_validation_SOC - mean_model_final_SOC)^2 / mean((validation_SOC - mean(validation_SOC))^2))
#       NSE <- 1 - (sum((difference)^2) / sum((validation_SOC - mean(validation_SOC))^2))
#       PBIAS <- 100 * difference / mean_validation_SOC
# 
#       new_metrics <- data.frame(
#         Metric = c("Mean Absolute Error (MAE)-wMEAN",
#                    "Root Mean Squared Error (RMSE)-wMEAN",
#                    "Bias-wMEAN",
#                    "Mean Squared Error (MSE)-wMEAN",
#                    "Coefficient of Determination (R²)-wMEAN",
#                    "Nash-Sutcliffe Efficiency (NSE)-wMEAN",
#                    "Percentage Bias (PBIAS)-wMEAN"),
#         Value = c(MAE, RMSE, Bias, MSE, R2, NSE, PBIAS)
#       )
# 
#       # Combine the original results with the new metrics
#       validation_results <- rbind(validation_results, new_metrics)
#       validation_results$stratum <- stratum
# 
#       validation <- rbind(validation, validation_results)
#     #}
#   #}
# 
# validation$Value<-round(validation$Value,2)
# print(validation)
# #write.csv(validation_results_wide, paste0(path,"validation-all_v3_500.csv"))
# 
# 
# 
# 
# #### SIMPLIFYED VALIDATION ####
# # 1. Select validation samples for Stratum 1
# df_validation_stratum1 <- SOCmeasurements[-calib_indices, ] %>%
#   filter(Strata == 1)
# 
# # 2. Compute mean observed SOC
# SOC_mean_val <- mean(df_validation_stratum1$SOC_ESM_0.50, na.rm = TRUE)
# 
# data_sample<-df_calibration$SOC_ESM_0.50[which(df_calibration$Strata==stratum)] ### use the same SOC than for calib to calculate FALLIOM
# SOC <- mean(data_sample)
# 
# clay = as.numeric(df_clay$Clay[df_clay$Strata==stratum])
# clay = mean(clay)
# 
# # 3. Use  C inputs and pool values from spin-up
# Cinput<- CPOOL$Cinp_s[which(CPOOL$stratum == stratum)]#from SPIN UP
# DPM<-CPOOL$DPMs[which(CPOOL$stratum == stratum)]#from SPIN UP
# RPM<-CPOOL$RPMs[which(CPOOL$stratum == stratum)]#from SPIN UP
# HUM<-CPOOL$HUMs[which(CPOOL$stratum == stratum)]#from SPIN UP
# BIO<-CPOOL$BIOs[which(CPOOL$stratum == stratum)]#from SPIN UP
# 
# # 4. Extract climate data for Stratum 1
# warmup_data = df_climate_warmup[df_climate_warmup$Strata==stratum,]
# Temp = df_climate_warmup[,c(1,5)]
# Temp$mean_TempERA<-as.numeric(Temp$mean_TempERA)
# Precip = df_climate_warmup[,c(1,2)]
# Precip$mean_Prec<-as.numeric(Precip$mean_Prec)
# Evp =  df_climate_warmup[,c(1,4)]
# Evp$mean_ET<-as.numeric(Evp$mean_ET)
# 
# # 5. Run RothC model for Stratum 1
# Pred_SOC <- Roth_C(
#   Cinputs = Cinput, years = seq(1/12, 5, 1/12), DPMptf = DPM, RPMptf = RPM, BIOptf = BIO,
#   HUMptf = HUM, FallIOM = FallIOM, Temp = Temp, Precip = Precip, Evp = Evp,
#   soil_cover = 1, soil.thick = 50, SOC = SOC, clay = Clay, bare1 = rep(F, 12), DR = 1.44
# ) %>% tail(1) %>% sum()  # Extract final SOC value
# 
# # 6. Compute validation metrics
# difference <- SOC_mean_val - Pred_SOC
# 
# print(difference)
# 
# 
# ### Validation without running spin up and warm up #### -->As it gives same result for all samples as for the only thing it is used the SOC is to calculate FallIOM that we use it fixed. 
# # Validation can be done for all samples but just comparing final soc obtained compared to every sample
# # Load validation dataset
# df_validation <- SOCmeasurements[-calib_indices, ]  # Exclude calibration samples
# 
# # Extract SOCeq from calibration run for each stratum
# SOCeq_calib <- as.data.frame(CPOOL$SOCmod_warm[which(CPOOL$stratum == stratum)])# Dataframe with SOC from warmup per stratum from calibration
# SOCeq_calib$Strata<-1
# colnames(SOCeq_calib)<-c("SOCmod_warm","Strata")
# 
# # Merge validation data with SOCeq from calibration
# validation_data <- merge(df_validation, SOCeq_calib, by = "Strata", all.x = TRUE)
# 
# # Compute differences
# differences <- validation_data$SOCmod_warm - validation_data$SOC_ESM_0.50# predicted - observed
# 
# # Compute Validation Metrics
# validation_metrics <- data.frame(
#   Metric = c("MAE", "RMSE", "Bias", "MSE", "NSE", "PBIAS"),
#   Value = c(
#     mean(abs(differences)),  # MAE
#     sqrt(mean(differences^2)),  # RMSE
#     mean(differences),  # Bias
#     mean(differences^2),  # MSE
#     1 - (sum(differences^2) / sum((validation_data$SOC_ESM_0.50 - mean(validation_data$SOC_ESM_0.50))^2)),  # NSE
#     100 * mean(differences) / mean(validation_data$SOC_ESM_0.50)  # PBIAS
#   )
# )
# validation_metrics$Value<-round(validation_metrics$Value,2)
# print(validation_metrics)
# #MAE: Measures the average absolute difference between observed SOC and SOCeq. In this case, on average, SOCeq deviates by 29.69 units from observed SOC.
# #RMSE: Similar to MAE but gives more weight to larger errors.
# #BIAS: Measures systematic overestimation (positive) or underestimation (negative) by the model.
# #NSE: Compares the model’s performance to simply using the mean of observed SOC. 1.0 = perfect match, 0.0 = as good as the mean of observations, <0.0 = worse than using the mean. Here, NSE = 0.00, meaning SOCeq is no better than just using the average SOC as a prediction.
# #PBIAS: Shows how much the model overestimates or underestimates SOC in percentage terms. Negative means slight underestimation, positive means overestimation. -1.81% means SOCeq is 1.81% lower than observed SOC on average, which is a small bias.

#### FOWARD RUN #####
# 2. Use carbon input from spin up and every 4 years cinput of Hemp
Tyears = 40
YEARS = seq(from=2024,to=2024+Tyears,by=1)
mean_values <- rep(Cinput, length(YEARS))
# Assign the Cinput value for all years
df_CINPUT_yearly_BASE <- data.frame(Year = YEARS, Mean = mean_values)
# Assign every 4th year starting from 2024 set to Cinputs from Rusty (hemp)
mean_values[seq(1, length(YEARS), by = 4)] <- 5.834
df_CINPUT_yearly_PRACT <- data.frame(Year = YEARS, Mean = mean_values)
print(df_CINPUT_yearly_PRACT)
# 3. Use pool values from warm-up
DPM <- CPOOL$DPMw[which(CPOOL$stratum == stratum)]
RPM <- CPOOL$RPMw[which(CPOOL$stratum == stratum)]
HUM <- CPOOL$HUMw[which(CPOOL$stratum == stratum)]
BIO <- CPOOL$BIOw[which(CPOOL$stratum == stratum)]


# 4. Extract climate data for Stratum 1
warmup_data = df_climate_warmup[df_climate_warmup$Strata==stratum,]
Temp = df_climate_warmup[,c(1,5)]
Temp$mean_TempERA<-as.numeric(Temp$mean_TempERA)
Precip = df_climate_warmup[,c(1,2)]
Precip$mean_Prec<-as.numeric(Precip$mean_Prec)
Evp =  df_climate_warmup[,c(1,4)]
Evp$mean_ET<-as.numeric(Evp$mean_ET)

# 5. Run RothC model for each year in the dataset and store results for Stratum 1 BASELINE
# Initialize SOC pools
pools <- c(DPM, RPM,BIO, HUM)
# Store results in a vector
Pred_SOC_by_yearBASE <- numeric(length(YEARS))

# Loop through each year and update SOC pools
for (i in seq_along(YEARS)) {
  y <- YEARS[i]
  
  # Get C input for the year
  Cinput_value <- df_CINPUT_yearly_BASE$Mean[df_CINPUT_yearly_BASE$Year == y]
  
  # Run RothC for this year using previous pools
  pools <- Roth_C(
    Cinputs = Cinput_value, years = seq(0, 1, length.out = 12),  # Simulate one year
    DPMptf = pools[1], RPMptf = pools[2], BIOptf = pools[3],
    HUMptf = pools[4], FallIOM = FallIOM, Temp = Temp, 
    Precip = Precip, Evp = Evp, soil_cover = 1, soil.thick = 50, 
    SOC = SOC, clay = Clay, bare1 = rep(F, 12), DR = 1.44) # Extract all month's pools
  
  # update pools with new year 
  pools <- pools%>% tail(1) # Extract last month's pools
  
  # Store total SOC for this year
  Pred_SOC_by_yearBASE[i] <- sum(pools)
}


# 6.  Run RothC model for each year in the dataset and store results for Stratum 1  FOR PRACTICED SCENARIO
# Loop through each year and update SOC pools
pools <- c(DPM, RPM, BIO, HUM)
Pred_SOC_by_yearPRACT <- numeric(length(YEARS))

for (i in seq_along(YEARS)) {
  y <- YEARS[i]
  
  # Get C input for the year
  Cinput_value <- df_CINPUT_yearly_PRACT$Mean[df_CINPUT_yearly_PRACT$Year == y]
  
  # Run RothC for this year using previous pools
  pools <- Roth_C(
    Cinputs = Cinput_value, years = seq(0, 1, length.out = 12),  # Simulate one year
    DPMptf = pools[1], RPMptf = pools[2], BIOptf = pools[3],
    HUMptf = pools[4], FallIOM = FallIOM, Temp = Temp, 
    Precip = Precip, Evp = Evp, soil_cover = 1, soil.thick = 50, 
    SOC = SOC, clay = Clay, bare1 = rep(F, 12), DR = 1.44) # Extract all month's pools
  
  # update pools with new year 
  pools <- pools%>% tail(1) # Extract last month's pools
  
  # Store total SOC for this year
  Pred_SOC_by_yearPRACT[i] <- sum(pools)
}

# 7. difference
Pred_SOC_seq<-Pred_SOC_by_yearPRACT-Pred_SOC_by_yearBASE

library(ggplot2)
library(tidyr)

# Create a dataframe for plotting
df <- data.frame(
  Year = YEARS,
  BASELINE = Pred_SOC_by_yearBASE,
  PRACTICED = Pred_SOC_by_yearPRACT
) %>%
  pivot_longer(cols = -Year, names_to = "Scenario", values_to = "SOC")

# Plot
ggplot(df, aes(x = Year, y = SOC, color = Scenario)) +
  geom_line(size = 1) +
  geom_point(size = 2) + 
  ylim(130,140)+
  labs(title = "SOC Baseline vs Practiced",
       x = "Year", y = "SOC (tons/ha)")+
  theme_minimal()

# Create a data frame with your values
df_soc <- data.frame(Year = YEARS, SOC_tons_ha = Pred_SOC_seq)

# Plot using ggplot2
ggplot(df_soc, aes(x = Year, y = SOC_tons_ha)) +
  geom_line(color = "blue", size = 1.2) +  # Line plot
  geom_point(color = "red", size = 2) +    # Points on the line
  #ylim(0,5)+
  labs(title = "Carbon Sequestration Over Time",
       x = "Year",
       y = "SOC (tons/ha)") +
  theme_minimal()

df_soc$total<-df_soc$SOC_tons_ha*192.82  # 192.82   ha
