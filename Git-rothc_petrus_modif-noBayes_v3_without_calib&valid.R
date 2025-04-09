
### LIBRARIES ###
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(gridExtra)
library(SoilR)
library(raster)



### GLOBAL CONSTANTS ###
pE = 1.0 #Evaporation coefficient .Since we use ET-timeseries
soil.thick = 50 #Soil thickness (organic layer topsoil), in cm
bare = FALSE
Tspinup = 500 #Years for spin-up runs

### DATA ###
path<-"C:/Users/a_fer/Documents/SoilWatch/EcoGaia/field_campaign/"
#Load data
SOCmeasurements <- read.csv("C:/Users/a_fer/Documents/SoilWatch/EcoGaia/field_campaign/soc_winsorized_0-50.csv")
SOCmeasurements$Strata<-1 # we are not using by stratum now, so all the same
df_NPP <-  read.csv2("C:/Users/a_fer/Documents/SoilWatch/EcoGaia/RothC/NPP_annualMean_SD.csv")
df_NPP$Strata<-1# we are not using by stratum now, so all the same
#consider only not harvested part
df_NPP$meanNPPresid<-as.numeric(df_NPP$mean_NPP)*(1-0.37)


climate<-read.csv2("C:/Users/a_fer/Documents/SoilWatch/EcoGaia/RothC/climate_monthly_covar.csv")
climate$Strata<-1# we are not using by stratum now, so all the same
df_climate_spinup <- climate[which(climate$period=="jan2001-dec2019"),]  
df_climate_warmup <- climate[which(climate$period=="jan2020-dec2024"),]
df_clay<-SOCmeasurements[,c("Clay","Strata")] #17


spinupyears=seq(1/12,Tspinup,by=1/12)
setData <- function(stratum) {
  #Function that loads strata specific data required for model
  #Input: Stratum number
  #Output: Data for RothC
  stratum = as.character(stratum)
  
  #Spin-up phase data
  #Data 2001-2019
  spinup_data = df_climate_spinup[df_climate_spinup$Strata==stratum,]
  TempSpinup = df_climate_spinup[,c(1,5)]
  TempSpinup$mean_TempERA<-as.numeric(TempSpinup$mean_TempERA)
  PrecipSpinup = df_climate_spinup[,c(1,2)]
  PrecipSpinup$mean_Prec<-as.numeric(PrecipSpinup$mean_Prec)
  EvpSpinup = df_climate_spinup[,c(1,4)]
  EvpSpinup$mean_ET<-as.numeric(EvpSpinup$mean_ET)
  #CinputSpinup = Annual C inputs to soil (in Mg/ha/yr) during spin-up period
  
  #Warm-up phase data
  warmup_period = c(2020,2025)
  warmup_data = df_climate_warmup[df_climate_warmup$Strata==stratum,]
  Temp = df_climate_warmup[,c(1,5)]
  Temp$mean_TempERA<-as.numeric(Temp$mean_TempERA)
  Precip = df_climate_warmup[,c(1,2)]
  Precip$mean_Prec<-as.numeric(Precip$mean_Prec)
  Evp =  df_climate_warmup[,c(1,4)]
  Evp$mean_ET<-as.numeric(Evp$mean_ET)
  #CinputWarmup = Annual C inputs to soil (in Mg/ha/yr) during warm-up period
  Twarmup = warmup_period[2]-warmup_period[1]   #Years for warm-up. Round up the number of years: E.g. if period 1/2022 - 3/2023, then Twarmup = 2
  warmupyears=seq(1/12,Twarmup,by=1/12) 
  
  #SOC data
  SOCdata = as.numeric(SOCmeasurements$SOC_ESM_0.50[SOCmeasurements$Strata==stratum]) #Soil organic carbon in Mg/ha 
  monthSOCsamples = 10 #Month (as a number 1-12) when SOC measurements were taken
  monthWARMUPstarts = 1 #Month (as a number 1-12) when warm-up period starts
  FallIOM=0.049*mean(SOCdata)^(1.139) #IOM using Falloon method #Note we use SOC mean instead of individual samples
  clay = as.numeric(df_clay$Clay[df_clay$Strata==stratum])#as.numeric(df_clay)
  clay = mean(clay) #mean clay of stratum
  
  return(list(TempSpinup=TempSpinup,PrecipSpinup=PrecipSpinup,EvpSpinup=EvpSpinup,Temp=Temp,Precip=Precip,
              Evp=Evp,Twarmup=Twarmup,warmupyears=warmupyears,SOCdata=SOCdata,
              monthSOCsamples=monthSOCsamples,monthWARMUPstarts=monthWARMUPstarts,FallIOM=FallIOM,clay=clay))
}



### ROTHC AUXILIARY FUNCTIONS ###

source("C:/Users/a_fer/Documents/SoilWatch/EcoGaia/RothC/rothcscriptscarbonolocalproject/modifiedFunctions.R")

runRothC <- function(a1,a2,b1,resp_pwp,CinputSpinup,CinputWarmup,DR,clay,FallIOM,
                     Temp,Precip,Evp,TempSpinup,PrecipSpinup,EvpSpinup,
                     warmupyears,Twarmup,monthWARMUPstarts,monthSOCsamples) {
  #Function that runs both spin-up and warm-up phases of RothC 
  
  #Spin-up phase of Tspinup years using average climate data
  Cinitial = c(0,0,0,0,FallIOM)
  fT=fT_modified.RothC(TempSpinup[,2],a1,a2) #Temperature effects per month
  fW=fW_modified.RothC(P=(PrecipSpinup[,2]),E=(EvpSpinup[,2]),S.Thick=soil.thick,pClay=clay,pE=pE,bare=bare,b1=b1,resp_pwp=resp_pwp)$b #Moisture effects per month
  xi=data.frame(spinupyears,rep(fT*fW,length.out=length(spinupyears)))
  spinupModel=RothCModel(t=spinupyears,C0=c(DPM=Cinitial[1], RPM=Cinitial[2], BIO=Cinitial[3], HUM=Cinitial[4], IOM=Cinitial[5]),
                         In=CinputSpinup, clay=clay, xi=xi, DR=DR)
  SOCeq_spinup = getC(spinupModel)
  SOCeq_spinup = SOCeq_spinup[Tspinup*12 - (12-monthWARMUPstarts),]
  #Warm-up phase of Twarmup years using climate data from Twarmup years period
  fT=fT_modified.RothC(Temp[,2],a1,a2) #Temperature effects per month
  fW=fW_modified.RothC(P=(Precip[,2]),E=(Evp[,2]),S.Thick=soil.thick,pClay=clay,pE=pE,bare=bare,b1=b1,resp_pwp=resp_pwp)$b #Moisture effects per month
  xi=data.frame(warmupyears,rep(fT*fW,length.out=length(warmupyears)))
  warmupModel=RothCModel(t=warmupyears,C0=c(DPM=SOCeq_spinup[1],RPM=SOCeq_spinup[2],BIO=SOCeq_spinup[3],HUM=SOCeq_spinup[4],IOM=Cinitial[5]),
                         In=CinputWarmup, clay=clay, xi=xi, DR=DR)
  SOCeq_warmup = getC(warmupModel)
  SOCeq_warmup = SOCeq_warmup[Twarmup*12 - (12-monthSOCsamples),] #return equilibrium SOC at the month when SOC was measured
  
  return(list(SOCeq_spinup,SOCeq_warmup))
}



### USING ROTH C function #### 
stratum=1
a1= 47.9
a2= 106
b1= 0.444
resp_pwp=0.19
DR=1.44
data <- setData(stratum)
Cinput_spinup = as.numeric(df_NPP[df_NPP$Strata==stratum & df_NPP$period=="jan2001-dec2019",'meanNPPresid'])
Cinput_warmup = as.numeric(df_NPP[df_NPP$Strata==stratum & df_NPP$period=="jan2020-dec2024",'meanNPPresid'])


#run RothC
pools = runRothC(a1,a2,b1,resp_pwp,Cinput_spinup,Cinput_warmup,DR,data$clay,data$FallIOM,
                 data$Temp,data$Precip,data$Evp,data$TempSpinup,data$PrecipSpinup,data$EvpSpinup,
                 data$warmupyears,data$Twarmup,data$monthWARMUPstarts,data$monthSOCsamples)
SOCequi <- sum(pools[[2]]) #SOC = sum of warm-up carbon pools #for 2024
pools_spinup<-pools[[1]]   #SOC = sum of spin-up carbon pools #for 2019





#### FOWARD RUN #####
### LOAD MODEL ###

SOCequipools <- pools_spinup #take spin-up SOC pools
pools_warmup<-pools[[2]]
SOCeq<- pools_warmup

source("C:/Users/a_fer/Documents/SoilWatch/EcoGaia/RothC/rothcscriptscarbonolocalproject/modifiedFunctions.R")

simulateRothC <- function(a1,a2,b1,resp_pwp,DR,clay,FallIOM,Tyears,SOCeq,Cinput) {
  
  #Assume climate data 2016-2024 by default
  warmup_data = df_climate_warmup[df_climate_warmup$Strata==stratum,]
  Temp = df_climate_warmup[,c(1,5)]
  Temp$mean_TempERA<-as.numeric(Temp$mean_TempERA)
  Precip = df_climate_warmup[,c(1,2)]
  Precip$mean_Prec<-as.numeric(Precip$mean_Prec)
  Evp =  df_climate_warmup[,c(1,4)]
  Evp$mean_ET<-as.numeric(Evp$mean_ET)
  
  years=seq(1/12,Tyears,by=1/12) 
  fT=fT_modified.RothC(Temp[,2],a1,a2) #Temperature effects per month
  fW=fW_modified.RothC(P=(Precip[,2]),E=(Evp[,2]),S.Thick=soil.thick,pClay=clay,pE=pE,bare=bare,b1=b1,resp_pwp=resp_pwp)$b #Moisture effects per month
  xi=data.frame(years,rep(fT*fW,length.out=length(years)))
  FYM = 0
  
  Model=RothCModel(t=years,C0=c(DPM=SOCeq[1],RPM=SOCeq[2],BIO=SOCeq[3],HUM=SOCeq[4],IOM=FallIOM),
                   In=Cinput, clay=clay, xi=xi, DR=DR, FYM=FYM)
  pools = getC(Model)
  pools = pools #return SOCs starting from January 2016 until December 2016+Tyears
  
  return(pools)
}


simulateSOC <- function(Tyears,stratum,NPP){
  soc <- matrix(nrow=0, ncol=(Tyears))
  data <- setData(stratum)
  pools <- SOCeq
  soc_yearly <- Tyears
  for (year in 1:Tyears){ #loop through years to allow yearly NPP values to be used
    Cinput <- NPP[year]
    pools <- simulateRothC(a1,a2,b1,resp_pwp,DR,data$clay,data$FallIOM,
                           1,pools,Cinput) #Average climate data from 2016-2024 assumed
    pools <- pools[nrow(pools),]  #take final time value i.e. December
    soc_yearly[year] <- sum(pools) #sum over Cpools
  }
  soc <- rbind(soc,soc_yearly) #sum over posterior samples
  rownames(soc) <- NULL
  return(soc)
}



###PRACTICED SCENARIO ##
#-	get monthly baresoil/cover  this is something it can be changed in the forward (for baseline some bare soil, and in the scenario without bare soil -if they are applying conservative tillage). 
#-	agb = yield 9.27 ton/ha + residue 1.23 ton/ha
#-	bgb= allom equations hemp (other crops?)  BGB:AGB ratio 0.183 (same as root:shoot)  bgb = (9.27+1.23) x 0.183  = 1.92 ton/ha 
# Cinputs = 12.42 x 0.47 = 5.837 tonnes C/ha

#So instead of 5.02 for spin up and 5.19 for warm up (from modis) and 5.837 for hemp 
#C spin up= 5.02 * (1-0.37)= 3.16, 
#C warm up= 5.19 * (1-0.37)= 3.27  ////////// OR Cmean 2020-2023 (before hemp in 2024) = 5.176366* (1-0.37)= 3.26 
#C hemp= 3.67 (using harvest index of 10%)

Cinput20_23=3.26

Tyears = 41
YEARS = seq(from=2024,to=2024+Tyears-1,by=1)
soc_tons_ha <- matrix(nrow=0,ncol=Tyears)
#soc_tons_ha_baseline<- matrix(nrow=0,ncol=Tyears)
soc_tons_ha_practiced <- matrix(nrow=0,ncol=Tyears)

mean_values <- rep(Cinput20_23, length(YEARS))#rep(Cinput_warmup, length(YEARS))
# Assign the Cinput_warmup mean value for all years
#df_CINPUT_yearly_BASE <- data.frame(Year = YEARS, Mean = mean_values)
# Assign the Cinput_warmup mean value, with every 3th year starting from 2024 set to Cinputs from Rusty (hemp)
mean_values[seq(1, length(YEARS), by = 3)] <- 3.67
df_CINPUT_yearly_PRACT <- data.frame(Year = YEARS, Mean = mean_values)
print(df_CINPUT_yearly_PRACT)


stratum=1
a1= 47.9
a2= 106
b1= 0.444
resp_pwp=0.19
DR=1.44


STRATA<-1
for (stratum in STRATA){
  #baseline <- simulateSOC(Tyears=Tyears,stratum=stratum,NPP=df_CINPUT_yearly_BASE$Mean)
  practiced <- simulateSOC(Tyears=Tyears,stratum=stratum,NPP=df_CINPUT_yearly_PRACT$Mean)
  #carbon_sequestration <- apply(practiced,2,mean) - apply(baseline,2,mean) #mean is over posterior samples
  #soc_tons_ha <- rbind(soc_tons_ha,carbon_sequestration)
  #soc_tons_ha_baseline <- rbind(soc_tons_ha_baseline,apply(baseline,2,mean))
  soc_tons_ha_practiced <- rbind(soc_tons_ha_practiced,apply(practiced,2,mean))
}
#rownames(soc_tons_ha) <- STRATA
#colnames(soc_tons_ha) <- YEARS
#rownames(soc_tons_ha_baseline) <- STRATA
#colnames(soc_tons_ha_baseline) <- YEARS
rownames(soc_tons_ha_practiced) <- STRATA
colnames(soc_tons_ha_practiced) <- YEARS
#print(soc_tons_ha)
#mean_annual <- apply(soc_tons_ha,2,mean)


library(ggplot2)
library(tidyr)  # For pivoting data

# Convert matrices to data frames
# df_baseline <- data.frame(Year = as.numeric(colnames(soc_tons_ha_baseline)), 
#                           SOC = as.vector(soc_tons_ha_baseline),
#                           Scenario = "Baseline")

df_baseline<- data.frame(Year = as.numeric(colnames(soc_tons_ha_practiced)), 
                         SOC = SOCequi, ### CHANGE TO SOC AT DEC 2023
                         Scenario = "Baseline")

df_practiced <- data.frame(Year = as.numeric(colnames(soc_tons_ha_practiced)), 
                           SOC = as.vector(soc_tons_ha_practiced),
                           Scenario = "Practiced")


# Combine both datasets
df_combined <- bind_rows(df_baseline, df_practiced)

# Plot
ggplot(df_combined, aes(x = Year, y = SOC, color = Scenario, group = Scenario)) +
  geom_line(size = 1.2) +
  ylim(70, 80) +
  geom_point(size = 2) +
  labs(title = "SOC Baseline vs. Practiced",
       x = "Year",
       y = "SOC (tons/ha)") +
  scale_color_manual(values = c("Baseline" = "blue", "Practiced" = "red")) +
  theme_minimal()

# Convert to dataframe
#df_soc <- as.data.frame(t(soc_tons_ha))  # Transpose matrix
#colnames(df_soc)<-"SOC_tons_ha"
#df_soc$Year <- as.numeric(rownames(df_soc))  # Extract years
df_soc<-df_practiced[,c(1:2)]
df_soc$SOC_tons_ha<-df_soc$SOC-SOCequi

# Plot carbon sequestration
ggplot(df_soc, aes(x = Year, y = SOC_tons_ha)) +
  geom_line(color = "blue", size = 1.2) +
  ylim(0, 4) +
  geom_point(color = "red", size = 2) +
  labs(title = "Carbon Sequestration Over Time",
       x = "Year",
       y = "SOC (tons/ha)") +
  theme_minimal()

# Define area values for 2024-2063 from PDD
areas <- c(192.82, 1192, 3192, 5192, 7192, 9192, 1192, 13192)  
areas <- c(areas, rep(15192, 2064 - 2032 + 1))  
print(areas)
df_soc$areas<-areas
df_soc$total<-df_soc$SOC_tons_ha*areas # 192.82   ha # total annual sequestration
df_soc2<-df_soc[,c(1,3,4,5)]#df_soc[,c(2,1,3)]
df_soc2$CO2eq_tons <- df_soc2$total * 3.67 # Convert to CO2 equivalents
df_soc2$CO2eq_cumsum <- cumsum(df_soc2$CO2eq_tons) # Accumulative CO2eq stock (tons)
df_soc2[,c(2:6)]<-round(df_soc2[,c(2:6)],2)
print(df_soc2)
accumulated_soc40years<-sum(df_soc2$SOC_tons_ha)
accumulated_soc40years_AOI<-sum(df_soc2$total)

df_soc2$CO2eq_tons_ha <- df_soc2$CO2eq_tons/df_soc2$areas
