library(tidyverse)
library(patchwork)

source("graduatePoisson.R")
source("PlottingFunctions.R")
source("helperFunctions.R")
source("ExtrapolateFit.R")

load("data/Canada_HMD_df.Rda")
load("data/Sweden_Df_Year.Rda")
load("data/CPP_QPP.Rda")

#### Prep ####
##### Data #####
LAND <- "Sweden"
GENDER <- "Female"
if(LAND == "Canada") {
  # Valid Cohorts for Canada 1890 to 1910
  cohortMin <- 1900
  cohortMax <- 1908
  COHORT <- "1890 to 1910"
  DF <- QPP_CPP_DF
} else if(LAND == "Sweden") {
  # Valid Cohort for Sweden 1900 to 1911
  cohortMin <- 1900
  cohortMax <- 1908
  COHORT <- "1900 to 1910"
  DF <- Sweden_Df
}

# In general Age above 70
minAge <- 70
maxAge <- 109
maxAgeExtra <- 109

hmdDfMaxAge104 <- filter(DF, Age >= minAge & Age <= maxAge & Year <= 2020) %>%
  mutate(Cohort = Year - Age)

hmdDfMaxAge109 <-
  filter(DF, Age >= minAge & Age <= maxAgeExtra & Year <= 2020) %>%
  mutate(Cohort = Year - Age)

# Find full Cohorts
findFullCohortYears(hmdDfMaxAge109)

#### NOTES ####
# Canada
# Female
# if (LAND == "Canada" & GENDER == "Female") {
#   X1HermiteII <- 112
#   thresholdAge <- 98
#   X1HermiteV <- 109
# }
# # Male
# if (LAND == "Canada" & GENDER == "Male") {
#   X1HermiteII <- 120
#   thresholdAge <- 98
#   X1HermiteV <- 113
# }
#
# # Sweden
# # Female
# if (LAND == "Sweden" & GENDER == "Female") {
#   X1HermiteII <- 116
#   thresholdAge <- 98
#   X1HermiteV <- 112
# }
# # Male
# if (LAND == "Sweden" & GENDER == "Male") {
#   X1HermiteII <- 120
#   thresholdAge <- 98
#   X1HermiteV <- 113
# }

#### Singel Cohort ####
##### CREATE PARAMETER DATA #####
MakehamBeardParameters <- data.frame(Model = character(0),
                                     Cohort = numeric(0),
                                     Parameter = character(0),
                                     Value = numeric(0),
                                     SD = numeric(0))
MakehamBeardResult <- data.frame(Model = character(0),
                                 AIC = numeric(0),
                                 Rho = numeric(0))

for(cohort in cohortMin:cohortMax){
  df104 <- hmdDfMaxAge104 %>%
    filter(Cohort == cohort & Gender == GENDER)
  MakehamBeardModel104 <-
    graduateMakehamBeard(df104, analysis = FALSE)
  MakehamBeardModel104Analysis <-
    createGradAnalysis(MakehamBeardModel104)
  
  MakehamBeardParameters <- MakehamBeardParameters %>%
    add_row(tibble_row(
      Model = "Makeham Beard",
      Parameter = "Alpha",
      Cohort = cohort,
      Value = MakehamBeardModel104Analysis$modelResults$Parameters$alpha,
      SD = MakehamBeardModel104Analysis$modelResults$SD$alpha))
  
  MakehamBeardParameters <- MakehamBeardParameters %>%
    add_row(tibble_row(
      Model = "Makeham Beard",
      Parameter = "Beta",
      Cohort = cohort,
      Value = MakehamBeardModel104Analysis$modelResults$Parameters$beta,
      SD = MakehamBeardModel104Analysis$modelResults$SD$beta))
  
  MakehamBeardParameters <- MakehamBeardParameters %>%
    add_row(tibble_row(
      Model = "Makeham Beard",
      Parameter = "Epsilon",
      Cohort = cohort,
      Value = MakehamBeardModel104Analysis$modelResults$Parameters$epsilon,
      SD = MakehamBeardModel104Analysis$modelResults$SD$epsilon))
  
  MakehamBeardParameters <- MakehamBeardParameters %>%
    add_row(tibble_row(
      Model = "Makeham Beard",
      Parameter = "Rho",
      Cohort = cohort,
      Value = MakehamBeardModel104Analysis$modelResults$Parameters$rho,
      SD = MakehamBeardModel104Analysis$modelResults$SD$rho))
  
  
  MakehamBeardResult <- MakehamBeardResult %>%
    add_row(tibble_row(
      Model = "Makeham Beard",
      AIC = MakehamBeardModel104Analysis$modelResults$AIC),
            Rho = MakehamBeardModel104Analysis$modelResults$Parameters$rho)
}


HermiteIIParameter <- data.frame(Model = character(0),
                                     Cohort = numeric(0),
                                     Parameter = character(0),
                                     Value = numeric(0),
                                     SD = numeric(0))
HermiteIIResult <- data.frame(Model = character(0),
                              AIC = numeric(0),
                              HermiteIIX1 = numeric(0))
# X1HermiteII <- 116

for(cohort in c(cohortMin:cohortMax)){
  df104 <- hmdDfMaxAge104 %>%
    filter(Cohort == cohort & Gender == GENDER)
  bestHermiteIIDf <- bestHermiteII(df104, range = c(110:150))
  bestX1 <- slice_min(bestHermiteIIDf, order_by = AIC, n = 1)$X1
  hermiteII104 <- graduateHermite(df104, "II",
                                  X1 = bestX1,
                                  analysis = FALSE)
  hermiteII104Analysis <- createGradAnalysis(hermiteII104)
  
  HermiteIIParameter <- HermiteIIParameter %>%
    add_row(tibble_row(
      Model = "Hermite II",
      Parameter = "hermiteAlpha",
      Cohort = cohort,
      Value = hermiteII104Analysis$modelResults$Parameters$hermiteAlpha,
      SD = hermiteII104Analysis$modelResults$SD$hermiteAlpha))
  
  HermiteIIParameter <- HermiteIIParameter %>%
    add_row(tibble_row(
      Model = "Hermite II",
      Parameter = "hermiteM0",
      Cohort = cohort,
      Value = hermiteII104Analysis$modelResults$Parameters$hermiteM0,
      SD = hermiteII104Analysis$modelResults$SD$hermiteM0))
  
  HermiteIIParameter <- HermiteIIParameter %>%
    add_row(tibble_row(
      Model = "Hermite II",
      Parameter = "hermiteOmega",
      Cohort = cohort,
      Value = hermiteII104Analysis$modelResults$Parameters$hermiteOmega,
      SD = hermiteII104Analysis$modelResults$SD$hermiteOmega))
  
  HermiteIIResult <- HermiteIIResult %>%
    add_row(tibble_row(
      Model = "Hermite II",
      AIC = hermiteII104Analysis$modelResults$AIC,
      HermiteIIX1 = bestX1))
}


MakehamGPDParameter <- data.frame(Model = character(0),
                                     Cohort = numeric(0),
                                     Parameter = character(0),
                                     Value = numeric(0),
                                     SD = numeric(0))
MakehamGPDResult <- data.frame(Model = character(0),
                              AIC = numeric(0),
                               N = numeric(0))
# thresholdAge <- 98

for(cohort in c(cohortMin:cohortMax)){
  df104 <- hmdDfMaxAge104 %>%
    filter(Cohort == cohort & Gender == GENDER)
  bestMakehamGPDDf <- bestMakehamGPDThresholdAge(df104)
  bestThresholdAge <- slice_min(bestMakehamGPDDf, order_by = AIC, n = 1)$tA
  MakehamGPD104 <- graduateMakehamGPD(df104, thresholdAge = bestThresholdAge, analysis = FALSE)
  MakehamGPD104Analysis <- createGradAnalysis(MakehamGPD104)
  
  MakehamGPDParameter <- MakehamGPDParameter %>%
    add_row(tibble_row(
      Model = "MakehamGPD",
      Parameter = "Alpha",
      Cohort = cohort,
      Value = MakehamGPD104Analysis$modelResults$Parameters$alpha,
      SD = MakehamGPD104Analysis$modelResults$SD$alpha))
  
  MakehamGPDParameter <- MakehamGPDParameter %>%
    add_row(tibble_row(
      Model = "MakehamGPD",
      Parameter = "Beta",
      Cohort = cohort,
      Value = MakehamGPD104Analysis$modelResults$Parameters$beta,
      SD = MakehamGPD104Analysis$modelResults$SD$beta))
  
  MakehamGPDParameter <- MakehamGPDParameter %>%
    add_row(tibble_row(
      Model = "MakehamGPD",
      Parameter = "Epsilon",
      Cohort = cohort,
      Value = MakehamGPD104Analysis$modelResults$Parameters$epsilon,
      SD = MakehamGPD104Analysis$modelResults$SD$epsilon))
  
  MakehamGPDParameter <- MakehamGPDParameter %>%
    add_row(tibble_row(
      Model = "MakehamGPD",
      Parameter = "Xi",
      Cohort = cohort,
      Value = MakehamGPD104Analysis$modelResults$Parameters$xi,
      SD = MakehamGPD104Analysis$modelResults$SD$xi))
  
  MakehamGPDResult <- MakehamGPDResult %>%
    add_row(tibble_row(
      Model = "MakehamGPD",
      AIC = MakehamGPD104Analysis$modelResults$AIC,
    N = bestThresholdAge))
}




HermiteVParameter <- data.frame(Model = character(0),
                                 Cohort = numeric(0),
                                 Parameter = character(0),
                                 Value = numeric(0),
                                 SD = numeric(0))
HermiteVResult <- data.frame(Model = character(0),
                             AIC = numeric(0),
                             HermiteVX1 = numeric(0))
# X1HermiteV <- 112

for(cohort in c(cohortMin:cohortMax)){
  df104 <- hmdDfMaxAge104 %>%
    filter(Cohort == cohort & Gender == GENDER)
  bestHermiteVDf <- bestHermiteV(df104, range = 105:125)
  bestX1 <- slice_min(bestHermiteVDf, order_by = AIC, n = 1)$X1
  hermiteV104 <- graduateHermite(df104, "V", X1 = bestX1, analysis = FALSE)
  hermiteV104Analysis <- createGradAnalysis(hermiteV104)
  
  HermiteVParameter <- HermiteVParameter %>%
    add_row(tibble_row(
      Model = "Hermite V",
      Parameter = "hermiteAlpha",
      Cohort = cohort,
      Value = hermiteV104Analysis$modelResults$Parameters$hermiteAlpha,
      SD = hermiteV104Analysis$modelResults$SD$hermiteAlpha))
  
  HermiteVParameter <- HermiteVParameter %>%
    add_row(tibble_row(
      Model = "Hermite V",
      Parameter = "hermiteM0",
      Cohort = cohort,
      Value = hermiteV104Analysis$modelResults$Parameters$hermiteM0,
      SD = hermiteV104Analysis$modelResults$SD$hermiteM0))
  
  HermiteVParameter <- HermiteVParameter %>%
    add_row(tibble_row(
      Model = "Hermite V",
      Parameter = "hermiteOmega",
      Cohort = cohort,
      Value = hermiteV104Analysis$modelResults$Parameters$hermiteOmega,
      SD = hermiteV104Analysis$modelResults$SD$hermiteOmega))
  
  HermiteVParameter <- HermiteVParameter %>%
    add_row(tibble_row(
      Model = "Hermite V",
      Parameter = "hermiteV",
      Cohort = cohort,
      Value = hermiteV104Analysis$modelResults$Parameters$hermiteV,
      SD = hermiteV104Analysis$modelResults$SD$hermiteV))
  
  HermiteVResult <- HermiteVResult %>%
    add_row(tibble_row(
      Model = "Hermite V",
      AIC = hermiteV104Analysis$modelResults$AIC,
      HermiteVX1 = bestX1))
}


fullParameterList <- bind_rows(
  MakehamBeardParameters,
  HermiteIIParameter,
  MakehamGPDParameter,
  HermiteVParameter  
)

fullResultList <- bind_rows(
  MakehamBeardResult,
  HermiteIIResult,
  MakehamGPDResult,
  HermiteVResult  
)

fullResultListTable <- data.frame(
  "MakehamBeard" = MakehamBeardResult$AIC,
  "Rho" = MakehamBeardResult$Rho,
  "MakehamGPD" = MakehamGPDResult$AIC,
  "N" = MakehamGPDResult$N,
  "HermiteII" = HermiteIIResult$AIC,
  "HermiteIIX1" = HermiteIIResult$HermiteIIX1,
  "HermiteV" = HermiteVResult$AIC,
  "HermiteVX1" = HermiteVResult$HermiteVX1)
fullResultListTable

# for(para in unique(fullParameterList$Parameter)) {
#   print(para)
#   temp <- filter(fullParameterList, Parameter == para)
#   print(paste(round(min(temp$Value), digits = 4), round(max(temp$Value), digits = 4)))
#   print(paste(round(min(temp$Value) - max(temp$SD[is.finite(temp$SD)]), digits = 4), round(max(temp$Value) +  max(temp$SD[is.finite(temp$SD)]), digits = 4)))
# }

#
# PARAMETER <- "Alpha"
# title <- paste(unique(MakehamBeardParameters$Model), PARAMETER, COHORT, LAND, GENDER)
# plot <- plotParameter(MakehamBeardParameters, PARAMETER, title = title)
# plot
#
# PARAMETER <- "Beta"
# title <- paste(unique(MakehamBeardParameters$Model), PARAMETER, COHORT, LAND, GENDER)
# plot <- plotParameter(MakehamBeardParameters, PARAMETER, title = title)
# plot
#
# PARAMETER <- "Epsilon"
# title <- paste(unique(MakehamBeardParameters$Model), PARAMETER, COHORT, LAND, GENDER)
# plot <- plotParameter(MakehamBeardParameters, PARAMETER, title = title)
# plot
#
# PARAMETER <- "Rho"
# title <- paste(unique(MakehamBeardParameters$Model), PARAMETER, COHORT, LAND, GENDER)
# plot <- plotParameter(MakehamBeardParameters, PARAMETER, title = title, errorBar = FALSE)
# plot
#
#
#
# PARAMETER <- "hermiteAlpha"
# title <- paste(unique(HermiteIIParameter$Model), PARAMETER, COHORT, LAND, GENDER)
# plot <- plotParameter(HermiteIIParameter, PARAMETER, title = title)
# plot
#
# PARAMETER <- "hermiteM0"
# title <- paste(unique(HermiteIIParameter$Model), PARAMETER, COHORT, LAND, GENDER)
# plot <- plotParameter(HermiteIIParameter, PARAMETER, title = title)
# plot
#
# PARAMETER <- "hermiteOmega"
# title <- paste(unique(HermiteIIParameter$Model), PARAMETER, COHORT, LAND, GENDER)
# plot <- plotParameter(HermiteIIParameter, PARAMETER, title = title)
# plot
#
#
#
# PARAMETER <- "Alpha"
# title <- paste(unique(MakehamGPDParameter$Model), PARAMETER, COHORT, LAND, GENDER)
# plot <- plotParameter(MakehamGPDParameter, PARAMETER, title = title)
# plot
#
# PARAMETER <- "Beta"
# title <- paste(unique(MakehamGPDParameter$Model), PARAMETER, COHORT, LAND, GENDER)
# plot <- plotParameter(MakehamGPDParameter, PARAMETER, title = title)
# plot
#
# PARAMETER <- "Epsilon"
# title <- paste(unique(MakehamGPDParameter$Model), PARAMETER, COHORT, LAND, GENDER)
# plot <- plotParameter(MakehamGPDParameter, PARAMETER, title = title)
# plot
#
# PARAMETER <- "Xi"
# title <- paste(unique(MakehamGPDParameter$Model), PARAMETER, COHORT, LAND, GENDER)
# plot <- plotParameter(MakehamGPDParameter, PARAMETER, title = title)
# plot
#
#
#
#
# PARAMETER <- "hermiteAlpha"
# title <- paste(unique(HermiteVParameter$Model), PARAMETER, COHORT, LAND, GENDER)
# plot <- plotParameter(HermiteVParameter, PARAMETER, title = title)
# plot
#
# PARAMETER <- "hermiteM0"
# title <- paste(unique(HermiteVParameter$Model), PARAMETER, COHORT, LAND, GENDER)
# plot <- plotParameter(HermiteVParameter, PARAMETER, title = title)
# plot
#
# PARAMETER <- "hermiteOmega"
# title <- paste(unique(HermiteVParameter$Model), PARAMETER, COHORT, LAND, GENDER)
# plot <- plotParameter(HermiteVParameter, PARAMETER, title = title)
# plot
#
# PARAMETER <- "hermiteV"
# title <- paste(unique(HermiteVParameter$Model), PARAMETER, COHORT, LAND, GENDER)
# plot <- plotParameter(HermiteVParameter, PARAMETER, title = title)
# plot
#
#
#
#
# ##### STLT AICs #####
# GompertzGPDResult <- data.frame(Model = character(0),
#                                AIC = numeric(0))
#
# for(cohort in c(cohortMin:cohortMax)){
#   df104 <- hmdDfMaxAge104 %>%
#     filter(Cohort == cohort & Gender == GENDER)
#   GompertzGPD104 <- graduateGompertzGPD(df104, thresholdAge = thresholdAge, analysis = FALSE)
#   GompertzGPD104Analysis <- createGradAnalysis(GompertzGPD104)
#
#   GompertzGPDResult <- GompertzGPDResult %>%
#     add_row(tibble_row(
#       Model = "STLT",
#       AIC = GompertzGPD104Analysis$modelResults$AIC))
# }
