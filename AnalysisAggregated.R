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
LAND <- "CanadaCQ"
GENDER <- "Female"
if(LAND == "Canada") {
  # Valid Cohorts for Canada 1890 to 1910
  cohortMin <- 1890
  cohortMax <- 1910
  COHORT <- "1890 to 1910"
  DF <- Canada_HMD_df
} else if(LAND == "Sweden") {
  # Valid Cohort for Sweden 1900 to 1911
  cohortMin <- 1900
  cohortMax <- 1910
  COHORT <- "1900 to 1910"
  DF <- Sweden_Df
} else if (LAND == "CanadaCQ") {
  # Valid Cohorts for Canada 1890 to 1910
  cohortMin <- 1890
  cohortMax <- 1910
  COHORT <- "1890 to 1910"
  DF <- QPP_CPP_DF
  LAND == "Canada"
}

# In general Age above 70
minAge <- 70
maxAge <- 109
maxAgeExtra <- 109
cohortLower <- 1900
cohortUpper <- 1908
COHORT <- paste(cohortLower, "to", cohortUpper)

if(cohortLower < cohortMin) {
  cohortLower <- cohortMin
}

if(cohortUpper > cohortMax) {
  cohortUpper <- cohortMax
}

hmdDfMaxAge104 <- filter(DF, Age >= minAge & Age <= maxAge & Year <= 2020) %>%
  mutate(Cohort = Year - Age)

hmdDfMaxAge109 <-
  filter(DF, Age >= minAge & Age <= maxAgeExtra & Year <= 2020) %>%
  mutate(Cohort = Year - Age)

# Find full Cohorts
findFullCohortYears(hmdDfMaxAge109)

for(cohort in unique(hmdDfMaxAge109$Cohort)) {
  print(cohort)
  c <- hmdDfMaxAge109 %>%
    filter(.,Cohort == cohort)

  print(cohort)
  print(min(c$Age))
  print(max(c$Age))
}

#### NOTES ####
# Canada
# Female
if (LAND == "Canada" & GENDER == "Female") {
  X1HermiteII <- 116
  thresholdAgeGomp <- 99
  thresholdAgeMake <- 98
  X1HermiteV <- 112
}
# Male
if (LAND == "Canada" & GENDER == "Male") {
  X1HermiteII <- 120
  thresholdAgeGomp <- 99
  thresholdAgeMake <- 98
  X1HermiteV <- 113
}

# Sweden
# Female
if (LAND == "Sweden" & GENDER == "Female") {
  X1HermiteII <- 116
  thresholdAgeGomp <- 99
  thresholdAgeMake <- 98
  X1HermiteV <- 112
}
# Male
if (LAND == "Sweden" & GENDER == "Male") {
  X1HermiteII <- 120
  thresholdAgeGomp <- 99
  thresholdAgeMake <- 98
  X1HermiteV <- 113
}

#### Aggregated Data 1890 to 1910 ####

df104 <- hmdDfMaxAge104 %>%
  filter(Cohort >= cohortLower & Cohort <= cohortUpper & Gender == GENDER) %>%
  group_by(Age) %>%
  summarise(Deaths = sum(Deaths), Exposure = sum(Exposure))

df109 <- hmdDfMaxAge109 %>%
  filter(Cohort >= cohortLower & Cohort <= cohortUpper & Gender == GENDER) %>%
  group_by(Age) %>%
  summarise(Deaths = sum(Deaths), Exposure = sum(Exposure))


####  Hermite II ####
bestHermiteIIDf <- bestHermiteII(df104, range = c(110:150))
#### Hermite V ####
bestHermiteVDf <- bestHermiteV(df104, range = c(105:125))
#### Makeham GPD ####
bestMakehamGPDDf <- bestMakehamGPDThresholdAge(df104)
####  Hermite II ####
slice_min(bestHermiteIIDf, order_by = AIC, n = 3)
#### Hermite V ####
slice_min(bestHermiteVDf, order_by = AIC, n = 3)
#### Makeham GPD ####
slice_min(bestMakehamGPDDf, order_by = AIC, n = 3)
X1HermiteII <- 112
thresholdAgeMake <- 96
X1HermiteV <- 107


#### MakehamBeard ####
MakehamBeardModel104 <-
  graduateMakehamBeard(df104, analysis = FALSE)
MakehamBeardModel104Analysis <-
  createGradAnalysis(MakehamBeardModel104)
MakehamBeardModel104Analysis$modelResults$AIC


# title <- paste("Makeham Beard", COHORT, GENDER)
# plotLogMortality(MakehamBeardModel104Analysis,
#                  title = title,
#                  errorBar = TRUE)
# plotResMortality(
#   MakehamBeardModel104Analysis,
#   ylim = 1,
#   title = title
# )

MakehamBeardModel109Analysis <-
  applyFittedCLAModelOnData(MakehamBeardModel104, df109)
MakehamBeardModel109Analysis <-
  extrapolateAnalysisFrom109(MakehamBeardModel109Analysis)
MakehamBeardModel109Analysis <-
  addFittedColumnTo109(MakehamBeardModel109Analysis, maxAge)


title <- paste("Makeham Beard", COHORT, GENDER, "and Extrapolation to", maxAgeExtra)
plotMakehamBeard109 <- plotLogMortalitySplitData(
  MakehamBeardModel109Analysis,
  title = title
)
plotMakehamBeard109
resMakehamBeard109 <- plotResMortality(
  MakehamBeardModel109Analysis,
  ylim = 2,
  title = title
)
resMakehamBeard109


####  Hermite II #### 
bestHermiteIIDf <- bestHermiteII(df104)
slice_min(bestHermiteIIDf, order_by = AIC, n = 3)
# X1HermiteII <- 120
hermiteII104 <-
  graduateHermite(df104, "II", X1 = X1HermiteII, analysis = FALSE)
hermiteII104Analysis <- createGradAnalysis(hermiteII104)
hermiteII104Analysis$modelResults$AIC

# title <- paste("Hermite II", COHORT, GENDER, "with X1", X1HermiteII)
# plotLogMortality(hermiteII104Analysis,
#                  title = title)
# plotResMortality(
#   hermiteII104Analysis,
#   ylim = 1,
#   title = title
# )

hermiteII109Analysis <-
  applyFittedHMTModelOnData(hermiteII104, df109)
hermiteII109Analysis <-
  extrapolateAnalysisFrom109(hermiteII109Analysis)
hermiteII109Analysis <- addFittedColumnTo109(hermiteII109Analysis, maxAge)

title <- paste(
  "Hermite II",
  COHORT,
  GENDER,
  "with X1",
  X1HermiteII,
  "and Extrapolation to",
  maxAgeExtra
)
plotHermiteII109 <- plotLogMortalitySplitData(
  hermiteII109Analysis,
  title = title
)
plotHermiteII109
resHermiteII109 <- plotResMortality(
  hermiteII109Analysis,
  ylim = 2,
  title = title
)
resHermiteII109


#### EVT #### 


bestMakehamGPDDf <- bestMakehamGPDThresholdAge(df104)
slice_min(bestMakehamGPDDf, order_by = AIC, n = 3)
# thresholdAgeMake <- 98
MakehamGPD104 <-
  graduateMakehamGPD(df104, thresholdAge = thresholdAgeMake, analysis = FALSE)
MakehamGPD104Analysis <- createGradAnalysis(MakehamGPD104)
MakehamGPD104Analysis$modelResults$AIC

# title <- paste(
#   "Makeham with GPD",
#   COHORT,
#   GENDER,
#   "with Threshold Age",
#   thresholdAgeGomp
# )
# plotLogMortality(
#   MakehamGPD104Analysis,
#   title = title
# )
# plotResMortality(
#   MakehamGPD104Analysis,
#   ylim = 1,
#   title = title
# )

MakehamGPD109Analysis <-
  applyFittedEVTModelOnData(MakehamGPD104, df109)
MakehamGPD109Analysis <-
  extrapolateAnalysisFrom109(MakehamGPD109Analysis)
MakehamGPD109Analysis <- addFittedColumnTo109(MakehamGPD109Analysis, maxAge)

title <- paste(
  "MakehamGPD",
  COHORT,
  GENDER,
  "with Threshold Age",
  thresholdAgeGomp,
  "and Extrapolation to",
  maxAgeExtra
)
plotMakehamGPD109 <- plotLogMortalitySplitData(
  MakehamGPD109Analysis,
  title = title
)
plotMakehamGPD109
resMakehamGPD109 <- plotResMortality(
  MakehamGPD109Analysis,
  ylim = 2,
  title = title
)
resMakehamGPD109

#### Hermite V ####
bestHermiteVDf <- bestHermiteV(df104)
slice_min(bestHermiteVDf, order_by = AIC, n = 3)
# X1HermiteV <- 110
hermiteV104 <-
  graduateHermite(df104, "V", X1 = X1HermiteV, analysis = FALSE)
hermiteV104Analysis <- createGradAnalysis(hermiteV104)
hermiteV104Analysis$modelResults$AIC

# title <- paste("Hermite V", COHORT, GENDER, "with X1", X1HermiteV)
# plotLogMortality(hermiteV104Analysis,
#                  title = title)
# plotResMortality(
#   hermiteV104Analysis,
#   ylim = 1,
#   title = title
# )

hermiteV109Analysis <- applyFittedHMTModelOnData(hermiteV104, df109)
hermiteV109Analysis <-
  extrapolateAnalysisFrom109(hermiteV109Analysis)
hermiteV109Analysis <- addFittedColumnTo109(hermiteV109Analysis, maxAge)

title <- paste(
  "Hermite V",
  COHORT,
  GENDER,
  "with X1",
  X1HermiteV,
  "and Extrapolation to",
  maxAgeExtra
)
plotHermiteV109 <- plotLogMortalitySplitData(
  hermiteV109Analysis,
  title = title
)
plotHermiteV109
resHermiteV109 <- plotResMortality(
  hermiteV109Analysis,
  ylim = 2,
  title = title
)
resHermiteV109


#### Mega Plot ####
megaPlot <- bind_rows(MakehamBeardModel109Analysis$plottingData,
                          hermiteII109Analysis$plottingData,
                          MakehamGPD109Analysis$plottingData,
                          hermiteV109Analysis$plottingData)
title <- paste(
  "Models", COHORT, GENDER,
  "and Extrapolation to",
  maxAgeExtra
)
plotLogMortalitySplitData(
  megaPlot,
  title = title,
  hideLegend = FALSE
)
megaPlotExtra <- bind_rows(MakehamBeardModel109Analysis$extrapolatedData,
                           hermiteII109Analysis$extrapolatedData,
                           MakehamGPD109Analysis$extrapolatedData,
                           hermiteV109Analysis$extrapolatedData)
maxAgeExtrapol <- 125
title <- paste(
  "Models", COHORT, LAND, GENDER,
  "and Extrapolation to",
  maxAgeExtrapol
)
plotLogMortalitySplitData(
  megaPlotExtra,
  xlim = c(70, 125),
  ylim = c(-4.5, 2.5),
  title = title,
  hideLegend = FALSE
)

AggregatedPatchwork <- plotMakehamBeard109 + plotHermiteII109 + plotMakehamGPD109 + plotHermiteV109
AggregatedPatchworkRes <- resMakehamBeard109 + resHermiteII109 + resMakehamGPD109 + resHermiteV109

title <- paste(LAND, "Aggregated", GENDER, "Patchwork.png", sep = "")
titleRes <- paste(LAND, "Aggregated", GENDER, "PatchworkRes.png", sep = "")
# ggsave(title, plot = AggregatedPatchwork, width = 30, height = 20, units = "cm", path = "plots")
# ggsave(titleRes, plot = AggregatedPatchworkRes, width = 30, height = 20, units = "cm", path = "plots")



MakehamBeardModel109Analysis$modelResults$AIC
hermiteII109Analysis$modelResults$AIC
MakehamGPD109Analysis$modelResults$AIC
hermiteV109Analysis$modelResults$AIC

