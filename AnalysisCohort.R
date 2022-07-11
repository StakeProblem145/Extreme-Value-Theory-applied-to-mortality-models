library(tidyverse)
library(patchwork)

source("graduatePoisson.R")
source("PlottingFunctions.R")
source("helperFunctions.R")

load("data/Canada_HMD_df.Rda")
load("data/Sweden_Df_Year.Rda")

#### Prep ####
##### Data #####
LAND <- "Canada"
GENDER <- "Female"
COHORT <- 1900

if(LAND == "Canada") {
  # Valid Cohorts for Canada 1890 to 1910
  DF <- Canada_HMD_df
} else if(LAND == "Sweden") {
  # Valid Cohort for Sweden 1900 to 1911
  DF <- Sweden_Df
}

# In general Age above 70
minAge <- 70
maxAge <- 104
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



#### Singel Cohort ####
##### Data #####
df104 <- hmdDfMaxAge104 %>%
  filter(Cohort == COHORT & Gender == GENDER)

df109 <- hmdDfMaxAge109 %>%
  filter(Cohort == COHORT  & Gender == GENDER)

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
  addFittedColumnTo109(MakehamBeardModel109Analysis)


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
hermiteII109Analysis <- addFittedColumnTo109(hermiteII109Analysis)

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
bestGompertzGPDDf <- bestGompertzGPDThresholdAge(df104)
slice_min(bestGompertzGPDDf, order_by = AIC, n = 3)
# thresholdAgeGomp <- 98
GompertzGPD104 <-
  graduateGompertzGPD(df104, thresholdAge = thresholdAgeGomp, analysis = FALSE)
GompertzGPD104Analysis <- createGradAnalysis(GompertzGPD104)
GompertzGPD104Analysis$modelResults$AIC

# title <- paste(
#   "Gompertz with GPD",
#   COHORT,
#   GENDER,
#   "with Threshold Age",
#   thresholdAgeGomp
# )
# plotLogMortality(
#   GompertzGPD104Analysis,
#   title = title
# )
# plotResMortality(
#   GompertzGPD104Analysis,
#   ylim = 1,
#   title = title
# )

GompertzGPD109Analysis <-
  applyFittedEVTModelOnData(GompertzGPD104, df109)
GompertzGPD109Analysis <-
  addFittedColumnTo109(GompertzGPD109Analysis)

title <- paste(
  "GompertzGPD",
  COHORT,
  GENDER,
  "with Threshold Age",
  thresholdAgeGomp,
  "and Extrapolation to",
  maxAgeExtra
)
plotGompertzGPD109 <- plotLogMortalitySplitData(
  GompertzGPD109Analysis,
  title = title
)
plotGompertzGPD109
resGompertzGPD109 <- plotResMortality(
  GompertzGPD109Analysis,
  ylim = 2,
  title = title
)
resGompertzGPD109


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
MakehamGPD109Analysis <- addFittedColumnTo109(MakehamGPD109Analysis)

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
hermiteV109Analysis <- addFittedColumnTo109(hermiteV109Analysis)

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

CohortPatchwork <- plotMakehamBeard109 + plotHermiteII109 + plotMakehamGPD109 + plotHermiteV109
CohortPatchworkRes <- resMakehamBeard109 + resHermiteII109 + resMakehamGPD109 + resHermiteV109
CohortPatchwork
CohortPatchworkRes

title <- paste(LAND, COHORT, GENDER, "Patchwork.png", sep = "")
titleRes <- paste(LAND, COHORT, GENDER, "PatchworkRes.png", sep = "")
# ggsave(title, plot = CohortPatchwork, width = 30, height = 20, units = "cm", path = "plots")
# ggsave(titleRes, plot = CohortPatchworkRes, width = 30, height = 20, units = "cm", path = "plots")



MakehamBeardModel109Analysis$modelResults$AIC
hermiteII109Analysis$modelResults$AIC
MakehamGPD109Analysis$modelResults$AIC
hermiteV109Analysis$modelResults$AIC










