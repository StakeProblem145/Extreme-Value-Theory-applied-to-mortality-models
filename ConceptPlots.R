library(tidyverse)
library(patchwork)

source("graduatePoisson.R")
source("PlottingFunctions.R")

load("data/Canada_HMD_df.Rda")

#### Prep ####
# Valid Cohorts for Canada 1890 to 1910
DF <- Canada_HMD_df

minAge <- 60
maxAge <- 104
maxAgeExtra <- 109

hmdDfMaxAge104 <- filter(DF, Age >= minAge & Age <= maxAge & Year <= 2020) %>%
  mutate(Cohort = Year - Age)

hmdDfMaxAge109 <-
  filter(DF, Age >= minAge & Age <= maxAgeExtra & Year <= 2020) %>%
  mutate(Cohort = Year - Age)

#### NOTES ####
# Canada
# Female
X1HermiteII <- 119
thresholdAgeGomp <- 99
thresholdAgeMake <- 98
X1HermiteV <- 111

#### Aggregated Data 1890 to 1910 ####
##### Data #####
LAND <- "Canada"
cohortMin <- 1905
cohortMax <- 1905
COHORT <- "1905"
GENDER <- "Female"

df104 <- hmdDfMaxAge104 %>%
  filter(Cohort >= cohortMin & Cohort <= cohortMax & Gender == GENDER) %>%
  group_by(Age) %>%
  summarise(Deaths = sum(Deaths), Exposure = sum(Exposure))

df109 <- hmdDfMaxAge109 %>%
  filter(Cohort >= cohortMin & Cohort <= cohortMax & Gender == GENDER) %>%
  group_by(Age) %>%
  summarise(Deaths = sum(Deaths), Exposure = sum(Exposure))



#### Mortality ####
plot <- ggplot(GompertzModelAna$plottingData, aes(x = Age)) +
  geom_point(aes(y = exp(obs)), alpha = 7/10) +
  xlim(60,110) +
  ylim(NA, 0.8) +
  ggtitle(paste("Mortality", COHORT, GENDER)) +
  #geom_col(width=0.8)+
  #ggtitle(paste(title)) +
  xlab("Age") +
  ylab("Mortality") +
  theme_clean() +
  MODEL_PLOT_THEME
plot
# ggsave("MortalityConceptPlot.png", plot = plot, width = 30, height = 20, units = "cm", path = "plots")


#### Gompertz ####
GompertzModel <-
  graduateGompertz(df109, analysis = FALSE)
GompertzModelAna <-
  createGradAnalysis(GompertzModel)
GompertzModelAna$modelResults$AIC


title <- paste("Gompertz", COHORT, GENDER)
GompertzPlot <- plotLogMortality(GompertzModelAna,
                 xlim = c(60,110),
                 ylim = c(-5, 0.5),
                 title = title,
                 errorBar = TRUE)
GompertzPlotRes <- plotResMortality(
  GompertzModelAna,
  xlim = c(60,110),
  ylim = 1,
  title = title
)

GompertzPlot
GompertzPlotRes
# ggsave("GompertzConceptPlot.png", plot = GompertzPlot, width = 30, height = 20, units = "cm", path = "plots")
# ggsave("GompertzConceptPlotRes.png", plot = GompertzPlotRes, width = 30, height = 20, units = "cm", path = "plots")


#### Makeham ####
MakehamModel <-
  graduateMakeham(df109, analysis = FALSE)
MakehamModelAna <-
  createGradAnalysis(MakehamModel)
MakehamModelAna$modelResults$AIC


title <- paste("Makeham", COHORT, GENDER)
MakehamPlot <- plotLogMortality(MakehamModelAna,
                 xlim = c(60,110),
                 ylim = c(NA, 0.5),
                 title = title,
                 errorBar = TRUE)
MakehamPlotRes <- plotResMortality(
  MakehamModelAna,
  xlim = c(60,110),
  ylim = 1,
  title = title
)


#### Perks ####
PerksModel <-
  graduatePerks(df109, analysis = FALSE)
PerksModelAna <-
  createGradAnalysis(PerksModel)
PerksModelAna$modelResults$AIC


title <- paste("Perks", COHORT, GENDER)
PerksPlot <- plotLogMortality(PerksModelAna,
                 xlim = c(60,110),
                 ylim = c(NA, 0.5),
                 title = title,
                 errorBar = TRUE)
PerksPlotRes <- plotResMortality(
  PerksModelAna,
  xlim = c(60,110),
  ylim = 1,
  title = title
)


#### Beard ####
BeardModel <-
  graduateBeard(df109, analysis = FALSE)
BeardModelAna <-
  createGradAnalysis(BeardModel)
BeardModelAna$modelResults$AIC


title <- paste("Beard", COHORT, GENDER)
BeardPlot <- plotLogMortality(BeardModelAna,
                              xlim = c(60,110),
                 ylim = c(NA, 0.5),
                 title = title,
                 errorBar = TRUE)
BeardPlotRes <- plotResMortality(
  BeardModelAna,
  xlim = c(60,110),
  ylim = 1,
  title = title
)


#### MakehamPerks ####
MakehamPerksModel <-
  graduateMakehamPerks(df109, analysis = FALSE)
MakehamPerksModelAna <-
  createGradAnalysis(MakehamPerksModel)
MakehamPerksModelAna$modelResults$AIC


title <- paste("Makeham Perks", COHORT, GENDER)
MakehamPerksPlot <- plotLogMortality(MakehamPerksModelAna,
                                     xlim = c(60,110),
                 ylim = c(NA, 0.5),
                 title = title,
                 errorBar = TRUE)
MakehamPerksPlotRes <- plotResMortality(
  MakehamPerksModelAna,
  xlim = c(60,110),
  ylim = 1,
  title = title
)



MakehamToMakehamPerksConceptPlot <- MakehamPlot + PerksPlot + BeardPlot + MakehamPerksPlot
# ggsave("MakehamToMakehamPerksConceptPlot.png", plot = MakehamToMakehamPerksConceptPlot, width = 30, height = 20, units = "cm", path = "plots")

MakehamToMakehamPerksConceptPlotRes <- MakehamPlotRes + PerksPlotRes + BeardPlotRes + MakehamPerksPlotRes
# ggsave("MakehamToMakehamPerksConceptPlotRes.png", plot = MakehamToMakehamPerksConceptPlotRes, width = 30, height = 20, units = "cm", path = "plots")




#### MakehamBeard ####
MakehamBeardModel <-
  graduateMakehamBeard(df109, analysis = FALSE)
MakehamBeardModelAna <-
  createGradAnalysis(MakehamBeardModel)
MakehamBeardModelAna$modelResults$AIC


title <- paste("Makeham Beard", COHORT, GENDER)
MakehamBeardPlot <- plotLogMortality(MakehamBeardModelAna,
                                     xlim = c(60,110),
                 ylim = c(NA, 0.5),
                 title = title,
                 errorBar = TRUE)
MakehamBeardPlotRes <- plotResMortality(
  MakehamBeardModelAna,
  xlim = c(60,110),
  ylim = 1,
  title = title
)
MakehamBeardPlot
MakehamBeardPlotRes
# ggsave("MakehamBeardConceptPlot.png", plot = MakehamBeardPlot, width = 30, height = 20, units = "cm", path = "plots")
# ggsave("MakehamBeardConceptPlotRes.png", plot = MakehamBeardPlotRes, width = 30, height = 20, units = "cm", path = "plots")



GompertzModelAna$modelResults$AIC
MakehamModelAna$modelResults$AIC
PerksModelAna$modelResults$AIC
BeardModelAna$modelResults$AIC
MakehamPerksModelAna$modelResults$AIC
MakehamBeardModelAna$modelResults$AIC

### Mega Plot
megaPlotDf <- bind_rows(MakehamModelAna$plottingData,
                        PerksModelAna$plottingData,
                        BeardModelAna$plottingData,
                        MakehamPerksModelAna$plottingData)
megaPlotDf$fitted <- "Fit"
title <- paste(
  "Models", COHORT, GENDER
)
megaPlot <- plotLogMortalitySplitData(
  megaPlotDf,
  title = title,
  xlim = c(60,110),
  ylim = c(NA, 0.5),
  hideLegend = FALSE
)
megaPlot
# ggsave("MakehamToMakehamPerksConceptPlot.png", plot = megaPlot, width = 30, height = 20, units = "cm", path = "plots")




#### HERMITE ####

scaleAge <- function(x, x0, x1) {
  (x - x0) / (x1 - x0)
}
h00x <- function(x, x0, x1) {
  t <- scaleAge(x, x0, x1)
  return((1.0 + 2 * t) * (1.0 - t) * (1.0 - t))
}
h10x <- function(x, x0, x1) {
  t <- scaleAge(x, x0, x1)
  return(t * (1.0 - t) * (1.0 - t))
}
h01x <- function(x, x0, x1) {
  t <- scaleAge(x, x0, x1)
  return(t * t * (3.0 - 2 * t))
}
h11x <- function(x, x0, x1) {
  t <- scaleAge(x, x0, x1)
  return(t * t * (t - 1.0))
}
hQuarticx <- function(x, x0, x1) {
  t <- scaleAge(x, x0, x1)
  return(16 * t ^ 2 * (1 - t) ^ 2)
}

x <- seq(0, 1, length.out = 101)

X0 = 0
X1 = 1
h00 <- h00x(x, X0, X1)
h01 <- h01x(x, X0, X1)
h10 <- h10x(x, X0, X1)
h11 <- h11x(x, X0, X1)
hQuartic <- hQuarticx(x, X0, X1)
  
HermitePolynoms <- ggplot() +
  geom_line(aes(x = x, y = h00, linetype = "h00")) +
  geom_line(aes(x = x, y = h01, linetype = "h01")) +
  geom_line(aes(x = x, y = h10, linetype = "h10")) +
  geom_line(aes(x = x, y = h11, linetype = "h11")) +
  geom_line(aes(x = x, y = hQuartic, linetype = "hQuartic")) +
  ggtitle("Hermite Polynoms") +
  labs(linetype='Hermite Polynom') +
  theme_clean() +
  MODEL_PLOT_THEME

# ggsave("HermitePolynoms.png", plot = HermitePolynoms, width = 30, height = 20, units = "cm", path = "plots")





#### Hermite I and III ####
HermiteIModel <-
  graduateHermite(df109, "I")
HermiteIIIModel <-
  graduateHermite(df109, "III")
HermiteIModel$modelResults$AIC
HermiteIIIModel$modelResults$AIC

#### Hermite II  ####
HermiteIIModel <-
  graduateHermite(df109, "II", X1 = 115, analysis = FALSE)
HermiteIIModelAna <-
  createGradAnalysis(HermiteIIModel)
HermiteIIModelAna$modelResults$AIC


title <- paste("Hermite II", COHORT, GENDER)
HermiteIIPlot <- plotLogMortality(HermiteIIModelAna,
                                  xlim = c(60,110),
                                     ylim = c(NA, 0.5),
                                     title = title,
                                     errorBar = TRUE)
HermiteIIPlotRes <- plotResMortality(
  HermiteIIModelAna,
  xlim = c(60,110),
  ylim = 1,
  title = title
)
HermiteIIPlot
HermiteIIPlotRes


#### Hermite IV  ####
HermiteIVModel <-
  graduateHermite(df109, "IV", analysis = FALSE)
HermiteIVModelAna <-
  createGradAnalysis(HermiteIVModel)
HermiteIVModelAna$modelResults$AIC


title <- paste("Hermite IV", COHORT, GENDER)
HermiteIVPlot <- plotLogMortality(HermiteIVModelAna,
                                  xlim = c(60,110),
                                  ylim = c(NA, 0.5),
                                  title = title,
                                  errorBar = TRUE)
HermiteIVPlotRes <- plotResMortality(
  HermiteIVModelAna,
  xlim = c(60,110),
  ylim = 1,
  title = title
)
HermiteIVPlot
HermiteIVPlotRes

HermiteIIHermiteIV <- HermiteIIPlot + HermiteIVPlot + HermiteIIPlotRes + HermiteIVPlotRes
# ggsave("HermiteIIHermiteIVConcept.png", plot = HermiteIIHermiteIV, width = 30, height = 20, units = "cm", path = "plots")



#### Hermite V  ####
HermiteVModel <-
  graduateHermite(df109, "V", X1=110, analysis = FALSE)
HermiteVModelAna <-
  createGradAnalysis(HermiteVModel)
HermiteVModelAna$modelResults$AIC


title <- paste("Hermite V", COHORT, GENDER)
HermiteVPlot <- plotLogMortality(HermiteVModelAna,
                                 xlim = c(60,110),
                                  ylim = c(NA, 0.5),
                                  title = title,
                                  errorBar = TRUE)
HermiteVPlotRes <- plotResMortality(
  HermiteVModelAna,
  xlim = c(60,110),
  ylim = 1,
  title = title
)
HermiteVPlot
HermiteVPlotRes
# ggsave("HermiteVConceptPlot.png", plot = HermiteVPlot, width = 30, height = 20, units = "cm", path = "plots")
# ggsave("HermiteVConceptPlotRes.png", plot = HermiteVPlotRes, width = 30, height = 20, units = "cm", path = "plots")


