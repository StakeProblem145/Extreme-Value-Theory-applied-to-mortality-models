library(tidyverse)

source("graduatePoisson.R")

load("HMD_df.Rda")
load("OAS_df.Rda")
load("CPP_QPP.Rda")



#### Prep ####

# In general Age above 60
# Add cohort Column

HMD_df <- filter(HMD_df, Age >= 60) %>%
  mutate(Cohort = Year - Age)
OAS_df <- filter(OAS_df, Age >= 60) %>%
  mutate(Cohort = Year - Age)
CPP_df <- filter(CPP_df, Age >= 60) %>%
  mutate(Cohort = Year - Age)
QPP_df <- filter(QPP_df, Age >= 60) %>%
  mutate(Cohort = Year - Age)


#### Aggregate All ####

HMDAggregated <- HMD_df %>%
  filter(Cohort == 1900) %>%
  group_by(Age) %>%
  summarise(Deaths = sum(Deaths), Exposure = sum(Exposure))


resultGomp <- graduateGompertz(HMDAggregated)
ggplot(data = resultGomp$PlottingData) +
  geom_point(aes(x = Age, y = obs)) +
  geom_line(aes(x = Age, y = mod))


result <- graduateHermite(HMDAggregated, "V")
ggplot(data = result$PlottingData) +
  geom_point(aes(x = Age, y = obs)) +
  geom_line(aes(x = Age, y = mod))

result <- graduateMakehamBeard(HMDAggregated)
ggplot(data = result$PlottingData) +
  geom_point(aes(x = Age, y = obs)) +
  geom_line(aes(x = Age, y = mod))


for(tA in c(85:95)){
  result <- graduateMakehamGPD(HMDAggregated, thresholdAge = tA)
  print(tA)
  print(result$ModelResults$AIC)
  plot <- ggplot(data = result$PlottingData) +
    geom_point(aes(x = Age, y = obs)) +
    geom_line(aes(x = Age, y = mod))
  print(plot)
  Sys.sleep(1)
}
