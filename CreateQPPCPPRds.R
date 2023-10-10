library(tidyverse)

Exposures <- read.table('data/ExposuresCPPundQPP_forHMD.txt',sep='\t',header=TRUE,fileEncoding = "UTF-8")
Deaths <- read.table('data/DeathratesCPPundQPP_forHMD.txt',sep='\t',header=TRUE,fileEncoding = "UTF-8")

Exposures <- Exposures %>%
  select(-Total)

ExposuresMale <- Exposures %>%
  select(-Female) %>%
  mutate(Gender = "Male") %>%
  rename(Exposure = Male)

ExposuresFemale <- Exposures %>%
  select(-Male) %>%
  mutate(Gender = "Female") %>%
  rename(Exposure = Female)

Exposures <- bind_rows(ExposuresMale, ExposuresFemale) %>%
  group_by(Year, Age, Gender, Exposure)

remove(ExposuresMale, ExposuresFemale)

Deaths <- Deaths %>%
  select(-Total)

DeathsMale <- Deaths %>%
  select(-Female) %>%
  mutate(Gender = "Male") %>%
  rename(DeathRates = Male)

DeathsFemale <- Deaths %>%
  select(-Male) %>%
  mutate(Gender = "Female") %>%
  rename(DeathRates = Female)

Deaths <- bind_rows(DeathsMale, DeathsFemale) %>%
  group_by(Year, Age, Gender, DeathRates)

remove(DeathsMale, DeathsFemale)

Combined <- inner_join(Exposures, Deaths) %>%
  filter(Exposure != "." | DeathRates != ".")

Combined <- Combined %>%
  mutate_at(c("Exposure", "DeathRates"), as.numeric) %>%
  mutate(Deaths = Exposure * (1 - exp(-DeathRates)))

QPP_CPP_DF <- Combined %>%
  mutate(Cohort = Year - Age) %>%
  select(-DeathRates)

save(QPP_CPP_DF, file = "data/CPP_QPP.Rda")









