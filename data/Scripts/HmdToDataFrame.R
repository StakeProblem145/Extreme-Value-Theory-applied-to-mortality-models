library(dplyr)
library(HMDHFDplus)

deaths <- readHMD("data/Canada_Deaths_1x1.txt")
# test <- deaths %>%
#   mutate(Cohort = Year - Age) %>%
#   filter(Age == 81 & Cohort >= 1851 & Cohort <= 1915)
# sum(test$Female)
exposure <- readHMD("data/Canada_Exposures_1x1.txt")
# test <- exposure %>%
#   mutate(Cohort = Year - Age) %>%
#   filter(Age == 81 & Cohort >= 1851 & Cohort <= 1915)
# sum(test$Female)

exposure_F <- exposure %>%
  select(Year, Age, Female) %>%
  rename(Exposure = Female) %>%
  mutate(Gender = "Female")
exposure_M <- exposure %>%
  select(Year, Age, Male) %>%
  rename(Exposure = Male) %>%
  mutate(Gender = "Male")
deaths_F <- deaths %>%
  select(Year, Age, Female) %>%
  rename(Deaths = Female) %>%
  mutate(Gender = "Female")
deaths_M <- deaths %>%
  select(Year, Age, Male) %>%
  rename(Deaths = Male) %>%
  mutate(Gender = "Male")
exposure <- full_join(exposure_M, exposure_F)
deaths <- full_join(deaths_M, deaths_F)
df <- full_join(exposure, deaths) %>%
  relocate(Age, Gender, Year) %>%
  filter(Age >= 13) %>%
  mutate(Cohort = Year - Age)
Canada_HMD_df <- df

save(Canada_HMD_df, file = "data/Canada_HMD_df.Rda")
