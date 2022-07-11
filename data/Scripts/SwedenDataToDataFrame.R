library(dplyr)


load("data/CrudeAgeCohort_F_annual.rda")
load("data/CrudeAgeCohort_M_annual.rda")


formatGroupedDF <- function(df){
  select(df, -hazardRates) %>% 
    rename("Age" = "time", "Exposure" = "sumExposure", "Deaths" = "sumDecrement", "Cohort" = "category")
}

data_F <- formatGroupedDF(CrudeAgeCohort_F_annual)
data_M <- formatGroupedDF(CrudeAgeCohort_M_annual)

data_F$Cohort <- as.double(data_F$Cohort)
data_F$Age <- as.double(data_F$Age)

data_M$Cohort <- as.double(data_M$Cohort)
data_M$Age <- as.double(data_M$Age)


data_F <- data_F %>%
  mutate(Gender = "Female", Year = Cohort + Age)
data_M <- data_M %>%
  mutate(Gender = "Male", Year = Cohort + Age)

Sweden_Df <- bind_rows(data_F, data_M) %>%
  select(Age, Gender, Year, Exposure, Deaths, Cohort)


Sweden_Df <- Sweden_Df %>%
  filter(Age <= 80 & Exposure != 0 | Age > 80)



save(Sweden_Df, file = "data/Sweden_Df_Year.Rda")
