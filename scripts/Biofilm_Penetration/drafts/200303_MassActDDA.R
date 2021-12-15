# Files ----
library(tidyverse)

#Loading in data ----
BBA_Info<- read_csv(here("data_in", "Biofilm_Penetration", "200302_BBA_StrainNames.csv")) %>% rename(Strain='A Number')
BA_MIC_activity_biomass <- read_csv(here("data_out", "Biofilm_Penetration", "210303BA_MIC_biomass_activity.csv")) %>% 
  merge(BBA_Info) %>% select(!Strain) %>%  rename(Strain='yAG Number')%>%
  mutate(Drug="BA") %>% filter(Enviro=="3.12")

write_csv(BA_MIC_activity_biomass,here("data_out", "Biofilm_Penetration", "210303BA_MIC_biomass_activity_edited.csv"))

FLC_MIC_activity_biomass <- read_csv(here("data_out", 
                                          "Biofilm_Penetration", 
                                          "210303FLC_MIC_biomass_activity.csv")) %>% 
  merge(BBA_Info) %>% select(!Strain) %>%
  rename(Strain='yAG Number') %>%
  mutate(Drug="FLC") %>% 
  filter(Enviro ==8)

write_csv(FLC_MIC_activity_biomass,here("data_out", "Biofilm_Penetration", "210303FLC_MIC_biomass_activity_edited.csv"))


MassAct<-rbind(BA_MIC_activity_biomass, FLC_MIC_activity_biomass)

BA_DDA <- read_csv(here("data_out", "DDA", "2020_DDAComb.csv"))


MassActDDA<-merge(MassAct, BA_DDA, by= c("Strain", "Drug")) %>% select(!Type)
write_csv(MassActDDA,here("data_out", "Biofilm_Penetration", "210303MassActDDA.csv"))




