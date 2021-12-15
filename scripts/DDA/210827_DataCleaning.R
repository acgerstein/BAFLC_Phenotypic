library(here)
library(tidyverse)
library(ggforce)
library(ggpubr)
library(knitr)
library(kableExtra)
library(webshot)
library(extrafont)
library(viridis)
library(plotly)


# Upload data ----

#* HSC001-105 strain information (species, anatomical site) ----

H1Info <- read_csv(here("data_in", "DDA","2014_HSC_info.csv")) %>% 
  select(!X1) %>% 
  rename("Strain" = "line") %>% mutate(Year="2013")

#* HSC001-105 biological replicate 1 (B1), Boric Acid (BA)
H1B1BA <- read_csv(here("data_in", "DDA","HSC001-105_48h_B1_df.csv")) %>%
  rename("Strain" = "line") %>% mutate(Drug="BA") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  merge(H1Info) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis", "C.tropicalis")) %>%group_by(Strain) %>% mutate(Rep=1:length(Strain)) 

# Check the quality of technical replicates

ggplotly(H1B1BA %>% ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
  scale_y_continuous(limits = c(0, 25)) + scale_x_discrete(label = seq(0, 100)))

ggplotly(H1B1BA %>% ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
  scale_y_continuous(limits = c(0, 1)) + scale_x_discrete(label = seq(0, 100)))

# Take average of technical replicates
H1B1BA <- read_csv(here("data_in", "DDA","HSC001-105_48h_B1_df.csv")) %>%
  rename("Strain" = "line") %>% mutate(Drug="BA") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain)) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) %>% 
  merge(H1Info) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis", "C.tropicalis")) 

#* HSC001-105 biological replicate 1 (B1), Fluconazole 
H1B1FLC <- read_csv(here("data_in", "DDA","dda_200309_Abdul_FLC_48h_B1_df.csv")) %>%
  rename("Strain" = "line") %>% mutate(Drug="FLC") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  subset(Strain%in%H1B1BA$Strain)%>% 
  merge(H1Info) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis",
                      "C.tropicalis")) %>% 
  group_by(Strain) %>% 
  mutate(Rep=1:length(Strain)) 

# Check the quality of technical Replicates
ggplotly(H1B1FLC %>% ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
                    scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,100)))
# RAD20 of HSC026 & HSC059 Rep 2 does not agree with the DDA image -> Remove

ggplotly(H1B1FLC %>% ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
  scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,100)))

# Take average of technical replicates
H1B1FLC <- read_csv(here("data_in", "DDA","dda_200309_Abdul_FLC_48h_B1_df.csv")) %>%
  rename("Strain" = "line") %>% mutate(Drug="FLC") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  subset(Strain%in%H1B1BA$Strain)%>% 
  merge(H1Info) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis",
                      "C.tropicalis")) %>% 
  group_by(Strain) %>% 
  mutate(Rep=1:length(Strain)) %>% 
  filter(!(Strain %in% c("HSC026", "HSC059") & Rep == 2)) %>% 
  group_by(Strain, Drug,Species, Origin, Year ) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) 

#* HSC001-105 biorep 1, HSC1-105, Boric Acid + Fluconazole ----
H1B1 <- rbind(H1B1BA, H1B1FLC) %>% mutate(Biorep="B1")


#* HSC001-105 biorep 2, Boric Acid & Fluconazole ----
H1B2 <- read_csv(here("data_in", "DDA","HSC001-105_48h_B2_df.csv")) %>%
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% filter(!Strain=="Blank") %>% 
  group_by(Strain, Drug) %>% 
  subset(Strain%in%H1B1$Strain)%>% 
  merge(H1Info) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata",
                      "C.parapsilosis", "C.tropicalis")) %>% 
  mutate(Biorep="B2") %>% 
  group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain)) 

# Check the quality of technical replicates
ggplotly(H1B2 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,100)))

ggplotly(H1B2 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,100)))
# RAD20 HSC091 Rep1 does not agree with the DDA image -> Remove

ggplotly(H1B2 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,100)))

ggplotly(H1B2 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,100)))

# Take average of technical replicates
H1B2 <- read_csv(here("data_in", "DDA","HSC001-105_48h_B2_df.csv")) %>%
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% filter(!Strain=="Blank") %>% 
  group_by(Strain, Drug) %>% 
  subset(Strain%in%H1B1$Strain)%>% 
  merge(H1Info) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata",
                      "C.parapsilosis", "C.tropicalis")) %>% 
  group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain)) %>%
  group_by(Strain, Drug,Species, Origin, Year ) %>% 
  filter(!(Strain == "HSC091" & Rep == 2)& Drug == "FLC") %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) %>% 
  mutate(Biorep="B2") 

# Only keep strains that have two bioreps 
H1B1 <- H1B1 %>% subset(Strain%in%H1B2$Strain)
H1B2 <- H1B2 %>% subset(Strain%in%H1B1$Strain)

#* All HSC001-105, two bioreps, BA & FLC ----
H1_df <- rbind(H1B1, H1B2)

##########################
#* HSC106-210 strain information (species, anatomical site) ----
H2Info <- read_csv(here("data_in", "DDA","Aleeza_HSC_Yeast stocks_2018.csv" )) %>% 
  select('HSC#', "SPECIMEN", "ORGANISM") %>% 
  rename("Strain"='HSC#', "Species"="ORGANISM", "Origin"="SPECIMEN") %>% mutate(Year="2018")

#* HSC106-150 biorep 1, FLC & BA ----
H2B1pt1 <- read_csv(here("data_in", "DDA","HSC106_159_48h_df_p1_B1.csv")) %>% 
  filter(!line=="Blank") %>%
  rename("Strain"="line", "Drug"="type") %>% 
  select(Strain, Drug, RAD20, FoG20, RAD50, FoG50) %>%
  group_by(Strain, Drug) %>% 
  merge(H2Info, by = "Strain") %>% 
  group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain)) 

# Check the quality of technical replicates
ggplotly(H2B1pt1 %>% filter(Drug == "BA") %>%  
           ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0, 25)) + scale_x_discrete(label = seq(0,100)))

ggplotly(H2B1pt1 %>% filter(Drug == "FLC") %>%  
           ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0, 25)) + scale_x_discrete(label = seq(0,100)))

ggplotly(H2B1pt1 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,100)))

ggplotly(H2B1pt1 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,100)))

# Take average of technical replicates
H2B1pt1 <- read_csv(here("data_in", "DDA","HSC106_159_48h_df_p1_B1.csv")) %>% 
  filter(!line=="Blank") %>%
  rename("Strain"="line", "Drug"="type") %>% 
  select(Strain, Drug, RAD20, FoG20, RAD50, FoG50) %>%
  group_by(Strain, Drug) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50), Biorep ="B1" ) %>% 
  merge(H2Info, by = "Strain") %>% 
  filter(!(Strain == "HSC113" & Drug == "BA"))

#* HSC151-210 biorep 1, FLC & BA ----
H2B1pt2 <- read_csv(here("data_in", "DDA","HSC106_159_48h_df_p2_B1.csv")) %>% 
  filter(!line=="Blank") %>%
  rename("Strain"="line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% 
  group_by(Strain, Drug) %>% 
  merge(H2Info, by = "Strain")%>% 
  group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain))

# Check the quality of the technical replicates
ggplotly(H2B1pt2 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,100)))

ggplotly(H2B1pt2 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,100)))
# RAD20 of HSC166 FLC Rep 1 does not agree with the DDA image -> Remove

ggplotly(H2B1pt2 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,100)))

ggplotly(H2B1pt2 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,100)))


# Take average of technical replicates
H2B1pt2 <- read_csv(here("data_in", "DDA","HSC106_159_48h_df_p2_B1.csv")) %>% 
  filter(!line=="Blank") %>%
  rename("Strain"="line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% 
  mutate(Rep=1:length(Strain)) %>% 
  group_by(Strain, Drug) %>% 
  filter(!(Strain == "HSC166" & Rep == 1 & Drug == "FLC")) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50), Biorep ="B1" ) %>% 
  merge(H2Info, by = "Strain")

# HSC106-210 biorep 1, FLC & BA
H2B1 <- rbind(H2B1pt1, H2B1pt2) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis", "C.tropicalis"))

#* HSC106-210 biorep 2, FLC & BA ----
H2B2 <- read_csv(here("data_in", "DDA","HSC106_205_48h_df_B2.csv")) %>% 
  filter(!line=="Blank") %>%
  rename("Strain"="line", "Drug"="type") %>% 
  select(Strain, Drug, RAD20, FoG20, RAD50, FoG50) %>% 
  group_by(Strain, Drug) %>% 
  merge(H2Info, by = "Strain") %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis", 
                      "C.tropicalis")) %>% 
  filter(!Strain == "HSC167") %>% 
  group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain))

# Check the quality of the technical Reps
ggplotly(H2B2 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,100)))
# RAD20 of HSC159 Rep 1 BA does not agree with the DDA image -> Remove

ggplotly(H2B2 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,100)))
# RAD20 of HSC131 Rep 2 FLC does not agree with the DDA image -> Remove

ggplotly(H2B2 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,100)))

ggplotly(H2B2 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,100)))

# Take average of technical replicates
H2B2 <- read_csv(here("data_in", "DDA","HSC106_205_48h_df_B2.csv")) %>% 
  filter(!line=="Blank") %>%
  rename("Strain"="line", "Drug"="type") %>% 
  select(Strain, Drug, RAD20, FoG20, RAD50, FoG50) %>% 
  group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain)) %>% 
  filter(!(Strain == "HSC131" &  Drug == "FLC" &Rep == 2)) %>% 
  filter(!(Strain == "HSC159" &  Drug == "FLC" &Rep == 1)) %>% 
  group_by(Strain, Drug) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50), Biorep ="B2" ) %>% 
  merge(H2Info, by = "Strain") %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis", 
                      "C.tropicalis"))

#* HSC106-210 biorep 1&2, FLC & BA ----
H2_df <- rbind(H2B1, H2B2) 
 
#* HSC001-210 biorep 1&2, FLC & BA ----
H_df <- rbind(H1_df, H2_df) 

#* Biorep 3 - HSC years 1 & 2
alb_glb_B03 <- read_csv(here("data_in","DDA", "DDA201218_48h_df.csv")) %>%
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% 
  filter(!Strain=="Blank") %>% 
  group_by(Strain) %>% 
  mutate(Rep=1:length(Strain)) 

# Check the quality of the technical replicates
ggplotly(alb_glb_B03 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,1000)))

ggplotly(alb_glb_B03 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,1000)))

ggplotly(alb_glb_B03 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,1000)))

ggplotly(alb_glb_B03 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,1000)))

# Take average of technical replicates
alb_glb_B03 <- read_csv(here("data_in","DDA", "DDA201218_48h_df.csv")) %>%
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% 
  filter(!Strain=="Blank") %>% 
  group_by(Strain) %>% 
  mutate(Rep=1:length(Strain)) %>% 
  group_by(Strain, Drug) %>%
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50),
            FoG50 = mean(FoG50)) %>% 
  filter(!(Strain == "HSC105")) %>% 
  filter(!(Strain == "HSC131"))
#HSC105 & HSC131 of this biorep looks different from the other 2 biological reps -> Remove

para_trop_B03 <- read_csv(here("data_in","DDA", "Para_Trop_B03_df.csv")) %>% 
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% filter(!Strain=="Blank") %>%  group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain)) 

# Check the quality of the technical replicates
ggplotly(para_trop_B03 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,1000)))

ggplotly(para_trop_B03 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,1000)))
# RAD20 of HSC184 Rep 2 does not agree with DDA image -> Remove
ggplotly(para_trop_B03 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,1000)))

ggplotly(para_trop_B03 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,1000)))

# Take average of technical replicates
para_trop_B03 <- read_csv(here("data_in","DDA", "Para_Trop_B03_df.csv")) %>% 
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% filter(!Strain=="Blank") %>%  group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain)) %>% 
  filter(!(Strain == "HSC184" & Rep == 2 & Drug == "FLC")) %>% 
group_by(Strain, Drug) %>%
summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50) )

#* HSC001-210 biorep 3, FLC & BA ----
Info <- rbind(H1Info, H2Info)
HB3 <- rbind(alb_glb_B03, para_trop_B03) %>% merge(Info, by = "Strain") %>% mutate(Biorep = "B3")


#* HSC001-210 biorep 1, 2 & 3, FLC & BA ----
H_wide <- rbind(H1B1, H2B1) %>% 
 merge(rbind(H1B2, H2B2), by = c("Strain","Drug","Species", "Origin", "Year")) %>% 
  rename("RAD20_B01" = "RAD20.x", "RAD20_B02" = "RAD20.y",
         "RAD50_B01" = "RAD50.x", "RAD50_B02" = "RAD50.y",
         "FoG20_B01" = "FoG20.x", 
         "FoG50_B01" = "FoG50.x", 
         "FoG20_B02" = "FoG20.y",
         "FoG50_B02" = "FoG50.y") %>% 
  select(!c(Biorep.x,Biorep.y)) %>% 
  merge(HB3, by = c("Strain","Drug")) %>% 
  rename("RAD20_B03" = "RAD20", "FoG20_B03" = "FoG20", "RAD50_B03" = "RAD50", "FoG50_B03" = "FoG50")

write_csv(H_wide, here("data_out", "DDA","210827_HSC_B1B2B3.csv"))


HSC <- rbind(H_df, HB3) %>% 
  group_by(Strain, Drug, Species) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50) )

# Check whether the biological replicates correlate
cor.test(H_wide$RAD20_B01, H_wide$RAD20_B02)
cor.test(H_wide$RAD20_B02, H_wide$RAD20_B03)
cor.test(H_wide$RAD20_B01, H_wide$RAD20_B03)
cor.test(H_wide$FoG20_B01, H_wide$FoG20_B02)
cor.test(H_wide$FoG20_B01, H_wide$FoG20_B03)

# Number os HSC strains 
HSC %>% filter(Drug=="BA") %>% 
  group_by(Species) %>% 
  summarise(n=n()) %>% 
  kable("html") %>% 
  kable_styling(latex_options = c("striped", "scale_down"))

#85 c.alb, 50 c. glb, 20 c.para, 7 c. trop

#* AG strains ----
AGInfo <- read_csv(here("data_in", "DDA","BermanStrains_GersteinInventory.csv")) %>% 
  select('AG strain', Species, 'Strain Name:') %>% 
  rename("Strain" = 'AG strain', "Strain_Name"='Strain Name:')

AGInfo$Strain <- sub("^", "yAG", AGInfo$Strain )
AGInfo$Species <- sub(" ", "", AGInfo$Species)

#* AG strains biorep 1 FLC & BA ----
AGB1 <- read_csv(here("data_in", "DDA","AG_201001_48h_B1_df.csv")) %>%
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain, Drug, RAD20, FoG20, RAD50, FoG50) %>% 
  filter(!Strain == "Blank") %>% 
  group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain))

# Check the quality of the technical replicates
ggplotly(AGB1 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,1000)))

ggplotly(AGB1 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,1000)))

ggplotly(AGB1 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,1000)))

ggplotly(AGB1 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,1000)))

#Strains that have been filtered do not match the image
# Take average of technical replicates
AGB1 <- read_csv(here("data_in", "DDA","AG_201001_48h_B1_df.csv")) %>%
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  filter(!Strain=="Blank") %>% 
  group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain)) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) %>% 
  merge(AGInfo) %>% mutate(Biorep = "B1") 
  

#* AG strains biorep 2 FLC & BA ----
AGB2 <- read_csv(here("data_in", "DDA","AG_201008_48h_B2_df.csv")) %>%
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  
  filter(!Strain=="Blank") %>% 
  group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain)) 

# Check the quality of the technical replicates
ggplotly(AGB2 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,1000)))

ggplotly(AGB2 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,1000)))
# RAD20 yAG97 Rep 2 does not agree with the DDA image -> Remove

ggplotly(AGB2 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,1000)))
# RAD20 yAG77 Rep 2 does not agree with the DDA image -> Remove

ggplotly(AGB2 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,1000)))

# Take average of technical replicates
AGB2 <- read_csv(here("data_in", "DDA","AG_201008_48h_B2_df.csv")) %>%
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  
  filter(!Strain=="Blank") %>% 
  group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain)) %>% 
  filter(!(Strain == "yAG97" & Rep == 2 & Drug == "FLC")) %>%
  filter(!(Strain == "yAG77" & Rep == 2 & Drug == "BA")) %>% 
  group_by(Strain, Drug) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20) , RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) %>% 
  merge(AGInfo) %>% mutate(Biorep = "B2") 

# Identify strains missing in the biological replicates
AGB1 %>% filter(!Strain %in% AGB2$Strain)
AGB2 %>% filter(!Strain %in% AGB1$Strain)
#yAG92 is not in B2
#yAG40 and 55 is missing in B1

# Filter out the missing strains
AGB1 <- AGB1 %>% filter(Strain %in% AGB2$Strain)
AGB2 <- AGB2 %>% filter(Strain %in% AGB1$Strain)

# Checkpoint
AGB1 %>% filter(!Strain %in% AGB2$Strain)
AGB2 %>% filter(!Strain %in% AGB1$Strain)

AG <- rbind(AGB1,AGB2)  

#* Missing strains biorep 1 & 2 FLC & BA ----
MissingStrains <- read_csv(here("data_in", "DDA","MissingStrains_df.csv")) %>% 
  filter(!line == "Blank") %>% 
  select(name, RAD20, FoG20, RAD50, FoG50) %>%
  separate(name, c("Strain", "Drug", "Rep", "Biorep"), sep = "_") %>% 
  group_by(Strain, Drug, Biorep) %>% 
  mutate(Rep=1:length(Strain))

# Check the quality of the technical replicates
ggplotly(MissingStrains %>% filter(Drug == "BA" & Biorep == "B01") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,1000)))

ggplotly(MissingStrains %>% filter(Drug == "BA" & Biorep == "B02") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,1000)))

ggplotly(MissingStrains %>% filter(Drug == "FLC"& Biorep == "B01") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,1000)))
# RAD20 yAG054 Rep 1 B01 does not agree with the DDA image -> Remove
ggplotly(MissingStrains %>% filter(Drug == "FLC"& Biorep == "B02") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD20, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,1000)))
# RAD20 yAG054 Rep 2 B02 does not agree with the DDA image -> Remove

ggplotly(MissingStrains %>% filter(Drug == "BA" & Biorep == "B01") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,1000)))

ggplotly(MissingStrains %>% filter(Drug == "BA" & Biorep == "B02") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,1000)))

ggplotly(MissingStrains %>% filter(Drug == "FLC" & Biorep == "B01") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,1000)))

ggplotly(MissingStrains %>% filter(Drug == "FLC" & Biorep == "B02") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG20, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,1000)))


# Take yAG003 out, the culture is too diluted 
# Take yAG046 out, the culture is too diluted 
# Take yAG096 out, the culture is too diluted 

# Take average of technical replicates
MissingStrains <- read_csv(here("data_in", "DDA","MissingStrains_df.csv")) %>% 
  filter(!line == "Blank") %>% 
  select(name, RAD20, FoG20, RAD50, FoG50) %>%
  separate(name, c("Strain", "Drug", "Rep", "Biorep"), sep = "_") %>% 
  group_by(Strain, Drug, Biorep) %>% 
  mutate(Rep=1:length(Strain)) %>% 
  filter(!(Strain %in% c("yAG003", "yAG046", "yAG096"))) %>%
  filter(!(Strain == "yAG054" & Rep == 1 & Biorep == "B01")) %>% 
  filter(!(Strain == "yAG054" & Rep == 2 & Biorep == "B02")) %>% 
  group_by(Strain, Drug) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20) , RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) %>% 
  mutate(strain = Strain) %>% 
  separate(strain, c(NA, "nAG"), sep = "yAG") 

# Add strain info
MissingStrains$nAG <- as.numeric(MissingStrains$nAG)
MissingStrains <- MissingStrains %>% select(!Strain) %>% 
  mutate(AG = "yAG") %>% unite(Strain, c("AG", "nAG"), sep = "" )
MissingStrains <- merge(MissingStrains, AGInfo, by = "Strain") %>% 
  mutate(Biorep = "B3")

AG_full <- rbind(AG, MissingStrains) %>% 
  group_by(Strain, Drug, Species) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50))

write_csv(AG_full, here("data_out", "DDA","210827_AG_B1B2.csv"))

# Number of AG strain
AG_full %>% filter(Drug=="BA") %>% 
  group_by(Species) %>% 
  summarise(n=n()) %>% 
  kable("html") %>% 
  kable_styling(latex_options = c("striped", "scale_down"))

#80 C.alb strains


AGHSC <- HSC %>% rbind(AG_full)
write_csv(AGHSC, here("data_out", "DDA","210827_HSCAG.csv"))

# Number of HSC & AG strain
AGHSC %>% filter(Drug=="BA") %>% 
  group_by(Species) %>% 
  summarise(n=n()) %>% 
  kable("html") %>% 
  kable_styling(latex_options = c("striped", "scale_down"))
# 165 C.alb, 50 C.glb, 20 C.para, 7 C.trop