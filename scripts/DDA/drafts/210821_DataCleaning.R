library(here)
library(tidyverse)
library(ggforce)
library(ggpubr)
library(knitr)
library(kableExtra)
library(webshot)
library(extrafont)
library(viridis)


# 1 Upload data ----

#* HSC001-105 strain information (species, anatomical site) ----

H1Info<-read_csv(here("data_in", "DDA","2014_HSC_info.csv")) %>% 
  select(!X1) %>% 
  rename("Strain" = "line") %>% mutate(Year="2013")


H1B1BA <- read_csv(here("data_in", "DDA","HSC001-105_48h_B1_df.csv")) %>%
  rename("Strain" = "line") %>% mutate(Drug="BA") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  merge(H1Info) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis", "C.tropicalis")) %>%group_by(Strain) %>% mutate(Rep=1:length(Strain)) 


H1B1BA %>% ggplot() + geom_point(aes(x = Strain, y = RAD50, color = Rep)) +
  scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,100))


H1B1BA %>% ggplot() + geom_point(aes(x = Strain, y = FoG50, color = Rep)) +
  scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,100))


H1B1BA <- read_csv(here("data_in", "DDA","HSC001-105_48h_B1_df.csv")) %>%
  rename("Strain" = "line") %>% mutate(Drug="BA") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain)) %>% 
  filter(!(Strain == "HSC003" & Rep == 1)) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) %>% 
  merge(H1Info) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis", "C.tropicalis")) 


H1B1FLC <- read_csv(here("data_in", "DDA","dda_200309_Abdul_FLC_48h_B1_df.csv")) %>%
  rename("Strain" = "line") %>% mutate(Drug="FLC") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  subset(Strain%in%H1B1BA$Strain)%>% 
  merge(H1Info) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis",
                      "C.tropicalis")) %>% 
  group_by(Strain) %>% 
  mutate(Rep=1:length(Strain)) 

ggplotly(H1B1FLC %>% ggplot() + geom_point(aes(x = Strain, y = RAD50, color = Rep)) +
                    scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,100)))


ggplotly(H1B1FLC %>% ggplot() + geom_point(aes(x = Strain, y = FoG50, color = Rep)) +
  scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,100)))

H1B1FLC <- read_csv(here("data_in", "DDA","dda_200309_Abdul_FLC_48h_B1_df.csv")) %>%
  rename("Strain" = "line") %>% mutate(Drug="FLC") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  subset(Strain%in%H1B1BA$Strain)%>% 
  merge(H1Info) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis",
                      "C.tropicalis")) %>% 
  group_by(Strain) %>% 
  mutate(Rep=1:length(Strain)) %>% 
  filter(!(Strain == "HSC059" & Rep == 2)) %>% 
  filter(!(Strain == "HSC026"& Rep ==1)) %>% 
  filter(!(Strain == "HSC086"& Rep ==1)) %>% 
  group_by(Strain, Drug,Species, Origin, Year ) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) 

#* HSC001-105 biorep 1, HSC1-105, Boric Acid + Fluconazole ----
H1B1 <- rbind(H1B1BA, H1B1FLC) %>% mutate(Biorep="B1")


#* HSC001-105 biorep 2, HSC1-105, Boric Acid & Fluconazole ----
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

ggplotly(H1B2 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD50, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,100)))

ggplotly(H1B2 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG50, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,100)))

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
  filter(!(Strain == "HSC007" & Drug == "FLC" & Rep == "a")) %>% 
  group_by(Strain, Drug,Species, Origin, Year ) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) %>% 
  mutate(Biorep="B2") 


# Only keep strains that have two bioreps 
H1B1 <- H1B1 %>% subset(Strain%in%H1B2$Strain)
H1B2 <- H1B2 %>% subset(Strain%in%H1B1$Strain)

#* All HSC1-105, two bioreps, BA & FLC ----
H1_df<-rbind(H1B1, H1B2)


H1_ag <-rbind(H1B1,H1B2) %>% group_by(Strain, Drug, Species, Origin, Year) %>% 
  summarise(RAD20 = mean(RAD20), 
            FoG20 = mean(FoG20),
    RAD50 = mean(RAD50), 
            FoG50 = mean(FoG50))

H1stat <- H1_ag %>% group_by(Species, Drug, Year) %>% 
  summarise(RAD50Med = median(RAD50, na.rm = TRUE),
            RAD50Mean = mean(RAD50, na.rm = TRUE),
            FoG50Med = median(FoG50, na.rm = TRUE),
            FoG50mean = mean(FoG50, na.rm = TRUE))
H1statBA<- H1stat %>% filter(Drug=="BA")
H1statFLC<- H1stat %>% filter(Drug=="FLC")

Tab_H1<-tibble(Strain=unique(H1_ag$Strain)) %>% merge(H1Info) %>% 
  group_by(Species) %>% 
  summarise(n=n()) 

#We have 41 C.alb, 17 C.glab, 11 C.para, 5 C.trop from H1



##########################
H2Info<- read_csv(here("data_in", "DDA","Aleeza_HSC_Yeast stocks_2018.csv" )) %>% 
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

ggplotly(H2B1pt1 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD50, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,100)))

ggplotly(H2B1pt1 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG50, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,100)))

H2B1pt1 <- read_csv(here("data_in", "DDA","HSC106_159_48h_df_p1_B1.csv")) %>% 
  filter(!line=="Blank") %>%
  rename("Strain"="line", "Drug"="type") %>% 
  select(Strain, Drug, RAD20, FoG20, RAD50, FoG50) %>%
  group_by(Strain, Drug) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50), Biorep ="B1" ) %>% 
  merge(H2Info, by = "Strain")


#* HSC151-210 biorep 1, FLC & BA ----
H2B1pt2<- read_csv(here("data_in", "DDA","HSC106_159_48h_df_p2_B1.csv")) %>% 
  filter(!line=="Blank") %>%
  rename("Strain"="line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% 
  group_by(Strain, Drug) %>% 
  merge(H2Info, by = "Strain")%>% 
  group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain))


ggplotly(H2B1pt2 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD50, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,100)))

ggplotly(H2B1pt2 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG50, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,100)))



H2B1pt2<- read_csv(here("data_in", "DDA","HSC106_159_48h_df_p2_B1.csv")) %>% 
  filter(!line=="Blank") %>%
  rename("Strain"="line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% 
  group_by(Strain, Drug) %>% 
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

ggplotly(H2B2 %>% filter(Drug == "BA") %>%  ggplot() + geom_point(aes(x = Strain, y = RAD50, color = Rep)) +
           scale_y_continuous(limits = c(0,25)) + scale_x_discrete(label = seq(0,100)))

ggplotly(H2B2 %>% filter(Drug == "FLC") %>%  ggplot() + geom_point(aes(x = Strain, y = FoG50, color = Rep)) +
           scale_y_continuous(limits = c(0,1)) + scale_x_discrete(label = seq(0,100)))

#Overwritting H2B2 
H2B2 <- read_csv(here("data_in", "DDA","HSC106_205_48h_df_B2.csv")) %>% 
  filter(!line=="Blank") %>%
  rename("Strain"="line", "Drug"="type") %>% 
  select(Strain, Drug, RAD20, FoG20, RAD50, FoG50) %>% 
  group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain)) %>% 
  filter(!(Strain == "HSC204" &  Drug == "FLC" &Rep == 1)) %>% 
  filter(!(Strain == "HSC138" &  Drug == "FLC" & Rep == 1)) %>% 
  filter(!(Strain == "HSC131" &  Drug == "FLC" )) %>% 
  group_by(Strain, Drug) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50), Biorep ="B2" ) %>% 
  merge(H2Info, by = "Strain") %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis", 
                      "C.tropicalis")) %>% 
  filter(!Strain == "HSC167") %>% 
  filter(!Strain == "HSC169")


#* HSC106-210 biorep 1&2, FLC & BA ----
H2_df <- rbind(H2B1, H2B2) 
 

H2_ag <- rbind(H2B1, H2B2) %>% group_by(Strain, Drug, Species, Origin, Year) %>% 
  summarise(RAD50=mean(RAD20), 
            FoG50=mean(FoG20),
    RAD50=mean(RAD50), 
            FoG50=mean(FoG50))

H2stat <- H2_ag %>% group_by(Species, Drug, Year) %>% 
  summarise(RAD50Med = median(RAD50, na.rm = TRUE),
            RAD50Mean = mean(RAD50, na.rm = TRUE),
            FoG50Med = median(FoG50, na.rm = TRUE),
            FoG50mean = mean(FoG50, na.rm = TRUE))
H2statBA<- H2stat %>% filter(Drug=="BA")
H2statFLC<- H2stat %>% filter(Drug=="FLC")

Tab_H2<-tibble(Strain=unique(H2_ag$Strain)) %>% merge(H2Info) %>% 
  group_by(Species) %>% 
  summarise(n=n()) 

#* Combining the HSC 2 datasets (batch1: 2012-2013, batch2: 2018-2019) ----
HSC <- rbind(H1_ag, H2_ag) 

Info<-rbind(H1Info, H2Info)

HSCstat<- HSC %>% group_by(Species, Drug ) %>% 
  summarise(RAD50Med = median(RAD50, na.rm = TRUE),
            RAD50Mean = mean(RAD50, na.rm = TRUE),
            FoG50Med = median(FoG50, na.rm = TRUE),
            FoG50mean = mean(FoG50, na.rm = TRUE))
HSCstatBA<- HSCstat %>% filter(Drug=="BA")
HSCstatFLC<-HSCstat %>% filter(Drug=="FLC")

HSC %>% filter(Drug=="BA") %>% 
  group_by(Species) %>% 
  summarise(n=n()) %>% 
  kable("html") %>% 
  kable_styling(latex_options = c("striped", "scale_down"))  

#* Biorep 3 - HSC years 1 & 2
alb_glb_B03 <- read_csv(here("data_in","DDA", "DDA201218_48h_df.csv")) %>%
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% filter(!Strain=="Blank") %>% 
  group_by(Strain, Drug) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50)) %>% 
  subset(Strain%in%HSC$Strain)%>% 
  merge(Info) %>% mutate( Biorep="B03") %>% 
  filter(!Strain == "HSC105" & Biorep == "B03")

para_trop_B03 <- read_csv(here("data_in","DDA", "Para_Trop_B03_df.csv")) %>% 
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% filter(!Strain=="Blank") %>% 
  group_by(Strain, Drug) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) %>% 
  subset(Strain%in%HSC$Strain)%>% 
  merge(Info) %>% mutate( Biorep="B03") 

HB3 <- rbind(alb_glb_B03, para_trop_B03)

H_df<- rbind(H1_df, H2_df) 

H_wide <- rbind(H1B1, H2B1) %>% 
 merge(rbind(H1B2, H2B2), by = c("Strain","Drug","Species", "Origin", "Year")) %>% 
  rename("RAD20_B01" = "RAD20.x", "RAD20_B02" = "RAD20.y",
         "RAD50_B01" = "RAD50.x", "RAD50_B02" = "RAD50.y",
         "FoG20_B01" = "FoG20.x", 
         "FoG50_B01" = "FoG50.x", 
         "FoG20_B02" = "FoG20.y",
         "FoG50_B02" = "FoG50.y") %>% 
  select(!c(Biorep.x,Biorep.y)) %>% 
  merge(HB3, by = c("Strain","Drug","Species", "Origin", "Year")) %>% 
  rename("RAD20_B03" = "RAD20", "FoG20_B03" = "FoG20", "RAD50_B03" = "RAD50", "FoG50_B03" = "FoG50") %>% 
  select(!Biorep)

write_csv(H_wide, here("data_out", "DDA","210821_HSC_B1B2B3.csv"))

cor.test(H_wide$RAD20_B01, H_wide$RAD50_B01)
cor.test(H_wide$RAD20_B02, H_wide$RAD50_B02)
cor.test(H_wide$RAD20_B03, H_wide$RAD50_B03)

cor.test(H_wide$FoG20_B01, H_wide$FoG50_B01)
cor.test(H_wide$FoG20_B02, H_wide$FoG50_B02)
cor.test(H_wide$FoG20_B03, H_wide$FoG50_B03)

cor.test(H_wide$RAD50_B01, H_wide$RAD50_B02)
cor.test(H_wide$RAD50_B02, H_wide$RAD50_B03)
cor.test(H_wide$RAD50_B01, H_wide$RAD50_B03)
cor.test(H_wide$FoG50_B01, H_wide$FoG50_B02)
cor.test(H_wide$FoG50_B01, H_wide$FoG50_B03)

# I am going to keep all of the bioreps since they correlated with each other and the correlation was significat 

H_all_stat <- H_wide %>% group_by(Strain, Drug, Species, Year) %>%
  summarise(RAD50_sd = sd(c(RAD50_B01, RAD50_B02, RAD50_B03)),
            FoG50_sd = sd(c(FoG50_B01, FoG50_B02, FoG50_B03)),
            RAD50 = (RAD50_B01 + RAD50_B02 + RAD50_B03)/3,
            FoG50 = (FoG50_B01 + FoG50_B02 + FoG50_B03)/3)

write_csv(H_all_stat, here("data_out", "DDA","210629_HSC_sd.csv"))

H_all_stat %>% filter(Drug=="BA") %>% 
  group_by(Species) %>% 
  summarise(n=n()) %>% 
  kable("html") %>% 
  kable_styling(latex_options = c("striped", "scale_down"))

#70 c.alb, 44 c. glb, 20 c.para, 7 c. trop



#* AG strains ----
AGInfo <- read_csv(here("data_in", "DDA","BermanStrains_GersteinInventory.csv")) %>% 
  select('AG strain', Species, 'Strain Name:') %>% 
  rename("Strain" = 'AG strain', "Strain_Name"='Strain Name:')

AGInfo$Strain <- sub("^", "yAG", AGInfo$Strain )
AGInfo$Species <- sub(" ", "", AGInfo$Species)

#* AG strains biorep 1 FLC & BA ----
#Strains that have been filtered do not match the image
AGB1 <- read_csv(here("data_in", "DDA","AG_201001_48h_B1_df.csv")) %>%
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  filter(!Strain=="Blank") %>% 
  group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain)) %>% 
  filter(!(Strain == "yAG16" & Drug == "FLC"& Rep == 1 )) %>% 
  filter(!(Strain == "yAG20" & Drug == "FLC"& Rep == 1 )) %>% 
  filter(!(Strain == "yAG81" & Drug == "FLC"& Rep == 1 )) %>% 
  filter(!(Strain == "yAG108" & Drug == "FLC"& Rep == 1 )) %>% 
  filter(!(Strain == "yAG11" & Drug == "FLC"& Rep == 2 )) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) %>% 
  merge(AGInfo) %>% mutate(Biorep = "B1") %>% 
  filter(!Strain== "yAG41") %>%
  filter(!Strain=="yAG92") %>% 
  filter(!Strain=="yAG109") %>% 
filter(!Strain=="yAG85")

#* AG strains biorep 2 FLC & BA ----
AGB2 <- read_csv(here("data_in", "DDA","AG_201008_48h_B2_df.csv")) %>%
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD20, FoG20, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  
  filter(!Strain=="Blank") %>% 
  group_by(Strain, Drug) %>% 
  mutate(Rep=1:length(Strain)) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20) , RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) %>% 
  merge(AGInfo) %>% mutate(Biorep = "B2") %>% 
  filter(!Strain== "yAG40") %>%
  filter(!Strain== "yAG41") %>%
  filter(!Strain== "yAG55") %>% 
  filter(!Strain=="yAG109") %>% 
  filter(!Strain== "yAG51")
#yAG92 is not in B2
#yAG40 and 55 is missing in B1

AGB1 <- AGB1 %>% filter(Strain %in% AGB2$Strain)
AGB2 <- AGB2 %>% filter(Strain %in% AGB1$Strain)

AGB1 %>% filter(!Strain %in% AGB2$Strain)
AGB2 %>% filter(!Strain %in% AGB1$Strain)
#yAG92 is not in B2
#yAG40 and 55 is missing in B1

AG <- rbind(AGB1,AGB2)  


AG_df <- AG %>% 
  group_by(Strain, Drug, Species, Type) %>% 
  summarise(RAD50_sd = sd(RAD50),
            FoG50_sd = sd(FoG50),
    RAD50 = mean(RAD50), 
            FoG50 = mean(FoG50))



MissingStrains <- read_csv(here("data_in", "DDA","MissingStrains_df.csv")) %>% 
  filter(!line == "Blank") %>% 
  select(name, RAD20, FoG20, RAD50, FoG50) %>%
  separate(name, c("Strain", "Drug", "Rep", "Biorep"), sep = "_") %>% 
  filter(!(Strain== "yAG046")) %>% 
  filter(!(Strain== "yAG117")) %>% 
  filter(!(Strain== "yAG058")) %>% 
  filter(!(Strain== "yAG003"& Drug == "FLC" & Rep == "a" & Biorep == "B02")) %>% 
  filter(!(Strain== "yAG054"& Drug == "FLC" & Rep == "a" & Biorep == "B01")) %>% 
  filter(!(Strain== "yAG98"& Drug == "FLC" & Rep == "a" & Biorep == "B02")) %>% 
  filter(!(Strain== "yAG109"& Drug == "FLC" & Rep == "b" & Biorep == "B01")) %>% 
  group_by(Strain, Drug, Biorep) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50)) %>% 
  mutate(strain = Strain) %>% 
  separate(strain, c(NA, "nAG"), sep = "yAG") 

MissingStrains$nAG <- as.numeric(MissingStrains$nAG)

MissingStrains <- MissingStrains %>% select(!Strain) %>% 
  mutate(AG = "yAG") %>% unite(Strain, c("AG", "nAG"), sep = "" )
 

MissingStrains <- merge(MissingStrains, AGInfo, by = "Strain") %>% 
  mutate(Biorep = "B3")


AG %>% filter(Strain %in% MissingStrains$Strain)

MissingStrains %>% filter(Strain %in% AG$Strain)

AG_full <- rbind(AG, MissingStrains) %>% 
  group_by(Strain, Drug, Species) %>% 
  summarise(RAD20 = mean(RAD20), FoG20 = mean(FoG20), RAD50 = mean(RAD50), FoG50 = mean(FoG50))

write_csv(AG_full, here("data_out", "DDA","210821_AG_B1B2.csv"))

AG_full %>% filter(Drug=="BA") %>% 
  group_by(Species) %>% 
  summarise(n=n()) %>% 
  kable("html") %>% 
  kable_styling(latex_options = c("striped", "scale_down"))

#78 C.alb strains
HSC$Origin <- NULL

AGHSC<- HSC %>% rbind(AG_full)
write_csv(AGHSC, here("data_out", "DDA","210821_HSCAG.csv"))

AGHSC %>% filter(Drug=="BA") %>% 
  group_by(Species) %>% 
  summarise(n=n()) %>% 
  kable("html") %>% 
  kable_styling(latex_options = c("striped", "scale_down"))
