library(here)
library(tidyverse)
library(ggforce)
library(ggpubr)
library(knitr)
library(kableExtra)
library(webshot)
library(extrafont)
library(plotly)


#upload the raw data
H1Info<-H1Info<-read_csv(here("data_in","2014_Isolate list_info.csv")) %>% 
  select(HSC, 'Isolate ID', SPECIMEN) %>% 
  rename("Strain" = "HSC", "Species"='Isolate ID', "Origin"="SPECIMEN") %>% 
  mutate(Year="2014")

H1B1BA <- read_csv(here("data_in","HSC001-105_48h_B1_df.csv")) %>%
  rename("Strain" = "line") %>% mutate(Drug="BA") %>% 
  select(Strain,Drug, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  summarise(RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) %>% 
  merge(H1Info) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis", "C.tropicalis"))


H1B1FLC <- read_csv(here("data_in","dda_200309_Abdul_FLC_48h_B1_df.csv")) %>%
  rename("Strain" = "line") %>% mutate(Drug="FLC") %>% 
  select(Strain,Drug, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  summarise(RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) %>% 
  subset(Strain%in%H1B1BA$Strain)%>% 
  merge(H1Info) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis",
                      "C.tropicalis"))


H1B1 <-rbind(H1B1BA,H1B1FLC) %>% mutate(Biorep="B1")

#HSC50 is missing from biorep2
#HSC59,61,71 is missing from Biorep 1

H1B2 <- read_csv(here("data_in","HSC001-105_48h_B2_df.csv")) %>%
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD50, FoG50) %>% filter(!Strain=="Blank") %>% 
  group_by(Strain, Drug) %>% 
  summarise(RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) %>% 
  subset(Strain%in%H1B1$Strain)%>% 
  merge(H1Info) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata",
                      "C.parapsilosis", "C.tropicalis")) %>% 
  mutate(Biorep="B2")

H1B1<- H1B1 %>% subset(Strain%in%H1B2$Strain)

H1<-rbind(H1B1,H1B2) %>% group_by(Strain, Drug, Species, Origin, Year) %>% 
  summarise(RAD50=mean(RAD50), 
            FoG50=mean(FoG50))

H1stat<-H1 %>% group_by(Species, Drug, Year) %>% 
  summarise(RAD50Med = median(RAD50, na.rm = TRUE),
            RAD50Mean = mean(RAD50, na.rm = TRUE),
            FoG50Med = median(FoG50, na.rm = TRUE),
            FoG50mean = mean(FoG50, na.rm = TRUE))
H1statBA<- H1stat %>% filter(Drug=="BA")
H1statFLC<- H1stat %>% filter(Drug=="FLC")

Tab_H1<-tibble(Strain=unique(H1$Strain)) %>% merge(H1Info) %>% 
  group_by(Species) %>% 
  summarise(n=n()) 

#We have 40 C.alb, 18 C.glab, 12 C.para, 5 C.trop from H1

#Plotting the first Batch 

#BA Resistance graph
H1p1<-H1 %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( H1statBA, mapping = aes(x=Species, y= RAD50Med),
              color="red",pch="_",size=10, alpha=0.8)+
  ylab("RAD50")+expand_limits(x = 0, y = c(0,28))+
  theme_bw()+ labs(subtitle = "BA")+
  scale_y_reverse()+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

#FLC Resistance graph
H1p2<-H1 %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( H1statFLC, mapping = aes(x=Species, y= RAD50Med),
              color="red",pch="_",size=10, alpha=0.8)+
  ylab("RAD50")+expand_limits(x = 0, y = c(0,28))+
  theme_bw()+ labs(subtitle = "FLC")+
  scale_y_reverse()+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

#BA Tolerance graph 
H1p3<-H1 %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG50), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( H1statBA, mapping = aes(x=Species, y= FoG50Med),
              color="red",pch="_",size=10, alpha=0.8)+
  ylab("FoG50")+expand_limits(x = 0, y = c(0,1))+
  theme_bw()+ labs(subtitle = "BA")+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

#FLC Tolerance graph 
H1p4<-H1 %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG50), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( H1statFLC, mapping = aes(x=Species, y= FoG50Med),
              color="red",pch="_",size=10, alpha=0.8)+
  ylab("FoG50") + expand_limits(x = 0, y = c(0,1))+
  theme_bw()+ labs(subtitle = "FLC")+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

ggarrange(H1p1, H1p2, H1p3, H1p4, align = "hv")

#Second Batch of HSC 
H2Info<- read_csv(here("data_in","Aleeza_HSC_Yeast stocks_2018.csv" )) %>% 
  select('HSC#', "SPECIMEN", "ORGANISM") %>% 
  rename("Strain"='HSC#', "Species"="ORGANISM", "Origin"="SPECIMEN") %>% 
  mutate(Year="2018")

H2B1pt1<- read_csv(here("data_in","HSC106_159_48h_df_p1_B1.csv")) %>% 
  filter(!line=="Blank") %>%
  rename("Strain"="line", "Drug"="type") %>% 
  select(Strain, Drug, RAD50, FoG50) %>%
  group_by(Strain, Drug) %>% 
  summarise(RAD50 = mean(RAD50), FoG50 = mean(FoG50), Biorep ="B1" ) %>% 
  merge(H2Info, by = "Strain")

H2B1pt2<- read_csv(here("data_in","HSC106_159_48h_df_p2_B1.csv")) %>% 
  filter(!line=="Blank") %>%
  rename("Strain"="line", "Drug"="type") %>% 
  select(Strain,Drug, RAD50, FoG50) %>% 
  group_by(Strain, Drug) %>% 
  summarise(RAD50 = mean(RAD50), FoG50 = mean(FoG50), Biorep ="B1" ) %>% 
  merge(H2Info, by = "Strain")

H2B1<- rbind(H2B1pt1, H2B1pt2) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis", "C.tropicalis"))

H2B2<- read_csv(here("data_in","HSC106_205_48h_df_B2.csv")) %>% 
  filter(!line=="Blank") %>%
  rename("Strain"="line", "Drug"="type") %>% 
  select(Strain, Drug, RAD50, FoG50) %>% 
  group_by(Strain, Drug) %>% 
  summarise(RAD50 = mean(RAD50), FoG50 = mean(FoG50), Biorep ="B2" ) %>% 
  merge(H2Info, by = "Strain") %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis", 
                      "C.tropicalis")) %>% 
  filter(!Strain == "HSC167")

#HSC167 is missing from Biorep 1

H2<-rbind(H2B1, H2B2) %>% group_by(Strain, Drug, Species, Origin, Year) %>% 
  summarise(RAD50=mean(RAD50), 
            FoG50=mean(FoG50))

H2stat<-H2 %>% group_by(Species, Drug,Year) %>% 
  summarise(RAD50Med = median(RAD50, na.rm = TRUE),
            RAD50Mean = mean(RAD50, na.rm = TRUE),
            FoG50Med = median(FoG50, na.rm = TRUE),
            FoG50mean = mean(FoG50, na.rm = TRUE))
H2statBA<- H2stat %>% filter(Drug=="BA")
H2statFLC<- H2stat %>% filter(Drug=="FLC")

Tab_H2<-tibble(Strain=unique(H2$Strain)) %>% merge(H2Info) %>% 
  group_by(Species) %>% 
  summarise(n=n()) 

#there are 38 C.alb, 30 C.26, 9 C.para and 2 C. trop 

#Plotting the second Batch 

#BA Resistance graph
H2p1<-H2 %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( H2statBA, mapping = aes(x=Species, y= RAD50Med),
              color="red",pch="_",size=10, alpha=0.8)+
  ylab("RAD50")+expand_limits(x = 0, y = c(0,28))+
  theme_bw()+ labs(subtitle = "BA")+
  scale_y_reverse()+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

#FLC Resistance graph
H2p2<-H2 %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( H2statFLC, mapping = aes(x=Species, y= RAD50Med),
              color="red",pch="_",size=10, alpha=0.8)+
  ylab("RAD50")+expand_limits(x = 0, y = c(0,28))+
  theme_bw()+ labs(subtitle = "FLC")+
  scale_y_reverse()+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

#BA Tolerance graph 
H2p3<-H2 %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG50), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( H2statBA, mapping = aes(x=Species, y= FoG50Med),
              color="red",pch="_",size=10, alpha=0.8)+
  ylab("FoG50")+expand_limits(x = 0, y = c(0,1))+
  theme_bw()+ labs(subtitle = "BA")+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

#FLC Tolerance graph 
H2p4<-H2 %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG50), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( H2statFLC, mapping = aes(x=Species, y= FoG50Med),
              color="red",pch="_",size=10, alpha=0.8)+
  ylab("FoG50") + expand_limits(x = 0, y = c(0,1))+
  theme_bw()+ labs(subtitle = "FLC")+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

H2P_all<-ggarrange(H2p1, H2p2, H2p3, H2p4, align = "hv")
ggsave("201202_HSC_Batch2.jpg", 
       path = here("figures_out"), 
       plot = H2P_all,
       width = 6.5,
       height = 3.5)

#combining the HSC 2 batches
HSC<-rbind(H1, H2) 

HSCstat<- HSC %>% group_by(Species, Drug, Year) %>% 
  summarise(RAD50Med = median(RAD50, na.rm = TRUE),
            RAD50Mean = mean(RAD50, na.rm = TRUE),
            FoG50Med = median(FoG50, na.rm = TRUE),
            FoG50mean = mean(FoG50, na.rm = TRUE))%>%
  filter(Species%in%c("C.albicans","C.glabrata"))
HSCstatBA<- HSCstat %>% filter(Drug=="BA")
HSCstatFLC<-HSCstat %>% filter(Drug=="FLC")

HSC %>% filter(Drug=="BA") %>% 
  group_by(Species, Year) %>% 
  summarise(n=n()) %>% 
  kable("html") %>% 
  kable_styling(latex_options = c("striped", "scale_down"))  


#Plotting the full data
HSCp1<-HSC %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50, color= Year), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( HSCstatBA, mapping = aes(x=Species, y= RAD50Med, group= Year),
              color="red",pch="_",size=6, alpha=0.8, 
              position= position_dodge(width = 0.5))+
  ylab("RAD50")+expand_limits(x = 0, y = c(0,28))+
  theme_bw()+ labs(subtitle = "BA")+
  scale_y_reverse()+
  scale_color_manual(values = c("#333333", "#660066"))+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

#FLC Resistance graph
HSCp2<-HSC %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50, color = Year), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( HSCstatFLC, mapping = aes(x=Species, y= RAD50Med,group= Year),
              color="red",pch="_",size=6, alpha=0.8, 
              position= position_dodge(width = 0.5)) +
  ylab("RAD50")+expand_limits(x = 0, y = c(0,28))+
  theme_bw()+ labs(subtitle = "FLC")+
  scale_y_reverse()+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#333333", "#660066"))

#BA Tolerance graph 
HSCp3<-HSC %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG50, color =Year), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( HSCstatBA, mapping = aes(x=Species, y= FoG50Med, group= Year),
              color="red",pch="_",size=6, alpha=0.8, 
              position= position_dodge(width = 0.5))+
  ylab("FoG50") + 
  expand_limits(x = 0, y = c(0,1))+
  theme_bw() + labs(subtitle = "BA") +
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#333333", "#660066"))

#FLC Tolerance graph 
HSCp4<-HSC %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG50, color = Year), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( HSCstatFLC, mapping = aes(x=Species,
                                        y= FoG50Med, group= Year),
              color="red",pch="_",size=6, alpha=0.8, 
              position= position_dodge(width = 0.5))+
  ylab("FoG50") + expand_limits(x = 0, y = c(0,1))+
  theme_bw()+ labs(subtitle = "FLC")+
  theme(legend.position="none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#333333", "#660066"))

ggarrange(HSCp1, HSCp2, HSCp3, HSCp4, align = "hv")

########################################################################
#upload the raw data


Info<-rbind(H1Info, H2Info)


#Comb1418= Combined 2014 and 2018 strains (2014 and 2018 strains processed in the same batch)
Comb1418 <- read_csv(here("data_in","DDA201218_48h_df.csv")) %>%
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD50, FoG50) %>% filter(!Strain=="Blank") %>% 
  group_by(Strain, Drug) %>% 
  summarise(RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) %>% 
  subset(Strain%in%HSC$Strain)%>% 
  merge(Info) %>% mutate(Type="Comb")
#Sep1418= 2014 and 2018 strains that were processed in seperate batches 
Sep1418<-HSC%>%subset(Strain%in%Comb1418$Strain) %>% mutate(Type="Sep")

CombSep<-rbind(Comb1418,Sep1418)
CS<-unite(CombSep, Year, c("Year","Type"), sep = "_")

CSstat<- CS %>% group_by(Species, Drug, Year) %>% 
  summarise(RAD50Med = median(RAD50, na.rm = TRUE),
            RAD50Mean = mean(RAD50, na.rm = TRUE),
            FoG50Med = median(FoG50, na.rm = TRUE),
            FoG50mean = mean(FoG50, na.rm = TRUE))
CSstatBA<- CSstat %>% filter(Drug=="BA")
CSstatFLC<-CSstat %>% filter(Drug=="FLC")


#Plotting the full data
CSp1<-CS %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50, color= Year), 
            maxwidth = 1,alpha = 0.6, size=2) + 
  geom_point( CSstatBA, mapping = aes(x=Species, y= RAD50Med, group= Year),
              color="red",pch="_",size=6, alpha=0.8, 
              position= position_dodge(width = 1))+
  ylab("RAD50")+expand_limits(x = 0, y = c(0,28))+
  theme_bw()+ labs(subtitle = "BA")+
  scale_y_reverse()+
  scale_color_manual(values = c("#333333","grey40","#660066", "maroon4"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

#FLC Resistance graph
CSp2<-CS%>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50, color = Year), 
            alpha = 0.6, size=2, maxwidth=1) + 
  geom_point( CSstatFLC, mapping = aes(x=Species, y= RAD50Med,group= Year),
              color="red",pch="_",size=6, alpha=0.8, 
              position= position_dodge(width = 1)) +
  ylab("RAD50")+expand_limits(x = 0, y = c(0,28))+
  theme_bw()+ labs(subtitle = "FLC")+
  scale_y_reverse()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#333333","grey40","#660066", "maroon4"))

#BA Tolerance graph 
CSp3<-CS %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG50, color =Year), 
            maxwidth = 1,alpha = 0.6, size=2) + 
  geom_point( CSstatBA, mapping = aes(x=Species, y= FoG50Med, group= Year),
              color="red",pch="_",size=6, alpha=0.8, 
              position= position_dodge(width = 1))+
  ylab("FoG50") + 
  expand_limits(x = 0, y = c(0,1))+
  theme_bw() + labs(subtitle = "BA") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#333333","grey40","#660066", "maroon4"))

#FLC Tolerance graph 
CSp4<-CS %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG50, color = Year), 
            maxwidth = 1,alpha = 0.6, size=2) + 
  geom_point( CSstatFLC, mapping = aes(x=Species,
                                          y= FoG50Med, group= Year),
              color="red",pch="_",size=6, alpha=0.8, 
              position= position_dodge(width = 1))+
  ylab("FoG50") + expand_limits(x = 0, y = c(0,1))+
  theme_bw()+ labs(subtitle = "FLC")+
  theme( panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         text=element_text(family="Times New Roman", 
                           face="bold", size = 8),
         axis.text.x = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#333333","grey40","#660066", "maroon4"))

CS_all<-ggarrange(CSp1, CSp2, CSp3, CSp4, align = "hv",
                     common.legend = TRUE)

ggsave("201228_HSC2013vs2018.jpg", 
       path = here("figures_out"), 
       plot = Pilot_all,
       width = 6.5,
       height = 3.5)


