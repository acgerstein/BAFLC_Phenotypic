#install.packages("here")
#install.packages("tidyverse")
#install.packages("ggforce")
#install.packages("ggpubr)
#install.packages("kableExtra")

library(here)
library(tidyverse)
library(ggforce)
library(ggpubr)
library(knitr)
library(kableExtra)
library(webshot)
library(extrafont)
library(viridis)


#upload the raw data
H1Info<-read_csv(here("data_in","DDA","2014_HSC_info.csv")) %>% 
  select(!X1) %>% 
  rename("Strain" = "line")

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

H1<-rbind(H1B1,H1B2) %>% group_by(Strain, Drug, Species, Origin) %>% 
  summarise(RAD50=mean(RAD50), 
            FoG50=mean(FoG50), Type= "HSC_2014")

H1stat<-H1 %>% group_by(Species, Drug) %>% 
  summarise(RAD50Med = median(RAD50, na.rm = TRUE),
            RAD50Mean = mean(RAD50, na.rm = TRUE),
            FoG50Med = median(FoG50, na.rm = TRUE),
            FoG50mean = mean(FoG50, na.rm = TRUE))
H1statBA<- H1stat %>% filter(Drug=="BA")
H1statFLC<- H1stat %>% filter(Drug=="FLC")

Tab_H1<-tibble(Strain=unique(H1$Strain)) %>% merge(H1Info) %>% 
  group_by(Species) %>% 
  summarise(n=n()) 

#We have 41 C.alb, 17 C.glab, 11 C.para, 5 C.trop from H1

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

H1P_all<-ggarrange(H1p1, H1p2, H1p3, H1p4, align = "hv")
ggsave("201202_HSC_Batch1.jpg", 
       path = here("figures_out"),
       plot = H1P_all, 
       width = 6.5,
       height = 3.5)


#Second Batch of HSC 
H2Info<- read_csv(here("data_in","Aleeza_HSC_Yeast stocks_2018.csv" )) %>% 
  select('HSC#', "SPECIMEN", "ORGANISM") %>% 
  rename("Strain"='HSC#', "Species"="ORGANISM", "Origin"="SPECIMEN")

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

H2<-rbind(H2B1, H2B2) %>% group_by(Strain, Drug, Species, Origin) %>% 
  summarise(RAD50=mean(RAD50), 
            FoG50=mean(FoG50), Type="HSC_2018")

H2stat<-H2 %>% group_by(Species, Drug) %>% 
  summarise(RAD50Med = median(RAD50, na.rm = TRUE),
            RAD50Mean = mean(RAD50, na.rm = TRUE),
            FoG50Med = median(FoG50, na.rm = TRUE),
            FoG50mean = mean(FoG50, na.rm = TRUE))
H2statBA<- H2stat %>% filter(Drug=="BA")
H2statFLC<- H2stat %>% filter(Drug=="FLC")

Tab_H2<-tibble(Strain=unique(H2$Strain)) %>% merge(H2Info) %>% 
  group_by(Species) %>% 
  summarise(n=n()) 

#there are 31 C.alb, 26 C.26, 4 C.para and 1 C. trop 

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
HSC<-rbind(H1, H2) %>% select(!Origin)

HSCstat<- HSC %>% group_by(Species, Drug, Type) %>% 
  summarise(RAD50Med = median(RAD50, na.rm = TRUE),
            RAD50Mean = mean(RAD50, na.rm = TRUE),
            FoG50Med = median(FoG50, na.rm = TRUE),
            FoG50mean = mean(FoG50, na.rm = TRUE))
HSCstatBA<- HSCstat %>% filter(Drug=="BA")
HSCstatFLC<-HSCstat %>% filter(Drug=="FLC")

HSC %>% filter(Drug=="BA") %>% 
  group_by(Species, Type) %>% 
  summarise(n=n()) %>% 
  kable("html") %>% 
  kable_styling(latex_options = c("striped", "scale_down"))  
  

#Plotting the full data
HSCp1<-HSC %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50, color= Type), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( HSCstatBA, mapping = aes(x=Species, y= RAD50Med, group= Type),
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
  geom_sina(mapping= aes(x=Species, y=RAD50, color = Type), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( HSCstatFLC, mapping = aes(x=Species, y= RAD50Med,group= Type),
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
  geom_sina(mapping= aes(x=Species, y=FoG50, color =Type), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( HSCstatBA, mapping = aes(x=Species, y= FoG50Med, group= Type),
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
  geom_sina(mapping= aes(x=Species, y=FoG50, color = Type), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( HSCstatFLC, mapping = aes(x=Species,
                                       y= FoG50Med, group= Type),
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

HSCP_all<-ggarrange(HSCp1, HSCp2, HSCp3, HSCp4, align = "hv")

ggsave("201202_HSC_comb.jpg", 
       path = here("figures_out"), 
       plot = HSCP_all,
       width = 6.5,
       height = 3.5)

HSCstatall<- HSC %>% group_by(Species, Drug) %>% 
  summarise(RAD50Med = median(RAD50, na.rm = TRUE),
            RAD50Mean = mean(RAD50, na.rm = TRUE),
            FoG50Med = median(FoG50, na.rm = TRUE),
            FoG50mean = mean(FoG50, na.rm = TRUE))
HSCstatBA_all<- HSCstatall %>% filter(Drug=="BA")
HSCstatFLC_all<-HSCstatall %>% filter(Drug=="FLC")

HSC %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( HSCstatBA_all, mapping = aes(x=Species, y= RAD50Med),
              color="red",pch="_",size=10, alpha=0.8, 
              position= position_dodge(width = 0.5))+
  ylab("RAD50")+expand_limits(x = 0, y = c(0,28))+
  theme_bw()+ labs(subtitle = "BA")+
  scale_y_reverse()+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

AGInfo<-read_csv(here("data_in","BermanStrains_GersteinInventory.csv")) %>% 
  select('AG strain', Species, 'Strain Name:') %>% 
  rename("Strain" = 'AG strain', "Strain_Name"='Strain Name:')

AGInfo$Strain <- sub("^", "yAG", AGInfo$Strain )
AGInfo$Species <- sub(" ", "", AGInfo$Species)

AGB1 <- read_csv(here("data_in","AG_201001_48h_B1_df.csv")) %>%
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  filter(!Strain=="Blank") %>% 
  summarise(RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) %>% 
  merge(AGInfo) %>% mutate(Biorep = "B1") %>% 
  filter(!Strain=="yAG92")

AGB2 <- read_csv(here("data_in","AG_201008_48h_B2_df.csv")) %>%
  rename("Strain" = "line", "Drug"="type") %>% 
  select(Strain,Drug, RAD50, FoG50) %>% group_by(Strain, Drug) %>% 
  filter(!Strain=="Blank") %>% 
  summarise(RAD50 = mean(RAD50), FoG50 = mean(FoG50) ) %>% 
  merge(AGInfo) %>% mutate(Biorep = "B2") %>% 
  filter(!Strain== "yAG40") %>% filter(!Strain== "yAG55")

AGB1 %>% filter(!Strain %in% AGB2$Strain)
AGB2 %>% filter(!Strain %in% AGB1$Strain)
#yAG92 is not in B2
#yAG40 and 55 is missing in B1

AG<-rbind(AGB1,AGB2) %>% mutate(Type ="yAG") %>% select(!Strain_Name)  
  
AG_df<- AG %>% 
  group_by(Strain, Drug, Species, Type) %>% 
  summarise(RAD50=mean(RAD50), 
            FoG50=mean(FoG50))

df<-rbind(AG_df, HSC)

dfp1<-df %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50, color= Type), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  ylab("RAD50")+expand_limits(x = 0, y = c(0,28))+
  theme_bw()+ labs(subtitle = "BA")+
  scale_y_reverse()+
  scale_color_manual(values = c("#333333", "#660066", "#000066"))+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

dfp2<-df %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50, color= Type), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  ylab("RAD50")+expand_limits(x = 0, y = c(0,28))+
  theme_bw()+ labs(subtitle = "FLC")+
  scale_y_reverse()+
  scale_color_manual(values = c("#333333", "#660066", "#000066"))+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

dfp3<-df %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG50, color= Type), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  ylab("FoG50")+expand_limits(x = 0, y = c(0,1))+
  theme_bw()+ labs(subtitle = "BA")+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#333333", "#660066", "#000066"))

dfp4<-df %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG50, color= Type), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  ylab("FoG50")+expand_limits(x = 0, y = c(0,1))+
  theme_bw()+ labs(subtitle = "FLC")+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#333333", "#660066", "#000066"))

df_all<-ggarrange(dfp1, dfp2, dfp3, dfp4, align = "hv")

ggsave("201203_HSC_AG.jpg", 
       path = here("figures_out"), 
       plot = df_all,
       width = 6.5,
       height = 3.5)

###################################################################
df_BAStat<-df %>% filter (Drug=="BA") %>%
  group_by(Species) %>% 
  summarise(RAD50Med = median(RAD50, na.rm = TRUE),
            RAD50Mean = mean(RAD50, na.rm = TRUE),
            FoG50Med = median(FoG50, na.rm = TRUE),
            FoG50mean = mean(FoG50, na.rm = TRUE))

df_FLCStat<-df %>% filter (Drug=="FLC") %>%
  group_by(Species) %>% 
  summarise(RAD50Med = median(RAD50, na.rm = TRUE),
            RAD50Mean = mean(RAD50, na.rm = TRUE),
            FoG50Med = median(FoG50, na.rm = TRUE),
            FoG50mean = mean(FoG50, na.rm = TRUE))

DFp1<-df %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50), 
            alpha = 0.6, size=2) +
  geom_point( df_BAStat, mapping = aes(x=Species, y= RAD50Med),
              color="red",pch="_",size=10, alpha=0.8)+
  ylab("RAD50")+expand_limits(x = 0, y = c(0,28))+
  theme_bw()+ labs(subtitle = "BA")+
  scale_y_reverse()+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

DFp2<-df %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50), 
            alpha = 0.6, size=2) + 
  geom_point( df_FLCStat, mapping = aes(x=Species, y= RAD50Med),
              color="red",pch="_",size=10, alpha=0.8)+
  ylab("RAD50")+expand_limits(x = 0, y = c(0,28))+
  theme_bw()+ labs(subtitle = "FLC")+
  scale_y_reverse()+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

DFp3<-df %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG50), 
            alpha = 0.6, size=2) + 
  geom_point( df_BAStat, mapping = aes(x=Species, y= FoG50Med),
              color="red",pch="_",size=10, alpha=0.8)+
  ylab("FoG50")+expand_limits(x = 0, y = c(0,1))+
  theme_bw()+ labs(subtitle = "BA")+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

DFp4<-df %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG50), 
           alpha = 0.6, size=2) + 
  geom_point( df_FLCStat, mapping = aes(x=Species, y= FoG50Med),
              color="red",pch="_",size=10, alpha=0.8)+
  ylab("FoG50")+expand_limits(x = 0, y = c(0,1))+
  theme_bw()+ labs(subtitle = "FLC")+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) 

DF_all<-ggarrange(DFp1, DFp2, DFp3, DFp4, align = "hv")

ggsave("201203_Candida_DDAFull.jpg", 
       path = here("figures_out"), 
       plot = DF_all,
       width = 6.5,
       height = 3.5)

df %>% 
  group_by(Species) %>% 
  summarise(n=n()) 


df_BA<-df %>% filter(Drug=="BA") %>% rename("BA_RAD50"="RAD50", "BA_FoG50"="FoG50") 
df_FLC<-df %>% filter(Drug=="FLC") %>% rename("FLC_RAD50"="RAD50", "FLC_FoG50"="FoG50") 
Res_Tol<-merge(df_BA, df_FLC, by= c("Strain", "Species", "Type"))%>% select(!c(Drug.x,Drug.y))

R_BAFLC<-Res_Tol %>% ggplot() +
  geom_jitter( mapping = aes(x=FLC_RAD50, y=BA_RAD50, color=Species), size=2, alpha=0.5)+
  geom_smooth(mapping = aes(x=FLC_RAD50, y=BA_RAD50,color=Species, group= Species),
              se=FALSE, method = lm, linetype="dashed",  size=1)+
  ylab("BA Resistance")+xlab("FLC Resistance")+
  scale_x_continuous(breaks = seq(0,30, by=10))+
  theme_bw()+ 
  scale_color_viridis(discrete = TRUE)+expand_limits(x = 0)+
  theme(legend.position="right",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", face="bold", size = 8),axis.text.x = element_text(hjust = 0.5))

T_BAFLC<-Res_Tol %>% ggplot() +
  geom_jitter(mapping = aes(x=FLC_FoG50, y=BA_FoG50, color=Species), size=2, alpha=0.5)+
  geom_smooth(mapping = aes(x=FLC_FoG50, y=BA_FoG50, color=Species, group =Species),
              se=FALSE, method = lm, linetype="dashed",  size=1)+
  xlim(0,1)+ ylim(0,0.5)+
  ylab(expression("BA "*FoG[50]))+xlab(expression("FLC "*FoG[50]))+
  theme_bw()+ 
  scale_color_viridis(discrete = TRUE)+expand_limits(x = 0)+
  theme(legend.position="right",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", face="bold", size = 8))
RT<-ggarrange(R_BAFLC,T_BAFLC, common.legend = TRUE)
ggsave("210301_ResistanceToleranceCor.jpg", 
       path = here("figures_out"), 
       plot = RT,
       width = 5,
       height = 3)
