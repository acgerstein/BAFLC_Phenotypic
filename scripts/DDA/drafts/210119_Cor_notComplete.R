install.packages("corr")
library(tidyverse)
library(corrr)

H1Info<-H1Info<-read_csv(here("data_in","2014_Isolate list_info.csv")) %>% 
  select(HSC, 'Isolate ID', SPECIMEN) %>% 
  rename("Strain" = "HSC", "Species"='Isolate ID', "Origin"="SPECIMEN") %>% 
  mutate(Year="2014")

H1B1BA <- read_csv(here("data_in","Updated_HSC_Analysis",
                        "210118_H1B1_BA.csv")) %>%
  separate(name, c("Strain", "Type","Rep")) %>% 
  mutate(Drug="BA") %>% 
  select(Strain,Drug,Rep, RAD50, FoG50) %>% 
  merge(H1Info) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis", "C.tropicalis"))


H1B1FLC <- read_csv(here("data_in","Updated_HSC_Analysis",
                         "210118_H1B1_FLC_AA.csv")) %>%
  separate(name, c("Strain","Drug", "Rep")) %>%
  select(Strain,Drug, Rep, RAD50, FoG50) %>% 
  subset(Strain%in%H1B1BA$Strain)%>% 
  merge(H1Info) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis",
                      "C.tropicalis"))

H1B1 <-rbind(H1B1BA,H1B1FLC) %>% mutate(Biorep="B1")

#HSC50 is missing from biorep2
#HSC59,61,71 is missing from Biorep 1

H1B2 <- read_csv(here("data_in","HSC001-105_48h_B2_df.csv")) %>%
  separate(name, c("Strain", "Drug","Rep")) %>%
  select(Strain,Drug, Rep, RAD50, FoG50) %>% filter(!Strain=="Blank") %>% 
  subset(Strain%in%H1B1$Strain)%>% 
  merge(H1Info) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata",
                      "C.parapsilosis", "C.tropicalis")) %>% 
  mutate(Biorep="B2")

H1B1<- H1B1 %>% subset(Strain%in%H1B2$Strain)


#Second Batch of HSC 
H2Info<- read_csv(here("data_in","Aleeza_HSC_Yeast stocks_2018.csv" )) %>% 
  select('HSC#', "SPECIMEN", "ORGANISM") %>% 
  rename("Strain"='HSC#', "Species"="ORGANISM", "Origin"="SPECIMEN") %>% 
  mutate(Year="2018")

H2B1pt1<- read_csv(here("data_in","HSC106_159_48h_df_p1_B1.csv")) %>% 
  filter(!line=="Blank") %>%
  separate(name, c("Strain", "Drug","Rep")) %>%
  select(Strain, Drug, Rep, RAD50, FoG50) %>%
  merge(H2Info, by = "Strain")

H2B1pt2<- read_csv(here("data_in","HSC106_159_48h_df_p2_B1.csv")) %>% 
  filter(!line=="Blank") %>%
  rename("Strain"="line", "Drug"="type") %>% 
  separate(name, c("Strain", "Drug","Rep")) %>%
  select(Strain, Drug, Rep, RAD50, FoG50) %>%
  merge(H2Info, by = "Strain")



H2B1<- rbind(H2B1pt1, H2B1pt2) %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis",
                      "C.tropicalis"))%>% mutate(Biorep="B1")

H2B2<- read_csv(here("data_in","Updated_HSC_Analysis",
                     "H2B2_48h_df.csv")) %>% 
  filter(!line=="Blank") %>%
  separate(name, c("Strain", "Drug","Rep")) %>%
  select(Strain, Drug, Rep, RAD50, FoG50) %>%

  merge(H2Info, by = "Strain") %>% 
  filter(Species%in%c("C.albicans", "C.glabrata", "C.parapsilosis", 
                      "C.tropicalis")) %>% 
  subset(Strain%in%H2B1$Strain)%>% mutate(Biorep="B2")

H2B1<-H2B1 %>%  subset(Strain%in%H2B2$Strain)


Info<-rbind(H1Info, H2Info)

Comb1418 <- read_csv(here("data_in",
                          "DDA201218_48h_df.csv")) %>%
  filter(!line=="Blank") %>%
  separate(name, c("Strain", "Drug","Rep")) %>%
  select(Strain,Drug,Rep, RAD50, FoG50) %>%  
  merge(Info) %>% mutate(Type="Comb")


H1Ball<-merge(H1B1,H1B2, by= c("Strain", "Drug","Rep", "Year", "Origin", "Species"))

H2Ball<-merge(H2B1, H2B2, by =c("Strain", "Drug","Rep", "Year", "Origin", "Species"))

HSCAllB<-rbind(H1Ball,H2Ball) %>% rename("RAD50.B1"="RAD50.x", "RAD50.B2"="RAD50.y", "FoG50.B1"="FoG50.x", "FoG50.B2"="FoG50.y") %>% select(Strain, Drug,Rep, Year, RAD50.B1, RAD50.B2, FoG50.B1, FoG50.B2)

HSCAllSave<-merge(HSCAllB, Comb1418, by = c("Strain", "Drug", "Rep", "Year"))%>%
  select(Strain, Drug,Species, Rep, Year, RAD50.B1, RAD50.B2, RAD50, FoG50.B1, FoG50.B2, FoG50) %>% rename("RAD50.Comb"="RAD50", "FoG50.Comb"="FoG50")

write_csv(HSCAllSave, "210118_HSCAllBatches.csv")


Hall<-HSCAllSave %>% group_by(Strain, Drug, Species, Year) %>% 
  summarise(RAD50.B1= mean(RAD50.B1), RAD50.B2=mean(RAD50.B2), RAD50.Comb=mean(RAD50.Comb),
            FoG50.B1=mean(FoG50.B1), FoG50.B2=mean(FoG50.B2),
            FoG50.Comb=mean(FoG50.Comb))

T<-Hall %>% group_by(Strain, Drug, Species, Year) %>% 
  summarise(RAD50.B= (RAD50.B1+RAD50.B1)/2, RAD50.Comb = RAD50.Comb,
            FoG50.B=(FoG50.B1+FoG50.B2)/2, FoG50.Comb=FoG50.Comb)

Hall %>% group_by(Species, Drug) %>% summarise(Correlation= cor(RAD50.B2, RAD50.Comb))


Hall %>% ggplot(aes(x= FoG50.B1, y=FoG50.B2))+geom_point()+
  expand_limits(x=c(0,1), y=c(0,1))+
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE)+
  facet_grid(Drug~Species)


Hall %>% ggplot(aes(x= RAD50.B2, y=RAD50.Comb))+geom_point()+
  expand_limits(x=c(0,1), y=c(0,1))+
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE)+
  facet_grid(Drug~Species*Year)

Hall %>% ggplot(aes(x= FoG50.B1, y=FoG50.Comb))+geom_point()+
  expand_limits(x=c(0,1), y=c(0,1))+
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE)+
  facet_grid(Drug~Species*Year)

Hall_BA<-Hall %>% filter(Drug=="BA") %>% filter(Species =="C.albicans")
Hall_FLC<-Hall %>% filter(Drug=="FLC" ) %>% filter(Species =="C.albicans")

cor.test(Hall_BA$RAD50.B1, Hall_BA$RAD50.B2)
cor.test(Hall_BA$RAD50.B2, Hall_BA$RAD50.Comb)
cor.test(Hall_BA$RAD50.B1, Hall_BA$RAD50.Comb)

cor.test(Hall_BA$FoG50.B1, Hall_BA$FoG50.B2)
cor.test(Hall_BA$FoG50.B2, Hall_BA$FoG50.Comb)
cor.test(Hall_BA$FoG50.B1, Hall_BA$FoG50.Comb)


cor.test(Hall_FLC$RAD50.B1, Hall_FLC$RAD50.B2)
cor.test(Hall_FLC$RAD50.B2, Hall_FLC$RAD50.Comb)
cor.test(Hall_FLC$RAD50.B1, Hall_FLC$RAD50.Comb)




cor.test(Hall_FLC$FoG50.B1, Hall_FLC$FoG50.B2)
cor.test(Hall_FLC$FoG50.B2, Hall_FLC$FoG50.Comb)
cor.test(Hall_FLC$FoG50.B1, Hall_FLC$FoG50.Comb)





H1Ball<-merge(H1B1,H1B2, by= c("Strain", "Drug","Rep", "Year", "Origin", "Species")) %>% 


H2Ball<-merge(H2B1, H2B2, by =c("Strain", "Drug","Rep", "Year", "Origin", "Species"))

HSCAllB<-rbind(H1Ball,H2Ball) %>% rename("RAD50.B1"="RAD50.x", "RAD50.B2"="RAD50.y", "FoG50.B1"="FoG50.x", "FoG50.B2"="FoG50.y") %>% select(Strain, Drug,Rep, RAD50.B1, RAD50.B2, FoG50.B1, FoG50.B2)

HSCAllSave<-merge(HSCAllB, Comb1418, by = c("Strain", "Drug", "Rep")) %>% select(Strain, Drug,Species, Rep, RAD50.B1, RAD50.B2, RAD50, FoG50.B1, FoG50.B2, FoG50) %>% rename("RAD50.Comb"="RAD50", "FoG50.Comb"="FoG50")

write_csv(HSCAllSave, "210118_HSCAllBatches.csv")


Hall<-HSCAllSave %>% group_by(Strain, Drug, Species) %>% 
  summarise(RAD50.B1= mean(RAD50.B1), RAD50.B2=mean(RAD50.B2), RAD50.Comb=mean(RAD50.Comb),
            FoG50.B1=mean(FoG50.B1), FoG50.B2=mean(FoG50.B2),
            FoG50.Comb=mean(FoG50.Comb))
