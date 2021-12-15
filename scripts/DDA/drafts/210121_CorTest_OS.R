library(here)
library(tidyverse)
library(broom)

#Loading in the cleaned data
d<-read_csv(here("210121_HSCAllBatches_cleaned.csv"))

#B1 -> biorep1, Comb -> biorep2, Comb -> (combined) HSC 2014 and 2018 strains
#processed in the same batch 

#combining the technical reps for each biorep 
df<-d %>% group_by(Strain, Drug, Species) %>% 
  summarise(RAD50.B1= mean(RAD50.B1, na.rm = TRUE),
            RAD50.B2=mean(RAD50.B2, na.rm = TRUE),
            RAD50.Comb=mean(RAD50.Comb, na.rm = TRUE),
            FoG50.B1=mean(FoG50.B1),
            FoG50.B2=mean(FoG50.B2, na.rm = TRUE),
            FoG50.Comb=mean(FoG50.Comb, na.rm = TRUE))

#Plotting Resistance 
RB1B2<-df %>% ggplot(aes(x= RAD50.B1, y=RAD50.B2))+geom_point()+
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE)+
  facet_grid(Drug~Species) + theme_bw()

RB1Comb<-df %>% ggplot(aes(x= RAD50.B1, y=RAD50.Comb))+geom_point()+
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE)+
  facet_grid(Drug~Species) + theme_bw()

RB2Comb<-df %>% ggplot(aes(x= RAD50.B2, y=RAD50.Comb))+geom_point()+
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE)+
  facet_grid(Drug~Species) + theme_bw()


#ploting tolernce 
TB1B2<-df %>% ggplot(aes(x= FoG50.B1, y=FoG50.B2))+geom_point()+
  expand_limits(x=c(0,1), y=c(0,1))+
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE)+
  facet_grid(Drug~Species) + theme_bw()

TB1Comb<-df %>% ggplot(aes(x= FoG50.B1, y=FoG50.Comb))+geom_point()+
  expand_limits(x=c(0,1), y=c(0,1))+
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE)+
  facet_grid(Drug~Species) + theme_bw()

TB2Comb<- df %>% ggplot(aes(x= FoG50.B2, y=FoG50.Comb))+geom_point()+
  expand_limits(x=c(0,1), y=c(0,1))+
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE)+
  facet_grid(Drug~Species) + theme_bw()


#Correlation test 
df_BA<-df%>% filter(Drug=="BA")
df_FLC<-df %>% filter(Drug=="FLC" )

#Boric Acid
#Biorep1 vs Biorep 2
BARB1B2_Calb<-tidy(cor.test(data = df_BA,
              ~ RAD50.B1+ RAD50.B2,
              subset=(Species=="C.albicans"))) %>%
  mutate(Species="C.albicans", Test ="B1B2 Resistance")

BARB1B2_Cglb<-tidy(cor.test(data = df_BA,
              ~ RAD50.B1+ RAD50.B2,
              subset=(Species=="C.glabrata")))%>%
  mutate(Species="C.glabrata", Test ="B1B2 Resistance")

BATB1B2_Calb<-tidy(cor.test(data = df_BA,
                          ~ FoG50.B1+ FoG50.B2,
                          subset=(Species=="C.albicans")))%>%
  mutate(Species="C.albicans", Test ="B1B2  Tolerance")

BATB1B2_Cglb<-tidy(cor.test(data = df_BA,
                          ~ FoG50.B1+ FoG50.B2,
                          subset=(Species=="C.glabrata")))%>%
  mutate(Species="C.glabrata", Test ="B1B2  Tolerance")

BAB1B2<-rbind(BARB1B2_Calb, BARB1B2_Cglb, BATB1B2_Calb, BATB1B2_Cglb )


#Biorep1 vs Comb
BARB1Comb_Calb<-tidy(cor.test(data = df_BA,
                         ~ RAD50.B1+ RAD50.Comb,
                         subset=(Species=="C.albicans")))%>%
  mutate(Species="C.albicans", Test ="B1Comb Resistance")

BARB1Comb_Cglb<-tidy(cor.test(data = df_BA,
                         ~ RAD50.B1+ RAD50.Comb,
                         subset=(Species=="C.glabrata")))%>%
  mutate(Species="C.albicans", Test ="B1Comb Resistance")

BATB1Comb_Calb<-tidy(cor.test(data = df_BA,
                          ~ FoG50.B1+ FoG50.Comb,
                          subset=(Species=="C.albicans")))%>%
  mutate(Species="C.albicans", Test ="B1Comb Tolerance")

BATB1Comb_Cglb<-tidy(cor.test(data = df_BA,
                          ~ FoG50.B1+ FoG50.Comb,
                          subset=(Species=="C.glabrata")))%>%
  mutate(Species="C.glabrata", Test ="B1Comb Tolerance")

BAB1Comb<-rbind(BARB1Comb_Calb, BARB1Comb_Cglb, BATB1Comb_Calb, BATB1Comb_Cglb )

#Biorep2 vs Comb
BARB2Comb_Calb<-tidy(cor.test(data = df_BA,
                         ~ RAD50.B2+ RAD50.Comb,
                         subset=(Species=="C.albicans")))%>%
  mutate(Species="C.albicans", Test ="B2Comb Resistance")

BARB2Comb_Cglb<-tidy(cor.test(data = df_BA,
                         ~ RAD50.B2+ RAD50.Comb,
                         subset=(Species=="C.glabrata")))%>%
  mutate(Species="C.albicans", Test ="B2Comb Resistance")

BATB2Comb_Calb<-tidy(cor.test(data = df_BA,
                            ~ FoG50.B2+ FoG50.Comb,
                            subset=(Species=="C.albicans")))%>%
  mutate(Species="C.albicans", Test ="B2Comb Tolerance")

BATB2Comb_Cglb<-tidy(cor.test(data = df_BA,
                            ~ FoG50.B2+ FoG50.Comb,
                            subset=(Species=="C.glabrata")))%>%
  mutate(Species="C.glabrata", Test ="B2Comb  Tolerance")

BAB2Comb<-rbind(BARB2Comb_Calb, BARB2Comb_Cglb, BATB2Comb_Calb, BATB2Comb_Cglb )



BAAllCor<-rbind(BAB1B2, BAB1Comb, BAB2Comb) %>% mutate(Drug= "BA")
colnames(BAAllCor)
BAAllCor<-BAAllCor[, c(10, 9, 11, 1, 2, 3, 4, 5, 6, 7, 8)]
#########################################################


#FLC
#Biorep1 vs Biorep 2
FLCRB1B2_Calb<-tidy(cor.test(data = df_FLC,
                            ~ RAD50.B1+ RAD50.B2,
                            subset=(Species=="C.albicans"))) %>%
  mutate(Species="C.albicans", Test ="B1B2 Resistance")

FLCRB1B2_Cglb<-tidy(cor.test(data = df_FLC,
                            ~ RAD50.B1+ RAD50.B2,
                            subset=(Species=="C.glabrata")))%>%
  mutate(Species="C.glabrata", Test ="B1B2 Resistance")

FLCTB1B2_Calb<-tidy(cor.test(data = df_FLC,
                            ~ FoG50.B1+ FoG50.B2,
                            subset=(Species=="C.albicans")))%>%
  mutate(Species="C.albicans", Test ="B1B2 Tolerance")

FLCTB1B2_Cglb<-tidy(cor.test(data = df_FLC,
                            ~ FoG50.B1+ FoG50.B2,
                            subset=(Species=="C.glabrata")))%>%
  mutate(Species="C.glabrata", Test ="B1B2 Tolerance")

FLCB1B2<-rbind(FLCRB1B2_Calb, FLCRB1B2_Cglb, FLCTB1B2_Calb, FLCTB1B2_Cglb )


#Biorep1 vs Comb
FLCRB1Comb_Calb<-tidy(cor.test(data = df_FLC,
                              ~ RAD50.B1+ RAD50.Comb,
                              subset=(Species=="C.albicans")))%>%
  mutate(Species="C.albicans", Test ="B1Comb Resistance")

FLCRB1Comb_Cglb<-tidy(cor.test(data = df_FLC,
                              ~ RAD50.B1+ RAD50.Comb,
                              subset=(Species=="C.glabrata")))%>%
  mutate(Species="C.albicans", Test ="B1Comb Resistance")

FLCTB1Comb_Calb<-tidy(cor.test(data = df_FLC,
                              ~ FoG50.B1+ FoG50.Comb,
                              subset=(Species=="C.albicans")))%>%
  mutate(Species="C.albicans", Test ="B1Comb Tolerance")

FLCTB1Comb_Cglb<-tidy(cor.test(data = df_FLC,
                              ~ FoG50.B1+ FoG50.Comb,
                              subset=(Species=="C.glabrata")))%>%
  mutate(Species="C.glabrata", Test ="B1Comb Tolerance")

FLCB1Comb<-rbind(FLCRB1Comb_Calb,
                 FLCRB1Comb_Cglb, FLCTB1Comb_Calb, FLCTB1Comb_Cglb )

#Biorep2 vs Comb
FLCRB2Comb_Calb<-tidy(cor.test(data = df_FLC,
                              ~ RAD50.B2+ RAD50.Comb,
                              subset=(Species=="C.albicans")))%>%
  mutate(Species="C.albicans", Test ="B2Comb Resistance")

FLCRB2Comb_Cglb<-tidy(cor.test(data = df_FLC,
                              ~ RAD50.B2+ RAD50.Comb,
                              subset=(Species=="C.glabrata")))%>%
  mutate(Species="C.albicans", Test ="B2Comb Resistance")

FLCTB2Comb_Calb<-tidy(cor.test(data = df_FLC,
                              ~ FoG50.B2+ FoG50.Comb,
                              subset=(Species=="C.albicans")))%>%
  mutate(Species="C.albicans", Test ="B2Comb Tolerance")

FLCTB2Comb_Cglb<-tidy(cor.test(data = df_FLC,
                              ~ FoG50.B2+ FoG50.Comb,
                              subset=(Species=="C.glabrata")))%>%
  mutate(Species="C.glabrata", Test ="B2Comb Tolerance")

FLCB2Comb<-rbind(FLCRB2Comb_Calb, FLCRB2Comb_Cglb, 
                 FLCTB2Comb_Calb, FLCTB2Comb_Cglb )



FLCAllCor<-rbind(FLCB1B2, FLCB1Comb, FLCB2Comb) %>% mutate(Drug= "FLC")
colnames(FLCAllCor)
FLCAllCor<-FLCAllCor[, c(10, 9, 11, 1, 2, 3, 4, 5, 6, 7, 8)]

#Combining BA and FLC results 

AllCor<- rbind(BAAllCor, FLCAllCor)

write_csv( AllCor, "210121_HSC_All_CorResults.csv")
