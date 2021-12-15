library(tidyverse)
library(here)
library(lme4)
library(broom)
library(reshape2)
library(extrafont)
library(viridis)
library(growthcurver)

dfB01 <- read_delim(here("data_in","Biofilm_Formation","210223_Evos_B01_editted2.csv"), delim=",")%>%select(c("Filename","Background","Fungi"))

dfB02 <- read_delim(here("data_in","Biofilm_Formation","210529_Evos_TL_B02.csv"), delim=",")%>%select(c("Filename","Background","Celltype 1")) %>% rename("Fungi" = "Celltype 1")

#B02 has the same plate layout of B01. Can use the info file from B01
Info<- read_delim(here("data_in","Biofilm_Formation","210218_Well_Info.csv"), delim=",")

#col<-c("A02"="#1C8678","A03"= "#ACDB1F","A04"= "#65CA44",
#       "A08" ="#28A26B","A10"= "#2B4279","A12"= "#371663",
#      "A17"= "#FEE51F","A18"= "#370755")

DB01<-dfB01%>%separate(Filename,c("Well","Time"),"_t")
DB02<-dfB02%>%separate(Filename,c("Well","Time"),"_t")


DFB01 <- DB01%>%
  separate(Time, c("Time","Image"),sep=".T")%>%
  select(!"Image")%>%select(!Background) %>%
  filter(!Well%in%c("A12", "B12", "C12",
                    "D12","E12", "F12",
                    "G12", "H12")) %>% 
  merge(Info, by ="Well") %>%  
  select(Strain, Time, Drug, Concentration, Fungi) %>% 
  mutate(Biorep="B01")%>% 
  na.omit()

DFB02<-DB02%>%
  separate(Well, c("Strain","Drug", "Concentration"),sep="_")%>%
  separate(Time, c("Time","Image"),sep=".j")%>%
  select(!"Image")%>%select(!Background) %>%
  select(Strain, Time, Drug, Concentration, Fungi) %>% mutate(Biorep="B02") %>% 
  na.omit()

DFB01$Time<-as.numeric(DFB01$Time)
DFB02$Time<-as.numeric(DFB02$Time)
DFB01$Concentration<-as.factor(DFB01$Concentration)
DFB02$Concentration<-as.factor(DFB02$Concentration)

DF <- rbind(DFB01, DFB02) %>%  
  group_by(Strain, Time, Drug, Concentration) %>%
  summarise( Fungi_sd= sd(Fungi, na.rm = TRUE), Fungi = mean(Fungi))

DFB01 %>% ggplot(aes(x=Time, y=Fungi, color=Concentration))+
  geom_point(aes(), alpha=0.5) + 
  scale_x_continuous(breaks=seq(0,24,4))+
  facet_grid(Strain ~ Drug) +
  scale_color_viridis(discrete = TRUE)+ theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))

DFB02 %>% ggplot(aes(x=Time, y=Fungi, color=Concentration))+
  geom_point(aes(), alpha=0.5) + 
  scale_x_continuous(breaks=seq(0,24,4))+
  facet_grid(Strain ~ Drug) +
  scale_color_viridis(discrete = TRUE)+ theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))




T<- DF %>% unite(Strain, Strain,Drug, Concentration, sep = "_" ) %>% 
  pivot_wider( names_from = "Strain", values_from = "Fungi") %>% 
  rename("time"="Time") 
S<-SummarizeGrowthByPlate(plate=T) 

#d.comb <- d %>% unite(Strain, c("Strain", "Drug", "Concentration"))


######################################################################################
#Applying the model to all of the strains 

summG <- function(x) {SummarizeGrowth(T$time,x)}
lapply(T[2:ncol(T)], summG)
models.all <- lapply(T[2:ncol(T)], summG)


df.predicted.plate <- data.frame(time = T$time)
for (i in names(T[2:ncol(T)])) 
{df.predicted.plate[[i]] <- stats:::predict(models.all[[i]]$model)}

melt1 <- melt(T, id.vars = "time", variable.name = "sample", value.name = "area")
melt2 <- melt(df.predicted.plate, id.vars = "time", variable.name = "sample", value.name = "pred.area")
df.final <- cbind(melt1, pred.area=melt2[,3]) %>% rename("Strain"="sample")
rm(melt1)
rm(melt2)

df_final<-df.final %>% separate("Strain", into= c("Strain", "Drug", "Concentration"), sep="_")




#Taking the BA data out
BA_final<-df_final%>%filter(Drug=="BA"|Drug=="ND")

BA_final$Concentration <- factor(BA_final$Concentration, levels = c("0","0.4","0.8","1.6","3.2","6.4"))



#I am going to plot timepoints 2-24 
BA_final%>%ggplot(aes(x=time, y=area, color=Concentration))+
  geom_point(aes(), alpha=0.5) + geom_line(aes(y=pred.area)) +
  scale_x_continuous(breaks=seq(0,24,4))+
  facet_wrap(~Strain, ncol = 2) +
  scale_color_viridis(discrete = TRUE)+ theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))

ggsave(plot=last_plot(),"210218_BA_LogFit_Timelapse.jpg", height = 7.2, width = 4.8)


#Working on FLC data
FLC_final<-df_final%>%filter(Drug=="FLC"|Drug=="ND")

FLC_final$Concentration<- as.factor(FLC_final$Concentration)

FLC_final%>%
  ggplot(mapping=aes(x=time, y=area, color=Concentration))+
  geom_point(aes(), alpha=0.5, size =2) + ylim(0,1)+
  geom_line(aes(y=pred.area)) +
  scale_x_continuous(breaks=seq(0,24,4))+
  facet_wrap(~Strain, ncol = 2) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))

ggsave(plot=last_plot(),"210218_FLC_LogFit_Timelapse.jpg", 
       height = 7.2, width = 4.8)

#K is the asymptote and r is the slope
S.all<-S %>% rename("Strain"="sample") %>%
  separate("Strain", into= c("Strain", "Drug", "Concentration"), sep="_")


#B<-BA_final %>% filter(Concentration >1.6) %>% 
group_by(Strain, Concentration) %>% 
  do(fit = tidy(lm(Concentration ~ area, data = .))) %>% 
  unnest(fit) %>% filter(term!="(Intercept)") %>% 
  select(Strain, Concentration, estimate) %>% 
  rename("r"="estimate")


BAbreak<-unique(BA_final$Concentration)
S.all$Concentration<-as.numeric(S.all$Concentration)


BA_k<-S.all%>%filter(Drug=="BA"|Drug=="ND")%>%
  ggplot(aes(x=Concentration, y=k, color=Strain, group= Strain))+
  geom_line()+
  
  geom_point(aes(), alpha=0.5, size=2)+ geom_line(size=1)+
  expand_limits(y= c(0,1))+
  scale_color_manual(values = col)+ ylab("Asymptote")+
  xlab("BA Concentration")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))

S.all%>%filter(Drug=="BA"|Drug=="ND")%>% 
  ggplot(aes(x=Concentration, y=r, color=Strain, group = Strain))+
  geom_point(aes(), alpha=0.5, size=2)+ geom_line(size=1)+
  expand_limits(y= c(0,1))+
  scale_color_viridis(discrete = TRUE)+ ylab("Asymptote")+
  xlab("BA Concentration")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))

ggsave(plot= BA_k,"210225_EvosB01_BA_k_Timelapse.jpg", height = 2, width = 3)

S.FLC<-S.all%>%filter(Drug=="FLC"|Drug=="ND")
S.FLC$Concentration <- factor(S.FLC$Concentration, levels = c("0","2","4","8","16","32"))

S.FLC%>%
  ggplot(aes(x=Concentration, y=k, color=Strain, group = Strain))+
  geom_point(aes(), alpha=0.5, size=2)+ 
  geom_line(size=1)+
  expand_limits( y=c(0,1))+
  scale_y_continuous(breaks = seq(0,1, 0.25))+
  
  scale_color_manual(values = col)+   ylab("Asymptote")+
  xlab("FLC Concentration")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text())
ggsave(plot= last_plot(),"210225_EvosB01_FLC_k_Timelapse.jpg", height = 2, width = 3)

S.all%>%filter(Drug=="FLC"|Drug=="ND")%>%
  ggplot(aes(x=Concentration, y=r, color=Strain, group = Strain))+
  geom_point(aes(), alpha=0.5, size=2)+ 
  geom_line()+
  expand_limits(y= c(0,1))+ 
  scale_color_viridis(discrete = TRUE)+ ylab("Slope")+
  xlab("FLC Concentration")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))

ggsave(plot=FLC_k,"210218_FLC_k_Timelapse.jpg", height = 3, width = 4)

