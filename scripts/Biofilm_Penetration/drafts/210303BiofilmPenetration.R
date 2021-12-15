library(here)
library(tidyverse)
library(reshape2)
library(viridis)
library(reshape2)
library(ggforce)
library(scales)
library(extrafont)


# Load data ----

PreBA<-read_delim(here("data_in", "Biofilm_Penetration", "BBA_PreDrug_FullDF.csv"), delim=",")
prebackground<-mean(PreBA$Background)
PreBA<-PreBA%>%select(!Background)
preBA<-melt(PreBA,c("Strain","Biorep"))
preBA$value= preBA$value- prebackground
preBA$value<-ifelse(preBA$value<0, 0,preBA$value)
preBA<-preBA%>%rename(enviro=variable,OD=value)
preBA$enviro=as.numeric(levels(preBA$enviro))[preBA$enviro]
preBA$enviro<-ifelse(preBA$enviro=="0", 0.001,preBA$enviro)

preBA<-preBA%>%group_by(Strain)%>%summarise(mean=mean(OD))

PostBA<-read_delim(here("data_in", "Biofilm_Penetration","BBA_PostDrug_FullDF.csv"), delim=",")
postbackground<-mean(PostBA$Background)
PostBA<-PostBA%>%select(!Background)

postBA<-melt(PostBA,c("Strain","Biorep"))
postBA$value= postBA$value- postbackground
postBA$value<-ifelse(postBA$value<0, 0,postBA$value)

postBA<-postBA%>%rename(enviro=variable,OD=value)
postBA$enviro=as.numeric(levels(postBA$enviro))[postBA$enviro]
postBA$enviro<-ifelse(postBA$enviro==0,0.390625/2 ,postBA$enviro)

#extracting the diffierent concentrations of BA. I have 7 different concentrations of BA
BREAKBA<-unique(postBA$enviro)
LABBA<-unique(postBA$enviro)%>%round(digits = 2)
LABBA[7]<-0

#Normalizing the BA data
BA.ag<-aggregate(data=postBA,postBA["OD"],postBA[c("Strain","enviro")],mean)

BA.norm<-merge(preBA,BA.ag,by="Strain")%>%
  group_by(Strain,enviro)%>%summarise(OD=OD/mean) %>% 
  rename(Enviro = enviro)

BA50 <-read_delim(here("data_in", "Planktonic_MIC50","200721_BA_MIC50_editted.csv"),delim = ",")

BA_MIC_biomass<-merge(BA.norm, BA50, by="Strain")
BA_MIC_biomass$Enviro <- round(BA_MIC_biomass$Enviro, 2) 

#save the normalized BA data 
write_csv(BA50norm, here("data_out", "Biofilm_Penetration", "210303BA_MIC_biomass.csv"))

#Importing FLC data
PreFLC<-read_delim(here("data_in", "Biofilm_Penetration","BFLC_PreDrug_FullDF.csv"), delim=",")
prebackground<-mean(PreFLC$Background)
PreFLC<-PreFLC%>%select(!Background)
preFLC<-melt(PreFLC,c("Strain","Biorep"))
preFLC$value= preFLC$value- prebackground
preFLC$value<-ifelse(preFLC$value<0, 0,preFLC$value)
preFLC<-preFLC%>%rename(enviro=variable,OD=value)
preFLC$enviro=as.numeric(levels(preFLC$enviro))[preFLC$enviro]
preFLC$enviro<-ifelse(preFLC$enviro=="0", 0.001,preFLC$enviro)

preFLC<-preFLC%>%group_by(Strain)%>%summarise(mean=mean(OD))

PostFLC<-read_delim(here("data_in", "Biofilm_Penetration","BFLC_PostDrug_FullDF.csv"), delim=",")
postbackground<-mean(PostFLC$Background)
PostFLC<-PostFLC%>%select(!Background)
postFLC<-melt(PostFLC,c("Strain","Biorep"))
postFLC$value= postFLC$value- postbackground
postFLC$value<-ifelse(postFLC$value<0, 0,postFLC$value)
postFLC<-postFLC%>%rename(enviro=variable,OD=value)
postFLC$enviro=as.numeric(levels(postFLC$enviro))[postFLC$enviro]
postFLC$enviro<-ifelse(postFLC$enviro==0,0.5 ,postFLC$enviro)


#aggregating and normalizing the FLC data 
postflc.ag<-aggregate(data=postFLC,postFLC["OD"],postFLC[c("Strain","Biorep","enviro")],mean)

FLC.ag<-aggregate(data=postflc.ag,postflc.ag["OD"],postflc.ag[c("Strain","enviro")],mean)

postFLC.norm<-merge(preFLC,FLC.ag,by="Strain")%>%filter(Strain!="A52")%>%
  group_by(Strain,enviro)%>%summarise(OD=mean(OD/mean))

BREAKFLC<-unique(postFLC$enviro)
LABFLC<-unique(postFLC$enviro)%>%round(digits = 2)
LABFLC[7]<-0

FLCMIC<- read_delim(here("data_in", "Planktonic_MIC50","200721_FLC_MIC50_editted.csv"),delim=",")
FLC_MIC_biomass <-merge(postFLC.norm,FLCMIC, by= "Strain")
FLC_MIC_biomass$Enviro <- round(FLC_MIC_biomass$Enviro, 2) 
names(FLC_MIC_biomass)[2] <- "Enviro"

#save the normalized BA data 
write_csv(FLC50norm, here("data_out", "Biofilm_Penetration", "210303FLC_MIC_biomass.csv"))

#adding XTT data 
XTTBA<-read_delim(here("data_in", "Biofilm_Penetration","BBA_XTT_FullDF.csv"),
                delim=",")
xttBAbackground<-mean(XTTBA$Background)
XTTBA<-XTTBA%>%select(!Background)
xttBA<-melt(XTTBA,c("Strain","Biorep"))
xttBA$value= xttBA$value- xttBAbackground
xttBA$value<-ifelse(xttBA$value<0, 0.0625,xttBA$value)
xttBA<-xttBA%>%rename(enviro=variable,activity=value)
xttBA$enviro=as.numeric(levels(xttBA$enviro))[xttBA$enviro]
xttBA$enviro<-ifelse(xttBA$enviro==0,0.39/2,xttBA$enviro)

#normalizing the BA XTT data
BA.ND <- xttBA%>%filter(enviro==min(enviro))
normBA <- merge(xttBA,BA.ND, by=c("Strain","Biorep"))%>%
  rename(activity=activity.x, min.env=enviro.y, max.activity=activity.y)%>%
  group_by(Strain, enviro.x)%>%
  summarise(Activity_norm = mean(activity/max.activity)) %>% 
  rename(Enviro = enviro.x)
normBA$Enviro <- round(normBA$Enviro, 2)  

#FLCXTT
XTTFLC <- read_delim(here("data_in", "Biofilm_Penetration","BFLC_XTT_FullDF.csv"), delim=",")
xttFLCbackground <- mean(XTTFLC$Background)
XTTFLC <- XTTFLC%>%select(!Background)
xttFLC <- melt(XTTFLC,c("Strain","Biorep"))
xttFLC$value = xttFLC$value- xttFLCbackground
xttFLC$value <-ifelse(xttFLC$value<0, 0.0625,xttFLC$value)
xttFLC <- xttFLC%>%rename(enviro=variable,activity=value)
xttFLC$enviro <- as.numeric(levels(xttFLC$enviro))[xttFLC$enviro]
xttFLC$enviro <- ifelse(xttFLC$enviro==0,0.5,xttFLC$enviro)

#Normalizing
FLC.ND<-xttFLC%>%filter(enviro==min(enviro))

normFLC <- merge(xttFLC,FLC.ND, by=c("Strain","Biorep"))%>%
  filter(Strain!="A61")%>%
  rename(activity=activity.x, min.env=enviro.y, max.activity=activity.y)%>%
  group_by(Strain, enviro.x)%>%
  summarise(Activity_norm=mean(activity/max.activity)) %>% 
  rename(Enviro = enviro.x) 
normFLC$Enviro <- round(normFLC$Enviro, 2)  


BA_MIC_activity <- merge(normBA,BAMIC, by="Strain")
FLC_MIC_activity <- merge(normFLC,FLCMIC, by="Strain")

write_csv(BA_MIC_activity,here("data_out", "Biofilm_Penetration", "210303BA_MIC_activity.csv"))
write_csv(FLC_MIC_activity,here("data_out", "Biofilm_Penetration", "210303FLC_MIC_activity.csv"))

# Files to look at ----

BA_MIC_activity_biomass <- merge(BA_MIC_biomass, BA_MIC_activity, by=c("Strain", "Enviro", "MIC50"))
FLC_MIC_activity_biomass <- merge(FLC_MIC_biomass, FLC_MIC_activity, by=c("Strain", "Enviro", "MIC50"))

lmtest <- lm(BA_MIC_activity_biomass$OD~BA_MIC_activity_biomass$Activity_norm*BA_MIC_activity_biomass$Enviro)

as.numeric(unique(BA_MIC_activity_biomass$Enviro))
# 0.20  0.39  0.78  1.56 12.50  3.12  6.25

sub0.20 <- subset(BA_MIC_activity_biomass, Enviro == "0.20")
plot(sub0.20$OD, sub0.20$Activity_norm)
t020 <- cor.test(sub0.20$OD, sub0.20$Activity_norm)

sub039 <- subset(BA_MIC_activity_biomass, Enviro == "0.39")
plot(sub39$OD, sub$Activity_norm)
t039 <- cor.test(sub039$OD, sub039$Activity_norm)

sub078 <- subset(BA_MIC_activity_biomass, Enviro == "0.78")
plot(sub078$OD, sub078$Activity_norm)
t078 <- cor.test(sub078$OD, sub078$Activity_norm)

sub1.56 <- subset(BA_MIC_activity_biomass, Enviro == "1.56")
plot(sub1.56$OD, sub1.56$Activity_norm)
t156 <- cor.test(sub1.56$OD, sub1.56$Activity_norm)

sub3.12 <- subset(BA_MIC_activity_biomass, Enviro == "3.12")
plot(sub3.12$OD, sub3.12$Activity_norm)
t312 <- cor.test(sub3.12$OD, sub3.12$Activity_norm)

sub6.25 <- subset(BA_MIC_activity_biomass, Enviro == "6.25")
plot(sub6.25$OD, sub6.25$Activity_norm)
t625 <- cor.test(sub6.25$OD, sub6.25$Activity_norm)

sub12.50 <- subset(BA_MIC_activity_biomass, Enviro == "0.20")
plot(sub0.20$OD, sub0.20$Activity_norm)
t020 <- cor.test(sub0.20$OD, sub0.20$Activity_norm)


#Ola to do : add in ALL strains!
#Run all drug levels
#Create table to show correlations

S1 <- subset(BA_MIC_activity_biomass, Strain == "A01" & Enviro != "0.20", Enviro != "12.50")
par(mfrow=c(2, 1), mar=c(1,1,1,1))
plot(S1$Enviro, S1$OD)
plot(S1$Enviro, S1$Activity_norm)

S312 <- subset(BA_MIC_activity_biomass, Enviro == "3.12")

#looking the the correlation coefficient of biomass and activity
postBA%>%aggrecor.test(postBA$OD,xttBA$activity)
cor.test(postFLC$OD,xttFLC$activity)


# Plot data ----
#Plotting normalized BA data 
BAn<-
  ggplot(BA_MIC_biomass,mapping= aes(x=enviro, y=OD, color=MIC50, group=Strain )) +
  geom_line(stat="smooth", method="loess", 
            mapping= aes(x=enviro, y=OD,color=MIC50),se=FALSE,alpha=0.5,
            size=1.5, span=0.8)+
  geom_hline(mapping= aes(yintercept =1), linetype= "dashed", color="Red")+
  scale_x_continuous( trans = 'log2', #trans=reverselog_trans(2),
                      breaks=BREAKBA,
                      labels = LABBA )+
  xlab("BA concentration (mg/ml)")+
  ylab("Biomass OD600")+expand_limits(x = 0.5, y = c(0,2.1))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   text=element_text(family="Times New Roman", face="bold", size = 12),
                   legend.position="none",axis.text.x = element_text( vjust = 1.5, size =12))+ 
  scale_color_viridis()


FLCn<-ggplot(FLC_MIC_biomass,mapping= aes(x=enviro, y=OD, group=Strain )) +
  geom_line(stat="smooth", method="loess", 
            mapping= aes(x=enviro, y=OD,color=MIC50),se=FALSE,alpha=0.5,
            size=1.5, span=0.8)+
  geom_hline(mapping= aes(yintercept =1), linetype= "dashed", color="Red")+
  scale_x_continuous( trans = 'log2', #trans=reverselog_trans(2),
                      breaks=BREAKFLC,
                      labels = LABFLC )+
  xlab("FLC concentration (ug/ml)")+
  ylab("Biomass OD600")+expand_limits(x = 0, y = c(0,2.1))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   text=element_text(family="Times New Roman", face="bold", size = 12),
                   legend.position="none",axis.text.x = element_text( vjust = 1.5, size =12))+ 
  scale_color_viridis(discrete = TRUE)
#saving the FLC normalized figure

ggsave("201118_FLC_Biomass.png", plot=FLCn, height = 4, width = 4)


ggplot(BA_MIC_activity,mapping= aes(x=enviro.x, y=norm,color=Strain )) +
  geom_smooth(se=FALSE, span=0.8)+
  scale_x_continuous( trans = 'log2', #trans=reverselog_trans(2),
                      breaks=BREAKBA,
                      labels = LABBA )+
  coord_cartesian(clip = "off")+
  xlab("BA concentration (mg/ml)")+
  ylab("Biofilm Activity")+
  theme_bw()+theme(axis.text.x = element_text(angle = -90),legend.position = "None")+ 
  scale_color_viridis(discrete = TRUE)

#ploting all the strains and making the lines more transparent

NormBA<-ggplot(ActMIC,mapping= aes(x=enviro.x, y=norm )) +
  geom_line(stat="smooth", method="loess", 
            mapping= aes(x=enviro.x, y=norm,color=MIC50,group=Strain),se=FALSE,alpha=0.5,
            size=1.5, span=0.8)+
  scale_x_continuous( trans = 'log2', #trans=reverselog_trans(2),
                      breaks=BREAKBA,
                      labels = LABBA )+
  xlab("BA concentration (mg/ml)")+
  ylab("Biofilm Activity")+expand_limits(x = 0, y=c(0.5,1.2))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   text=element_text(family="Times New Roman", face="bold", size = 12),
                   legend.position="none",axis.text.x = element_text( vjust = 1.5, size =12))+ 
  scale_color_viridis(discrete = TRUE )

ggsave("201118_BBA_norm_XTT.png", plot=NormBA, height = 4, width = 4)

NormFLC<-ggplot(FMICact,mapping= aes(x=enviro.x, y=norm,group=Strain )) +
  geom_line(stat="smooth", method="loess", 
            mapping= aes(x=enviro.x, y=norm,color=MIC50),se=FALSE,alpha=0.5,
            size=1.5, span=0.8)+
  scale_x_continuous( trans = 'log2', #trans=reverselog_trans(2),
                      breaks=BREAKFLC,
                      labels = LABFLC )+
  xlab("FLC concentration (ug/ml)")+
  ylab("Biofilm Activity")+expand_limits(x = 0, y=c(0.5,1.2))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   text=element_text(family="Times New Roman", face="bold", size = 12),
                   legend.position="none",axis.text.x = element_text( vjust = 1.5, size =12))+ 
  scale_color_viridis(discrete = TRUE)

ggsave("201118_BFLC_norm_XTT.png", plot=NormFLC, height = 4, width = 4)
