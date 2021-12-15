library(here)
library(tidyverse)
library(reshape2)
library(viridis)
library(reshape2)
library(ggforce)
library(scales)
library(extrafont)


# Load data ----
PreA02A18BA <- read_delim(here("data_in", "Biofilm_Penetration", "A02A18_BBA_PreDrug.csv"), delim=",") %>% select(!Drug)


PreBA <- read_delim(here("data_in", "Biofilm_Penetration", "BBA_PreDrug_FullDF.csv"), delim=",") %>% rbind(PreA02A18BA)

prebackground <- mean(PreBA$Background)
PreBA <- PreBA%>%select(!Background)
preBA <- melt(PreBA,c("Strain","Biorep"))
preBA$value = preBA$value - prebackground
preBA$value <- ifelse(preBA$value<0, 0, preBA$value)
preBA <- preBA %>% rename(enviro = variable, OD = value)
preBA$enviro = as.numeric(levels(preBA$enviro))[preBA$enviro]
preBA$enviro <- ifelse(preBA$enviro=="0", 0.001,preBA$enviro)


preBA <- preBA%>% group_by(Strain) %>% 
  summarise(mean=mean(OD))

############################################
PostA02A18BA <-  read_delim(here("data_in", "Biofilm_Penetration", 
                               "A02A18_BBA_PostDrug.csv"), delim=",") %>% select(!Drug)

PostBA <- read_delim(here("data_in",
                          "Biofilm_Penetration","BBA_PostDrug_FullDF.csv"), delim=",") %>%
  rbind(PostA02A18BA)


postbackground <- mean(PostBA$Background)
PostBA <- PostBA %>% select(!Background)

postBA <- melt(PostBA, c("Strain", "Biorep"))
postBA$value = postBA$value - postbackground
postBA$value <-ifelse(postBA$value < 0, 0, postBA$value)



postBA <- postBA %>% rename(enviro = variable, OD = value)
postBA$enviro = as.numeric(levels(postBA$enviro))[postBA$enviro]

#Which strains did not form biofilms 
BANoBiofilm <- postBA %>% filter(enviro == 0 & OD < 0.1) %>% select(Strain)

postBA <- postBA %>% filter(!Strain%in%BANoBiofilm$Strain)

postBA$enviro <- ifelse(postBA$enviro == 0, 0.390625/2, postBA$enviro)




#extracting the different concentrations of BA. I have 7 different concentrations of BA
BREAKBA <- unique(postBA$enviro)
LABBA <- c(12.8, 6.4, 3.2, 1.6, 0.8, 0.4, 0)


#Normalizing the BA data 
BA.ag <- postBA %>% group_by(Strain, enviro) %>% 
  summarise(OD = mean(OD))
  

BA.norm <- merge(preBA, BA.ag, by = "Strain") %>%
  group_by(Strain, enviro) %>% summarise(OD = OD-mean) %>% 
  rename(Enviro = enviro)


BAA02A18MIC <- read_delim(here("data_in", "Planktonic_MIC50", "A02A18_BA_MIC50.csv"),delim = ",") %>% select(!OD)

BAMIC <- read_delim(here("data_in", "Planktonic_MIC50", "200721_BA_MIC50_editted.csv"),delim = ",") %>% rbind(BAA02A18MIC) %>% mutate(MIC50 = round(MIC50, 1)) 

BA_MIC_biomass <- merge(BA.norm, BAMIC, by="Strain")
BA_MIC_biomass$Enviro <- round(BA_MIC_biomass$Enviro, 2) 



#Importing FLC data
PreA02A18FLC <- read_delim(here("data_in", "Biofilm_Penetration","A02A18_BFLC_PreDrug.csv"), delim=",") %>% select(!Drug)

PreFLC <- read_delim(here("data_in", "Biofilm_Penetration","BFLC_PreDrug_FullDF.csv"), delim=",") %>% rbind(PreA02A18FLC)

prebackground <- mean(PreFLC$Background)
PreFLC <- PreFLC%>%select(!Background)
preFLC <- melt(PreFLC,c("Strain","Biorep"))
preFLC$value = preFLC$value - prebackground
preFLC$value <-ifelse(preFLC$value < 0, 0,preFLC$value)
preFLC <- preFLC %>% rename(enviro = variable, OD = value)
preFLC$enviro = as.numeric(levels(preFLC$enviro))[preFLC$enviro]
preFLC$enviro <- ifelse(preFLC$enviro == "0", 0.001, preFLC$enviro)

preFLC <- preFLC %>% group_by(Strain) %>% summarise(mean = mean(OD))
##################################

PostA02A18FLC <-  read_delim(here("data_in", "Biofilm_Penetration","A02A18_BFLC_PostDrug.csv"), delim=",") %>% select(!Drug)

PostFLC <- read_delim(here("data_in", "Biofilm_Penetration","BFLC_PostDrug_FullDF.csv"), delim=",") %>% rbind(PostA02A18FLC)


postbackground <- mean(PostFLC$Background)
PostFLC <- PostFLC%>%select(!Background)
postFLC <- melt(PostFLC, c("Strain", "Biorep"))
postFLC$value = postFLC$value - postbackground
postFLC$value <- ifelse(postFLC$value < 0, 0, postFLC$value)
postFLC <- postFLC %>% rename(enviro = variable, OD = value)
postFLC$enviro=as.numeric(levels(postFLC$enviro))[postFLC$enviro]

FLCNoBiofilm <- postFLC %>% filter(enviro == 0 & OD < 0.1) %>% select(Strain)

postFLC <- postFLC %>% filter(!Strain %in% FLCNoBiofilm$Strain)

postFLC$enviro <- ifelse(postFLC$enviro == 0,0.5 ,postFLC$enviro)


#aggregating and normalizing the FLC data 
FLC.ag <- postFLC %>% group_by(Strain, enviro) %>% 
  summarise(OD = mean(OD)) 
  



# postFLC.norm <- merge(preFLC, FLC.ag, by = "Strain") %>% 
  # group_by(Strain, enviro) %>% summarise(Mass = OD, OD = mean(OD/mean))

postFLC.norm <- merge(preFLC, FLC.ag, by = "Strain") %>% 
  group_by(Strain, enviro) %>% summarise(Mass = OD, OD = mean(OD-mean))

BREAKFLC <- unique(postFLC$enviro)
LABFLC <- unique(postFLC$enviro)%>%round(digits = 2)
LABFLC[7] <- 0
##########################
FLCA02A18MIC <- read_delim(here("data_in", "Planktonic_MIC50","A02A18_FLC_MIC50.csv"),delim=",") 

FLCMIC <- read_delim(here("data_in", "Planktonic_MIC50","200721_FLC_MIC50_editted.csv"),delim=",") %>% rbind(FLCA02A18MIC) %>% mutate(MIC50 = replace(MIC50, MIC50 >= 4, "≥ 4")) %>% mutate(MIC50= replace(MIC50, MIC50 != "≥ 4", ' ≤ 2')) 

FLC_MIC_biomass <- merge(postFLC.norm, FLCMIC, by= "Strain")
FLC_MIC_biomass$enviro <- round(FLC_MIC_biomass$enviro, 2) 


#adding XTT data 
A02A18XTTBA <- read_delim(here("data_in", "Biofilm_Penetration","A02A18_BBA_XTT.csv"),
                          delim=",") %>% select(!Drug)

XTTBA <- read_delim(here("data_in", "Biofilm_Penetration","BBA_XTT_FullDF.csv"),
                delim=",") %>% rbind(A02A18XTTBA)

xttBAbackground <- mean(XTTBA$Background)
XTTBA <- XTTBA %>% select(!Background)
xttBA <- melt(XTTBA, c("Strain", "Biorep"))
xttBA$value = xttBA$value - xttBAbackground
xttBA$value <- ifelse(xttBA$value < 0, 0.0625, xttBA$value)
xttBA <- xttBA %>% rename(enviro = variable, activity = value)
xttBA$enviro = as.numeric(levels(xttBA$enviro))[xttBA$enviro]
xttBA$enviro <- ifelse(xttBA$enviro == 0, 0.39/2, xttBA$enviro)
xttBA <- xttBA %>% filter(!Strain %in% BANoBiofilm$Strain)
#normalizing the BA XTT data
BA.ND <- xttBA%>%filter(enviro==min(enviro))

normBA <- merge(xttBA,BA.ND, by = c("Strain","Biorep"))%>%
  rename(activity = activity.x, min.env = enviro.y, max.activity = activity.y)%>%
  group_by(Strain, enviro.x)%>%
  summarise(Activity = mean(activity), Activity_norm = mean(activity-max.activity)) %>% 
  rename(Enviro = enviro.x)
normBA$Enviro <- round(normBA$Enviro, 2)  

#FLCXTT
A02A18XTTFLC <- read_delim(here("data_in", "Biofilm_Penetration","A02A18_BFLC_XTT.csv"),
                       delim=",") %>% select(!Drug)

XTTFLC <- read_delim(here("data_in", "Biofilm_Penetration","BFLC_XTT_FullDF.csv"), delim=",") %>% rbind(A02A18XTTFLC)

xttFLCbackground <- mean(XTTFLC$Background)
XTTFLC <- XTTFLC %>% select(!Background)
xttFLC <- melt(XTTFLC, c("Strain", "Biorep"))
xttFLC$value = xttFLC$value - xttFLCbackground
xttFLC$value <- ifelse(xttFLC$value < 0, 0.0625, xttFLC$value)
xttFLC <- xttFLC %>% rename(enviro = variable, activity = value)
xttFLC$enviro <- as.numeric(levels(xttFLC$enviro))[xttFLC$enviro]
xttFLC$enviro <- ifelse(xttFLC$enviro == 0, 0.5, xttFLC$enviro)
xttFLC <- xttFLC %>% filter(!Strain %in% FLCNoBiofilm$Strain)
#Normalizing
FLC.ND <- xttFLC %>% filter(enviro == min(enviro))

normFLC <- merge(xttFLC, FLC.ND, by = c("Strain","Biorep"))%>%
  rename(activity = activity.x, min.env = enviro.y, max.activity = activity.y)%>%
  group_by(Strain, enviro.x)%>%
  summarise(Activity = mean(activity), Activity_norm = mean(activity-max.activity)) %>% 
  rename(Enviro = enviro.x) 
normFLC$Enviro <- round(normFLC$Enviro, 2)  


BA_MIC_activity <- merge(normBA,BAMIC, by="Strain")
FLC_MIC_activity <- merge(normFLC,FLCMIC, by="Strain")


# Combine the dataframes ----

BA_MIC_activity_biomass <- merge(BA_MIC_biomass, BA_MIC_activity, by=c("Strain", "Enviro", "MIC50"))

FLC_MIC_biomass <- FLC_MIC_biomass %>% rename(Enviro = enviro)
FLC_MIC_activity_biomass <- merge(FLC_MIC_biomass, FLC_MIC_activity, by=c("Strain", "Enviro", "MIC50"))



# Plot data ----
#Plotting normalized BA data 
BA_MIC_biomass$MIC50 <- as.factor(BA_MIC_biomass$MIC50)

BAn <- BA_MIC_biomass %>%
  ggplot(mapping= aes(x=Enviro, y=OD, color=Strain)) +
  geom_point(alpha=0.5) +
  geom_line(aes(group = Strain) , alpha = 0.5,
            size=1.5) +
  geom_smooth(aes(x = Enviro, y = OD), size=1.5, se = FALSE, color = "black") +
  geom_hline(mapping= aes(yintercept =0), linetype= "dashed", color="Red")+
  scale_x_continuous( trans = 'log2', #trans=reverselog_trans(2),
                      breaks=BREAKBA,
                      labels = LABBA )+ ylim(-0.4, 0.4)+
  xlab("BA concentration (mg/ml)")+
  ylab("Biomass OD600")+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   text=element_text(family="Times New Roman", face="bold", size = 12), axis.text.x = element_text( vjust = 1.5, size =12))+ scale_color_viridis(discrete = TRUE)


FLCn <-  FLC_MIC_biomass %>%
  ggplot(mapping= aes(x=Enviro, y=OD, color=Strain)) +
  geom_point(alpha=0.5) +
  geom_line(aes(group = Strain) , alpha = 0.5,
            size=1.5) +
  geom_smooth(aes(x = Enviro, y = OD), size=1.5, se = FALSE, color = "black") +
  geom_hline(mapping= aes(yintercept =0), linetype= "dashed", color="Red")+
  scale_x_continuous( trans = 'log2', #trans=reverselog_trans(2),
                      breaks=BREAKBA,
                      labels = LABBA )+ ylim(-0.4, 0.4)+
  xlab("FLC concentration (μg/ml)")+
  ylab("Biomass OD600")+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   text=element_text(family="Times New Roman", face="bold", size = 12), axis.text.x = element_text( vjust = 1.5, size =12))+ scale_color_viridis(discrete = TRUE)
 



BA_act <- ggplot(BA_MIC_activity,
       mapping= aes(x = Enviro, y = Activity_norm,
                    color = Strain)) +
  geom_point(alpha = 0.5) + 
  geom_line(aes(group = Strain), alpha = 0.5, size = 1.5) +
  geom_smooth(aes(x = Enviro, y = Activity_norm), size=1.5, se = FALSE, color = "black") +
  scale_x_continuous( trans = 'log2', #trans=reverselog_trans(2),
                      breaks=BREAKBA,
                      labels = LABBA )+ ylim(-1.2, 1.2) +
  
  xlab("BA concentration (mg/ml)")+
  ylab("Biofilm Activity") +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   text=element_text(family="Times New Roman", face="bold", size = 12),
                   axis.text.x = element_text( vjust = 1.5, size =12)) + 
  scale_color_viridis(discrete = TRUE)


FLC_act <- ggplot(FLC_MIC_activity,
       mapping= aes(x=Enviro, y=Activity_norm, color = Strain)) +
  geom_point(alpha = 0.5) + 
  geom_line(aes (group = Strain), alpha = 0.5, size = 1.5) +
  geom_smooth(aes(x = Enviro, y = Activity_norm), size=1.5, se = FALSE, color = "black") +
  scale_x_continuous( trans = 'log2', #trans=reverselog_trans(2),
                      breaks=BREAKFLC,
                      labels = LABFLC )+ ylim(-1.2, 1.2) +
  xlab("FLC concentration (μg/ml)")+
  ylab("Biofilm Activity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Times New Roman", face="bold", size = 12), legend.position="right",axis.text.x = element_text( vjust = 1.5, size =12)) + scale_color_viridis(discrete = TRUE)


FLCMassAct <- ggarrange(FLCn, FLC_act, common.legend = TRUE, ncol = 1, legend = "none")
BAMassAct <- ggarrange(BAn, BA_act, common.legend = TRUE, ncol = 1, legend = "none")

FLCBAMassAct <- ggarrange(FLCMassAct,BAMassAct, align = "hv")

ggsave(here("figures_out","Biofilm_Penetration", "210831_MassAct.jpg"), 
             FLCBAMassAct,
                          height = 6, 
                          width = 6)


# Adding in the Strain Info and Exporting the data as .CSV
Info <- read_csv(here("data_in", "Biofilm_Penetration", "200302_BBA_StrainNames.csv"))
BA.ag <- BA.ag %>% rename(Enviro = enviro) %>% mutate(Enviro = round(Enviro, digits = 2)) 

BA_MIC_MassAct <- BA_MIC_activity_biomass %>%  rename(Mass_norm = OD) %>% merge(BA.ag, by = c("Strain", "Enviro")) %>% rename(Mass = OD) %>% merge(Info, by = "Strain") %>% merge(preBA, by = "Strain") %>% rename(PreMass = mean) %>% 
  select(!Strain) %>% 
  rename(Strain = StrainAGNo) %>% relocate(Strain, Enviro, MIC50, Mass, Mass_norm) 

FLC_MIC_MassAct <- merge(FLC_MIC_activity_biomass, Info, by = "Strain") %>% select(!Strain) %>% 
  rename(Strain = StrainAGNo, Mass_norm = OD) %>% relocate(Strain) 



BA_MIC_MassAct$Enviro <- ifelse(BA_MIC_MassAct$Enviro == "0.2","0", BA_MIC_MassAct$Enviro)


FLC_MIC_MassAct$Enviro <- ifelse(FLC_MIC_MassAct$Enviro == "0.5", "0", FLC_MIC_MassAct$Enviro)


write_csv(BA_MIC_MassAct, here("data_out", "Biofilm_Penetration", "210705BA_MIC_MassAct.csv"))
write_csv(FLC_MIC_MassAct, here("data_out", "Biofilm_Penetration", "210705FLC_MIC_MassAct.csv"))

