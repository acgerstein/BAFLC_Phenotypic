library(here)
library(tidyverse)
library(reshape2)
library(viridis)
library(reshape2)
library(ggforce)
library(scales)
library(extrafont)


BA_MassAct <- read_csv(here("data_out",
                            "Biofilm_Penetration",
                            "210303BA_MIC_biomass_activity.csv"))

FLC_MassAct <- read_csv(here("data_out",
                            "Biofilm_Penetration",
                            "210303FLC_MIC_biomass_activity.csv")) %>% 
  filter(Strain%in% BA_MassAct$Strain)

BA_MassAct <- BA_MassAct %>% filter(Strain%in% FLC_MassAct$Strain)
BA_MassAct$MIC50 <- as.factor(BA_MassAct$MIC50)

FLCBreaks <- unique(FLC_MassAct$Enviro)
BABreaks <- c(0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8)

# Plot data ----
#Plotting normalized BA data 
BAn <- BA_MassAct %>% 
  ggplot(mapping= aes(x = Enviro, y = OD, color = MIC50, group = Strain )) +
  
  geom_line(stat="smooth", method="loess", 
            mapping= aes(x=Enviro, y=OD,color=MIC50),se=FALSE,alpha=0.5,
            size=1.5, span=0.8)+
  geom_hline(mapping= aes(yintercept =1), linetype= "dashed", color="Red")+
  scale_x_continuous( trans = 'log2', 
                      breaks = BABreaks) +
  xlab("BA concentration (mg/ml)")+
  ylab("Biomass OD600")+expand_limits(x = 0.5, y = c(0,2.1))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   text=element_text(family="Times New Roman", face="bold", size = 12),
                   legend.position="right",axis.text.x = element_text( vjust = 1.5, size =12))+ 
scale_color_viridis(discrete = TRUE, label = c("0.40", "0.80"))


FLCn <- ggplot(FLC_MassAct, mapping = aes(x=Enviro, y=OD, group=Strain )) +
  geom_line(stat="smooth", method="loess", 
            mapping= aes(x=Enviro, y=OD,color=MIC50),se=FALSE,alpha=0.5,
            size=1.5, span=0.8)+
  geom_hline(mapping= aes(yintercept =1), linetype= "dashed", color="Red")+
  scale_x_continuous( trans = 'log2', #trans=reverselog_trans(2),
                      breaks=FLCBreaks)+
  xlab("FLC concentration (ug/ml)")+
  ylab("Biomass OD600")+expand_limits(x = 0, y = c(0,2.1))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   text=element_text(family="Times New Roman", face="bold", size = 12),
                   legend.position="right",axis.text.x = element_text( vjust = 1.5, size =12))+ 
  scale_color_manual(values = c("#FF0000","#9900CC", "#6666FF","#66FFFF"  ))
#saving the FLC normalized figure




#ploting all the strains and making the lines more transparent

NormBA <- ggplot(BA_MassAct, mapping= aes(x = Enviro, y = Activity_norm )) +
  geom_line(stat="smooth", method="loess", 
            mapping= aes(color=MIC50,group=Strain), se=FALSE, alpha=0.5,
            size=1.5, span=0.8)+
  scale_x_continuous( trans = 'log2', #trans=reverselog_trans(2),
                      breaks= BABreaks)+
  xlab("BA concentration (mg/ml)")+
  ylab("Biofilm Activity")+expand_limits(x = 0, y=c(0.5,1.3))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   text=element_text(family="Times New Roman", face="bold", size = 12),
                   legend.position="right",axis.text.x = element_text( vjust = 1.5, size =12))+ 
  scale_color_viridis(discrete = TRUE, label = c("0.40", "0.80") )



NormFLC <- ggplot(FLC_MassAct,mapping= aes(x=enviro.x, y=norm,group=Strain )) +
  geom_line(stat="smooth", method="loess", 
            mapping= aes(x=Enviro, y=Activity_norm,color=MIC50),se=FALSE,alpha=0.5,
            size=1.5, span=0.8)+
  scale_x_continuous( trans = 'log2', #trans=reverselog_trans(2),
                      breaks= FLCBreaks)+
  xlab("FLC concentration (ug/ml)")+
  ylab("Biofilm Activity")+expand_limits(x = 0, y=c(0.5,1.3))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   text=element_text(family="Times New Roman", face="bold", size = 12),
                   legend.position="right",axis.text.x = element_text( vjust = 1.5, size =12))+  scale_color_manual(values = c("#FF0000","#9900CC", "#6666FF","#66FFFF"  ))
  


BA <- ggarrange(BAn, NormBA, common.legend = TRUE,  nrow = 2)
FLC <- ggarrange( FLCn,NormFLC,  common.legend = TRUE,  nrow = 2)

FLCBA_MassAct <- ggarrange(FLC, BA)

ggsave(here("figures_out","Biofilm_Penetration","210704_FLCBA_MassAct.jpg"), plot = FLCBA_MassAct, height = 6, width = 6)

tidy(aov(Activity_norm~MIC50*Enviro, data = BA_MassAct))
tidy(aov(OD~MIC50*Enviro, data = BA_MassAct))

tidy(aov(Activity_norm~MIC50*Enviro, data = FLC_MassAct))
tidy(aov(OD~MIC50*Enviro, data = FLC_MassAct))
