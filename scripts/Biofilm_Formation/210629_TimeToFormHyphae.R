library(here)
library(tidyverse)
library(reshape2)
library(viridis)
library(ggsci)
library(extrafont)
library(ggpubr)

Strain_Names <- tibble( Strain = c("A02","A03","A04","A08","A10","A12","A17","A18"), StrainName = c("P87", "GC75","P78048", "P75016", "P76055", "T101","SC5314", "FH1"))

df <- read_delim(here("data_in","Biofilm_Formation","210608_TimeToHyphalFormatiom.csv"), delim=",") %>% group_by(Strain, Drug, Concentration) %>% 
  pivot_longer(cols= starts_with("T_"), names_to = "BioRep", values_to = "Time") %>% 
  summarise(Time_sd= sd(Time, na.rm = TRUE), Time = mean(Time)) %>% 
  merge(Strain_Names, by = "Strain")


df$Concentration<-as.factor(df$Concentration)




#BA
df%>%filter(Drug=="BA"|Drug=="ND")%>% ggplot(aes(color=StrainName)) +

  geom_point(mapping= aes(x=Concentration, y=Time),size=3,
             position=position_dodge(width=0.5), alpha=0.8) +
  geom_errorbar(aes(ymin=Time-Time_sd, ymax = Time+Time_sd, x= Concentration, group = StrainName), width = 0, position=position_dodge(width=0.5))+
  annotate("text",x = 1, y = 25, fontface =2,  color = "red", size = 4,  label = "No hyphae") +
  scale_y_continuous(breaks = seq(0,24,4))+ 
  geom_hline(yintercept = 24, color="red", linetype="dashed")+
  coord_cartesian(clip = "off")+ expand_limits(y=c(0,28))+
  xlab("BA Concentration (mg/mL)")+  theme_bw()+ 
  theme(panel.grid.minor=element_blank())+
  ylab("Time to First Hyphae (hours)")+
  labs(color = "Strain") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", face="bold")) +
  scale_colour_brewer(palette = "Set1")

BA <- df%>%filter(Drug=="BA"|Drug=="ND")%>% ggplot(aes(color=StrainName)) +
  geom_point(mapping= aes(x=Concentration, y=Time),size=3,
             position=position_dodge(width=0.5), alpha=0.8) +
  geom_errorbar(aes(ymin=Time-Time_sd, ymax = Time+Time_sd, x= Concentration, group = StrainName), width = 0, position=position_dodge(width=0.5))+
  annotate("text",x = 1.5, y = 25.5, fontface =2,  color = "black", size = 4,  label = "No hyphae") +
  scale_y_continuous(breaks = seq(0,24,4))+ 
  geom_hline(yintercept = 24, color="grey", linetype="dashed")+
  coord_cartesian(clip = "off")+ expand_limits(y=c(0,28))+
  xlab("BA Concentration (mg/mL)")+  theme_bw()+ 
  theme(panel.grid.minor=element_blank())+
  ylab("Time to First Hyphae (h)")+
  labs(color = "Strain") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", face="bold", size = 16), 
        axis.text = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  scale_colour_brewer(palette = "Set1")


FLC <- df%>%filter(Drug=="FLC"|Drug=="ND")%>%ggplot(aes(color=StrainName)) +
  geom_point(mapping= aes(x=Concentration, y=Time),size=3,
             position=position_dodge(width=0.5), alpha=0.8) +
  geom_errorbar(aes(ymin=Time-Time_sd, ymax = Time + Time_sd, x= Concentration, group = StrainName), width = 0, position=position_dodge(width=0.5))+
  geom_hline(yintercept = 24, color="grey", linetype="dashed")+
  annotate("text",x = 1.5, y = 25.5, fontface =2,  color = "black", size = 4,  label = "No hyphae") +
  scale_y_continuous(breaks = seq(0,24,4))+
  expand_limits(y=c(0,28))+
  coord_cartesian(clip = "off")+ 
  xlab("FLC Concentration (μg/mL)")+  theme_bw() + 
  ylab("Time to First Hyphae (h)")+
  labs(color = "Strain") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", face="bold", size = 16),
        axis.text = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))  +
  scale_colour_brewer(palette = "Set1")


BAFLC <- ggarrange(FLC, BA, common.legend = TRUE, legend = "right")
BAFLC
ggsave(here("figures_out","Biofilm_Formation","210821_Time_to_Hyphae.jpg"), BAFLC,
       width = 8, height = 3.5)
