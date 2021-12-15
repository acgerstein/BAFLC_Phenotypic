library(here)
library(tidyverse)
library(cowplot)
library(ggforce)
library(gridExtra)
library(ggpubr)


df <- read_csv(here("data_out", "DDA","210629_HSCAG.csv")) 
df_stat <- df %>% group_by(Species, Drug) %>% 
  summarise(RAD50Med = median(RAD50, na.rm = TRUE), 
            FoG50Med = median(FoG50, na.rm = TRUE))
df_stat_BA <- df_stat %>% filter(Drug == "BA")
df_stat_FLC <- df_stat %>% filter(Drug == "FLC")
#Plotting the full data
R_BA <- df %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( df_stat_BA, mapping = aes(x=Species, y= RAD50Med),
              color="red",pch="_",size=10, alpha=0.8, 
              position= position_dodge(width = 0.5))+
  scale_x_discrete(label = c("C. albicans", "C. glabrata", "C. parapsilosis", "C. tropicalis")) +
  ylab("RAD50")+expand_limits(x = 0, y = c(0,28))+
  theme_bw()+ labs(subtitle = "BA")+
  scale_y_reverse()+
  scale_color_manual(values = c("#333333", "#660066"))+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 12)) 

#FLC Resistance graph
R_FLC <- df %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( df_stat_FLC, mapping = aes(x=Species, y= RAD50Med),
              color="red",pch="_",size=10, alpha=0.8, 
              position= position_dodge(width = 0.5))+
  scale_x_discrete(label = c("C. albicans", "C. glabrata", "C. parapsilosis", "C. tropicalis")) +
  ylab("RAD50")+expand_limits(x = 0, y = c(0,28))+
  theme_bw()+ labs(subtitle = "FLC")+
  scale_y_reverse()+
  scale_color_manual(values = c("#333333", "#660066"))+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

#BA Tolerance graph 
T_BA <- df %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG50), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( df_stat_BA, mapping = aes(x=Species, y= FoG50Med),
              color="red",pch="_",size=10, alpha=0.8, 
              position= position_dodge(width = 0.5))+
  scale_x_discrete(label = c("C. albicans", "C. glabrata", "C. parapsilosis", "C. tropicalis")) +
  ylab("FoG50") + 
  expand_limits(x = 0, y = c(0,1))+
  theme_bw() + 
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 12),
        axis.text.x = element_text(hjust = 0.5, face = "bold.italic"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

#FLC Tolerance graph 
T_FLC <- df %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG50), 
            maxwidth = 0.50,alpha = 0.6, size=2) + 
  geom_point( df_stat_FLC, mapping = aes(x=Species, y= FoG50Med),
              color="red",pch="_",size=10, alpha=0.8, 
              position= position_dodge(width = 0.5))+
  scale_x_discrete(label = c("C. albicans", "C. glabrata", "C. parapsilosis", "C. tropicalis"))+
  ylab("FoG50") + 
  expand_limits(x = 0, y = c(0,1))+
  theme_bw() + 
  theme(legend.position="none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 12),
        axis.text.x = element_text(hjust = 0.5, face = "bold.italic"))


p_all<- plot_grid(R_FLC + theme(axis.title.x = element_blank(),
                        axis.text.x = element_blank()),
          R_BA + theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank(), 
                       axis.text.y = element_blank(),
                       axis.title.y = element_blank()),
          T_FLC +theme(), 
          T_BA + theme(axis.text.y = element_blank(),
            axis.title.y = element_blank()),
          align = "v")


ggsave(plot = p_all, here("figures_out", "DDA",
                                "210629_HSCAG.jpg"))
