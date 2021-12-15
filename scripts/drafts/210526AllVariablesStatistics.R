#install.packages("corrr") 
library(corrr)
library(tidyverse)
library(here)
library(ggpubr)
library(extrafont)
library(broom)

# Files ----
d <- read_csv(here("data_out", "Combined", "210303MassActDDA.csv"))
d_BA<-d %>% filter(Drug=="BA") %>% filter(!Enviro==0.20 & FoG50 < 0.25)
d_FLC<-d %>% filter(Drug=="FLC")





require(plyr)
func <- function(d_BA)
{
  return(data.frame(tidy(cor.test(d_BA$FoG50, d_BA$Activity, method = "spearman"))))
}

d_BA_cor <- ddply(d_BA, .(Enviro), func)

require(plyr)
func <- function(d_FLC)
{
  return(data.frame(tidy(cor.test(d_FLC$FoG50, d_FLC$Activity, method = "spearman"))))
}

d_FLC_cor <- ddply(d_FLC, .(Enviro), func)

names(d)[5:6] <- c("Biomass", "Activity")




BA_TolAct_p <- merge(d_BA, d_BA_cor) %>% filter(FoG50 < 0.25) %>% 
  ggplot(aes(x = FoG50, y = Activity)) +
  geom_point() +  geom_smooth(aes(linetype = (p.value>0.05)),
                              method = "lm", se = FALSE, color = "black", size = 0.5) +
  stat_cor( method = "spearman",
           label.x = 0.15,
           label.y = 1.2,
           digits = 2,
           label.sep = ",", 
           size = 2) +
  xlab("Planktonic Tolerance") +
  ylab("Biofilm Activity") +
  xlim(0, 0.3) +
  ylim(0.45, 1.2) +
  facet_wrap(~Enviro) +
  theme_bw() +
  theme(legend.position = "None",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                            face = "bold", size = 10),
        axis.text.x = element_text(hjust = 0.5))

ggsave("210526_BATolActCor.jpg",BA_TolAct_p, 
       height = 3.5, 
       width = 5)

BA_ResAct_p <- merge(d_BA, d_BA_cor) %>% 
  ggplot(aes(x = RAD50, y = Activity)) +
  geom_point() +
  stat_cor(method = "spearman",
           label.x = 0.15,
           label.y = 1.2,
           digits = 2,
           label.sep = ",", 
           size = 2) +
  geom_smooth(linetype = "dashed",
              method = "lm", se = FALSE, color = "black", size = 0.5) +
  xlab("Planktonic Resistance") +
  ylab("Biofilm Activity") +
  ylim(0.45, 1.2) +
  facet_wrap(~Enviro) +
  theme_bw() +
  theme(legend.position = "None",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                            face = "bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))


ggsave("210526_BAResActCor.jpg",BA_ResAct_p, 
       height = 3.5, 
       width = 5)

FLC_TolAct_p <- merge(d_FLC, d_FLC_cor) %>% filter(Enviro > 0.5) %>% 
  ggplot(aes(x = FoG50, y = Activity)) +
  geom_point() +
  stat_cor(method = "spearman",
           label.x = 0.15,
           label.y = 1.25,
           digits = 2,
           label.sep = ",", 
           size = 2) +
  geom_smooth(linetype = "dashed",
              method = "lm", se = FALSE, color = "black", size = 0.5) +
  xlab("Planktonic Tolerance") +
  ylab("Biofilm Activity") +
  facet_wrap(~Enviro) +
  theme_bw() +
  theme(legend.position = "None",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                            face = "bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))

ggsave("210526_FLCTolActCor.jpg",FLC_TolAct_p, 
       height = 3.5, 
       width = 5)


FLC_ResAct_p <- merge(d_FLC, d_FLC_cor) %>% filter(Enviro > 0.5) %>% 
  ggplot(aes(x = RAD50, y = Activity)) +
  geom_point() +
  stat_cor(method = "spearman",
           label.x = 0.15,
           label.y = 1.25,
           digits = 2,
           label.sep = ",", 
           size = 2) +
  geom_smooth(linetype = "dashed",
              method = "lm", se = FALSE, color = "black", size = 0.5) +
  xlab("Planktonic Resistance") +
  ylab("Biofilm Activity") +
  facet_wrap(~Enviro) +
  theme_bw() +
  theme(legend.position = "None",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                            face = "bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))

ggsave("210519_FLCResActCor.jpg",FLC_ResAct_p, 
       height = 3.5, 
       width = 5)
