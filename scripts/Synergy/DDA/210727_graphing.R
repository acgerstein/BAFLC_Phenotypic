library(here)
library(tidyverse)
library(ggpubr)
library(ggforce)
library(extrafont)
library(gt)

Info <-  read_csv(here("FLCBA_Synergy_Checkerboard" ,"DDA", "HSC_Info.csv"))


df_alb <- read_csv(here("FLCBA_Synergy_Checkerboard" ,"DDA",
                        "SynergyDDA_Albicans_df.csv"))%>% 
  select(line, type, RAD50, FoG50) %>% 
  rename("Strain" = "line",
         "Drug" = "type") %>% 
  filter(!Strain == "Blank") %>% 
  group_by(Strain, Drug) %>% 
  summarise(RAD50 = mean(RAD50), 
            FoG50 = mean(FoG50)) %>% 
  mutate(Species = "C.albicans", Year = 2014)

df_nonalb <- read_csv(here("FLCBA_Synergy_Checkerboard", "DDA",
              "Non_alb_FLCBA_df.csv")) %>% 
  select(line, type, RAD50, FoG50) %>% 
  rename("Strain" = "line",
         "Drug" = "type") %>% 
  filter(!Strain == "Blank") %>% 
  group_by(Strain, Drug) %>% 
  summarise(RAD50 = mean(RAD50), 
            FoG50 = mean(FoG50)) %>% 
  merge(Info, by = "Strain")

df <- rbind(df_alb, df_nonalb)
write_csv(df, here("FLCBA_Synergy_Checkerboard" ,"DDA",
                   "210804_SynergyDDA_df.csv"))
df %>% group_by(Species) %>% summarise(Strain = unique(Strain)) %>% summarise(n = length(Strain))


#transform(df, Drug=factor(Drug,levels=c("FLC","BA","FLC+BA")))
ResistanceProfile <- ggplot(df) + 
  geom_histogram(aes(x = RAD50), binwidth= 1, fill = "purple", alpha = 0.6) + 
  theme_bw() + 
  scale_x_reverse() +
  theme(legend.position = "none",                                                                     
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                            face = "bold", size = 12), axis.text.x = element_text(hjust = 0.5)) +
  facet_wrap(~Drug, nrow = 3)

ToleranceProfile <- df %>%  ggplot() + 
  geom_histogram(aes(x = FoG50), binwidth= 0.05, fill = "purple", alpha = 0.6) + 
  theme_bw() + 
  theme(legend.position = "none",                                                                     
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                            face = "bold", size = 12), axis.text.x = element_text(hjust = 0.5)) +
  facet_wrap(~Drug, nrow = 3)

ResistanceToleranceProfile <- ggarrange(ResistanceProfile, ToleranceProfile)

ggsave(here("FLCBA_Synergy_Checkerboard",
            "DDA",
            "figures_out", 
            "210811_ResistanceToleranceProfile.jpg"), 
       plot = ResistanceToleranceProfile,
       width = 7, 
       height = 4)

df_Diff <- df %>% group_by(Strain) %>% 
  summarise(RAD_BA_Diff = RAD50[Drug == "FLC+BA"] - RAD50[Drug == "BA"],
            RAD_FLC_Diff = RAD50[Drug == "FLC+BA"] - RAD50[Drug == "FLC"],
            FoG_BA_Diff = FoG50[Drug == "FLC+BA"] - FoG50[Drug == "BA"],
            FoG_FLC_Diff = FoG50[Drug == "FLC+BA"] - FoG50[Drug == "FLC"]) 

df_Diff %>% 
  ggplot() +
  geom_histogram(aes(x = RAD_BA_Diff),bins = 10, fill = "purple", alpha = 0.6)

df_Diff %>% 
  ggplot() +
  geom_histogram(aes(x = RAD_FLC_Diff),bins = 30, fill = "purple", alpha = 0.6) + xlim(-5, 30)



df_Diff %>% 
  ggplot() +
  geom_histogram(aes(x = FoG_BA_Diff),bins = 30, fill = "purple", alpha = 0.6)

df_Diff %>% 
  ggplot() +
  geom_histogram(aes(x = FoG_FLC_Diff),bins = 30, fill = "purple", alpha = 0.6)

#write_csv(df, here("DDA","Pilot_DDA","210531_SynergyDDA_AllSpecies_df.csv"))


df_summary <- df %>% 
  group_by(Strain) %>% 
  mutate(Tol_Synergy = ifelse(Drug == "FLC+BA" &
                                RAD50[Drug == "FLC+BA"] >= RAD50[Drug == "FLC"] &
                                RAD50[Drug == "FLC+BA"] >= RAD50[Drug == "BA"] &
                                 0.95*FoG50[Drug == "FLC+BA"] <= FoG50[Drug == "FLC"] & 
                                0.95*FoG50[Drug == "FLC+BA"] <= FoG50[Drug == "BA"],
                              "TRUE", "FALSE" ), 
         Tol_Ant = ifelse(Drug == "FLC+BA" &
                                RAD50[Drug == "FLC+BA"] >= RAD50[Drug == "FLC"] &
                                RAD50[Drug == "FLC+BA"] >= RAD50[Drug == "BA"] &
                                0.95*FoG50[Drug == "FLC+BA"] >= FoG50[Drug == "FLC"] & 
                                0.95*FoG50[Drug == "FLC+BA"] >= FoG50[Drug == "BA"],
                              "TRUE", "FALSE" ), 
         Res_ant = ifelse(Drug == "FLC+BA" &
                                RAD50[Drug == "FLC"] - RAD50[Drug == "FLC+BA"] >= 2 & 
                                RAD50[Drug == "BA"] - RAD50[Drug == "FLC+BA"] >= 2, "TRUE", "FALSE"), 
         Res_Synergy = ifelse(Drug == "FLC+BA" &
                                RAD50[Drug == "FLC+BA"] > RAD50[Drug == "FLC"]&
                                RAD50[Drug == "FLC+BA"] > RAD50[Drug == "BA"], "TRUE", "FALSE"), 
         
         Res_more_two = ifelse(Drug == "FLC+BA" &
                                 RAD50[Drug == "FLC+BA"] - RAD50[Drug == "FLC"] >=2 &
                                 RAD50[Drug == "FLC+BA"] - RAD50[Drug == "BA"] >=2, "TRUE", "FALSE"))
  
df_summary %>% group_by(Strain) %>% 
  filter( Res_ant == "TRUE")

df_tb <- df_summary %>% group_by(Species) %>% 
  mutate( n = length(unique(Strain))) %>% 
  summarise(n = mean(n), 
            ToleranceSynergy = length(which(Tol_Synergy == "TRUE" )), 
            ResAnt = length(which(Res_ant == "TRUE")),
            ResSynergy = length(which(Res_Synergy == "TRUE")),
            ResTolSynergy = length(which(Res_more_two == "TRUE" & Tol_Synergy == "TRUE" )),
            MoreThan2 = length(which(Res_more_two == "TRUE"))) %>%
  mutate(TolSynergyPercent = ToleranceSynergy/n*100, ResSynergyPercent = MoreThan2/n*100)


DecR <- df_summary %>% group_by(Strain ,Species) %>% 
  filter(Res_more_two == "TRUE") 

DecT <- df_summary %>% group_by(Strain ,Species) %>% 
  filter(Tol_Synergy == "TRUE") 

S_list <- rbind(DecR, DecT) %>% select(!Res_ant)

write_csv(df_tb, file =  here("DDA",
                              "Pilot_DDA","210519_Synergy_dataSum.csv"))
df_stat %>%  df %>% 
  group_by(Species, Drug) %>%
  summarise(RAD50_med = mean(RAD50), 
            FoG50_med = mean(FoG50)) 

p2 <- df %>% ggplot() +
  geom_sina(mapping = aes(x = Drug, y = FoG50, color = Tol_Synergy, group = Drug), 
            alpha = 0.6, size=2) + 
  geom_point(df_stat, mapping = aes(x = Drug, y = FoG50_med),
              color = "red",pch = "_",size = 10, alpha = 0.8) +
  ylab("Tolerance") + 
  expand_limits(x = 0, y = c(0,1)) +
  theme_bw() + 
  theme(legend.position = "none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                          face = "bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))




p1<- df %>% ggplot() +
  geom_sina(mapping = aes(x = Drug, y = RAD50, color = Res_Synergy, group = Drug), 
            alpha = 0.6, size=2) + 
  geom_point(df_stat, mapping = aes(x = Drug, y = RAD50_med),
             color = "red",pch = "_",size = 10, alpha = 0.8) +
  scale_y_reverse() +
  ylab("Resistance") + 
  expand_limits(x = 0, y = c(25,0)) +
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                            face = "bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))

RT<-ggarrange(p1,p2)
ggsave("210511_TR_Syn.jpg",
       RT, 
       path = here("DDA",
            "Pilot_DDA",
            "210507_SynergyFLC+BA",
            "48h"), height = 3, width = 5)

p3 <- df %>% ggplot() +
  geom_jitter(mapping = aes(x = Strain, y = RAD50, color = Drug, shape = Drug), 
             alpha = 0.6, size=2, width = 0.05, height = 0.01) + 
  scale_x_discrete(breaks = unique(df$Strain), labels = seq(1, 78)) +
  scale_y_reverse() +
  ylab("Resistance") + 
  expand_limits(x = 0, y = c(30,0)) +
  theme_bw() + 
  theme(legend.position = "top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                            face = "bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("blue", "red", "purple4" ))


ggsave("210512_Syn_Resistance.jpg", p3, 
       path = here("DDA",
                   "Pilot_DDA",
                   "210507_SynergyFLC+BA",
                   "48h"), 
       width = 7, 
       height = 4)
p4 <- df %>% ggplot() +
  geom_point(mapping = aes(x = Strain, y = FoG50, color = Drug, shape = Drug), 
             alpha = 0.6, size=2) + 
  scale_x_discrete(breaks = unique(df$Strain), labels = seq(1, 78)) +
  ylab("Tolerance") + 
  expand_limits(x = 0, y = c(0,1)) +
  theme_bw() + 
  theme(legend.position = "top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                            face = "bold", size = 8),
        axis.text.x = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("blue", "red", "purple4" ))

ggsave("210512_Syn_Tolerance.jpg", p4, 
       path = here("DDA",
                   "Pilot_DDA",
                   "210507_SynergyFLC+BA",
                   "48h"), 
       width = 7, 
       height = 4)
