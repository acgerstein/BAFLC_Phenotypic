#install.packages("lmerTest") 
#library(corrr)
library(tidyverse)
library(here)
library(ggpubr)
library(extrafont)
library(broom)
library(lmerTest)

# Files ----
DDA <- read_csv(here("data_out", "DDA", "210827_HSCAG.csv"))
DDA_BA <- DDA %>% filter(Drug == "BA")
DDA_FLC <- DDA %>% filter(Drug == "FLC")

MassAct_BA <- read_csv(here("data_out", "Biofilm_Penetration", "210705BA_MIC_MassAct.csv"))

MassAct_BA$Enviro[MassAct_BA$Enviro == "0.39"] <- "0.4"
MassAct_BA$Enviro[MassAct_BA$Enviro == "0.78"] <- "0.8"
MassAct_BA$Enviro[MassAct_BA$Enviro == "1.56"] <- "1.6"
MassAct_BA$Enviro[MassAct_BA$Enviro == "3.12"] <- "3.2"
MassAct_BA$Enviro[MassAct_BA$Enviro == "6.25"] <- "6.4"
MassAct_BA$Enviro[MassAct_BA$Enviro == "12.5"] <- "12.8"
  

MassAct_FLC <- read_csv(here("data_out", "Biofilm_Penetration", "210705FLC_MIC_MassAct.csv"))

DDAMassAct_BA <- merge(DDA_BA, MassAct_BA, by = "Strain")
DDAMassAct_FLC <- merge(DDA_FLC, MassAct_FLC, by = "Strain")

DDAMassAct_BA <- subset(DDAMassAct_BA, Enviro != "0")
DDAMassAct_FLC <- subset(DDAMassAct_FLC, Enviro != "0")

# Are resistance (RAD50, MIC50) and biomass (Mass) correlated?
# There are 6 drug environments and a no drug

BAlmer_mass <- lmer(Mass ~ RAD50 + Enviro + (1|Strain), DDAMassAct_BA, REML = FALSE)
anova(BAlmer_mass, type = 3)

# Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq Mean Sq NumDF DenDF  F value Pr(>F)    
# RAD50  0.0020 0.00202     1    55   0.8253 0.3676    
# Enviro 3.4483 0.57472     6   330 235.1149 <2e-16 ***

BAlmer_mass_FoG <- lmer(Mass ~ FoG50 + Enviro + (1|Strain), DDAMassAct_BA, REML = FALSE)
anova(BAlmer_mass_FoG, type = 3)

# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF DenDF  F value Pr(>F)    
# FoG50  0.00178 0.00178     1    55   0.8429 0.3626    
# Enviro 2.53589 0.50718     5   275 239.9903 <2e-16 ***

BAlmer_activity <- lmer(Activity ~ RAD50 + Enviro + (1|Strain), DDAMassAct_BA, REML = FALSE)
anova(BAlmer_activity, type = 3)

# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF DenDF  F value Pr(>F)    
# RAD50   0.0326 0.03259     1    55   1.3457  0.251    
# Enviro 17.7306 2.95510     6   330 122.0242 <2e-16 ***

BAlmer_activity_FoG <- lmer(Activity ~ FoG50 + Enviro + (1|Strain), DDAMassAct_BA, REML = FALSE)
anova(BAlmer_activity_FoG, type = 3)
# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF DenDF  F value Pr(>F)    
# FoG50   0.0005 0.00048     1    55   0.0185 0.8923    
# Enviro 14.8223 2.96445     5   275 113.2594 <2e-16 ***

FLClmer_mass <- lmer(Mass ~ RAD50 + Enviro + (1|Strain), DDAMassAct_FLC, REML = FALSE)
anova(FLClmer_mass, type = 3)

# Type III Analysis of Variance Table with Satterthwaite's method
#           Sum Sq   Mean Sq NumDF DenDF F value   Pr(>F)   
# RAD50  0.0191682 0.0191682     1    55  8.6393 0.004802 **
# Enviro 0.0061533 0.0061533     1   330  2.7733 0.096795 . 

FLClmer_mass_FoG <- lmer(Mass ~ FoG50 + Enviro + (1|Strain), DDAMassAct_FLC, REML = FALSE)
anova(FLClmer_mass_FoG, type = 3)
# Type III Analysis of Variance Table with Satterthwaite's method
#           Sum Sq   Mean Sq NumDF DenDF F value Pr(>F)
# FoG50  0.0037572 0.0037572     1    55  2.3436 0.1315
# Enviro 0.0032324 0.0032324     1   275  2.0162 0.1568

FLClmer_activity <- lmer(Activity ~ RAD50 + Enviro + (1|Strain), DDAMassAct_FLC, REML = FALSE)
anova(FLClmer_activity, type = 3)
# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
# RAD50  0.04682 0.04682     1    55  7.7675  0.007287 ** 
# Enviro 0.50971 0.50971     1   275 84.5628 < 2.2e-16 ***


FLClmer_activity_FoG <- lmer(Activity ~ FoG50 + Enviro + (1|Strain), DDAMassAct_FLC, REML = FALSE)
anova(FLClmer_activity_FoG, type = 3)
# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
# FoG50  0.04992 0.04992     1    55  8.2815  0.005692 ** 
# Enviro 0.50971 0.50971     1   275 84.5628 < 2.2e-16 ***

plot(DDAMassAct_FLC$RAD50, DDAMassAct_FLC$Activity, col=DDAMassAct_FLC$Enviro)
abline(lm(DDAMassAct_FLC$Activity~DDAMassAct_FLC$RAD50))

plot(DDAMassAct_FLC$FoG50, DDAMassAct_FLC$Activity, col=DDAMassAct_FLC$Enviro)
abline(lm(DDAMassAct_FLC$Activity~DDAMassAct_FLC$FoG50))

plot(DDAMassAct_FLC$RAD50, DDAMassAct_FLC$Mass, col=DDAMassAct_FLC$Enviro)
abline(lm(DDAMassAct_FLC$Mass~DDAMassAct_FLC$RAD50))

SpearmanCor <- function (x, y) {
  cor <- tidy(cor.test(x, y, method = "spearman", na.rm = TRUE))}



# BA Resistance and Biomass Correlation 
DDAMassAct_BA %>% group_by(Enviro) %>% 
  summarise(cor = SpearmanCor(RAD50, Mass))  
# BA Resistance and Activity Correlation 
DDAMassAct_BA %>% group_by(Enviro) %>% 
  summarise(cor = SpearmanCor(RAD50, Activity)) 

# BA Tolerance and Biomass Correlation 
DDAMassAct_BA %>% group_by(Enviro) %>% 
  summarise(cor = SpearmanCor(FoG50, Mass)) 
# BA Tolerance and Activity Correlation 
DDAMassAct_BA %>% group_by(Enviro) %>% 
  summarise(cor = SpearmanCor(FoG50, Activity)) 


# FLC Resistance and Biomass Correlation
DDAMassAct_FLC %>% group_by(Enviro) %>% 
  summarise(cor = SpearmanCor(RAD50, Mass))  
# FLC Resistance and Activity Correlation
DDAMassAct_FLC %>% group_by(Enviro) %>% 
  summarise(cor = SpearmanCor(RAD50, Activity)) 
# FLC Tolerance and Biomass Correlation
DDAMassAct_FLC %>% group_by(Enviro) %>% 
  summarise(cor = SpearmanCor(FoG50, Mass)) 
# FLC Tolerance and Activity Correlation
DDAMassAct_FLC %>% group_by(Enviro) %>% 
  summarise(cor = SpearmanCor(FoG50, Activity)) 


BATolAct <- DDAMassAct_BA %>% 
  ggplot(aes(x = FoG50, y = Activity)) +
  geom_point() +  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5, linetype = "dashed") +
  stat_cor( aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
            method = "spearman",
            label.x = 0.2,
            label.y = 3,
            digits = 2,
            label.sep = ",", 
            size = 2) +
  xlab("BA Planktonic Tolerance") +
  ylab("BA Biofilm Activity") +
  facet_wrap(~Enviro) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white",
                                        color = "black",
                                        size = 1),
        strip.text = element_text(colour = "black"), legend.position = "None",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                            face = "bold", size = 14),
        axis.text.x = element_text(hjust = 0.5))

ggsave(here("figures_out", "AllStats", "210707_BATolAct.jpg"), plot = BATolAct)

FLCTolAct <- DDAMassAct_FLC %>% 
  ggplot(aes(x = FoG50, y = Activity)) +
  geom_point() +  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5, linetype = "dashed") +
  stat_cor( aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),method = "spearman",
            label.x = 0.1,
            label.y = 3.4,
            digits = 2,
            label.sep = ",", 
            size = 2) +
  xlab("FLC Planktonic Tolerance") +
  ylab("FLC Biofilm Activity") +
  facet_wrap(~Enviro) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white",
                                        color = "black",
                                        size = 1),
        strip.text = element_text(colour = "black"),legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                            face = "bold", size = 14),
        axis.text.x = element_text(hjust = 0.5),
        legend.text = element_text(family = "Times New Roman", 
                                   face = "bold",
                                   size = 14))
ggsave(here("figures_out", "AllStats", "210707_FLCTolAct.jpg"), plot = FLCTolAct)

BAResMass <- DDAMassAct_BA %>% 
  ggplot(aes(x = RAD50, y = Mass)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5,
              linetype = "dashed") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), method = "spearman",
            digits = 2,
            label.sep = ",", 
            label.x = 14,
            size = 2) +
  xlab("BA Planktonic Resistance") +
  ylab("BA Biofilm Biomass") +
  facet_wrap(~Enviro) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white",
                                        color = "black",
                                        size = 1),
        strip.text = element_text(colour = "black"), legend.position = "None",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                            face = "bold", size = 14),
        axis.text.x = element_text(hjust = 0.5))

ggsave(here("figures_out", "AllStats", "210707_BAResMass.jpg"), plot = BAResMass)

FLCResMass <- DDAMassAct_FLC %>% 
  ggplot(aes(x = RAD50, y = Mass)) +
  geom_point() +  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5, linetype = "dashed") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
            method = "spearman",
            label.x = 0.1,
            label.y = 3,
            digits = 2,
            label.sep = ",", 
            size = 2) +
  xlab("FLC Planktonic Resistance") +
  ylab("FLC Biofilm Biomass") +
  facet_wrap(~Enviro) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white",
                                        color = "black",
                                        size = 1),
        strip.text = element_text(colour = "black") ,legend.position = "None",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                            face = "bold", size = 14),
        axis.text.x = element_text(hjust = 0.5))


ggsave(here("figures_out", "AllStats", "210707_FLCResMass.jpg"), plot = FLCResMass)


