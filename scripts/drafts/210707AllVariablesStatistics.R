#install.packages("corrr") 
library(corrr)
library(tidyverse)
library(here)
library(ggpubr)
library(extrafont)
library(broom)

# Files ----
DDA <- read_csv(here("data_out", "DDA", "210827_HSCAG.csv")) %>% select(!c(RAD50, FoG50))
DDA_BA <- DDA %>% filter(Drug == "BA")
DDA_FLC <- DDA %>% filter(Drug == "FLC")

MassAct_BA <- read_csv(here("data_out", "Biofilm_Penetration", "210705BA_MIC_MassAct.csv"))

MassAct_BA$Enviro[MassAct_BA$Enviro == "0.39"] <- "0.4"
MassAct_BA$Enviro[MassAct_BA$Enviro == "0.78"] <- "0.8"
MassAct_BA$Enviro[MassAct_BA$Enviro == "1.56"] <- "1.6"
MassAct_BA$Enviro[MassAct_BA$Enviro == "3.12"] <- "3.2"
MassAct_BA$Enviro[MassAct_BA$Enviro == "6.25"] <- "6.4"
MassAct_BA$Enviro[MassAct_BA$Enviro == "12.5"] <- "12.8"

# What is the percent of isolates that had increased activity at 0.4 mg/mL

MassAct_BA %>% filter(Enviro %in% c("0", "0.4")) %>% group_by(Strain) %>% filter(Activity[Enviro == "0.4"]> Activity[Enviro == "0"]) 

MassAct_BA %>% filter(Enviro %in% c("0", "0.4")) %>% group_by(Strain) %>% filter(Mass[Enviro == "0.4"]> Mass[Enviro == "0"]) 
  

MassAct_FLC <- read_csv(here("data_out", "Biofilm_Penetration", "210705FLC_MIC_MassAct.csv"))

DDAMassAct_BA <- merge(DDA_BA, MassAct_BA, by = "Strain")
DDAMassAct_FLC <- merge(DDA_FLC, MassAct_FLC, by = "Strain")


# Are resistance (RAD20) and biomass (Mass) correlated?
# There are 6 drug environments and a no drug


FLClmer_mass <- lmer(Mass ~ RAD20 + Enviro + (1|Strain), DDAMassAct_FLC, REML = FALSE)
anova(FLClmer_mass, type = 3)

# Type III Analysis of Variance Table with Satterthwaite's method
#           Sum Sq   Mean Sq NumDF DenDF F value  Pr(>F)  
# RAD20  0.0035209 0.0035209     1    57  1.5922 0.21216  
# Enviro 0.0065919 0.0065919     1   342  2.9809 0.08515 .

FLClmer_activity <- lmer(Activity ~ RAD20 + Enviro + (1|Strain), DDAMassAct_FLC, REML = FALSE)
anova(FLClmer_activity, type = 3)
# Type III Analysis of Variance Table with Satterthwaite's method
#          Sum Sq  Mean Sq NumDF DenDF F value    Pr(>F)    
# RAD20  0.021143 0.021143     1    57  2.5104    0.1186    
# Enviro 0.273328 0.273328     1   342 32.4546 2.633e-08 ***

FLClmer_mass_FoG <- lmer(Mass ~ FoG20 + Enviro + (1|Strain), DDAMassAct_FLC, REML = FALSE)
anova(FLClmer_mass_FoG, type = 3)
# Type III Analysis of Variance Table with Satterthwaite's method
#           Sum Sq   Mean Sq NumDF DenDF F value Pr(>F)
# FoG20  0.0028951 0.0028951     1    54  1.2774 0.2634
# Enviro 0.0042536 0.0042536     1   324  1.8768 0.1716

FLClmer_Act_FoG <- lmer(Activity ~ FoG20 + Enviro + (1|Strain), DDAMassAct_FLC, REML = FALSE)
anova(FLClmer_Act_FoG, type = 3)
# Type III Analysis of Variance Table with Satterthwaite's method
#          Sum Sq  Mean Sq NumDF DenDF F value    Pr(>F)    
# FoG20  0.050703 0.050703     1    54  5.7646   0.01983 *  
# Enviro 0.264077 0.264077     1   324 30.0238 8.586e-08 ***

FLClmer_Act_FoG <- lmer(Activity ~ RAD20 + FoG20 + (1|Enviro) + (1|Strain), DDAMassAct_FLC, REML = FALSE)
anova(FLClmer_Act_FoG, type = 3)
# Type III Analysis of Variance Table with Satterthwaite's method
#          Sum Sq  Mean Sq NumDF DenDF F value    Pr(>F)    
# RAD20  0.017410 0.017410     1    54  1.9794   0.16518    
# FoG20  0.061731 0.061731     1    54  7.0184   0.01056 *  
# Enviro 0.264077 0.264077     1   324 30.0238 8.586e-08 ***
# ---


BAlmer_mass <- lmer(Mass ~ RAD20 + Enviro + (1|Strain), DDAMassAct_BA, REML = FALSE)
anova(BAlmer_mass, type = 3)

# Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq Mean Sq NumDF DenDF  F value Pr(>F)    
# RAD20  0.0025  0.0025     1    57   1.0418 0.3117    
# Enviro 3.5946  0.5991     6   342 249.9375 <2e-16 ***

BAlmer_mass_FoG <- lmer(Mass ~ FoG20 + Enviro + (1|Strain), DDAMassAct_BA, REML = FALSE)
anova(BAlmer_mass_FoG, type = 3)

# Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq Mean Sq NumDF DenDF  F value Pr(>F)    
# FoG20  0.0000  0.0000     1    57   0.0004 0.9837    
# Enviro 3.5946  0.5991     6   342 249.9375 <2e-16 ***



BAlmer_activity <- lmer(Activity ~ RAD20 + Enviro + (1|Strain), DDAMassAct_BA, REML = FALSE)
anova(BAlmer_activity, type = 3)

# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)    
# RAD20   0.1345 0.13454     1    57   5.548 0.02197 *  
# Enviro 17.8613 2.97688     6   342 122.753 < 2e-16 ***

BAlmer_activity_FoG <- lmer(Activity ~ FoG20 + Enviro + (1|Strain), DDAMassAct_BA, REML = FALSE)
anova(BAlmer_activity_FoG, type = 3)
# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF DenDF  F value Pr(>F)    
# FoG20   0.0017 0.00169     1    57   0.0695  0.793    
# Enviro 17.8613 2.97688     6   342 122.7539 <2e-16 ***


anova(lmer(Activity ~  Enviro + RAD20 + (1|Strain), DDAMassAct_FLC, REML = FALSE))
# Type III Analysis of Variance Table with Satterthwaite's method
#          Sum Sq  Mean Sq NumDF DenDF F value    Pr(>F)    
# RAD20  0.021143 0.021143     1    57  2.5104    0.1186    
# Enviro 0.273328 0.273328     1   342 32.4546 2.633e-08 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


DDAMassAct_FLC_zero <- DDAMassAct_FLC %>% filter(!Enviro == 0) 

anova(lmer(Mass ~ log(Enviro,2) + (1|Strain), DDAMassAct_FLC_zero, REML = FALSE))

anova(lmer(Activity ~ log(Enviro,2) + (1|Strain), DDAMassAct_FLC_zero, REML = FALSE))
# Type III Analysis of Variance Table with Satterthwaite's method
#           Sum Sq   Mean Sq NumDF DenDF F value Pr(>F)
# RAD20  0.0189472 0.0189472     1    57  2.4583 0.1224
# Enviro 0.0007224 0.0007224     1   285  0.0937 0.7597

anova(lmer(Activity ~ RAD20 + Enviro + (1|Strain), DDAMassAct_FLC %>% filter(Enviro != 32), REML = FALSE))
# Type III Analysis of Variance Table with Satterthwaite's method
#           Sum Sq   Mean Sq NumDF DenDF F value Pr(>F)
# RAD20  0.0189472 0.0189472     1    57  2.4583 0.1224
# Enviro 0.0007224 0.0007224     1   285  0.0937 0.7597

FLClmer_activity_FoG <- lmer(Activity ~ FoG50 + Enviro + (1|Strain), DDAMassAct_FLC, REML = FALSE)
anova(FLClmer_activity_FoG, type = 3)
# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
# FoG50  0.04992 0.04992     1    55  8.2815  0.005692 ** 
# Enviro 0.50971 0.50971     1   275 84.5628 < 2.2e-16 ***




SpearmanCor <- function (x, y) {
  cor <- tidy(cor.test(x, y, method = "spearman", na.rm = TRUE))}

# BA Resistance and Biomass Correlation 
DDAMassAct_BA %>% group_by(Enviro) %>% 
  summarise(cor = SpearmanCor(RAD20, Mass))  
# BA Resistance and Activity Correlation 
DDAMassAct_BA %>% group_by(Enviro) %>%
  summarise(cor = SpearmanCor(RAD20, Activity)) 

# BA Tolerance and Biomass Correlation 
DDAMassAct_BA %>% group_by(Enviro) %>% 
  summarise(cor = SpearmanCor(FoG20, Mass)) 
# BA Tolerance and Activity Correlation 
DDAMassAct_BA %>% group_by(Enviro) %>% 
  summarise(cor = SpearmanCor(FoG20, Activity)) 


# FLC Resistance and Biomass Correlation
DDAMassAct_FLC %>% group_by(Enviro) %>% 
  summarise(cor = SpearmanCor(RAD20, Mass))  
# FLC Resistance and Activity Correlation
DDAMassAct_FLC %>% group_by(Enviro) %>% 
  summarise(cor = SpearmanCor(RAD20, Activity)) 
# FLC Tolerance and Biomass Correlation
DDAMassAct_FLC %>% group_by(Enviro) %>% 
  summarise(cor = SpearmanCor(FoG20, Mass)) 
# FLC Tolerance and Activity Correlation
DDAMassAct_FLC %>% group_by(Enviro) %>% 
  summarise(cor = SpearmanCor(FoG20, Activity)) 


FLCResMass <- DDAMassAct_FLC %>% 
  ggplot(aes(x = RAD20, y = Mass)) +
  geom_point() +  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5, linetype = "dashed") +
  xlab("FLC Planktonic Susceptibility") +
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



BAResMass <- DDAMassAct_BA %>% 
  ggplot(aes(x = RAD20, y = Mass)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5,
              linetype = "dashed") +
  xlab("BA Planktonic Susceptibility") +
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


ResMass <- ggarrange(FLCResMass, BAResMass, ncol = 1, align = "hv")

ggsave(here("figures_out", "AllStats", "211027_ResMass.jpg"),
       plot = ResMass,
       height = 12, 
       width = 7)


FLCResAct <- DDAMassAct_FLC %>% 
  ggplot(aes(x = RAD20, y = Activity)) +
  geom_point() +  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5, linetype = "dashed") +
  xlab("FLC Planktonic Susceptibility") +
  ylab("FLC Biofilm Activity") +
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



BAResAct <- DDAMassAct_BA %>% 
  ggplot(aes(x = RAD20, y = Activity)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5,
              linetype = "dashed") +
  xlab("BA Planktonic Susceptibility") +
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


ResAct <- ggarrange(FLCResAct, BAResAct, ncol = 1, align = "hv")

ggsave(here("figures_out", "AllStats", "211027_ResAct.jpg"),
       plot = ResAct,
       height = 12, 
       width = 7)



FLCTolMass <- DDAMassAct_FLC %>% 
  ggplot(aes(x = FoG20, y = Mass)) +
  geom_point() +  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5, linetype = "dashed") +
  xlab("FLC Planktonic Tolerance") +
  ylab("FLC Biofilm Biomass") +
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

BATolMass <- DDAMassAct_BA %>% 
  ggplot(aes(x = FoG20, y = Mass)) +
  geom_point() +  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5, linetype = "dashed") +
  xlab("BA Planktonic Tolerance") +
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


TolMass <- ggarrange(FLCTolMass, BATolMass, ncol = 1, align = "hv")

ggsave(here("figures_out", "AllStats", "211027_TolMass.jpg"),
       plot = TolMass,
       height = 12, 
       width = 7)


FLCTolAct <- DDAMassAct_FLC %>% 
  ggplot(aes(x = FoG20, y = Activity)) +
  geom_point() +  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5, linetype = "dashed") +
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

BATolAct <- DDAMassAct_BA %>% 
  ggplot(aes(x = FoG20, y = Activity)) +
  geom_point() +  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5, linetype = "dashed") +
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


TolAct <- ggarrange(FLCTolAct, BATolAct, ncol = 1, align = "hv")

ggsave(here("figures_out", "AllStats", "211027_TolAct.jpg"),
       plot = TolAct,
       height = 12, 
       width = 7)


