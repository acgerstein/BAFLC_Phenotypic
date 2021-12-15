library(here)
library(tidyverse)
library(cowplot)
library(ggforce)
library(gridExtra)
library(broom)
library(ggpubr)

# Load the data 
df <- read_csv(here("data_out", "DDA","210629_HSCAG.csv")) 

######################################################################
# Determine whether the data is normally distributed 
shapiro.test(df[df$Species == "C.albicans" & df$Drug == "BA",]$RAD50)
# W = 0.96099, p-value = 0.0002217
shapiro.test(df[df$Species == "C.albicans" & df$Drug == "FLC",]$RAD50)
# W =0.85506, p-value = 4.167e-11
shapiro.test(df[df$Species == "C.glabrata" & df$Drug == "BA",]$RAD50)
# W = 0.96914, p-value = 0.2456
shapiro.test(df[df$Species == "C.glabrata" & df$Drug == "FLC",]$RAD50)
# W = 0.96949, p-value = 0.2532
shapiro.test(df[df$Species == "C.parapsilosis" & df$Drug == "BA",]$RAD50)
# W = 0.89192, p-value = 0.02917
shapiro.test(df[df$Species == "C.parapsilosis" & df$Drug == "FLC",]$RAD50)
# W = 0.76384, p-value = 0.0002618
shapiro.test(df[df$Species == "C.tropicalis" & df$Drug == "BA",]$RAD50)
# W = 0.8688, p-value = 0.1812
shapiro.test(df[df$Species == "C.tropicalis" & df$Drug == "FLC",]$RAD50)
# W = 0.9692, p-value = 0.8926

# C. alb, C. para are not normally distributed for BA and FLC

######################################################################
# The Fligner-Killeen’s test is one of the many tests for homogeneity of variances which is most robust against departures from normality.

# Resistance variance homogeneity
fligner.test(RAD50 ~ Drug, data = df[df$Species == "C.albicans",]) 
# Fligner-Killeen:med chi-squared = 55.573, df = 1, p-value = 9.003e-14
fligner.test(RAD50 ~ Drug, data = df[df$Species == "C.glabrata",])
# Fligner-Killeen:med chi-squared = 6.951, df = 1, p-value = 0.008377
fligner.test(RAD50 ~ Drug, data = df[df$Species == "C.parapsilosis",])
# Fligner-Killeen:med chi-squared = 0.85179, df = 1, p-value = 0.356
fligner.test(RAD50 ~ Drug, data = df[df$Species == "C.tropicalis",])
# Fligner-Killeen:med chi-squared = 0.77419, df = 1, p-value = 0.3789

# Tolerance variance homogeneity
fligner.test(FoG50 ~ Drug, data = df[df$Species == "C.albicans",])
# Fligner-Killeen:med chi-squared = 127.88, df = 1, p-value < 2.2e-16
fligner.test(FoG50 ~ Drug, data = df[df$Species == "C.glabrata",])
# Fligner-Killeen:med chi-squared = 4.4466, df = 1, p-value = 0.03497
fligner.test(FoG50 ~ Drug, data = df[df$Species == "C.parapsilosis",])
# Fligner-Killeen:med chi-squared = 0.01026, df = 1, p-value = 0.9193
fligner.test(FoG50 ~ Drug, data = df[df$Species == "C.tropicalis",])
# Fligner-Killeen:med chi-squared = 3.2457, df = 1, p-value = 0.07161
######################################################################
# var.test(RAD50 ~ Drug, data = df[df$Species == "C.albicans",])
# var.test(RAD50 ~ Drug, data = df[df$Species == "C.glabrata",])
# var.test(RAD50 ~ Drug, data = df[df$Species == "C.parapsilosis",])
# var.test(RAD50 ~ Drug, data = df[df$Species == "C.tropicalis",])
# 
# var.test(FoG50 ~ Drug, data = df[df$Species == "C.albicans",])
# var.test(FoG50 ~ Drug, data = df[df$Species == "C.glabrata",])
# var.test(FoG50 ~ Drug, data = df[df$Species == "C.parapsilosis",])
# var.test(FoG50 ~ Drug, data = df[df$Species == "C.tropicalis",])
######################################################################

OneWayAnovaBA_RAD <- aov(RAD50 ~ Species, data = df %>% filter(Drug == "BA"))
summary(OneWayAnovaBA_RAD) 
#              Df Sum Sq Mean Sq F value Pr(>F)    
#Species       3  377.1  125.70   35.12 <2e-16 ***
#Residuals   226  809.0    3.58                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# BA
df_BA <- df %>% filter(Drug == "BA")
kruskal.test(RAD50 ~ Species, data = df_BA)
# Kruskal-Wallis chi-squared = 80.543, df = 3, p-value < 2.2e-16
kruskal.test(FoG50 ~ Species, data = df_BA)
#Kruskal-Wallis chi-squared = 105.84, df = 3, p-value < 2.2e-16

pairwise.wilcox.test(df_BA$RAD50, df_BA$Species,
                     p.adjust.method = "BH")
#                C.albicans C.glabrata C.parapsilosis
#   C.glabrata     4.6e-16    -          -             
#   C.parapsilosis 0.03313    0.00172    -             
#   C.tropicalis   0.00017    0.90756    0.03313     

pairwise.wilcox.test(df_BA$FoG50, df_BA$Species,
                     p.adjust.method = "BH")
#                C.albicans C.glabrata C.parapsilosis
#  C.glabrata     < 2e-16    -          -             
#  C.parapsilosis 0.1124     2.2e-06    -             
#  C.tropicalis   3.4e-05    0.0518     0.0066     


# FLC
df_FLC <- df %>% filter(Drug == "FLC")

kruskal.test(RAD50 ~ Species, data = df_FLC)

# Kruskal-Wallis chi-squared = 45.054, df = 3, p-value = 9.012e-10
kruskal.test(FoG50 ~ Species, data = df_FLC)
# Kruskal-Wallis chi-squared = 91.483, df = 3, p-value < 2.2e-16

pairwise.wilcox.test(df_FLC$RAD50, df_FLC$Species,
                     p.adjust.method = "BH")
#                C.albicans C.glabrata C.parapsilosis
#  C.glabrata     5.2e-06    -          -             
#  C.parapsilosis 3.6e-05    2.0e-06    -             
#  C.tropicalis   0.62030    0.00038    0.02396  

pairwise.wilcox.test(df_FLC$FoG50, df_FLC$Species,
                     p.adjust.method = "BH")
#               C.albicans C.glabrata C.parapsilosis
#  C.glabrata     2.4e-12    -          -             
#  C.parapsilosis 1.6e-11    1.3e-06    -             
#  C.tropicalis   0.29801    0.01074    0.00018    

tidy(TukeyHSD(OneWayAnovaBA_RAD ))

OneWayAnovaFLC_RAD <- aov(RAD50 ~ Species, data = df %>% filter(Drug == "FLC"))
summary(OneWayAnovaFLC_RAD) 

tidy(TukeyHSD(OneWayAnovaFLC_RAD))

df_stat <- df %>% group_by(Species, Drug) %>% 
  summarise(RAD50Med = median(RAD50, na.rm = TRUE), 
            FoG50Med = median(FoG50, na.rm = TRUE))

df %>% group_by(Species, Drug) %>% summarise(RAD_min = min(RAD50), FoG_min = min(FoG50, na.rm = TRUE), RAD_max = max(RAD50), FoG_max = max(FoG50, na.rm = TRUE))

df %>% separate(Strain, c("type", "number"), sep = "AG") %>% 
  filter(type == "y") %>% 
  group_by(Species) %>% filter(Drug == "BA") %>% 
  summarise(n = length(type))

df %>% filter(Drug == "BA") %>% 
  summarise(n = length(Strain))

df_stat_BA <- df_stat %>% filter(Drug == "BA")
df_stat_FLC <- df_stat %>% filter(Drug == "FLC")
#Plotting the full data
R_BA <- df %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50), 
            maxwidth = 0.50,alpha = 0.6, size=4) + 
  annotate(geom = "text", x= c(1, 2, 3, 4), y= -2, label = c("a","b","c","b"),
           color="black", size = 8) +
  geom_point( df_stat_BA, mapping = aes(x=Species, y= RAD50Med),
              color="red",pch="_",size=10, alpha=0.8, 
              position= position_dodge(width = 0.5)) +
  scale_x_discrete(label = c("C. albicans", "C. glabrata", "C. parapsilosis", "C. tropicalis")) +
  ylab("RAD50") + expand_limits(y = c(-2,28))+
  theme_bw()+ labs(subtitle = "BA")+
  scale_y_reverse()+
  scale_color_manual(values = c("#333333", "#660066"))+
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 18)) 

#FLC Resistance graph
R_FLC <- df %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD50), 
            maxwidth = 0.50,alpha = 0.6, size=4) + 
  annotate(geom = "text", x= c(1, 2, 3, 4), y= -2, label = c("a","b","c","a"),
           color="black", size = 8) +
  geom_point( df_stat_FLC, mapping = aes(x=Species, y= RAD50Med),
              color="red",pch="_",size=8, alpha=0.8, 
              position= position_dodge(width = 0.5))+
  scale_x_discrete(label = c("C. albicans", "C. glabrata", "C. parapsilosis", "C. tropicalis")) +
  ylab("RAD50") + expand_limits(y = c(-2,28))+
  theme_bw()+ labs(subtitle = "FLC")+
  scale_y_reverse()+
  scale_color_manual(values = c("#333333", "#660066"))+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 18))

#BA Tolerance graph 
T_BA <- df %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG50), 
            maxwidth = 0.50,alpha = 0.6, size=4) + 
  annotate(geom = "text", x= c(1, 2, 3, 4), y = 1.05, label = c("a","b","a","b"),
           color="black", size = 8) +
  geom_point( df_stat_BA, mapping = aes(x=Species, y= FoG50Med),
              color="red",pch="_",size=10, alpha=0.8, 
              position= position_dodge(width = 0.5))+
  scale_x_discrete(label = c("C. albicans", "C. glabrata", "C. parapsilosis", "C. tropicalis")) +
  ylab("FoG50") + 
  expand_limits( y = c(0,1))+
  theme_bw() + 
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 18),
        axis.text.x = element_text(hjust = 0.5, face = "bold.italic"))

#FLC Tolerance graph 
T_FLC <- df %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG50), 
            maxwidth = 0.50,alpha = 0.6, size=4) + 
  annotate(geom = "text", x= c(1, 2, 3, 4), y = 1.05, label = c("a","b","c","a"),
           color="black", size = 8) +
  geom_point( df_stat_FLC, mapping = aes(x=Species, y= FoG50Med),
              color="red",pch="_",size=10, alpha=0.8, 
              position= position_dodge(width = 0.5))+
  scale_x_discrete(label = c("C. albicans", "C. glabrata", "C. parapsilosis", "C. tropicalis"))+
  ylab("FoG50") + 
  expand_limits( y = c(0,1))+
  theme_bw() + 
  theme(legend.position="none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 18),
        axis.text.x = element_text(hjust = 0.5, face = "bold.italic"))


p_all<- ggarrange(R_FLC,
          R_BA,
          T_FLC,
          T_BA, 
          align = "hv")


ggsave(plot = p_all, here("figures_out", "DDA",
                                "210809_HSCAG.jpg"), 
                          width = 12, 
                          height = 8)

rsq <- function (x, y) cor(x, y, method = "spearman") ^ 2

r_pvalue<- function (x, y) {
  cor <- tidy(cor.test(x, y, method = "spearman"))
paste(cor$p.value)}

SpearmanCor <- function (x, y) {
  cor <- tidy(cor.test(x, y, method = "spearman", na.rm = TRUE))}
  

df %>% group_by(Species, Drug) %>% 
  summarise(cor = SpearmanCor(RAD50, FoG50)) 

FLCBA_Cor<- df %>% filter(Drug == "BA") %>% rename("BA_RAD50" = "RAD50", "BA_FoG50" = "FoG50") %>% 
  select(!c(Drug, Year)) %>% 
  merge(df %>% filter(Drug == "FLC"), by = c("Strain", "Species")) %>% 
  rename("FLC_RAD50" = "RAD50", "FLC_FoG50" = "FoG50") %>% 
  select(!c(Drug, Year))
  
FLCBA_Cor %>% group_by(Species) %>% 
  summarise(cor = SpearmanCor(BA_RAD50, FLC_RAD50)) 
FLCBA_Cor %>% group_by(Species) %>% 
  summarise(cor = SpearmanCor(BA_FoG50, FLC_FoG50)) 

Correlation <- df %>% group_by(Species, Drug) %>% 
  na.omit() %>% 
  summarise(cor = cor(RAD50, FoG50),
    r_squared = rsq(RAD50, FoG50),
            p_value = r_pvalue(RAD50, FoG50))

Correlation$p_value <- as.numeric(Correlation$p_value)

 BART <- df %>% merge(Correlation, by = c("Species", "Drug")) %>% 
  filter(Drug == "BA") %>% 
  ggplot(aes(x = RAD50, y = FoG50, color = Species)) +
  geom_point(size = 3, alpha = 0.6) + 
 
  stat_cor( aes(label = paste(..rr.label..,
                              if_else(readr::parse_number(..p.label..) < 0.001, 
                                      "p<0.001", ..p.label..), sep = "~`,   `~")),
            method = "spearman",
            label.x = 0,
            digits = 2,
            label.y = c(0.98, 0.93, 0.88, 0.83),
            label.sep = ",",
            size = 3) +
  geom_smooth(aes(linetype = p_value < 0.05 & r_squared > 0.4, group = Species ),
              method = "lm", se = FALSE, 
              size = 0.5) +
  xlab("BA Resistance") +
  ylab("BA Tolerance") +
  ylim(seq(0,1.08))+
  theme_bw() +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "white",
                                        color = "black",
                                        size = 1),
        strip.text = element_text(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                            face = "bold",
                            size = 18),
        axis.text.x = element_text(hjust = 0.5),
        legend.text = element_text(family = "Times New Roman", 
                                     face = "italic",
                                     size = 18)) +
  scale_linetype_manual(values = c('dashed', 'solid'),guide = "none") +
  scale_color_manual(values = c("#882255", "#0072B2", "#009E73", "#999933"), 
                     label = c("C. albicans", "C. glabrata", "C. parapsilosis",
                               "C. tropicalis"))


 FLCRT <- df %>% merge(Correlation, by = c("Species", "Drug")) %>% 
   filter(Drug == "FLC") %>% 
   ggplot(aes(x = RAD50, y = FoG50, color = Species)) +
   geom_point(size = 3, alpha = 0.6) + 
   
   stat_cor( aes(label = paste(..rr.label..,
                               if_else(readr::parse_number(..p.label..) < 0.001, 
                                       "p<0.001", ..p.label..), sep = "~`,   `~")),
             method = "spearman",
             label.x = 0,
             digits = 2,
             label.y = c(0.15, 0.10, 0.05, 0.0),
             label.sep = ",",
             size = 3) +
   geom_smooth(aes(linetype = p_value < 0.05 & r_squared > 0.4, group = Species ),
               method = "lm", se = FALSE, 
               size = 0.5) +
   xlab("FLC Resistance") +
   ylab("FLC Tolerance") +
   ylim(seq(0,1.08))+
   theme_bw() +
   theme(legend.position = "right",
         strip.background = element_rect(fill = "white",
                                         color = "black",
                                         size = 1),
         strip.text = element_text(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         text = element_text(family = "Times New Roman", 
                             face = "bold",
                             size = 18),
         axis.text.x = element_text(hjust = 0.5),
         legend.text = element_text(family = "Times New Roman", 
                                    face = "italic",
                                    size = 18)) +
   scale_linetype_manual(values = c('dashed', 'solid'),guide = "none") +
   scale_color_manual(values = c("#882255", "#0072B2", "#009E73", "#999933"), 
                      label = c("C. albicans", "C. glabrata", "C. parapsilosis",
                                "C. tropicalis"))
 





FLCBARADCor <- df %>% filter(Drug == "BA") %>% rename("BA_RAD50" = "RAD50", "BA_FoG50" = "FoG50") %>% 
  select(!c(Drug, Year)) %>% 
  merge(df %>% filter(Drug == "FLC"), by = c("Strain", "Species")) %>% 
  rename("FLC_RAD50" = "RAD50", "FLC_FoG50" = "FoG50") %>% 
  select(!c(Drug, Year)) %>% 
  ggplot(aes(x = BA_RAD50, y = FLC_RAD50, color = Species)) +
  geom_point(size = 3, alpha = 0.6) + 
  geom_smooth(method = "lm", linetype = "longdash", se = FALSE)+
  stat_cor( aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman",
            label.sep = ",", 
            size = 3, 
           label.x = 3,
           label.y = c(5, 3.8, 2.6, 1.4),
           digits = 2,
            position = "identity") +
  xlab("BA Resistance") +
  ylab("FLC Resistance") +
  theme_bw() +
  guides(linetype = FALSE) +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                            face = "bold", size = 18),
        axis.text.x = element_text(hjust = 0.5),
        legend.text = element_text(family = "Times New Roman", 
                                   face = "italic",
                                   size = 18)) +
  scale_color_manual(values = c("#882255", "#0072B2", "#009E73", "#999933"), 
                     label = c("C. albicans", "C. glabrata", "C. parapsilosis",
                               "C. tropicalis"))

FLCBAFoGCor <- df %>% filter(Drug == "BA") %>% rename("BA_RAD50" = "RAD50", "BA_FoG50" = "FoG50") %>% 
  select(!c(Drug, Year)) %>% 
  merge(df %>% filter(Drug == "FLC"), by = c("Strain", "Species")) %>% 
  rename("FLC_RAD50" = "RAD50", "FLC_FoG50" = "FoG50") %>% 
  select(!c(Drug, Year)) %>% 
  ggplot(aes(x = BA_FoG50, y = FLC_FoG50, color = Species)) +
  geom_point(size = 3, alpha = 0.6) + 
  geom_smooth(method = "lm", linetype = "longdash", se = FALSE)+
  stat_cor( aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
            method = "spearman",
            label.sep = ",", 
            size = 3, 
            label.x = 0.35,
            digits = 2,
            label.y = c(1, 0.95, 0.90, 0.85),
            position = "identity") +
  ylim(0,1) +
  xlab("BA Tolerance") +
  ylab("FLC Tolerance") +
  theme_bw() +
  guides(linetype = FALSE) +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", 
                            face = "bold", size = 18),
        axis.text.x = element_text(hjust = 0.5),
        legend.text = element_text(family = "Times New Roman", 
                                   face = "italic",
                                   size = 18)) +
  scale_color_manual(values = c("#882255", "#0072B2", "#009E73", "#999933"), 
                     label = c("C. albicans", "C. glabrata", "C. parapsilosis",
                               "C. tropicalis"))

FLCBACor <- ggarrange(BART, FLCRT, FLCBAFoGCor,  FLCBARADCor, common.legend = TRUE,
                      legend = "right")

ggsave(plot = FLCBACor, here("figures_out", "DDA",
                          "210701_FLCBACor.jpg"), 
       width = 12,
       height = 8)


#Species impact drug susceptibility significantly 
df_BA <- df %>% filter(Drug == "BA")%>% na.omit()
df_FLC <- df %>% filter(Drug == "FLC") %>% na.omit()




tidy(aov(RAD50~Species, data = df_BA))
tidy(aov(FoG50~Species, data = df_BA))

tidy(aov(RAD50~Species, data = df_FLC))
tidy(aov(FoG50~Species, data = df_FLC))


TukeyHSD(aov(RAD50~Species, data = df_BA))
TukeyHSD(aov(FoG50~Species, data = df_BA))

TukeyHSD(aov(RAD50~Species, data = df_FLC))
TukeyHSD(aov(FoG50~Species, data = df_FLC))
