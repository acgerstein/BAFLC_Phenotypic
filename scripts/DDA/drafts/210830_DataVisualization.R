library(here)
library(tidyverse)
library(cowplot)
library(ggforce)
library(gridExtra)
library(broom)
library(ggpubr)

# Load the data 
df <- read_csv(here("data_out", "DDA","210827_HSCAG.csv")) 


######################################################################
# Determine whether the RAD20 is normally distributed 
shapiro.test(df[df$Species == "C.albicans" & df$Drug == "BA",]$RAD20)
# W = 0.94871, p-value = 1.037e-05
shapiro.test(df[df$Species == "C.albicans" & df$Drug == "FLC",]$RAD20)
# W = 0.83292, p-value = 1.863e-12
shapiro.test(df[df$Species == "C.glabrata" & df$Drug == "BA",]$RAD20)
# 0.92869, p-value = 0.004931
shapiro.test(df[df$Species == "C.glabrata" & df$Drug == "FLC",]$RAD20)
# W = 0.96418, p-value = 0.133
shapiro.test(df[df$Species == "C.parapsilosis" & df$Drug == "BA",]$RAD20)
# W = 0.91261, p-value = 0.0715
shapiro.test(df[df$Species == "C.parapsilosis" & df$Drug == "FLC",]$RAD20)
# W = 0.64966, p-value = 1.008e-05
shapiro.test(df[df$Species == "C.tropicalis" & df$Drug == "BA",]$RAD20)
# W = 0.91517, p-value = 0.4328
shapiro.test(df[df$Species == "C.tropicalis" & df$Drug == "FLC",]$RAD20)
# W = 0.93295, p-value = 0.5763

# C. glb FLC, C. para BA, & C.trop BA and FLC RAD20 are not normally distributed 

# Determine whether the FoG20 is normally distributed 
shapiro.test(df[df$Species == "C.albicans" & df$Drug == "BA",]$FoG20)
# W = 0.71998, p-value = 2.26e-16
shapiro.test(df[df$Species == "C.albicans" & df$Drug == "FLC",]$FoG20)
# W = 0.97913, p-value = 0.01883
shapiro.test(df[df$Species == "C.glabrata" & df$Drug == "BA",]$FoG20)
# 0.97131, p-value = 0.2614
shapiro.test(df[df$Species == "C.glabrata" & df$Drug == "FLC",]$FoG20)
# W = 0.92042, p-value = 0.002739
shapiro.test(df[df$Species == "C.parapsilosis" & df$Drug == "BA",]$FoG20)
# W = 0.67906, p-value = 2.19e-05
shapiro.test(df[df$Species == "C.parapsilosis" & df$Drug == "FLC",]$FoG20)
# W = 0.82621, p-value = 0.002184
shapiro.test(df[df$Species == "C.tropicalis" & df$Drug == "BA",]$FoG20)
# W = 0.937, p-value = 0.6119
shapiro.test(df[df$Species == "C.tropicalis" & df$Drug == "FLC",]$FoG20)
# W = 0.93819, p-value = 0.6446

# C. glb BA & C.trop BA and FLC FoG20 are not normally distributed 
######################################################################
# The Fligner-Killeen’s test is one of the many tests for homogeneity of variances which is most robust against departures from normality.

# Resistance variance homogeneity
fligner.test(RAD20 ~ Drug, data = df[df$Species == "C.albicans",]) 
# Fligner-Killeen:med chi-squared = 23.368, df = 1, p-value = 1.338e-06
fligner.test(RAD20 ~ Drug, data = df[df$Species == "C.glabrata",])
# Fligner-Killeen:med chi-squared = 7.7979, df = 1, p-value = 0.005231
fligner.test(RAD20 ~ Drug, data = df[df$Species == "C.parapsilosis",])
# Fligner-Killeen:med chi-squared =  0.022405, df = 1, p-value = 0.881
fligner.test(RAD20 ~ Drug, data = df[df$Species == "C.tropicalis",])
# Fligner-Killeen:med chi-squared = 0.97458, df = 1, p-value = 0.3235

# Tolerance variance homogeneity
fligner.test(FoG20 ~ Drug, data = df[df$Species == "C.albicans",])
# Fligner-Killeen:med chi-squared = 116.59, df = 1, p-value < 2.2e-16
fligner.test(FoG20 ~ Drug, data = df[df$Species == "C.glabrata",])
# Fligner-Killeen:med chi-squared = 4.942, df = 1, p-value = 0.02621
fligner.test(FoG20 ~ Drug, data = df[df$Species == "C.parapsilosis",])
# Fligner-Killeen:med chi-squared = 0.88105, df = 1, p-value = 0.3479
fligner.test(FoG20 ~ Drug, data = df[df$Species == "C.tropicalis",])
# Fligner-Killeen:med chi-squared = 1.6003, df = 1, p-value = 0.2059
######################################################################



# BA
df_BA <- df %>% filter(Drug == "BA")
kruskal.test(RAD20 ~ Species, data = df_BA)
# Kruskal-Wallis chi-squared = 50.9, df = 3, p-value = 5.138e-11

pairwise.wilcox.test(df_BA$RAD20, df_BA$Species,
                     p.adjust.method = "BH")
#              C.albicans C.glabrata C.parapsilosis
#   C.glabrata     3e-08      -          -             
#   C.parapsilosis 0.00197    0.98442    -             
#   C.tropicalis   0.00025    0.06125    0.14512     

kruskal.test(FoG20 ~ Species, data = df_BA)
#Kruskal-Wallis chi-squared = 109.61, df = 3, p-value < 2.2e-16
pairwise.wilcox.test(df_BA$FoG20, df_BA$Species,
                     p.adjust.method = "BH")
#                C.albicans C.glabrata C.parapsilosis
#   C.glabrata     < 2e-16    -          -             
#   C.parapsilosis 0.0454     1.6e-06    -             
#   C.tropicalis   3.1e-05    0.3492     0.0018   
# 

# FLC
df_FLC <- df %>% filter(Drug == "FLC")

kruskal.test(RAD20 ~ Species, data = df_FLC)
# Kruskal-Wallis chi-squared = 65.47, df = 3, p-value = 3.979e-14

pairwise.wilcox.test(df_FLC$RAD20, df_FLC$Species,
                     p.adjust.method = "BH")
#                C.albicans C.glabrata C.parapsilosis
#   C.glabrata     1.2e-10    -          -             
#   C.parapsilosis 3.0e-05    1.8e-07    -             
#   C.tropicalis   0.35024    0.00083    0.00399    

kruskal.test(FoG20 ~ Species, data = df_FLC)
# Kruskal-Wallis chi-squared = 82.649, df = 3, p-value < 2.2e-16
pairwise.wilcox.test(df_FLC$FoG20, df_FLC$Species,
                     p.adjust.method = "BH")
#               C.albicans C.glabrata C.parapsilosis
#   C.glabrata     4.3e-11    -          -             
#   C.parapsilosis 4.3e-11    2.5e-05    -             
#   C.tropicalis   0.60160    0.03568    0.00088    


df_stat <- df %>% group_by(Species, Drug) %>% 
  summarise(RAD20Med = median(RAD20, na.rm = TRUE), 
            FoG20Med = median(FoG20, na.rm = TRUE),
            RAD50Med = median(RAD50, na.rm = TRUE), 
            FoG50Med = median(FoG50, na.rm = TRUE))

RADFoG_MaxMin <- df %>% group_by(Species, Drug) %>% summarise(RAD_min = min(RAD20), FoG_min = min(FoG20, na.rm = TRUE), RAD_max = max(RAD20), FoG_max = max(FoG20, na.rm = TRUE))



df_stat_BA <- df_stat %>% filter(Drug == "BA")
df_stat_FLC <- df_stat %>% filter(Drug == "FLC")

#Plotting the full data
R_BA <- df %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=RAD20), 
            maxwidth = 0.50, alpha = 0.6, size=4) + 

  geom_point( df_stat_BA, mapping = aes(x=Species, y= RAD20Med),
              color="red", pch="_", alpha=0.8, lwd = 20.5,
              position= position_dodge(width = 0.4)) +
  annotate(y = RADFoG_MaxMin[RADFoG_MaxMin$Drug == "BA",]$RAD_min - 2, 
           geom = "text", x= c(1, 2, 3, 4), label = c("b","a","a","a"),
           color="black", size = 6) +
  scale_x_discrete(label = c("C. albicans", 
                             "C. glabrata", "C. parapsilosis",
                             "C. tropicalis")) +
  ylab("Susceptibility") + expand_limits(y = c(-4,28))+
  theme_bw()+ labs(subtitle = "Boric Acid")+
  scale_y_reverse()+
  scale_color_manual(values = c("#333333", "#660066"))+
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 18), 
        axis.text.x = element_text(hjust = 0.5, face = "bold.italic")) 


#FLC Resistance graph
R_FLC <- df %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x = Species, y = RAD20), 
            maxwidth = 0.50,alpha = 0.6, size = 4) + 
  geom_point( df_stat_FLC, mapping = aes(x = Species, y = RAD20Med), 
              color = "red", 
              pch="_", alpha=0.8, lwd = 20.5,
              position= position_dodge(width = 0.4)) +
  annotate(y = RADFoG_MaxMin[RADFoG_MaxMin$Drug == "FLC",]$RAD_min - 2, 
           geom = "text", x= c(1, 2, 3, 4), label = c("b","a","c","b"),
           color ="black", size = 6) +
  scale_x_discrete(label = c("C. albicans", "C. glabrata", "C. parapsilosis", "C. tropicalis")) +
  ylab("Susceptibility") + expand_limits(y = c(-4,28))+
  theme_bw()+ labs(subtitle = "Fluconazole")+
  scale_y_reverse()+
  scale_color_manual(values = c("#333333", "#660066"))+
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family ="Times New Roman", 
                          face ="bold", size = 18), 
        axis.text.x = element_text(hjust = 0.5, face = "bold.italic"))

#BA Tolerance graph 
T_BA <- df %>% filter(Drug=="BA") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG20), 
            maxwidth = 0.50,alpha = 0.6, size=4) + 
  annotate(y = RADFoG_MaxMin[RADFoG_MaxMin$Drug == "BA",]$FoG_max + 0.1, 
           geom = "text", x= c(1, 2, 3, 4), label = c("b","a","c","a"),
           color ="black", size = 6) +
  geom_point( df_stat_BA, mapping = aes(x=Species, y= FoG20Med),
              color = "red", 
              pch="_", alpha=0.8, lwd = 20.5,
              position= position_dodge(width = 0.5))+
  scale_x_discrete(label = c("C. albicans", "C. glabrata", "C. parapsilosis", "C. tropicalis")) +
  scale_y_continuous(breaks = seq(0,1, 0.25)) +
  ylab("Tolerance") + 
  expand_limits( y = c(0,1.1))+
  theme_bw() + 
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 18),
        axis.text.x = element_text(hjust = 0.5, face = "bold.italic"))

#FLC Tolerance graph 
T_FLC <- df %>% filter(Drug=="FLC") %>% 
  ggplot() +
  geom_sina(mapping= aes(x=Species, y=FoG20), 
            maxwidth = 0.50,alpha = 0.6, size=4) + 
  annotate(y = RADFoG_MaxMin[RADFoG_MaxMin$Drug == "FLC",]$FoG_max + 0.1, 
           geom = "text", x= c(1, 2, 3, 4), label = c("a","b","c","a"),
           color ="black", size = 6) +
  geom_point(df_stat_FLC, mapping = aes(x=Species, y= FoG20Med),
             color = "red", 
             pch="_", alpha=0.8, lwd = 20.5,
             position= position_dodge(width = 0.5))+
  scale_x_discrete(label = c("C. albicans", "C. glabrata", "C. parapsilosis", "C. tropicalis"))+
  scale_y_continuous(breaks = seq(0,1, 0.25)) +
  ylab("Tolerance") + 
  expand_limits( y = c(0,1.1))+
  theme_bw() + 
  theme(legend.position="none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 18),
        axis.text.x = element_text(hjust = 0.5, face = "bold.italic"))



p_all <- ggarrange(R_FLC,
          R_BA,
          T_FLC,
          T_BA, 
          align = "hv")


ggsave(plot = p_all, here("figures_out", "DDA",
                                "210830_HSCAG.jpg"), 
                          width = 12, 
                          height = 9)
# Susceptibility and Tolerance correlation

# r_pvalue <- function (x, y) {
#   cor_func <- tidy(cor.test(x, y, method = "spearman", na.rm = TRUE))
# paste(cor_func$p.value)}

SpearmanCor <- function (x, y) {
  cor <- tidy(cor.test(x, y, method = "spearman", na.rm = TRUE))}

# cor.test(df[df$Species == "C.albicans" & df$Drug == "BA",]$RAD20, 
#          df[df$Species == "C.albicans" & df$Drug == "BA",]$FoG20, 
#          method = "spearman", na.rm = TRUE)
# 
# cor.test(df[df$Species == "C.albicans" & df$Drug == "FLC",]$RAD20, 
#          df[df$Species == "C.albicans" & df$Drug == "FLC",]$FoG20, 
#          method = "spearman", na.rm = TRUE)

df %>% group_by(Species, Drug) %>% 
  summarise(cor = SpearmanCor(RAD20, FoG20)) 

FLCBA_Cor <- df %>% filter(Drug == "BA") %>% rename("BA_RAD20" = "RAD20", "BA_FoG20" = "FoG20", "BA_RAD50" = "RAD50", "BA_FoG50" = "FoG50") %>% 
  select(!c(Drug)) %>% 
  merge(df %>% filter(Drug == "FLC"), by = c("Strain", "Species")) %>% 
  rename("FLC_RAD50" = "RAD50", "FLC_FoG50" = "FoG50", 
         "FLC_RAD20" = "RAD20", "FLC_FoG20" = "FoG20") %>% 
  select(!c(Drug))
  
FLCBA_Cor %>% group_by(Species) %>% 
  summarise(cor = SpearmanCor(FLC_RAD20, FLC_FoG20)) 

FLCBA_Cor %>% group_by(Species) %>% 
  summarise(cor = SpearmanCor(BA_RAD20, BA_FoG20)) 

FLCBA_Cor %>% group_by(Species) %>% 
  summarise(cor = SpearmanCor(BA_RAD20, FLC_RAD20)) 
FLCBA_Cor %>% group_by(Species) %>% 
  summarise(cor = SpearmanCor(BA_FoG20, FLC_FoG20)) 

Correlation <- df %>% group_by(Species, Drug) %>% 
  na.omit() %>% 
  summarise(rho = cor(RAD20, FoG20, method = "spearman"),
            p_value = r_pvalue(RAD20, FoG20))

Correlation$p_value <- as.numeric(Correlation$p_value)

 BART <- df %>% merge(Correlation, by = c("Species", "Drug")) %>% 
  filter(Drug == "BA") %>% 
  ggplot(aes(x = FoG20, y = RAD20, color = Species)) +
  geom_point(size = 3, alpha = 0.6) + 
  geom_smooth(aes(linetype = p_value < 0.05 , group = Species ),
              method = "lm", se = FALSE, 
              size = 0.8) +
   scale_y_reverse(limits = c(25, 0)) +
  xlab("BA Tolerance") +
  ylab("BA Susceptibility") +
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
  scale_color_manual(values = c("#CC79A7", "#42007B", "#009E73", "#DD5500"), 
                     label = c("C. albicans", "C. glabrata", "C. parapsilosis",
                               "C. tropicalis"))

 BART
 
 FLCRT <- df %>% 
   merge(Correlation, by = c("Species", "Drug")) %>% 
   filter(Drug == "FLC") %>% 
   ggplot(aes(x = FoG20, y = RAD20, color = Species)) +
   geom_point(size = 3, alpha = 0.6) + 
   geom_smooth(aes(linetype = p_value < 0.05, group = Species ),
               method = "lm", se = FALSE, 
               size = 0.8) +
   scale_y_reverse(limits = c(25, 0)) +
   xlim(0, 1) +
   xlab("FLC Tolerance") +
   ylab("FLC Susceptibility") +
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
   scale_color_manual(values = c("#CC79A7", "#42007B", "#009E73", "#DD5500"), label = c("C. albicans", "C. glabrata", "C. parapsilosis",
                                "C. tropicalis"))
 
 
 FLCRT




FLCBARADCor <- df %>% filter(Drug == "BA") %>% rename("BA_RAD20" = "RAD20", "BA_FoG20" = "FoG20") %>% 
  select(!c(Drug)) %>% 
  merge(df %>% filter(Drug == "FLC"), by = c("Strain", "Species")) %>% 
  rename("FLC_RAD20" = "RAD20", "FLC_FoG20" = "FoG20") %>% 
  select(!c(Drug)) %>% 
  ggplot(aes(x = BA_RAD20, y = FLC_RAD20, color = Species)) +
  geom_point(size = 3, alpha = 0.6) + 
  geom_smooth(method = "lm", linetype = "longdash", se = FALSE)+
  xlab("BA Susceptibility") +
  ylab("FLC Susceptibility") +
  scale_x_reverse() +
  scale_y_reverse() +
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
  scale_color_manual(values = c("#CC79A7", "#42007B", "#009E73", "#DD5500"), 
                     label = c("C. albicans", "C. glabrata", "C. parapsilosis",
                               "C. tropicalis"))

FLCBARADCor

FLCBAFoGCor <- df %>% filter(Drug == "BA") %>% rename("BA_RAD20" = "RAD20", "BA_FoG20" = "FoG20") %>% 
  select(!c(Drug)) %>% 
  merge(df %>% filter(Drug == "FLC"), by = c("Strain", "Species")) %>% 
  rename("FLC_RAD20" = "RAD20", "FLC_FoG20" = "FoG20") %>% 
  select(!c(Drug)) %>% 
  ggplot(aes(x = BA_FoG20, y = FLC_FoG20, color = Species)) +
  geom_point(size = 3, alpha = 0.6) + 
  geom_smooth(method = "lm", linetype = "longdash", se = FALSE)+
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
  scale_color_manual(values = c("#CC79A7", "#42007B", "#009E73", "#DD5500"), 
                     label = c("C. albicans", "C. glabrata", "C. parapsilosis",
                               "C. tropicalis"))

FLCBAFoGCor
FLCBACor <- ggarrange(FLCRT, BART, FLCBARADCor, FLCBAFoGCor,  common.legend = TRUE,
                      legend = "right", align = "hv" )

FLCBACor

ggsave(plot = FLCBACor, here("figures_out", "DDA",
                          "210831_FLCBACor.jpg"), 
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
