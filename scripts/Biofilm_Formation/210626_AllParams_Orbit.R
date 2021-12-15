library(tidyverse)
library(here)
library(lme4)
library(broom)
library(reshape2)
library(extrafont)
library(ggpubr)
library(growthcurver)

d <- read_csv(here("data_in", "Biofilm_Formation","210623_OrbitOut_EvosB01B02.csv"))

Strain_Names <- tibble( Strain = c("A02","A03","A04","A08","A10","A12","A17","A18"), StrainName = c("P87", "GC75","P78048", "P75016", "P76055", "T101","SC5314", "FH1"))

DF <- d %>% separate(Filename,c("Well","Time"),"_t") %>%
  filter(!Well%in%c("A12", "B12", "C12",
                    "D12","E12", "F12",
                    "G12", "H12")) %>% 
  separate(Time, c("Time","Image"),sep=".j") %>% separate(Time, c("Time", "Biorep"), sep = "_") %>% mutate(Biorep = replace_na(Biorep, "B02")) %>% separate(Well, c("Strain", "Drug", "Concentration"), sep = "_") %>% select(Strain, Time, Drug, Concentration, Biorep, Fungi)

DF$Time<-as.numeric(DF$Time)
DF$Concentration<-as.factor(DF$Concentration)


DF %>% filter(Biorep== "B01") %>% ggplot(aes(x=Time, y=Fungi, color=Concentration))+
  geom_point(aes(), alpha=0.5) + 
  scale_x_continuous(breaks=seq(0,24,4))+
  facet_grid(Strain ~ Drug) +
  scale_color_viridis(discrete = TRUE)+ theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))

DF %>% filter(Biorep== "B02") %>% ggplot(aes(x=Time, y=Fungi, color=Concentration))+
  geom_point(aes(), alpha=0.5) + 
  scale_x_continuous(breaks=seq(0,24,4))+
  facet_grid(Strain ~ Drug) +
  scale_color_viridis(discrete = TRUE)+ theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))


DF_All <- DF %>% group_by(Strain, Time, Drug, Concentration) %>%
  summarise(Fungi_sd= sd(Fungi, na.rm = TRUE),
            Fungi = round (mean(Fungi, na.rm = TRUE), digits = 2)) %>% 
  na.omit()

DF_Wide <- DF_All %>% 
  unite(Strain, Strain, Drug, Concentration, sep = "_" ) %>% 
  select(!Fungi_sd) %>% 
  pivot_wider( names_from = "Strain", values_from = "Fungi") %>% 
  rename("time"="Time") %>% 
  na.omit()

DF_B1B2 <- DF %>% unite(Strain, c("Strain", "Biorep"))


DF_All %>% ggplot(aes(x=Time, y=Fungi, color=Concentration))+
  geom_point(aes(), alpha=0.5) + 
  scale_x_continuous(breaks=seq(0,24,4))+
  facet_grid(Strain ~ Drug) +
  scale_color_viridis(discrete = TRUE)+ theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))


DF_B1B2_wide<- DF_B1B2 %>% 
  unite(Strain, Strain, Drug, Concentration, sep = "_" ) %>% 
  pivot_wider( names_from = "Strain", values_from = "Fungi") %>% 
  rename("time"="Time") %>% 
  na.omit()

logfit_params<- SummarizeGrowthByPlate(plate = DF_Wide) 
B1B2_logfit_params <- SummarizeGrowthByPlate(plate = DF_B1B2_wide) 

# Create a new dataframe, DF_Wide, that contains the predicted % area of the logistic fit 
# This will be used to plot the predicted fit in Basic Protocol 6. 
models_all <- lapply(DF_Wide[2:ncol(DF_Wide)], 
                     function(x) SummarizeGrowth(DF_Wide$time, x)
)

df_predicted_plate <- data.frame(time = DF_Wide$time)

for (i in names(DF_Wide[2:ncol(DF_Wide)])){
  df_predicted_plate[[i]] <- 
    stats:::predict(models_all[[i]]$model)
}


# Create a new dataframe, df_logfit, that contains Orbit % area and logistic fit predicted % area.
df_logfit <- pivot_longer(df_predicted_plate, !time,
                           names_to = "Strain",
                           values_to = "Pred_Area") %>% 
  rename("Time" = "time") %>% 
  separate(Strain, c("Strain", "Drug", "Concentration"), sep = "_") %>% 
  merge(DF_All, by = c("Time", "Strain", "Drug", "Concentration"))
  



df_logfit %>% filter(Drug == "BA") %>% 
  ggplot(aes(x = Time, y = Fungi, color = Concentration)) +
  geom_point(aes(), alpha=0.5) + 
 geom_line(aes(y = Pred_Area)) +
  geom_errorbar(aes(ymin = Fungi-Fungi_sd, ymax = Fungi + Fungi_sd), width = 0) +
  scale_x_continuous(breaks=seq(0,24,4))+
  facet_wrap(~Strain, ncol = 2) +
  scale_color_viridis(discrete = TRUE)+ theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))

df_logfit %>% filter(Drug == "FLC") %>% 
  ggplot(aes(x = Time, y = Fungi, color = Concentration)) +
  geom_point(aes(), alpha=0.5) + 
  geom_line(aes(y = Pred_Area)) +
  geom_errorbar(aes(ymin = Fungi-Fungi_sd, ymax = Fungi + Fungi_sd), width = 0) +
  scale_x_continuous(breaks=seq(0,24,4))+
  facet_wrap(~Strain, ncol = 2) +
  scale_color_viridis(discrete = TRUE)+ theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", 
                          face="bold", size = 8),
        axis.text.x = element_text(hjust = 0.5))


write_csv(df_logfit, 
          file = here("data_out", "Biofilm_Formation", "210624_LogFitOut_B01B02.csv"))



logfit_corrected <- logfit_params %>%
  select(sample, k, r) %>% 
  separate(sample, c("Strain", "Drug", "Concentration"), sep = "_") %>% 
  merge(df_logfit, by = c("Strain", "Drug",
                         "Concentration")) %>%
  mutate(k_corrected = replace(k, round(k, digits = 2) > 1,
                               Pred_Area[which(Time == "24" )]), # Change 24 to your period of incubation
         r_corrected = replace(r, k < 0.15, 0),  # Change 0.15 to the k cutoff where you consider your cells not growing
         k = round(k, digits = 2),
         Pred_Area = round(Pred_Area, digits = 2)) %>% 
  group_by(Strain, Concentration, Drug) %>% 
  summarise(t_asym = Time[which.min(abs (Pred_Area - k))],
            k = mean(k),
            k_corrected = round(mean(k_corrected), digits = 2),
            r = mean(r),
            r_corrected = mean(r_corrected)) %>% 
  mutate(t_asym = replace(t_asym, r_corrected == 0, 28)) %>% 
  merge(Strain_Names, by = "Strain")

 
logfit_corrected_sd <- B1B2_logfit_params %>%
  select(sample, k, r) %>% 
  separate(sample, c("Strain", "Biorep", "Drug", "Concentration"), sep = "_") %>% 
  merge(df_logfit, by = c("Strain", "Drug",
                          "Concentration")) %>%
  group_by(Strain, Concentration, Drug, Biorep) %>% 
  mutate(k_corrected = replace(k, round(k, digits = 2) > 1,
                               Pred_Area[which(Time == "24" )]), # Change 24 to your period of incubation
         r_corrected = replace(r, k < 0.15, 0),  # Change 0.15 to the k cutoff where you consider your cells not growing
         k = round(k, digits = 2),
         Pred_Area = round(Pred_Area, digits = 2),
         t_asym = Time[which.min(abs (Pred_Area - k))]) %>% 
  group_by(Strain, Concentration, Drug) %>% 
  summarise(
            t_asym_sd = sd(t_asym),
            k_sd = sd(k_corrected),
            r_sd = sd(r_corrected),
            t_asym = mean(t_asym),
            k = mean(k),
            k_corrected = round(mean(k_corrected), digits = 2),
            r = mean(r),
            r_corrected = mean(r_corrected)) %>% 
  mutate(t_asym = replace(t_asym, r_corrected == 0, 28), 
         t_asym_sd = replace(t_asym_sd, t_asym == 28, 0)) %>% 
  merge(Strain_Names, by = "Strain")


write_csv(logfit_corrected, 
          file = here("data_out", "Biofilm_Formation", "210624_LogFitParams_B01B02.csv"))

logfit_corrected$Concentration <- as.factor(logfit_corrected$Concentration)




 
##########################################
# Plot the change in asymptote(k) as the concentration of drug increases
##########################################
k_BA <- logfit_corrected %>% filter(Drug %in% c("BA", "ND")) %>% 
  ggplot(aes(x = Concentration, y = k_corrected,
             group = StrainName, color = StrainName)) +
  geom_point(size = 2) +
  geom_line() +
  ylim(0,1) +
  xlab("BA Concentration (mg/mL)") +
  ylab("% Area at Asymptote") +
  labs(color = "Strain") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", face="bold", size = 16),
        axis.text = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  scale_colour_brewer(palette = "Set1")  


k_FLC <- logfit_corrected %>% filter(Drug %in% c("FLC", "ND")) %>% 
  mutate(Concentration = factor(Concentration, levels = c("0","2","4","8","16","32"))) %>% 
  ggplot(aes(x = Concentration, y = k_corrected,
             group = StrainName, color = StrainName)) +
  geom_point(size = 2) +
  geom_line() +
  ylim(0,1) +
  xlab("FLC Concentration (μg/mL)") +
  ylab("% Area at Asymptote") +
  labs(color = "Strain") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", face="bold", size = 16),
        axis.text = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) + 
  scale_colour_brewer(palette = "Set1")

k <- ggarrange(k_FLC, k_BA, 
          legend = "right", 
          common.legend = TRUE )
##########################################
# Plot the change in  slope (r) as the concentration of drug increases 
##########################################
r_BA <- logfit_corrected %>% filter(Drug %in% c("BA", "ND")) %>% 
  ggplot(aes(x = Concentration, y = r_corrected,
             group = StrainName, color = StrainName)) +
  geom_point(size = 2) +
  geom_line() +
  ylim(0,1) +
  xlab("BA Concentration (mg/mL)") +
  ylab("Growth Rate") +
  labs(color = "Strain") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", face="bold", size = 16),
        axis.text = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  scale_colour_brewer(palette = "Set1")

r_FLC <- logfit_corrected %>% filter(Drug %in% c("FLC", "ND")) %>% 
  mutate(Concentration = factor(Concentration, levels = c("0","2","4","8","16","32"))) %>% 
  ggplot(aes(x = Concentration, y = r_corrected,
             group = StrainName, color = StrainName)) +
  geom_point(size = 2) +
  geom_line() +
  ylim(0,1) +
  xlab("FLC Concentration (μg/mL)") +
  ylab("Growth Rate") +
  labs(color = "Strain") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", face="bold", size = 16),
        axis.text = element_text(size = 14, hjust = 0.5), 
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  scale_colour_brewer(palette = "Set1")

r <- ggarrange(r_FLC, r_BA, 
                    legend = "right", 
                    common.legend = TRUE )

##########################################
# Plot the change in  time to reach the asymptote as the concentration of drug increases
##########################################
# Since in this example the images were only captured for 24 h; any well that took more than 24 h to reach the asymptote will be indicated with an open circle 
t_asym_BA <- logfit_corrected %>% filter(Drug %in% c("BA", "ND")) %>% 
  ggplot(aes(x = Concentration, y = t_asym, 
             group = StrainName, color = StrainName)) +
  geom_point(aes(shape = t_asym > 24), size = 2) + # Change 24 to your period of incubation
  geom_line() +
  scale_y_continuous(breaks = seq(0, 28, 4), # Change 28 to your period of incubation + 4.
                     labels = c(seq(0 , 24, 4), "> 24"), # Change 24 to your period of incubation
                     limits = c(0, 28)) + # Change 28 to your period of incubation + 4
  scale_shape_manual(values = c(16, 1), guide = "none") +
  xlab("BA Concentration (mg/mL)") +
  ylab("Time to Asymptote (h)") +
  labs(color = "Strain") +
  theme_bw() +
  theme(legend.position = "right",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", face="bold", size = 16),
        axis.text = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  scale_colour_brewer(palette = "Set1")


t_asym_FLC <- logfit_corrected %>% filter(Drug %in% c("FLC", "ND")) %>% 
  mutate(Concentration = factor(Concentration, levels = c("0","2","4","8","16","32"))) %>% 
  ggplot(aes(x = Concentration, y = t_asym, 
             group = StrainName, color = StrainName)) +
  geom_point(aes(shape = t_asym > 24), size = 2) + # Change 24 to your period of incubation
  geom_line() +
  scale_y_continuous(breaks = seq(0, 28, 4), # Change 28 to your period of incubation + 4.
                     labels = c(seq(0 , 24, 4), "> 24"), # Change 24 to your period of incubation
                     limits = c(0, 28)) + # Change 28 to your period of incubation + 4
  scale_shape_manual(values = c(16, 1), guide = "none") +
  xlab("FLC Concentration (μg/mL)") +
  ylab("Time to Asymptote (h)") +
  labs(color = "Strain") +
  theme_bw() +
  theme(legend.position = "right",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="Times New Roman", face="bold", size = 16),
        axis.text = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  scale_colour_brewer(palette = "Set1")

t_asym <- ggarrange(t_asym_FLC, t_asym_BA, 
          legend = "right", 
          common.legend = TRUE )

param_all <- ggarrange(k_BA, r_BA, t_asym_BA, k_FLC, r_FLC, 
                       t_asym_FLC, common.legend = TRUE,
                       legend = "right", 
                       vjust = 2)



param_all_long <- ggarrange(k_FLC,k_BA, r_FLC, r_BA, t_asym_FLC,  
          t_asym_BA, common.legend = TRUE,
          ncol = 2,
          nrow = 3,
          legend = "right")


ggarrange(k_BA,k_FLC, r_BA, r_FLC, t_asym_BA,  
          t_asym_FLC, BA, FLC,  common.legend = TRUE,
          ncol = 2,
          nrow = 4,
          legend = "right")

ggsave(here("figures_out", "Biofilm_Formation", "210822_AllParams_B01B02.jpg"), 
       param_all, 
       width = 12, height = 6)

ggsave(here("figures_out", "Biofilm_Formation", "210831_AllParams_long_B01B02.jpg"), 
       param_all_long, 
       width = 9, height = 10)
