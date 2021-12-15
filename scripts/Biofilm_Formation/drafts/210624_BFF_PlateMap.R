#only install it once. select the "all" option
#install.packages("remotes")
#remotes::install_github("briandconnelly/ggplot2bdc")

#load libraries
library(tidyverse)
library(ggplot2)
library(ggplot2bdc)
library(here)
#load files

row <- tibble(row = c("A", "B", "C", "D", "E", "F", "G", "H"), strain =c("P87", "GC75","P78048", "P75016", "P76055", "T101","SC5314", "FH1")) 

column <- tibble(column =1:12, drug = c("No drug","FLC", "BA", "FLC", "BA","FLC", "BA","FLC", "BA","FLC" ,"BA", "Blank"), concentration = c(0, 2, 0.4, 4, 0.8, 8, 1.6, 16, 3.2, 32, 6.4, 0 ), 
                 concentration_level = c(0, 2, 2, 4, 4, 8, 8, 16, 16, 32, 32, 0 ))

platemap <- crossing(row, column) %>% 
  mutate(well = paste(row,column, sep = "")) %>% 
  mutate(strain = replace(strain, column==12, NA))


#platemap

platemap$row <- as.factor(platemap$row)
platemap$row <- as.numeric(platemap$row)


ggplot(data = platemap, aes(x=column, y = row)) +
  geom_point(aes(alpha= concentration_level, fill= drug),shape=21, size=15 ) +
  coord_fixed(ratio=(13/12)/(9/8), xlim=c(0.5, 12.5), ylim=c(8.5,0.5)) +
  scale_y_reverse(breaks=seq(1, 8), labels=c("A","B","C","D","E","F","G","H")) +
  scale_x_continuous(breaks= seq(1,12,1),  position = "top") +
  geom_text(data= platemap, label=platemap$strain, size= 3, color= "Black", family = "Times New Roman", fontface = "bold")+
  theme_bdc_microtiter()+theme(legend.position = "none", 
    legend.title = element_text( family="Times New Roman", color = "Black",face="bold", size = 8),
    legend.text = element_text(family="Times New Roman", color = "Black",face="bold", size = 8),
    axis.text=element_text(family="Times New Roman", size=12), 
    plot.margin = margin(t = 20, b= 20, unit = "pt"))+  guides(alpha = FALSE) +
  scale_fill_manual(values = c("dodgerblue3", "white", "firebrick1", "#999991"))
 
ggsave(here("figures_out", "Biofilm_Formation", "210624_PlateMap_B01B02.jpg"), 
       last_plot(), 
       width = 9, height = 5)
