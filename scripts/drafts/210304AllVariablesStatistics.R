library(tidyverse)
library(here)

# Files ----
d <- read_csv(here("data_out", "Combined", "210303MassActDDA.csv"))
names(d)[5:6] <- c("Biomass", "Activity")

d20 <- subset(d, Enviro=="0.2")
d39 <- subset(d, Enviro=="0.39")
d78 <- subset(d, Enviro=="0.78")
d312 <- subset(d, Enviro=="3.12")
d625 <- subset(d, Enviro=="6.25")
d1250 <- subset(d, Enviro=="12.5")

d20_BA <- subset(d20, Drug == "BA")
d39_BA <- subset(d39, Drug == "BA")
d78_BA <- subset(d78, Drug == "BA")
d312_BA <- subset(d312, Drug == "BA")
d625_BA <- subset(d625, Drug == "BA")
d1250_BA <- subset(d1250, Drug == "BA")

#Test for normality - planktonic variables are not normal so should be using non-parametric stats = spearman rank correlation instead of the default (pearson)
shapiro.test(d312_BA$RAD50)
shapiro.test(d312_BA$FoG50)
shapiro.test(d625_BA$RAD50)
shapiro.test(d625_BA$FoG50)


# BA planktonic & penetration @ 0.78
par(mfrow=c(2, 2), mar=c(1, 1, 1,1))
plot(d78_BA$Biomass, d78_BA$RAD50)
cor.test(d78_BA$Biomass, d78_BA$RAD50, method = "spearman")
abline(lm(d78_BA$RAD50~d78_BA$Biomass), lty=2)

plot(d78_BA$Activity, d78_BA$RAD50)
cor.test(d78_BA$Activity, d78_BA$RAD50, method = "spearman")
abline(lm(d78_BA$RAD50~d78_BA$Activity), lty=2)

plot(d78_BA$Biomass, d78_BA$FoG50)
cor.test(d78_BA$Biomass, d78_BA$FoG50, method = "spearman")
abline(lm(d78_BA$FoG50~d78_BA$Biomass), lty=2)

plot(d78_BA$Activity, d78_BA$FoG50)
cor.test(d78_BA$Activity, d78_BA$FoG50, method = "spearman")
abline(lm(d78_BA$FoG50~d78_BA$Activity), lty=2)

# BA planktonic & penetration @ 3.12
par(mfrow=c(2, 2), mar=c(1, 1, 1,1))
plot(d312_BA$Biomass, d312_BA$RAD50)
cor.test(d312_BA$Biomass, d312_BA$RAD50, method = "spearman")
abline(lm(d312_BA$RAD50~d312_BA$Biomass), lty=2)

plot(d312_BA$Activity, d312_BA$RAD50)
cor.test(d312_BA$Activity, d312_BA$RAD50, method = "spearman")
abline(lm(d312_BA$RAD50~d312_BA$Activity), lty=2)

plot(d312_BA$Biomass, d312_BA$FoG50)
cor.test(d312_BA$Biomass, d312_BA$FoG50, method = "spearman")
abline(lm(d312_BA$FoG50~d312_BA$Biomass), lty=2)

plot(d312_BA$Activity, d312_BA$FoG50)
cor.test(d312_BA$Activity, d312_BA$FoG50, method = "spearman")
abline(lm(d312_BA$FoG50~d312_BA$Activity), lty=2)

# BA planktonic & penetration @ 6.25
par(mfrow=c(2, 2), mar=c(2, 2, 1,1))
plot(d625_BA$Biomass, d625_BA$RAD50)
cor.test(d625_BA$Biomass, d625_BA$RAD50, method = "spearman")
abline(lm(d625_BA$RAD50~d625_BA$Biomass), lty=2)

plot(d625_BA$Activity, d625_BA$RAD50)
cor.test(d625_BA$Activity, d625_BA$RAD50, method = "spearman")
abline(lm(d625_BA$RAD50~d625_BA$Activity), lty=2)

plot(d625_BA$Biomass, d625_BA$FoG50)
cor.test(d625_BA$Biomass, d625_BA$FoG50, method = "spearman")
abline(lm(d625_BA$FoG50~d625_BA$Biomass), lty=2)

plot(d625_BA$Activity, d625_BA$FoG50)
cor.test(d625_BA$Activity, d625_BA$FoG50, method = "spearman")
abline(lm(d625_BA$FoG50~d625_BA$Activity))


# BA planktonic & penetration @ 12.50
par(mfrow=c(2, 2), mar=c(2, 2, 1,1))
plot(d1250_BA$Biomass, d1250_BA$RAD50)
cor.test(d1250_BA$Biomass, d1250_BA$RAD50, method = "spearman")
abline(lm(d1250_BA$RAD50~d1250_BA$Biomass), lty=2)

plot(d1250_BA$Activity, d1250_BA$RAD50)
cor.test(d1250_BA$Activity, d1250_BA$RAD50, method = "spearman")
abline(lm(d1250_BA$RAD50~d1250_BA$Activity), lty=2)

plot(d1250_BA$Biomass, d1250_BA$FoG50)
cor.test(d1250_BA$Biomass, d1250_BA$FoG50, method = "spearman")
abline(lm(d1250_BA$FoG50~d1250_BA$Biomass), lty=2)

plot(d1250_BA$Activity, d1250_BA$FoG50)
cor.test(d1250_BA$Activity, d1250_BA$FoG50, method = "spearman")
abline(lm(d1250_BA$FoG50~d1250_BA$Activity))

cor.test(d20_BA$Activity, d20_BA$FoG50, method = "spearman")
cor.test(d39_BA$Activity, d39_BA$FoG50, method = "spearman")