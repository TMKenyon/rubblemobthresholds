# 3.1  Mobilisation in flume

# For Kenyon et al. (2023) Biogeosciences (in review)

# Libraries ----

library(plyr)
library(tidyverse)
library(dplyr)
library(broom)
library(emmeans)
library(MuMIn)
library(gridExtra)
library(MASS)
library(DHARMa)

# Data ----

# Dataset with inertia calculations included in response to comments from Reviewer 2 Biogeosciences submission

datnew <- read.csv("/Users/taniakenyon/Dropbox/PhD/Data/Flume_Hydro_Lab/Wave_Flume/Flume_Analysis/CleanMovement_Both_Collated_01042020-nored2-InertiaChecks-git.csv")

names(datnew)
range(datnew$axial_length_mm)

dat2new <- datnew %>% 
  dplyr::select(locking, substrate, size, unique, axial_length_mm, branched, wet_weight,
                wave_period, wave_height_m, water_depth_m, E_value, calculated.velocity_corrected, linear_wave_theory_velocity,
                KC, FItoFD, X24overKC2, InertiaSig, MaximumFrom,
                rock, walk_slide, flip, transport_after_slide_or_flip,
                movement_at_all, transport, greatest_movement_type.no_transport,
                greatest_movement_type.no_transport_walk_slide_combined,
                greatest_movement_type.transport, Waves, Remove_for_flip_analysis_but_not_movement_analysis)

dat2new <- dat2new %>%  filter(locking == "Free", 
                               Waves == "Regular", 
                               size != "Boulder", 
                               branched != "Plate/Massive") %>% droplevels()

dat2new <- dat2new %>% 
  mutate(size=factor(size,
                     levels = c('Small', 'Medium','Large', 'XL', 'XXL'),
                     labels = c('4-8cm', '9-15cm', '16-23cm', "24-29cm", "30-39cm")),
         branched = factor(branched, levels = c("Unbranched", "Branched"))) 


datcombnew <- transform(dat2new, size=revalue(size,c("30-39cm" = "24-29cm"))) # Combining the 2 largest groups into a category 24-39 cm.

datcombnew <- datcombnew %>% mutate(size=factor(size,
                                                levels = c('4-8cm', '9-15cm', '16-23cm', '24-29cm'),
                                                labels = c('4-8 cm', '9-15 cm', '16-23 cm', '24-39 cm')))

toremove <- datcombnew %>% group_by(calculated.velocity_corrected, size, substrate) %>% 
  summarise(movement_at_all = n()) %>% as.data.frame()

toremove %>% summarise(minmove = min(movement_at_all),
                       meanmove = mean(movement_at_all),
                       maxmove = max(movement_at_all)) # min of 3 replicates (ensuring there are enough replicate runs for each velocity, size and substrate combination)

toremove2 <- datcombnew %>% group_by(calculated.velocity_corrected, size, branched) %>% 
  summarise(movement_at_all = n()) %>% as.data.frame()

toremove2 %>% summarise(minmove = min(movement_at_all),
                        meanmove = mean(movement_at_all),
                        maxmove = max(movement_at_all)) # min of 3 replicates

datcombnew$substrate <- as.factor(datcombnew$substrate)
datcombnew$InertiaSig <- as.factor(datcombnew$InertiaSig)
datcombnew$MaximumFrom <- as.factor(datcombnew$MaximumFrom)
datcombnew$wave_period <- as.factor(datcombnew$wave_period)

dim(datcombnew)

datcombnew %>% group_by(InertiaSig) %>% summarise(no_rows = length(InertiaSig))


# 3.1.1 Mobilisation Thresholds

# OVERALL THRESHOLDS ----

# Table S3 ----

# transport thresholds with just velocity as explanatory
# To compare to field for thresholds, we want to consider only pieces 4-23 cm in size.

datcombnewNOXL <- datcombnew %>% filter(size != "24-39 cm")

fl.t.thresh <- glm(transport ~ 
                     calculated.velocity_corrected,
                   data = datcombnewNOXL,
                   family = binomial(link='logit'))

car::Anova(fl.t.thresh)
write.csv(car::Anova(fl.t.thresh), file = "car::Anova(fl.t.thresh).csv")
r.squaredGLMM(fl.t.thresh)

plot(simulateResiduals(fittedModel = fl.t.thresh)) # good

ld <- dose.p(fl.t.thresh, p=c(0.1,0.5,0.9))

ld.SE = attr(ld, "SE")
ld = data.frame(LD = attr(ld, 'p'),
                Dose = as.vector(ld),
                SE = ld.SE) %>% 
  mutate(lower=Dose - SE*qnorm(0.975),
         upper=Dose + SE*qnorm(0.975)) #95% CI are ~1.96 times the SE.

ld
write.csv(ld, file = "fl.t.thresh.thresholds.csv")

# Table S4 ----

# flipping thresholds with just velocity as explanatory

fl.f.threshNOXL <- glm(flip ~ 
                         calculated.velocity_corrected,
                       data = datcombnewNOXL,
                       family = binomial(link='logit'))

car::Anova(fl.f.threshNOXL)

write.csv(car::Anova(fl.f.threshNOXL), file = "car::Anova(fl.f.threshNOXL).csv")
r.squaredGLMM(fl.f.threshNOXL)

plot(simulateResiduals(fittedModel = fl.f.threshNOXL))

ld <- dose.p(fl.f.threshNOXL, p=c(0.1,0.5,0.9))

ld.SE = attr(ld, "SE")
ld = data.frame(LD = attr(ld, 'p'),
                Dose = as.vector(ld),
                SE = ld.SE) %>% 
  mutate(lower=Dose - SE*qnorm(0.975),
         upper=Dose + SE*qnorm(0.975)) #95% CI are round about 1.96 times the SE.

ld
write.csv(ld, file = "thresholds_for_fl.f.threshNOXL.csv")

# 3.1.2	Rubble and substrate effects on mobilsation

# ROCKING -----

summary(datcombnew$greatest_movement_type.no_transport_walk_slide_combined)

rockdat <- datcombnew %>% filter(greatest_movement_type.no_transport_walk_slide_combined != "Flip") %>% droplevels()
rockdat <- rockdat %>% filter(greatest_movement_type.no_transport_walk_slide_combined != "Walk_Slide") %>% droplevels()

rockdat$greatest_movement_type.no_transport_walk_slide_combined <- as.factor(rockdat$greatest_movement_type.no_transport_walk_slide_combined)

summary(rockdat$greatest_movement_type.no_transport_walk_slide_combined)

names(rockdat)

# Table S5 ----

fl.r.global <- glm(rock ~ 
                     calculated.velocity_corrected +
                     size +
                     substrate +
                     branched +
                     calculated.velocity_corrected:size +
                     calculated.velocity_corrected:substrate +
                     calculated.velocity_corrected:branched +
                     size:substrate +
                     size:branched +
                     substrate:branched +
                     calculated.velocity_corrected:size:branched +
                     calculated.velocity_corrected:size:substrate,
                   data = rockdat,
                   family = binomial(link='logit'))

car::Anova(fl.r.global)

write.csv(car::Anova(fl.r.global), file = "car::Anova(fl.r.global).csv")

# validation

1-pchisq(fl.r.global$deviance, fl.r.global$null)

plot(simulateResiduals(fittedModel = fl.r.global))

r.squaredGLMM(fl.r.global)

# Diagnostics
plot(fl.r.global) 

# Rocking thresholds (again, removing XL size)

rockdatNOXL <- rockdat %>% filter(size != "24-39 cm")

fl.r.thresh <- glm(rock ~ 
                     calculated.velocity_corrected,
                   data = rockdatNOXL,
                   family = binomial(link='logit'))
write.csv(car::Anova(fl.r.thresh), file = "car::Anova(fl.r.thresh).csv")

plot(simulateResiduals(fittedModel = fl.r.thresh))

ld <- dose.p(fl.r.thresh, p=c(0.1,0.5,0.9))

ld.SE = attr(ld, "SE")
ld = data.frame(LD = attr(ld, 'p'),
                Dose = as.vector(ld),
                SE = ld.SE) %>% 
  mutate(lower=Dose - SE*qnorm(0.975),
         upper=Dose + SE*qnorm(0.975)) # 95% CI are ~1.96 times the SE.

ld
write.csv(ld, file = "thresholds_for_fl.r.threshNOXL.csv")
r.squaredGLMM(fl.r.thresh) #20%

# results

car::Anova(fl.r.global)

write.csv(car::Anova(fl.r.global), file = "car::Anova(fl.r.global).csv")

summary(fl.r.global)

# Table S6 ----

fl.r.global.range <- with(rockdat, 
                          list(calculated.velocity_corrected = c(0.1, 0.2, 0.3, 0.4),
                               branched = levels(branched),
                               size = levels(size),
                               substrate = levels(substrate)))

fl.r.global.interac2 <- emmeans(fl.r.global, pairwise ~ branched | size | calculated.velocity_corrected, 
                                at = fl.r.global.range, type = "response")

fl.r.global.interac2

write.csv(fl.r.global.interac2$contrasts, file = "fl.r.global.interac2-2.csv")

# Table S7 -----

fl.r.global.interac3 <- emmeans(fl.r.global, pairwise ~ substrate | size |calculated.velocity_corrected,
                                at = fl.r.global.range, type = "response")

fl.r.global.interac3

write.csv(fl.r.global.interac3$contrasts, file = "fl.r.global.interac3.csv")

# Table S8 -----

fl.r.global.interac1 <- emmeans(fl.r.global, pairwise ~ size | branched | calculated.velocity_corrected, 
                                at = fl.r.global.range, type = "response") %>% as.data.frame()

fl.r.global.interac1
write.csv(fl.r.global.interac1, file = "fl.r.global.interac1.csv")

# Figure 3a ----

fl.r.global.grid <- with(datcombnew, list(calculated.velocity_corrected = seq(min(calculated.velocity_corrected), max(calculated.velocity_corrected ), len = 100),
                                          size = levels(size),
                                          branched = levels(branched),
                                          substrate = levels(substrate)))

fl.r.global.preds <- emmeans(fl.r.global, ~ calculated.velocity_corrected | branched | size | substrate,
                             at = fl.r.global.grid, type = "response") %>% as.data.frame()
head(fl.r.global.preds)

cols <- c("4-8 cm" = "#FFC1C1", "9-15 cm" = "#FFEC8B", "16-23 cm" = "#87CEEB", "24-39 cm" = "#66CDAA")
fils <- c("4-8 cm" = "#FFC1C1", "9-15 cm" = "#FFEC8B", "16-23 cm" = "#87CEEB", "24-39 cm" = "#66CDAA")
lins <- c("4-8 cm" = 3, "9-15 cm" = 2, "16-23 cm" = 1, "24-39 cm" = 4)

fl.r.global.plot <- fl.r.global.preds %>%
  ggplot(aes(x = calculated.velocity_corrected, y = prob, colour = size)) +
  geom_ribbon(aes(ymin=asymp.LCL, ymax=asymp.UCL, fill=size), alpha =0.3) +
  geom_line(aes(color = size, linetype = size)) +
  theme_classic() + 
  facet_grid(substrate~branched, scales = "free") +
  scale_color_manual(name = "Rubble length", values= cols) +
  scale_fill_manual(name = "Rubble length", values= fils) +
  scale_linetype_manual(name = "Rubble length", values= lins) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), 
                     limits = c(0,1)) +
  scale_x_continuous(limits = c(0,0.43)) +
  xlab(bquote('Bottom orbital velocity ('*m/s*')')) +
  ylab(bquote('Probability of rocking')) +
  theme_classic() +
  theme(text = element_text(size=13, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position=c(0.6,0.865))

fl.r.global.plot

# TRANSPORT -----

# Table S9 -----

fl.t.global <- glm(transport ~ 
                     calculated.velocity_corrected +
                     size +
                     substrate +
                     branched +
                     calculated.velocity_corrected:size +
                     calculated.velocity_corrected:substrate +
                     calculated.velocity_corrected:branched +
                     size:substrate +
                     size:branched +
                     substrate:branched +
                     calculated.velocity_corrected:size:branched +
                     calculated.velocity_corrected:size:substrate,
                   data = datcombnew,
                   family = binomial(link='logit'))

car::Anova(fl.t.global)
r.squaredGLMM(fl.t.global)

plot(simulateResiduals(fittedModel = fl.t.global))

# Table S10 ----

fl.t.global.range <- with(datcombnew, 
                          list(calculated.velocity_corrected = c(0.1, mean(calculated.velocity_corrected), 0.3, 0.4),
                               branched = levels(branched),
                               size = levels(size),
                               substrate = levels(substrate)))

fl.t.global.interac2 <- emmeans(fl.t.global, pairwise ~ branched | size | calculated.velocity_corrected, 
                                at = fl.t.global.range, type = "response")
write.csv(fl.t.global.interac2, file = "fl.t.global.interac2.csv")
fl.t.global.interac2

# Table S11 ----


fl.t.global.interac1 <- emmeans(fl.t.global, pairwise ~ size | branched | calculated.velocity_corrected, 
                               at = fl.t.global.range, type = "response")
write.csv(fl.t.global.interac1, file = "fl.t.global.interac1.csv")
fl.t.global.interac1

# Table S12 -----

fl.t.global.interac4 <- emmeans(fl.t.global, pairwise ~ substrate | size | calculated.velocity_corrected, 
                                at = fl.t.global.range, type = "response")
write.csv(fl.t.global.interac4, file = "fl.t.global.interac4.csv")
fl.t.global.interac4

# Figure 3b ----

fl.t.global.grid2 <- with(datcombnew, list(calculated.velocity_corrected = seq(min(calculated.velocity_corrected), 0.6, len = 200),
                                           size = levels(size),
                                           branched = levels(branched),
                                           substrate = levels(substrate)))

fl.t.global.preds2 <- emmeans(fl.t.global, ~ calculated.velocity_corrected | branched | size | substrate,
                              at = fl.t.global.grid2, type = "response") %>% as.data.frame()

write.csv(fl.t.global.preds, file = "fl.t.global.preds-for-Table-in-Mob-Ch.csv")
write.csv(fl.t.global.preds2, file = "fl.t.global.preds2-for-Table-in-Mob-Ch.csv")
head(fl.t.global.preds)

cols <- c("4-8 cm" = "#FFC1C1", "9-15 cm" = "#FFEC8B", "16-23 cm" = "#87CEEB", "24-39 cm" = "#66CDAA")
fils <- c("4-8 cm" = "#FFC1C1", "9-15 cm" = "#FFEC8B", "16-23 cm" = "#87CEEB", "24-39 cm" = "#66CDAA")
lins <- c("4-8 cm" = 3, "9-15 cm" = 2, "16-23 cm" = 1, "24-39 cm" = 4)

fl.t.global.plot <- fl.t.global.preds %>%
  ggplot(aes(x = calculated.velocity_corrected, y = prob, colour = size)) +
  geom_ribbon(aes(ymin=asymp.LCL, ymax=asymp.UCL, fill=size), alpha =0.3) +
  geom_line(aes(color = size, linetype = size)) +
  theme_classic() + 
  facet_grid(substrate~branched, scales = "free") +
  scale_color_manual(name = "Rubble length", values= cols) +
  scale_fill_manual(name = "Rubble length", values= fils) +
  scale_linetype_manual(name = "Rubble length", values= lins) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), 
                     limits = c(0,1)) +
  scale_x_continuous(limits = c(0,0.43)) +
  xlab(bquote('Bottom orbital velocity ('*m/s*')')) +
  ylab(bquote('Probability of transport (walk/slide/flip)')) +
  theme_classic() +
  theme(text = element_text(size=13, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position=c(0.6,0.865))

fl.t.global.plot

# FLIPPING ----

# Table S13 ----

fl.f.global2 <- glm(flip ~ 
                      calculated.velocity_corrected +
                      size +
                      substrate +
                      branched +
                      calculated.velocity_corrected:size +
                      calculated.velocity_corrected:substrate +
                      calculated.velocity_corrected:branched +
                      size:substrate +
                      size:branched +
                      substrate:branched +
                      calculated.velocity_corrected:size:branched +
                      calculated.velocity_corrected:size:substrate,
                    data = datcombnew,
                    family = binomial(link='logit'))

# validation

r.squaredGLMM(fl.f.global2) #59%

plot(simulateResiduals(fittedModel = fl.f.global2))

# Model fit
fl.f.global.resid2 <- sum(resid(fl.f.global2, type = "pearson")^2)
fl.f.global.resid2
1-pchisq(fl.f.global.resid2, fl.f.global2$df.residual) # good model fit
1-pchisq(fl.f.global2$deviance, fl.f.global2$df.residual)  # good model fit

# Diagnostics
plot(fl.f.global2) 

# results

car::Anova(fl.f.global2)

summary(fl.f.global2)

# Table S14

# interactions

fl.f.global.range <- with(datcombnew, 
                          list(calculated.velocity_corrected = c(0.1, mean(calculated.velocity_corrected), 0.4),
                               branched = levels(branched),
                               size = levels(size),
                               substrate = levels(substrate)))
fl.f.global.range

fl.f.global.interac1.2 <- emmeans(fl.f.global2, pairwise ~ branched | size | calculated.velocity_corrected, 
                                  at = fl.f.global.range, type = "response")
fl.f.global.interac1.2

write.csv(fl.f.global.interac1.2, file = "fl.f.global.interac1.2.csv")

# Table S15 ----

fl.f.global.interac2.2 <- emmeans(fl.f.global2, pairwise ~ size | branched | calculated.velocity_corrected, 
                                  at = fl.f.global.range, type = "response")
fl.f.global.interac2.2

write.csv(fl.f.global.interac2.2, file = "fl.f.global.interac2.2.csv")


# Table S16 ----

fl.f.global.interac3.3 <- emmeans(fl.f.global2, pairwise ~ substrate | branched | calculated.velocity_corrected, 
                                  at = fl.f.global.range, type = "response")
fl.f.global.interac3.3

# Figure 3c ----

fl.f.global.preds2 <- emmeans(fl.f.global2, ~ calculated.velocity_corrected | branched | size | substrate,
                              at = fl.f.global.grid, type = "response") %>% as.data.frame()
head(fl.f.global.preds2)

cols <- c("4-8 cm" = "#FFC1C1", "9-15 cm" = "#FFEC8B", "16-23 cm" = "#87CEEB", "24-39 cm" = "#66CDAA")
fils <- c("4-8 cm" = "#FFC1C1", "9-15 cm" = "#FFEC8B", "16-23 cm" = "#87CEEB", "24-39 cm" = "#66CDAA")
lins <- c("4-8 cm" = 3, "9-15 cm" = 2, "16-23 cm" = 1, "24-39 cm" = 4)

fl.f.global.plot2 <- fl.f.global.preds2 %>%
  ggplot(aes(x = calculated.velocity_corrected, y = prob, colour = size)) +
  geom_ribbon(aes(ymin=asymp.LCL, ymax=asymp.UCL, fill=size), alpha =0.3) +
  geom_line(aes(color = size, linetype = size)) +
  theme_classic() + 
  facet_grid(substrate~branched, scales = "free") +
  scale_color_manual(name = "Rubble length", values= cols) +
  scale_fill_manual(name = "Rubble length", values= fils) +
  scale_linetype_manual(name = "Rubble length", values= lins) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), 
                     limits = c(0,1)) +
  scale_x_continuous(limits = c(0,0.43)) +
  xlab(bquote('Bottom orbital velocity ('*m/s*')')) +
  ylab(bquote('Probability of flipping')) +
  theme(plot.margin = unit(c(9, 9, 9, 9),"points"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=10,family = "sans"),
        legend.text = element_text(size=9,family = "sans"),
        legend.position=c(0.51,0.5), #0,0 is for bottom,left
        legend.justification = c(0.01,1), #1,1 is for top,right
        axis.text.x.bottom = element_text(size=11,vjust = -0.05,colour = "black"),
        axis.text.y.left = element_text(size=11,colour = "black"),
        axis.title.x=element_text(size=14,family = "sans", color="black", vjust = -1),  
        axis.title.y=element_text(size=14,family = "sans", color="black", vjust = 3),
        strip.text.x= element_text(size = 12, family = "sans"))

fl.f.global.plot2

fl.f.global.plot2 <- fl.f.global.preds2 %>%
  ggplot(aes(x = calculated.velocity_corrected, y = prob, colour = size)) +
  geom_ribbon(aes(ymin=asymp.LCL, ymax=asymp.UCL, fill=size), alpha =0.3) +
  geom_line(aes(color = size, linetype = size)) +
  theme_classic() + 
  facet_grid(substrate~branched, scales = "free") +
  scale_color_manual(name = "Rubble length", values= cols) +
  scale_fill_manual(name = "Rubble length", values= fils) +
  scale_linetype_manual(name = "Rubble length", values= lins) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), 
                     limits = c(0,1)) +
  scale_x_continuous(limits = c(0,0.43)) +
  xlab(bquote('Bottom orbital velocity ('*m/s*')')) +
  ylab(bquote('Probability of flipping')) +
  theme_classic() +
  theme(text = element_text(size=13, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position=c(0.6,0.865))


fl.f.global.plot2

# INTERLOCKED ----

dat <- read.csv("/Users/taniakenyon/Dropbox/PhD/Data/Flume_Hydro_Lab/Wave_Flume/Flume_Analysis/CleanMovement_Both_Collated_01042020-nored2-InertiaChecks-git.csv")

names(dat)

dat2 <- dat %>% 
  dplyr::select(locking, substrate, size, unique, axial_length_mm, branched, wet_weight,
                wave_period, wave_height_m, water_depth_m, E_value, calculated.velocity_corrected, 
                 rock, walk_slide, flip, transport_after_slide_or_flip,
                movement_at_all, greatest_movement_type.no_transport,
                greatest_movement_type.no_transport_walk_slide_combined,
                greatest_movement_type.transport, Waves, Remove_for_flip_analysis_but_not_movement_analysis)

dat3 <- dat2 %>% 
  mutate(size=factor(size,
                     levels = c('Small', 'Medium','Large', 'XL', 'XXL', 'Boulder'),
                     labels = c('4-8cm', '9-15cm', '16-23cm', "24-29cm", "30-39cm", "Boulder")))

datI <- dat3 %>%  filter(locking == "Interlocking", # looking at interlocking ONLY
                         Waves == "Regular") %>% droplevels()

levels(datI$size) #So, I am only looked at interlocking of small and medium sized pieces, no others were trialed.
# I could only interlock small, branched pieces, not larger or unbranched pieces.

# Fit model

lockmove <- glm(movement_at_all ~ calculated.velocity_corrected,
                family = binomial(link='logit'),
                data = datI)
summary(lockmove)  

plot(simulateResiduals(fittedModel = lockmove))

# No relationship between velocity and the probability of movement for interlocked rubble pieces

# Table S17 ----

# Fit model for any movement

lockmove2 <- glm(movement_at_all ~ calculated.velocity_corrected * size,
                 family = binomial(link='logit'),
                 data = datI)
summary(lockmove2) 

car::Anova(lockmove2)
write.csv(car::Anova(lockmove2), file = "car::Anova(lockmove2).csv")

plot(simulateResiduals(fittedModel = lockmove2))

# Fit model for flipping (second part of Table S16)

lockflip <- glm(flip ~ calculated.velocity_corrected * size,
                family = binomial(link='logit'),
                data = datI)
summary(lockflip) 

car::Anova(lockflip)

plot(simulateResiduals(fittedModel = lockflip))

write.csv(car::Anova(lockflip), file = "car::Anova(lockflip).csv")

# Comparisons of different kinds of movement:

# Add a column called 'moved' which places a 1 if the column "greatest_movement_type.no_transport_walk_slide_combined" is 'Flip', 'Rock' or 'Walk_Slide'

datInew <- datI %>% 
  mutate(Flip = ifelse(greatest_movement_type.no_transport_walk_slide_combined == "Flip", 1, 0),
         Rock = ifelse(greatest_movement_type.no_transport_walk_slide_combined == "Rock", 1, 0),
         Transport = ifelse(greatest_movement_type.no_transport_walk_slide_combined == "Walk_Slide", 1, 0)) %>% 
  pivot_longer(cols = c("Flip", "Rock", "Transport"),
               names_to = "MovementType",
               values_to = "Moved")

head(datInew)

datInew$MovementType <- as.factor(datInew$MovementType)

datInew <- datInew %>% mutate(MovementType = factor(MovementType, 
                                                    levels = c("Rock", 
                                                        "Transport", "Flip")))

levels(datInew$MovementType)

# use greatest_movement_type.no_transport_walk_slide_combined

probmove.m <- glm(Moved ~ calculated.velocity_corrected + MovementType + size,
                  family = binomial(link='logit'),
                  data = datInew)

car::Anova(probmove.m) # no interactions so removed them one by one

# and no effect of velocity so removed it

# Table S49 -----

probmove <- glm(Moved ~ MovementType,
                family = binomial(link='logit'),
                data = datInew)

car::Anova(probmove)

write.csv(car::Anova(probmove), "probmove.csv")

emmeans(probmove, pairwise ~ MovementType, type = "response")

# rock vs transport: z-ratio = 3.671, P < 0.001
# rock vs flip: z-ratio = -3.671, P < 0.001

# Diagnostics

simulateResiduals(probmove, plot = T)

move.range <- with(datInew, list(MovementType = levels(MovementType)))

probmove.preds <- emmeans(probmove, pairwise ~ MovementType, 
                          at = move.range, type = "response")

probmove.preds

# Figure S9 -----

probmove.preds2 <- emmeans(probmove, ~ MovementType,
                           at = move.range, type = "response") %>% as.data.frame()

cols <- c("Flip" = "red", "Rock" = "gold", "Transport" = "darkorange1")
fils <- c("Flip" = "red", "Rock" = "gold", "Transport" = "darkorange1")
lins <- c("Flip" = 2, "Rock" = 1, "Transport" = 3)

pd <- position_dodge(0.5)

probmove.preds2.plot <- probmove.preds2 %>%
  ggplot(aes(x = MovementType, y = prob, colour = MovementType)) +
  geom_point(position=pd, size =3) + 
  geom_errorbar(aes(ymin=prob-SE, ymax=prob+SE), width=.2, position = pd) +
  #geom_pointrange(aes(ymin=meanGrowth-seGrowth, ymax=meanGrowth+seGrowth), width=.3, position = pd), data=dat.sum)+ 
  scale_color_manual(name = "Movement type", values= cols) +
  scale_fill_manual(name = "Movement type", values= fils) +
  scale_linetype_manual(name = "Movement type", values= lins) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), 
                     limits = c(0,1)) + 
  xlab(bquote('Movement type')) +
  ylab(bquote('Probability of movement')) +
  theme_classic()+
  theme(text = element_text(size=13, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position=c(0.85,0.8)) +
  ylim(0,0.2)

probmove.preds2.plot


# Figure S2 -----

# Comparison of Linear wave theory vs Soulsbys Cosine Approximation

library(tidyverse)

lsdat <- read.csv("/Users/taniakenyon/Dropbox/Manuscripts/In_Review/Ch.3-Kenyon_et_al-Rubble-Mobilisation/Biogeosciences/Data_Sets_and_Code-Github/linear_soulsbys.csv")
names(lsdat)

ggplot(lsdat, aes(x = Soulsbys_Cosine_Approximation, y = Linear_wave_theory)) + 
  geom_point() + stat_smooth(method = "lm", se= TRUE, level = 0.95, formula = y ~ x) +
  theme_classic() +
  xlab("Soulsby Cosine Approximation (m/s)") +
  ylab("Linear wave theory (m/s)") +
  ggtitle("Bottom orbital velocity as estimated from wave height and period using two methods")

reg1 <- lm(Linear_wave_theory ~ Soulsbys_Cosine_Approximation, data = lsdat)

summary(reg1) # R squared for this relationship is 0.999

rangeS <- with(lsdat, list(Soulsbys_Cosine_Approximation = seq(0, 0.45, len = 100)))

preds <- emmeans(reg1, ~ Soulsbys_Cosine_Approximation, 
                 at = rangeS, type = "response") %>% as.data.frame()

plot <- ggplot(data = preds, aes(x = Soulsbys_Cosine_Approximation, y = emmean)) +
  geom_line(colour = "red") +
  geom_ribbon(aes(ymin=lower.CL, ymax=upper.CL, fill = "red")) +
  geom_point(data = lsdat, aes(x = Soulsbys_Cosine_Approximation, y = Linear_wave_theory)) +
  theme_classic() +
  xlab("Soulsby Cosine Approximation (m/s)") +
  ylab("Linear wave theory (m/s)")
plot

# Figure S5 ----

# Plot of velocity against FI/FD

datnew <- read.csv("/Users/taniakenyon/Dropbox/PhD/Data/Flume_Hydro_Lab/Wave_Flume/Flume_Analysis/CleanMovement_Both_Collated_01042020-nored2-InertiaChecks-git.csv")

dim(datnew)

names(datnew)
range(datnew$axial_length_mm)

dat2new <- datnew %>% 
  dplyr::select(Run, locking, substrate, size, unique, axial_length_mm, branched, wet_weight,
                wave_period, wave_height_m, water_depth_m, E_value, calculated.velocity_corrected, linear_wave_theory_velocity,
                KC, FItoFD, X24overKC2, InertiaSig, MaximumFrom,
                 rock, walk_slide, flip, transport_after_slide_or_flip,
                movement_at_all, transport, greatest_movement_type.no_transport,
                greatest_movement_type.no_transport_walk_slide_combined,
                greatest_movement_type.transport, Waves, Remove_for_flip_analysis_but_not_movement_analysis)

dat2newFI <- dat2new %>%  filter(Waves != "Jonswap_Waves",
                                 size != "Boulder",
                                 branched != "Plate/Massive",
                                 locking != "Not-interlocking_on_steel_mesh_base") %>% droplevels()  # considering Free AND Interlocking for this investigation

dim(dat2newFI)

dat2newFI <- dat2newFI %>% 
  mutate(size=factor(size,
                     levels = c('Small', 'Medium','Large', 'XL', 'XXL'),
                     labels = c('4-8cm', '9-15cm', '16-23cm', "24-29cm", "30-39cm")),
         branched = factor(branched, levels = c("Unbranched", "Branched"))) 


datcombnewFI <- transform(dat2newFI, size=revalue(size,c("30-39cm" = "24-29cm"))) 

datcombnewFI <- datcombnewFI %>% mutate(size=factor(size,
                                                    levels = c('4-8cm', '9-15cm', '16-23cm', '24-29cm'),
                                                    labels = c('4-8 cm', '9-15 cm', '16-23 cm', '24-39 cm')))

dim(datcombnewFI)

datcombnewFI$substrate <- as.factor(datcombnewFI$substrate)
datcombnewFI$InertiaSig <- as.factor(datcombnewFI$InertiaSig)
datcombnewFI$MaximumFrom <- as.factor(datcombnewFI$MaximumFrom)
datcombnewFI$movement_at_all <- as.factor(datcombnewFI$movement_at_all)
datcombnewFI$Run <- as.factor(datcombnewFI$Run)
datcombnewFI$locking <- as.factor(datcombnewFI$locking)

datcombnewFI %>% group_by(InertiaSig) %>% summarise(no_rows = length(InertiaSig))
datcombnewFI %>% group_by(InertiaSig, greatest_movement_type.transport) %>% summarise(no_rows = length(InertiaSig))

datcombnewFI$movement_at_all <- as.factor(datcombnewFI$movement_at_all)

plot <- ggplot(data = datcombnewFI, aes(x = calculated.velocity_corrected, y = FItoFD, colour = movement_at_all)) +
  geom_point() +
  theme_classic() +
  xlab("Bottom orbital velocity (m/s)") +
  ylab("FI / FD ratio") +
  labs(colour = "Movement (1) or not (0)") +
  theme(legend.position = "bottom")
plot

# Figure S6 ----

# Another plot as above but only for inertia-significant runs, and colour by greatest_movement

datcombnewFI3 <- datcombnewFI %>%  filter(InertiaSig !=  "No",
                                          movement_at_all != "0") %>% droplevels()  # considering only Inertia sig, movement

names(datcombnewFI3)

dim(datcombnewFI3)

ggplot(data = datcombnewFI3, aes(x = calculated.velocity_corrected, y = FItoFD, colour = greatest_movement_type.transport)) +
  geom_point() +
  theme_classic() +
  xlab("Bottom orbital velocity (m/s)") +
  ylab("FI / FD ratio") +
  labs(colour = "Movement type") +
  theme(legend.position = "bottom")

# Filter for FI / FD <2 total of above dataset, then plot the inertial component

datcombnewFI3.3 <- datcombnewFI3 %>%  filter(FItoFD < 2)

ggplot(data = datcombnewFI3.3, aes(x = calculated.velocity_corrected, y = X24overKC2, colour = greatest_movement_type.transport)) +
  geom_point() +
  theme_classic() +
  xlab("Bottom orbital velocity (m/s)") +
  ylab("Inertia force as a proportion of drag force") +
  labs(colour = "Movement (1) or not (0)") +
  theme(legend.position = "bottom")

# Filter for FI / FD <2 total dataset, then plot the inertial component

datcombnewFI3.2 <- datcombnewFI %>%  filter(FItoFD < 2)

ggplot(data = datcombnewFI3.2, aes(x = calculated.velocity_corrected, y = X24overKC2, colour = movement_at_all)) +
  geom_point() +
  theme_classic() +
  xlab("Bottom orbital velocity (m/s)") +
  ylab("Inertia force as a proportion of drag force") +
  labs(colour = "Movement (1) or not (0)") +
  theme(legend.position = "bottom")

ggplot(data = datcombnewFI3.2, aes(x = calculated.velocity_corrected, y = X24overKC2, colour = greatest_movement_type.transport)) +
  geom_point() +
  theme_classic() +
  xlab("Bottom orbital velocity (m/s)") +
  ylab("Inertia force as a proportion of drag force") +
  labs(colour = "Movement (1) or not (0)") +
  theme(legend.position = "bottom")

# Filter further for Figure AA2 in Attachment A of response to comments, of just the 18 cases identified as having movement (not just rocking) under conditions that had the potential to be inertia-dominant

datcombnewFI2 <- datcombnewFI %>%  filter(InertiaSig !=  "No",
                                          movement_at_all != "0",
                                          flip != "0") %>% droplevels()  # considering only Inertia sig, movement

dim(datcombnewFI2) # This is the 18 cases shown in Table AA2 of Attachment A

(plot <- ggplot(data = datcombnewFI2, aes(x = Run, y = X24overKC2, colour = size, fill = size)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    xlab("Case number") +
    ylab("Inertia force as a proportion of drag force")) +
  ylim(0,1.0) +
  labs(colour = "Rubble length", fill = "Rubble length")

# Figure S7 ----

ggplot(data = datcombnewFI2, aes(x = calculated.velocity_corrected, y = X24overKC2, colour = size, fill = size)) +
  geom_point() +
  theme_classic() +
  xlab("Bottom orbital velocity (m/s)") +
  ylab("Inertia force as a proportion of drag force") +
  ylim(0,1.0) +
  labs(colour = "Rubble length", fill = "Rubble length") + 
  theme(legend.position = "bottom")


# Figure S8 ----

datcombnewFI3.2 <- datcombnewFI %>%  filter(FItoFD < 2)

ggplot(data = datcombnewFI3.2, aes(x = calculated.velocity_corrected, y = X24overKC2, colour = movement_at_all)) +
  geom_point() +
  theme_classic() +
  xlab("Bottom orbital velocity (m/s)") +
  ylab("Inertia force as a proportion of drag force") +
  labs(colour = "Movement (1) or not (0)") +
  theme(legend.position = "bottom")

