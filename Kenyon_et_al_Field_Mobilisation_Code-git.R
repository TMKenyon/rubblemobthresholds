# Mobilisation in the field

# For Kenyon et al. (2023) Biogeosciences (in review)

# Libraries ----

library(plyr)
library(ggplot2)
library(tidyverse)
library(car) #for scatterplotMatrix
library(nlme)
library(lme4) #for glmer
library(glmmTMB) #for glmmTMB
library(emmeans) #for emmeans
library(MuMIn) #for AICc
library(gridExtra)
library(DHARMa)

# An assessment and analysis of the wave orbital velocities measured in the field

# Data ----

# This dataset has been corrected in terms of how the wave number has been calculated, in response to comment from Reviewer 2 on submission to Biogeosciences

overall <- read.csv("~/Dropbox/PhD/Data/Maldives/Rubble_Movement/Shortterm_Movement/Analysis/Rubble_Transport_Shortterm_COL_v3-Sheet4-corrected_u-git.csv")

overall$aspect<-factor(overall$aspect)
overall$site<-factor(overall$site)
overall$siteUn <- factor(overall$siteUn)
overall$depth<-factor(overall$depth)
overall$year<-factor(overall$year)
overall$trial_id_each_depth<-factor(overall$trial_id_each_depth)
overall$trial_id_combined_depths<-factor(overall$trial_id_combined_depths)
overall$season<-factor(overall$season)
overall$old_id<-factor(overall$old_id)
overall$new_id<-factor(overall$new_id)
overall$trial<-factor(overall$trial)
overall$branched<-factor(overall$branched)
overall$size_cat_adjusted_to_flume_exp<-factor(overall$size_cat_adjusted_to_flume_exp)
overall$branches_opposing<-factor(overall$branches_opposing)
overall$day<-factor(overall$day)
overall$starting_substrate_new<-factor(overall$starting_substrate_new)
overall$substrate_stability<-factor(overall$substrate_stability)
overall$buoy_wt<-as.numeric(overall$buoy_wt)
overall$branch3_l<-as.numeric(overall$branch3_l)
overall$dry_wt<-as.numeric(overall$dry_wt)
overall$wet_wt<-as.numeric(overall$wet_wt)
overall$vol<-as.numeric(overall$vol)
overall$flipped_on_this_day_tidying <- as.character(overall$flipped_on_this_day_tidying)

# reorder levels 

levels(overall$aspect)
levels(overall$depth)
levels(overall$size_cat_adjusted_to_flume_exp)
levels(overall$branched)
levels(overall$starting_substrate_new)

names(overall)

overall = overall %>% 
  mutate(aspect=factor(aspect, levels = c('Lagoon', 'SE', 'West')), #re-leveling
         depth=factor(depth, levels = c('2m', '2-3m','6-7m')),
         size_cat_adjusted_to_flume_exp=factor(size_cat_adjusted_to_flume_exp,
                                               levels = c('Small', 'Medium','Large '),
                                               labels = c('4-8 cm', '9-15 cm', '16-23 cm')),
         CORRECTED_fastest_peak_U_cm = CORRECTED_fastest_peak_U * 100,
         pt_avg_u_cm = pt_avg_u * 100,
         CORRECTED_fastest_peak_U_cm = CORRECTED_fastest_peak_U * 100,
         pt_daily_avg_peak_cm = pt_daily_avg_peak * 100,
         current_peak_u_cm = current_peak_u * 100,
         current_avg_u_cm = current_avg_u * 100,
         current_daily_avg_peak_cm = current_daily_avg_peak * 100,
         branched = factor(branched, levels = c('Yes', 'No'),
                           labels = c('Branches', 'No_branches')))

overall <- transform(overall, 
                     starting_substrate_new=revalue(starting_substrate_new,
                                                    c("rubble " = "rubble", 
                                                      "sand " = "sand")))
levels(overall$starting_substrate_new)

#Make a new variable (column) called depthv2 that is a direct copy of 'depth' variable.
#Then, rename the Reef flat/= (Lagoon) to 2-3m, not 2m.

overall <- overall %>% mutate(depth2 = depth)
levels(overall$depth2)

# recoding/revalueing

overall <- transform(overall, depth2=revalue(depth2,c("2m"="2-3m")))
levels(overall$depth)
levels(overall$depth2)

# Do the same for day

overall <- overall %>% mutate(day2 = day)
overall <- transform(overall, day2=revalue(day2,c("1"="One", "2"="Two", "3"="Three")))
levels(overall$day2)

# DepthHabitat

# Create a new variable, combining two columns into one.

overall$aspectDepth <- paste(overall$aspect, overall$depth2, sep="_")

overall$aspectDepth <- factor(overall$aspectDepth,
                              levels = c('Lagoon_2-3m', 'SE_2-3m','SE_6-7m',
                                         'West_2-3m','West_6-7m'), 
                              labels = c('Lag_Shal', 'Shelt_Shal', 
                                         'Shelt_Deep', 'Exp_Shal', 'Exp_Deep'))
levels(overall$aspectDepth)

overall$season <- factor(overall$season,
                         levels = c('NE_Monsoon','W_Monsoon'), 
                         labels = c('NE Monsoon', 'W Monsoon'))

# Figure 4 ----

fig4maxCORRECTED <- ggplot() +
  geom_boxplot(data=overall,aes(x = aspectDepth,
                                y = CORRECTED_fastest_peak_U, fill = season),  size=.5,
               position=position_dodge(.83)) + 
  scale_fill_manual(name = "Season", values= c("#6495ED","#EE6363")) +
  ylab(bquote('Peak near-bed wave orbital velocity ('*m/s*')')) +
  scale_x_discrete(name = "Habitat and depth") +
  scale_y_continuous(breaks=c(0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45),  limits = c(0.0,0.45)) +
  theme_classic() +
  theme(plot.margin = unit(c(9, 9, 9, 9),"points"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=12,family = "sans"),
        legend.text = element_text(size = 10,family = "sans"),
        legend.position=c(0.05,1), 
        legend.justification = c(0,1.1),
        axis.text.x.bottom = element_text(size=12,
                                          vjust = -0.05),
        axis.text.y.left = element_text(size=12),
        axis.title.x=element_text(size=14,family = "sans",
                                  color="#000000",vjust = -1),  
        axis.title.y=element_text(size=14,family = "sans",
                                  color="#000000", vjust = 3))

fig4maxCORRECTED

# Figure S2 ----

names(overall)

fig4avgCORRECTED <- ggplot() +
  geom_boxplot(data=overall,aes(x = aspectDepth,
                                y = CORRECTED_avg_peak_U, fill = season),  size=.5,
               position=position_dodge(.83)) + 
  scale_fill_manual(name = "Season", values= c("#6495ED","#EE6363")) +
  ylab(bquote('Avg near-bed wave orbital velocity ('*m/s*')')) +
  scale_x_discrete(name = "Habitat and depth") +
  scale_y_continuous(breaks=c(0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45),  limits = c(0.0,0.45)) +
  theme_classic() +
  theme(plot.margin = unit(c(9, 9, 9, 9),"points"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=12,family = "sans"),
        legend.text = element_text(size = 10,family = "sans"),
        legend.position=c(0.05,1), 
        legend.justification = c(0,1.1),
        axis.text.x.bottom = element_text(size=12,
                                          vjust = -0.05),
        axis.text.y.left = element_text(size=12),
        axis.title.x=element_text(size=14,family = "sans",
                                  color="#000000",vjust = -1),  
        axis.title.y=element_text(size=14,family = "sans",
                                  color="#000000", vjust = 3))

fig4avgCORRECTED

grid.arrange(fig4max, fig4maxCORRECTED, ncol = 2)


# Table S18 ----

names(overall)
velsumm <- overall %>% dplyr::group_by(season, aspectDepth) %>%
  dplyr::summarize(minvel = min(CORRECTED_fastest_peak_U),
            medvel = median(CORRECTED_fastest_peak_U),
            maxvel = max(CORRECTED_fastest_peak_U))

velsumm


# Table S19 -----

gamma1 <- glmmTMB(CORRECTED_fastest_peak_U_cm ~ aspectDepth * season +
                    (1|trial_id_combined_depths/siteUn),
                  data = overall, family=Gamma (link ='log'))

write.csv(car::Anova(gamma1), file = "car::Anova(gamma1).csv")
car::Anova(gamma1)
r.squaredGLMM(gamma1) 

plot(simulateResiduals(fittedModel = gamma1))

# Table S20 -----

preds <- emmeans(gamma1, pairwise ~ season | aspectDepth, type = "response")
write.csv(preds$contrasts,
          file = "emmeans(gamma1, pairwise ~ season | aspectDepth).csv")

# Table S21 -----

preds <- emmeans(gamma1, pairwise ~ aspectDepth | season, type = "response")
write.csv(preds$contrasts,
          file = "emmeans(gamma1, pairwise ~ aspectDepth | season).csv")

# ALL DAYS 1-3 -----

# Transport in the Western Monsoon

# Table S20 ----

rm(list=ls())

overall <- read.csv("~/Dropbox/PhD/Data/Maldives/Rubble_Movement/Shortterm_Movement/Analysis/Rubble_Transport_Shortterm_COL_v3-Sheet4-corrected_u-git.csv")

overall$aspect<-factor(overall$aspect)
overall$site<-factor(overall$site)
overall$siteUn <- factor(overall$siteUn)
overall$depth<-factor(overall$depth)
overall$year<-factor(overall$year)
overall$trial_id_each_depth<-factor(overall$trial_id_each_depth)
overall$trial_id_combined_depths<-factor(overall$trial_id_combined_depths)
overall$season<-factor(overall$season)
overall$new_id<-factor(overall$new_id)
overall$trial<-factor(overall$trial)
overall$branched<-factor(overall$branched)
overall$size_cat_adjusted<-factor(overall$size_cat_adjusted)
overall$branches_opposing<-factor(overall$branches_opposing)
overall$day<-factor(overall$day)
overall$starting_substrate_new<-factor(overall$starting_substrate_new)
overall$substrate_stability<-factor(overall$substrate_stability)
overall$buoy_wt<-as.numeric(overall$buoy_wt)
overall$branch3_l<-as.numeric(overall$branch3_l)
overall$dry_wt<-as.numeric(overall$dry_wt)
overall$wet_wt<-as.numeric(overall$wet_wt)
overall$vol<-as.numeric(overall$vol)
overall$flipped_on_this_day_tidying <- as.character(overall$flipped_on_this_day_tidying)

overall <- transform(overall, 
                     starting_substrate_new=revalue(starting_substrate_new,
                                                    c("rubble " = "rubble", 
                                                      "sand " = "sand")))

overall = overall %>% 
  mutate(aspect=factor(aspect, levels = c('Lagoon', 'SE', 'West')), #releveling
         depth=factor(depth, levels = c('2m', '2-3m','6-7m')),
         size_cat_adjusted_to_flume_exp=factor(size_cat_adjusted_to_flume_exp,
                                               levels = c('Small', 'Medium','Large '), 
                                               labels = c('4-8 cm', '9-15 cm', '16-23 cm')),
         CORRECTED_fastest_peak_U_cm = CORRECTED_fastest_peak_U * 100,
         pt_avg_u_cm = pt_avg_u * 100,
         pt_daily_avg_peak_cm = pt_daily_avg_peak * 100,
         current_peak_u_cm = current_peak_u * 100,
         current_avg_u_cm = current_avg_u * 100,
         current_daily_avg_peak_cm = current_daily_avg_peak * 100,
         starting_substrate_new=factor(starting_substrate_new,
                                       levels = c('hard_carb', 'rubble', 'sand'), 
                                       labels = c('Hard carbonate', 'Rubble', 'Sand')),
         branched = factor(branched, levels = c('No', 'Yes'),
                           labels = c('Unbranched', 'Branched')))

overall <- overall %>% mutate(depth2 = depth)

overall <- transform(overall, depth2=revalue(depth2,c("2m"="2-3m")))
overall <- overall %>% mutate(day2 = day)
overall <- transform(overall, day2=revalue(day2,c("1"="One", "2"="Two", "3"="Three")))
overall <- filter(overall, trial == "Yellow_Rubble") # cut out the control rubble trial (not painted yellow)

overall$aspectdepth <- as.factor(paste(overall$aspect, overall$depth, sep = "_"))

overall %>% group_by(season,aspectdepth) %>%
  summarise(min = min(pt_avg_u),
            median = median(pt_avg_u),
            max = max(pt_avg_u),
            min2 = min(CORRECTED_fastest_peak_U),
            median2 = median(CORRECTED_fastest_peak_U),
            max2 = max(CORRECTED_fastest_peak_U))

overall$flipped_from_previous_day_bin <- ifelse(overall$flipped_from_previous_day_y_or_n == "Y", 1, 0)

# Model fit

wonly <- overall %>% filter(season == "W_Monsoon") %>% droplevels()

# Table S22 ----

wmod = glmmTMB(moved_at_all_over_1_cm ~
                 CORRECTED_fastest_peak_U * day2 + (1|trial_id_combined_depths/siteUn/new_id),
               data = wonly, family = binomial(link='logit'))
r.squaredGLMM(wmod) 

# Results

car::Anova(wmod) # interaction
plot(simulateResiduals(fittedModel = wmod)) # looks good
write.csv(car::Anova(wmod), file = "car::Anova(wmod).csv")

# Table S24 ----

a1.grid <- with(wonly, list(CORRECTED_fastest_peak_U = c(0.1, 0.2, 0.3, 0.4),
                            day2=levels(day2)))
a1.grid 

a1.probs <- emmeans(wmod, ~ CORRECTED_fastest_peak_U|day2, 
                    at=a1.grid, type = 'response') %>% as.data.frame

preds <- emmeans(wmod, pairwise ~ CORRECTED_fastest_peak_U|day2, at = a1.grid, type = "response")

preds

write.csv(preds$emmeans, file = "emmeans(wmod, pairwise ~ CORRECTED_fastest_peak_U|day2).csv")

# To plot:

a1.grid <- with(wonly, list(CORRECTED_fastest_peak_U = seq(min(CORRECTED_fastest_peak_U), max(CORRECTED_fastest_peak_U), len = 100),
                            day2=levels(day2)))
a1.grid #the new grid of pt_peak values

a1.probs <- emmeans(wmod, ~ CORRECTED_fastest_peak_U|day2, 
                    at=a1.grid, type = 'response') %>% as.data.frame

# Transport Northeast Monsoon

neonly <- overall %>% filter(season == "NE_Monsoon") %>% droplevels()

# Table S28 ----

nemod = glmmTMB(moved_at_all_over_1_cm ~
                  CORRECTED_fastest_peak_U * day2 + (1|trial_id_combined_depths/siteUn/new_id),
                data = neonly, family = binomial(link='logit'))
car::Anova(nemod) # interaction not significant
nemod2 <- update(nemod, ~ . - CORRECTED_fastest_peak_U:day2)
anova(nemod, nemod2) 
AICc(nemod, nemod2) 

r.squaredGLMM(nemod2) 

car::Anova(nemod2) # effect of day but not of velocity
write.csv(car::Anova(nemod2), file = "car::Anova(nemod2).csv")
plot(simulateResiduals(fittedModel = nemod2))

# Table S30 ----

preds <- emmeans(nemod2, pairwise ~ day2 , type = "response")

write.csv(preds$emmeans, file = "emmeans(nemod2, pairwise ~ day2-1.csv")
write.csv(preds$contrasts, file = "emmeans(nemod2, pairwise ~ day2-2.csv")

# Flipping Western Monsoon

# Table S23 ----

wmod.fl = glmmTMB(flipped_from_previous_day_bin ~
                    CORRECTED_fastest_peak_U * day2 + (1|trial_id_combined_depths/siteUn/new_id),
                  data = wonly, family = binomial(link='logit'))

car::Anova(wmod.fl) 

write.csv(car::Anova(wmod.fl), file = "car::Anova(wmod.fl).csv")

r.squaredGLMM(wmod.fl) 

plot(simulateResiduals(fittedModel = wmod.fl))

# Table S25 ----

grid <- with(wonly, list(CORRECTED_fastest_peak_U = c(0.1,0.2,0.3, 0.4),
                         day2 = levels(day2)))

preds1 <- emmeans(wmod.fl, pairwise ~   CORRECTED_fastest_peak_U | day2 , at = grid, type = "response")
write.csv(preds1$emmeans, file = "emmeans(wmod.fl,pairwise~CORRECTED_fastest_peak_U|day2-1.csv")

# Flipping Northeast Monsoon

# Table S29 ----

neonly <- overall %>% filter(season == "NE_Monsoon") %>% droplevels()

nemod.fl = glmmTMB(flipped_from_previous_day_bin ~
                     CORRECTED_fastest_peak_U * day2 + (1|trial_id_combined_depths/siteUn/new_id),
                   data = neonly, family = binomial(link='logit'))

car::Anova(nemod.fl) # no interaction
r.squaredGLMM(nemod.fl)

nemod.fl2 <- update(nemod.fl, ~ . - CORRECTED_fastest_peak_U:day2)

write.csv(Anova(nemod.fl2), file = "car::Anova(nemod.fl2).csv")

car::Anova(nemod.fl2) # effect of day, marginal effect of velocity

plot(simulateResiduals(fittedModel = nemod.fl2))

# Table S31 ----

mean(neonly$CORRECTED_fastest_peak_U)

preds <- emmeans(nemod.fl2, pairwise ~ day2, type = "response")

write.csv(preds$emmeans, file = "emmeans(nemod.fl2, pairwise ~ day2-1.csv")
write.csv(preds$contrasts, file = "emmeans(nemod.fl2, pairwise ~ day2-2.csv")

# Distance Western Monsoon

overall2 <- filter(overall, diag_dist > 1) # only interested in those that moved 1cm or more as buffer

levels(overall2$day2)

hist(log(overall2$diag_dist))

wonly2 <- overall2 %>% filter(season == "W_Monsoon") %>% droplevels()

# Table S26 ----

wmod.dist = lme(log(diag_dist) ~ CORRECTED_fastest_peak_U * day2,
                random = ~1|siteUn, 
                data = wonly2, method = 'ML')

Anova(wmod.dist) 

wmod.dist2 <- update(wmod.dist, ~ . - CORRECTED_fastest_peak_U:day2)

car::Anova(wmod.dist2)
write.csv(Anova(wmod.dist2), file = "Anova(wmod.dist2).csv")

car::Anova(wmod.dist2)

plot(wmod.dist2) 
qqnorm(wmod.dist2) 

plot(simulateResiduals(fittedModel = wmod.dist2))

# Table S27 ----

preds <- emmeans(wmod.dist2, pairwise ~ day2, type = "response")
write.csv(preds$emmeans, file = "emmeans(wmod.dist2, pairwise ~ day2-1.csv")
write.csv(preds$contrasts, file = "emmeans(wmod.dist2, pairwise ~ day2-2.csv")

# Distance Northeast Monsoon

neonly2 <- overall2 %>% filter(season == "NE_Monsoon") %>% droplevels()

hist(log(neonly2$diag_dist))

# Table S32 ----

nemod.dist = lme(log(diag_dist) ~ CORRECTED_fastest_peak_U * day2,
                 random = ~1|siteUn, 
                 data = neonly2, method = 'ML')

car::Anova(nemod.dist) 

nemod.dist2 <- update(nemod.dist, ~ . - CORRECTED_fastest_peak_U:day2)

car::Anova(nemod.dist2)

write.csv(Anova(nemod.dist2), file = "car::Anova(nemod.dist2).csv")

plot(simulateResiduals(fittedModel = nemod.dist))

# Table S33 ----

preds <- emmeans(nemod.dist2, pairwise ~ day2, type = "response")

write.csv(preds$emmeans, file = "emmeans(nemod.dist2, pairwise ~ day2-1.csv")
write.csv(preds$contrasts, file = "emmeans(nemod.dist2, pairwise ~ day2-2.csv")

# Figure 5a ----

# Transport:

cols <- c("One" = "#FFC1C1", "Two" = "#FFEC8B", "Three" = "#636363")
fils <- c("One" = "#FFC1C1", "Two" = "#FFEC8B", "Three" = "#636363")
lins <- c("One" = 1, "Two" = 2, "Three" = 3)

a1.grid <- with(wonly, list(CORRECTED_fastest_peak_U = seq(min(CORRECTED_fastest_peak_U), max(CORRECTED_fastest_peak_U), len = 100),
                           day2=levels(day2),
                           aspectdepth = levels(aspectdepth)))
a1.grid #the new grid of pt_peak values

a1.probs <- emmeans(wmod, ~ CORRECTED_fastest_peak_U|day2, 
                    at=a1.grid, type = 'response') %>% as.data.frame

head(a1.probs)
wonly %>% dplyr::group_by(day2) %>%
  dplyr::summarise(max = max(CORRECTED_fastest_peak_U), 
            min = min(CORRECTED_fastest_peak_U))

levels(wonly$day2)

head(a1.probs)
chunk1 = a1.probs%>%
  group_by(day2)%>%
  filter((day2=="One" & CORRECTED_fastest_peak_U < 0.44)) # highest measured max velocity on day 1 was 0.43m/s
chunk1
chunk2 = a1.probs%>%
  group_by(day2)%>%
  filter((day2=="Two" & CORRECTED_fastest_peak_U < 0.36)) # highest measured max velocity on day 2 was 0.37m/s

chunk3 = a1.probs%>%
  group_by(day2)%>%
  filter((day2=="Three" & CORRECTED_fastest_peak_U < 0.33)) # highest measured max velocity on day 2 was 0.34m/s

a1.probs

neonly %>% dplyr::group_by(aspectdepth) %>%
  summarise(max = max(CORRECTED_fastest_peak_U_cm),
            mean(CORRECTED_fastest_peak_U_cm),
            min = min(CORRECTED_fastest_peak_U_cm))

wonly %>% dplyr::group_by(aspectdepth) %>%
  summarise(max = max(CORRECTED_fastest_peak_U_cm), 
            mean(CORRECTED_fastest_peak_U_cm),
            median(CORRECTED_fastest_peak_U_cm),
            min = min(CORRECTED_fastest_peak_U_cm))

wonly.probs.move.chunk <- ggplot() +
  #geom_point(data=wonly, aes(y=moved_at_all_over_1_cm, x=CORRECTED_fastest_peak_U, colour = day2), alpha = 0.3, show.legend = FALSE) + 
  geom_ribbon(data=chunk1, aes(x=CORRECTED_fastest_peak_U, ymax=asymp.UCL, ymin=asymp.LCL, colour = day2,
                               fill=day2),  alpha=0.3, show.legend = TRUE) + 
  geom_line(data=chunk1, aes(x=CORRECTED_fastest_peak_U, y=prob,color=day2, linetype = day2), show.legend = TRUE) +
  geom_ribbon(data=chunk2, aes(x=CORRECTED_fastest_peak_U, ymax=asymp.UCL, ymin=asymp.LCL, colour = day2,
                               fill=day2),  alpha=0.3, show.legend = FALSE) + 
  geom_line(data=chunk2, aes(x=CORRECTED_fastest_peak_U, y=prob,color=day2, linetype = day2), show.legend = FALSE) +
  geom_ribbon(data=chunk3, aes(x=CORRECTED_fastest_peak_U, ymax=asymp.UCL, ymin=asymp.LCL, colour = day2,
                               fill=day2),  alpha=0.3, show.legend = FALSE) + 
  geom_line(data=chunk3, aes(x=CORRECTED_fastest_peak_U, y=prob,color=day2, linetype = day2), show.legend = FALSE) +
  scale_color_manual(name = "Days after deployment", values= cols) +
  scale_fill_manual(name = "Days after deployment", values= fils) +
  scale_linetype_manual(name = "Days after deployment", values= lins) +
  labs(colour = "Day", linetype = "Day", fill = "Day") +
  xlab("Peak bottom orbital velocity (m/s)") +
  ylab("(a) Transport >1cm") + 
  ylim(0,1) +
  theme_classic() +
  theme(text = element_text(size=25, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position=c(0.3,0.81)) +
  theme(legend.text=element_text(size=20)) +
  theme(axis.title.x=element_text(size=25,face = "bold",color="black")) +
  theme(axis.title.y=element_text(size=25,face = "bold",color="black"))
  #scale_x_continuous(breaks=(c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)), limits = c(0.1, 0.46))
wonly.probs.move.chunk

# Figure 5b ----

# Flipping:

a1.grid <- with(wonly, list(CORRECTED_fastest_peak_U = seq(min(CORRECTED_fastest_peak_U), max(CORRECTED_fastest_peak_U), len = 100),
                            day2=levels(day2),
                            aspectdepth = levels(aspectdepth)))
a1.grid #the new grid of pt_peak values

a1.probs <- emmeans(wmod.fl, ~ CORRECTED_fastest_peak_U|day2, 
                    at=a1.grid, type = 'response') %>% as.data.frame

chunk1 = a1.probs%>%
  group_by(day2)%>%
  filter((day2=="One" & CORRECTED_fastest_peak_U < 0.44))

chunk2 = a1.probs%>%
  group_by(day2)%>%
  filter((day2=="Two" & CORRECTED_fastest_peak_U < 0.36))

chunk3 = a1.probs%>%
  group_by(day2)%>%
  filter((day2=="Three" & CORRECTED_fastest_peak_U < 0.33))

a1.probs

wonly.probs.flip.chunk <- ggplot() +
  #geom_point(data=wonly, aes(y=flipped_from_previous_day_bin, x=CORRECTED_fastest_peak_U, colour = day2), alpha = 0.3, show.legend = FALSE) + 
  geom_ribbon(data=chunk1, aes(x=CORRECTED_fastest_peak_U, ymin=asymp.LCL, ymax=asymp.UCL, colour = day2, 
                               fill=day2),  alpha=0.3, show.legend = FALSE) +
  geom_line(data=chunk1, aes(x=CORRECTED_fastest_peak_U, y=prob,color=day2, linetype = day2), show.legend = FALSE) +
  geom_ribbon(data=chunk2, aes(x=CORRECTED_fastest_peak_U, ymin=asymp.LCL, ymax=asymp.UCL, colour = day2,
                               fill=day2),  alpha=0.3, show.legend = FALSE) + 
  geom_line(data=chunk2, aes(x=CORRECTED_fastest_peak_U, y=prob,color=day2, linetype = day2), show.legend = FALSE) +
  geom_ribbon(data=chunk3, aes(x=CORRECTED_fastest_peak_U, ymin=asymp.LCL, ymax=asymp.UCL, colour = day2,
                               fill=day2),  alpha=0.3, show.legend = FALSE) + 
  geom_line(data=chunk3, aes(x=CORRECTED_fastest_peak_U, y=prob,color=day2, linetype = day2), show.legend = FALSE) +
  labs(colour = "Day", linetype = "Day", fill = "Day") +
  xlab("Peak bottom orbital velocity (m/s)") +
  ylab("(b) Flipping") + 
  scale_color_manual(name = "Days after deployment", values= cols) +
  scale_fill_manual(name = "Days after deployment", values= fils) +
  scale_linetype_manual(name = "Days after deployment", values= lins) +
  ylim(0,1) +
  theme_classic() +
  theme(text = element_text(size=25, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position=c(0.27,0.81)) +
  theme(legend.text=element_text(size=20)) +
  theme(axis.title.x=element_text(size=25,face = "bold",color="black")) +
  theme(axis.title.y=element_text(size=25,face = "bold",color="black"))
wonly.probs.flip.chunk

# Figure 5c ----

# Distance:

a1.grid <- with(wonly2, list(CORRECTED_fastest_peak_U = seq(min(CORRECTED_fastest_peak_U), max(CORRECTED_fastest_peak_U), len = 100),
                             day2=levels(day2)))

a1.grid #the new grid of pt_peak values

a1.probs <- emmeans(wmod.dist2, ~ CORRECTED_fastest_peak_U|day2, 
                    at=a1.grid, type = 'response') %>% as.data.frame

chunk1 = a1.probs%>%
  group_by(day2)%>%
  filter((day2=="One" & CORRECTED_fastest_peak_U < 0.44))

chunk2 = a1.probs%>%
  group_by(day2)%>%
  filter((day2=="Two" & CORRECTED_fastest_peak_U < 0.36))

chunk3 = a1.probs%>%
  group_by(day2)%>%
  filter((day2=="Three" & CORRECTED_fastest_peak_U < 0.33))

head(a1.probs)

wonly.dist.chunk <- ggplot() +
  #geom_point(data=wonly2, aes(y=diag_dist, x=CORRECTED_fastest_peak_U, colour = day2), alpha = 0.3, show.legend = FALSE) + 
  geom_ribbon(data=chunk1, aes(x=CORRECTED_fastest_peak_U, ymin=lower.CL, ymax=upper.CL, colour = day2, 
                               fill=day2),  alpha=0.3, show.legend = FALSE) +
  geom_line(data=chunk1, aes(x=CORRECTED_fastest_peak_U, y=response,color=day2, linetype = day2), show.legend = FALSE) +
  geom_ribbon(data=chunk2, aes(x=CORRECTED_fastest_peak_U, ymin=lower.CL, ymax=upper.CL, colour = day2,
                               fill=day2),  alpha=0.3, show.legend = FALSE) + 
  geom_line(data=chunk2, aes(x=CORRECTED_fastest_peak_U, y=response,color=day2, linetype = day2), show.legend = FALSE) +
  geom_ribbon(data=chunk3, aes(x=CORRECTED_fastest_peak_U, ymin=lower.CL, ymax=upper.CL, colour = day2,
                               fill=day2),  alpha=0.3, show.legend = TRUE) + 
  geom_line(data=chunk3, aes(x=CORRECTED_fastest_peak_U, y=response,color=day2, linetype = day2), show.legend = TRUE) +
  labs(colour = "Day", linetype = "Day", fill = "Day") +
  xlab("Peak bottom orbital velocity (m/s)") +
  ylab("(c) Distance (cm)") +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(2,11)) +
  scale_color_manual(name = "Days after deployment", values= cols) +
  scale_fill_manual(name = "Days after deployment", values= fils) +
  scale_linetype_manual(name = "Days after deployment", values= lins) +
  theme_classic() +
  theme(text = element_text(size=25, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position="none")+
  theme(legend.text=element_text(size=20)) +
  theme(axis.title.x=element_text(size=25,face = "bold",color="black")) +
  theme(axis.title.y=element_text(size=25,face = "bold",color="black"))
wonly.dist.chunk

grid.arrange(wonly.probs.move.chunk, wonly.probs.flip.chunk, wonly.dist.chunk)

# DAY 1 ONLY ----

# Rubble & substrate characteristics: Transport ----

rm(list=ls())

overall <- read.csv("~/Dropbox/PhD/Data/Maldives/Rubble_Movement/Shortterm_Movement/Analysis/Rubble_Transport_Shortterm_COL_v3-Sheet4-corrected_u-git.csv")

overall$aspect<-factor(overall$aspect)
overall$site<-factor(overall$site)
overall$siteUn <- factor(overall$siteUn)
overall$depth<-factor(overall$depth)
overall$year<-factor(overall$year)
overall$trial_id_each_depth<-factor(overall$trial_id_each_depth)
overall$trial_id_combined_depths<-factor(overall$trial_id_combined_depths)
overall$season<-factor(overall$season)
overall$new_id<-factor(overall$new_id)
overall$trial<-factor(overall$trial)
overall$branched<-factor(overall$branched)
overall$size_cat_adjusted<-factor(overall$size_cat_adjusted)
overall$branches_opposing<-factor(overall$branches_opposing)
overall$day<-factor(overall$day)
overall$starting_substrate_new<-factor(overall$starting_substrate_new)
overall$substrate_stability<-factor(overall$substrate_stability)
overall$buoy_wt<-as.numeric(overall$buoy_wt)
overall$branch3_l<-as.numeric(overall$branch3_l)
overall$dry_wt<-as.numeric(overall$dry_wt)
overall$wet_wt<-as.numeric(overall$wet_wt)
overall$vol<-as.numeric(overall$vol)
overall$flipped_on_this_day_tidying <- as.character(overall$flipped_on_this_day_tidying)

glimpse(overall)

# reorder levels

overall <- transform(overall, 
                     starting_substrate_new=revalue(starting_substrate_new,
                                                    c("rubble " = "rubble", 
                                                      "sand " = "sand")))

overall = overall %>% 
  mutate(aspect=factor(aspect, levels = c('Lagoon', 'SE', 'West')), #re-leveling
         depth=factor(depth, levels = c('2m', '2-3m','6-7m')),
         size_cat_adjusted_to_flume_exp=factor(size_cat_adjusted_to_flume_exp,
                                               levels = c('Small', 'Medium','Large '), 
                                               labels = c('4-8 cm', '9-15 cm', '16-23 cm')),
         CORRECTED_fastest_peak_U_cm = CORRECTED_fastest_peak_U * 100,
         pt_avg_u_cm = pt_avg_u * 100,
         pt_daily_avg_peak_cm = pt_daily_avg_peak * 100,
         current_peak_u_cm = current_peak_u * 100,
         current_avg_u_cm = current_avg_u * 100,
         current_daily_avg_peak_cm = current_daily_avg_peak * 100,
         starting_substrate_new=factor(starting_substrate_new,
                                       levels = c('hard_carb', 'rubble', 'sand'), 
                                       labels = c('Hard carbonate', 'Rubble', 'Sand')),
         branched = factor(branched, levels = c('No', 'Yes'),
                           labels = c('Unbranched', 'Branched')))

#Make a new variable (column) called depthv2 that is a direct copy of 'depth' variable.
#Then, rename the Reef flat/Lagoon ones to 2-3m, not 2m.

overall <- overall %>% mutate(depth2 = depth)

# recoding/revalueing

overall <- transform(overall, depth2=revalue(depth2,c("2m"="2-3m")))

#Do the same for day 

overall <- overall %>% mutate(day2 = day)
overall <- transform(overall, day2=revalue(day2,c("1"="One", "2"="Two", "3"="Three")))

overall <- filter(overall, trial == "Yellow_Rubble") # cut out the control rubble trial (where rubble was not yellow, but rather unpainted)

over.1 <- filter(overall, day == "1")

# Other investigations ----

# What zone had most movement on Day 1, when there was a velocity relationship?

over.1$aspectDepth <- paste(over.1$aspect, over.1$depth2, sep="_")

over.1$aspectDepth <- factor(over.1$aspectDepth,
                              levels = c('Lagoon_2-3m', 'SE_2-3m','SE_6-7m',
                                         'West_2-3m','West_6-7m'), 
                              labels = c('Lag_Shal', 'Shelt_Shal', 
                                         'Shelt_Deep', 'Exp_Shal', 'Exp_Deep'))

names(over.1)
movesumm <- over.1 %>% dplyr::group_by(season, aspectDepth) %>%
  dplyr::summarize(minmove = min(moved_at_all_over_1_cm),
                   meanmove = median(moved_at_all_over_1_cm),
                   maxmove = max(moved_at_all_over_1_cm))

ggplot(over.1, aes(x = aspectDepth, y = moved_at_all_over_1_cm)) +
  geom_bar(stat = "identity") + facet_wrap(~season)

# more movement in W monsoon and at shallow sites in west monsoon.
# most movement Exp deep NE monsoon - perhaps because steeper than SE deep

ggplot(over.1, aes(x = aspectDepth, y = (mean(diag_dist)))) +
  geom_bar(stat = "identity") + facet_wrap(~season)

overall$aspectDepth <- paste(overall$aspect, overall$depth2, sep="_")

overall$aspectDepth <- factor(overall$aspectDepth,
                             levels = c('Lagoon_2-3m', 'SE_2-3m','SE_6-7m',
                                        'West_2-3m','West_6-7m'), 
                             labels = c('Lag_Shal', 'Shelt_Shal', 
                                        'Shelt_Deep', 'Exp_Shal', 'Exp_Deep'))

ggplot(overall, aes(x = aspectDepth, y = moved_at_all_over_1_cm)) +
  geom_bar(stat = "identity") + facet_wrap(~season) +
  theme_classic() +
  ylab("Total instances where rubble was transported >1 cm") +
  xlab("Habitat & Depth") +
  theme(text = element_text(size=25, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position=c(0.3,0.81)) +
  theme(legend.text=element_text(size=20)) +
  theme(axis.title.x=element_text(size=25,face = "bold",color="black")) +
  theme(axis.title.y=element_text(size=25,face = "bold",color="black"))

hist(overall$moved_at_all_over_1_cm)

totalmod <- glmmTMB(moved_at_all_over_1_cm ~ 
                          season * aspectDepth + 
                          (1|siteUn),
                        data = overall,
                        family = binomial(link='logit'))

car::Anova(totalmod) # interaction between season and aspectDepth

plot(simulateResiduals(totalmod)) # good

emmeans(totalmod, pairwise ~ aspectDepth | season, type = "response")
# hardly any significant differences between aspectDepth comparisons

emmeans(totalmod, pairwise ~  season | aspectDepth, type = "response")

# exposed shallow sites had more movement in western than NE monsoon.
# sheltered shallow is a trend toward this, but non-sig.

totalmod2 <- glmmTMB(moved_at_all_over_1_cm ~ 
                      season * aspectDepth + 
                      (1|siteUn),
                    data = over.1,
                    family = binomial(link='logit'))

car::Anova(totalmod2) # interaction between season and aspectDepth

plot(simulateResiduals(totalmod2)) # good

emmeans(totalmod2, pairwise ~ aspectDepth | season, type = "response")
# hardly any significant differences between aspectDepths

emmeans(totalmod2, pairwise ~  season | aspectDepth, type = "response")


# Table S35 ----

fi.m.senglob <- glmmTMB(moved_at_all_over_1_cm ~ 
                          CORRECTED_fastest_peak_U + 
                          size_cat_adjusted_to_flume_exp + 
                          branched +
                          CORRECTED_fastest_peak_U:size_cat_adjusted_to_flume_exp +
                          CORRECTED_fastest_peak_U:branched +
                          size_cat_adjusted_to_flume_exp:starting_substrate_new +
                          size_cat_adjusted_to_flume_exp:branched +
                          starting_substrate_new:branched +
                          starting_substrate_new + 
                          (1|siteUn),
                        data = over.1,
                        family = binomial(link='logit'))

car::Anova(fi.m.senglob)

fi.m.senglob2 <- update(fi.m.senglob, ~ . - branched:starting_substrate_new)
car::Anova(fi.m.senglob2)
AICc(fi.m.senglob, fi.m.senglob2) # lower AIC

fi.m.senglob3 <- update(fi.m.senglob2, ~ . - size_cat_adjusted_to_flume_exp:branched)
car::Anova(fi.m.senglob3)
AICc(fi.m.senglob, fi.m.senglob3)

fi.m.senglob4 <- update(fi.m.senglob3, ~ . - size_cat_adjusted_to_flume_exp:starting_substrate_new)
car::Anova(fi.m.senglob4)
AICc(fi.m.senglob, fi.m.senglob4)

fi.m.senglob5 <- update(fi.m.senglob4, ~ . - CORRECTED_fastest_peak_U:branched)
car::Anova(fi.m.senglob5)
AICc(fi.m.senglob, fi.m.senglob5)

write.csv(car::Anova(fi.m.senglob5), file = "car::Anova(fi.m.senglob5).csv")

r.squaredGLMM(fi.m.senglob5)

plot(simulateResiduals(fittedModel = fi.m.senglob5)) # good

# Table S36 ----

fi.m.senglob5.grid <- with(over.1, list(CORRECTED_fastest_peak_U = c(min(CORRECTED_fastest_peak_U), 
                                                                     0.1, mean(CORRECTED_fastest_peak_U), 
                                                                     0.19, 0.25, 0.3, 0.4, max(CORRECTED_fastest_peak_U)),
                           size_cat_adjusted_to_flume_exp = levels(size_cat_adjusted_to_flume_exp)))
fi.m.senglob5.grid

# On day 1, mean is 0.15 m/s, max is 0.43 m/s (rather than 0.55 m/s which it was before the change to wave number, and the min is 0.01 m/s)

fi.m.senglob5.preds <- emmeans(fi.m.senglob5, pairwise ~ size_cat_adjusted_to_flume_exp | CORRECTED_fastest_peak_U, 
                               at=fi.m.senglob5.grid, type = 'response')
fi.m.senglob5.preds

write.csv(fi.m.senglob5.preds$emmeans, file = "emmeans(fi.m.senglob5, pairwise ~ size_cat_adjusted_to_flume_exp-1.csv") # estimates
write.csv(fi.m.senglob5.preds$contrasts, file = "emmeans(fi.m.senglob5, pairwise ~ size_cat_adjusted_to_flume_exp-2.csv") # contrasts

# Table S37 ----

preds <- emmeans(fi.m.senglob5, pairwise ~ branched, type = 'response')

write.csv(preds$emmeans, file = "emmeans(fi.m.senglob5, pairwise ~ branched.csv")
write.csv(preds$contrasts, file = "emmeans(fi.m.senglob5, pairwise ~ branched-2.csv")

# Slope

# Table S38 ----

compnewSLONLY <- glmmTMB(moved_at_all_over_1_cm ~ 
                           CORRECTED_fastest_peak_U * 
                           avg_slope_angle_transect.1 +
                           (1|siteUn),
                         data = over.1,
                         family = binomial(link='logit'))

r.squaredGLMM(compnewSLONLY)

simulateResiduals(compnewSLONLY, plot = T)  # good

summary(compnewSLONLY) 
write.csv(car::Anova(compnewSLONLY), file = "compnewSLONLY.csv")

car::Anova(compnewSLONLY)

# Table S39 ----

compnewSLONLY.grid3 <- with(over.1, list(CORRECTED_fastest_peak_U = c(0.1, 0.2, 0.3, 0.4),
                                         avg_slope_angle_transect.1 = c(3, 13, 22)))

compnewSLONLY.grid3.preds <- emmeans(compnewSLONLY, pairwise ~ avg_slope_angle_transect.1 | CORRECTED_fastest_peak_U, 
                                     at=compnewSLONLY.grid3, type = 'response')

compnewSLONLY.grid3.preds <- emmeans(compnewSLONLY, pairwise ~ CORRECTED_fastest_peak_U | avg_slope_angle_transect.1, 
                                     at=compnewSLONLY.grid3, type = 'response')

compnewSLONLY.grid3.preds

write.csv(compnewSLONLY.grid3.preds$emmeans, file = "emmeans(compnewSLONLY, pairwise ~ avg_slope_angle_transect.1 | CORRECTED_fastest_peak_U-1.csv")
write.csv(compnewSLONLY.grid3.preds$contrasts, file = "emmeans(compnewSLONLY, pairwise ~ avg_slope_angle_transect.1 | CORRECTED_fastest_peak_U-2.csv")

# Figure 6a ----

#(a) Transport

cols <- c("4-8 cm" = "#AAADAA", "9-15 cm" = "#64CC7A", "16-23 cm" = "#7DB9EF")
fils <- c("4-8 cm" = "#AAADAA", "9-15 cm" = "#64CC7A", "16-23 cm" = "#7DB9EF")
lins <- c("4-8 cm" = 3, "9-15 cm" = 2, "16-23 cm" = 1)

fi.m.senglob5.grid <- with(over.1, list(CORRECTED_fastest_peak_U = seq(min(CORRECTED_fastest_peak_U), 0.75, len = 100),
                                        size_cat_adjusted_to_flume_exp = levels(size_cat_adjusted_to_flume_exp),
                                        branched=levels(branched),
                                        starting_substrate_new = levels(starting_substrate_new)))
fi.m.senglob5.grid

fi.m.senglob5.preds <- emmeans(fi.m.senglob5, ~ CORRECTED_fastest_peak_U | size_cat_adjusted_to_flume_exp | branched, 
                               at=fi.m.senglob5.grid, type = 'response') %>% as.data.frame()

View(fi.m.senglob5.preds)

fi.m.senglob5.plot <- fi.m.senglob5.preds %>% ggplot(aes(x = CORRECTED_fastest_peak_U, y = prob, 
                                                         colour = size_cat_adjusted_to_flume_exp)) +
  geom_ribbon(aes(ymin=asymp.LCL, ymax=asymp.UCL, 
                  fill=size_cat_adjusted_to_flume_exp), alpha =0.3) +
  geom_line(aes(color = size_cat_adjusted_to_flume_exp, 
                linetype = size_cat_adjusted_to_flume_exp)) +
  theme_classic() + 
  facet_wrap(~branched) +
  scale_color_manual(name = "Rubble length", values= cols) +
  scale_fill_manual(name = "Rubble length", values= fils) +
  scale_linetype_manual(name = "Rubble length", values= lins) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), 
                     limits = c(0,1)) +
  xlab(bquote('Peak bottom orbital velocity ('*m/s*')')) +
  ylab(bquote('(a) Transport <1 cm')) +
  theme_classic() +
  theme(text = element_text(size=25, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position="right") +
  theme(legend.text=element_text(size=20)) +
  theme(axis.title.x=element_text(size=25,face = "bold",color="black")) +
  theme(axis.title.y=element_text(size=25,face = "bold",color="black")) +
  theme(panel.background =element_rect(colour = "black", linewidth=1)) + 
  theme(axis.ticks.length=unit(.2,"cm"))
fi.m.senglob5.plot

ggsave("fi.m.senglob5.plot.pdf", 
       fi.m.senglob5.plot, width = 1170/90, 
       height = 536/90, dpi = 900)

# (b) Flipping 

# see after line 1203 for Figure 6b


# Figure 6c ----

#(c) Heat Map

compnewSLONLY.grid2 <- with(over.1, list(CORRECTED_fastest_peak_U = seq(min(CORRECTED_fastest_peak_U), max(CORRECTED_fastest_peak_U), len = 100),
                                         avg_slope_angle_transect.1 = seq(min(avg_slope_angle_transect.1), max(avg_slope_angle_transect.1), len = 100)))
compnewSLONLY.grid2

compnewSLONLY.pred <- emmeans(compnewSLONLY, ~ avg_slope_angle_transect.1 | CORRECTED_fastest_peak_U, 
                              at=compnewSLONLY.grid2, type = 'response') %>% as.data.frame()
head(compnewSLONLY.pred)

summary(compnewSLONLY.pred)
b <- c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

slope1 <- ggplot(compnewSLONLY.pred) + 
  aes(x = CORRECTED_fastest_peak_U, y = avg_slope_angle_transect.1, z = prob, fill = prob) + 
  geom_tile() + 
  geom_contour(color = "white", breaks = b, alpha = 0.5) + 
  scale_fill_gradientn(name = "Probability", limits = c(0.0,1.0),
                       colours=c("white", "grey", "lightskyblue", "royalblue", "black"),
                       breaks=b, labels=format(b)) +
  # low = "lightblue", high = "red") + 
  #facet_grid(branched~substrate) + 
  xlab(bquote('Bottom orbital velocity ('*m/s*')')) +
  ylab(bquote('Slope angle (degrees)')) +
  theme_classic() +
  theme(text = element_text(size=25, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", linewidth=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position="top") +
  theme(legend.text=element_text(size=11)) +
  theme(legend.title=element_text(size=15))
  #theme(axis.title.x=element_text(size=25,face = "bold",color="black")) +
  #theme(axis.title.y=element_text(size=25,face = "bold",color="black"))

slope2 <- slope1 + guides(fill = guide_colourbar(barwidth = 15, barheight = 1))

slope2

#http://127.0.0.1:30971/graphics/plot_zoom_png?width=590&height=621

ggsave("slope3.pdf", 
       slope2, width = 536/90, 
       height = 536/90, dpi = 900)

# Table S34 (Overall threshold) ----

# The 50% and 90% mobilisation thresholds for transport (> 1 cm) in the field
# averaged across all rubble sizes (4–23 cm), branchiness and substrate characteristics

ldmod <- glmmTMB(moved_at_all_over_1_cm ~ 
                   CORRECTED_fastest_peak_U +
                   (1|siteUn),
                 data = over.1,
                 family = binomial(link='logit'))

car::Anova(ldmod)

r.squaredGLMM(ldmod)

plot(simulateResiduals(ldmod)) # good

write.csv(car::Anova(ldmod), file = "carAnova-ldmod.csv")

# Getting an ld50 and ld90 for rubble pieces from size 4cm to 23cm, on rubble, sand or hard carbonate, 
# on gentle or steep slope, branched or unbranched

summary(ldmod) # get the intercept and slope from your model.
# Intercept is -1.4418 (a)
# slope is 4.8526 (b)

#log(p/(1-p)) = a + b*x (y is the probability)

#For the LD50 we seek the value of x when p = 0.5. Plugging this p into the equation we get:

# y = a + bx
# log(0.5/(1-0.5)) = a + b*x (where log(0.5/(1-0.5)) = y)

# 0 = a + b*x

# implying that x = -a/b for 50% threshold

# For LD90, substitute p=0.9 for 0.5. 

# So to work out the threshold:

log(0.5/(1-0.5)) # = 0         
# thus 0 = -a + bx, 
# thus 0 = -1.441828 + 4.8526x
# thus x = --1.441828/4.8526

log(0.9/(1-0.9)) # = 2.197225,
# thus 2.197225 = -1.4418 + 4.8526x
# thus 2.197225 + 1.4418 = -1.4418 + 4.8526x + 1.4418
# thus 3.639025 = 4.8526x
# thus 3.639025 / 4.8526 = 4.8526x / 4.8526

xmove50 <- --1.441828/4.8526 # this is the velocity at which 50% of the pieces will move

xmove90 <- 3.639025 / 4.8526 # this is the velocity at which 90% of the pieces will move

xmove50 # 50% threshold  0.297 m/s (0.3m/s)
xmove90 # 90% threshold 0.75 m/s

# Now calculate the standard errors:

# For 50% threshold

summary(ldmod)
# Intercept is -1.441828 
# slope is 4.8526
a <- -1.441828
b <- 4.8526

V <- vcov(ldmod)
V
row1 <- c(0.03558478,-0.1097155) # intercept and slope of vcov for row1
row2 <- c(-0.10971555,0.7035909) # intercept and slope of vcov for row2

V <- rbind(row1,row2)
pds <- matrix(c(-1/b, a/(b^2)),nrow=1)
pds
pds2 <- matrix(c(-1/b, a/(b^2)),ncol=1)
pds2
matrix(V)

#%*% pds2

se50 <- sqrt (pds %*% V %*% pds2)

se50 # So, the se for the 50% threshold is 0.037 m/s 

# For90% threshold:

pds90 <- matrix(c(-1, (-2.2/(b^2))),nrow=1)
pds90
pds290 <- matrix(c(-1, (-2.2/(b^2))),ncol=1)
pds290

se90 <- sqrt (pds90 %*% V %*% pds290)
se90 # so the se for the 90% threshold is 0.14 m/s

# Plot to see if the 50% and 90% thresholds from the above calcs are right

grid <- with(over.1, list(CORRECTED_fastest_peak_U = seq(min(CORRECTED_fastest_peak_U),0.76, len = 100)))

grid

ldmodpreds <- emmeans(ldmod, ~ CORRECTED_fastest_peak_U, 
                              at=grid, type = 'response') %>% as.data.frame()
head(ldmodpreds)

ggplot(ldmodpreds, aes(x = CORRECTED_fastest_peak_U, y = prob)) + ylim(0,1) +
  geom_line() + geom_ribbon(aes(ymin=asymp.LCL, ymax=asymp.UCL), alpha =0.3) +
  xlim(0,0.8)

# 50% threshold is bang on.
# It is predicting above the max velocity measured for the 90% threshold, so this should be considered with caution.


# Rubble & substrate characteristics: Flipping ----

rm(list=ls())

overall <- read.csv("~/Dropbox/PhD/Data/Maldives/Rubble_Movement/Shortterm_Movement/Analysis/Rubble_Transport_Shortterm_COL_v3-Sheet4-corrected_u-git.csv")
names(overall)
overall$aspect<-factor(overall$aspect)
overall$site<-factor(overall$site)
overall$siteUn <- factor(overall$siteUn)
overall$depth<-factor(overall$depth)
overall$year<-factor(overall$year)
overall$trial_id_each_depth<-factor(overall$trial_id_each_depth)
overall$trial_id_combined_depths<-factor(overall$trial_id_combined_depths)
overall$season<-factor(overall$season)
overall$new_id<-factor(overall$new_id)
overall$trial<-factor(overall$trial)
overall$branched<-factor(overall$branched)
overall$size_cat_adjusted<-factor(overall$size_cat_adjusted)
overall$branches_opposing<-factor(overall$branches_opposing)
overall$day<-factor(overall$day)
overall$starting_substrate_new<-factor(overall$starting_substrate_new)
overall$substrate_stability<-factor(overall$substrate_stability)
overall$buoy_wt<-as.numeric(overall$buoy_wt)
overall$branch3_l<-as.numeric(overall$branch3_l)
overall$dry_wt<-as.numeric(overall$dry_wt)
overall$wet_wt<-as.numeric(overall$wet_wt)
overall$vol<-as.numeric(overall$vol)
overall$flipped_on_this_day_tidying <- as.character(overall$flipped_on_this_day_tidying)

# reorder levels

overall <- transform(overall, 
                     starting_substrate_new=revalue(starting_substrate_new,
                                                    c("rubble " = "rubble", 
                                                      "sand " = "sand")))

overall = overall %>% 
  mutate(aspect=factor(aspect, levels = c('Lagoon', 'SE', 'West')), #releveling
         depth=factor(depth, levels = c('2m', '2-3m','6-7m')),
         size_cat_adjusted_to_flume_exp=factor(size_cat_adjusted_to_flume_exp,
                                               levels = c('Small', 'Medium','Large '),
                                               labels = c('4-8 cm', '9-15 cm', '16-23 cm')),
         CORRECTED_fastest_peak_U_cm = CORRECTED_fastest_peak_U * 100,
         pt_avg_u_cm = pt_avg_u * 100,
         pt_daily_avg_peak_cm = pt_daily_avg_peak * 100,
         current_peak_u_cm = current_peak_u * 100,
         current_avg_u_cm = current_avg_u * 100,
         current_daily_avg_peak_cm = current_daily_avg_peak * 100,
         starting_substrate_new=factor(starting_substrate_new,
                                       levels = c('hard_carb', 'rubble', 'sand'), 
                                       labels = c('Hard carbonate', 'Rubble', 'Sand')),
         branched = factor(branched, levels = c('No', 'Yes'),
                           labels = c('Unbranched', 'Branched')))

#Make a new variable (column) called depthv2 that is a direct copy of 'depth' variable.
#Then, rename the Reef flat/Lagoon ones to 2-3m, not 2m.

overall <- overall %>% mutate(depth2 = depth)

# recoding/revalueing

overall <- transform(overall, depth2=revalue(depth2,c("2m"="2-3m")))

overall <- overall %>% mutate(day2 = day)

overall <- transform(overall, day2=revalue(day2,c("1"="One", "2"="Two", "3"="Three")))

overall <- filter(overall, trial == "Yellow_Rubble")

overall$flipped_from_previous_day_bin <- ifelse(overall$flipped_from_previous_day_y_or_n == "Y", 1, 0)

over.1 <- filter(overall, day == "1")

# Table S40 ----

fi.f.senglob <- glmmTMB(flipped_from_previous_day_bin ~ 
                          CORRECTED_fastest_peak_U + 
                          size_cat_adjusted_to_flume_exp + 
                          branched +
                          CORRECTED_fastest_peak_U:size_cat_adjusted_to_flume_exp +
                          CORRECTED_fastest_peak_U:branched +
                          size_cat_adjusted_to_flume_exp:starting_substrate_new +
                          size_cat_adjusted_to_flume_exp:branched +
                          starting_substrate_new:branched +
                          starting_substrate_new + 
                          (1|siteUn),
                        data = over.1,
                        family = binomial(link='logit'))

r.squaredGLMM(fi.f.senglob)

ggplot() +
  geom_point(data=NULL, aes(y=resid(fi.f.senglob), x = fitted(fi.f.senglob))) 

library(sjPlot)
plot_model(fi.f.senglob, type = "diag") %>% plot_grid()

car::Anova(fi.f.senglob)

summary(fi.f.senglob)

fi.f.senglob2 <- update(fi.f.senglob, ~ . - branched:starting_substrate_new)
AICc(fi.f.senglob, fi.f.senglob2) # lower AIC
car::Anova(fi.f.senglob2)

fi.f.senglob3 <- update(fi.f.senglob2, ~ . - size_cat_adjusted_to_flume_exp:starting_substrate_new)
AICc(fi.f.senglob, fi.f.senglob3) # lower AIC
car::Anova(fi.f.senglob3)

fi.f.senglob4 <- update(fi.f.senglob3, ~ . - CORRECTED_fastest_peak_U:branched)
AICc(fi.f.senglob, fi.f.senglob4) 
car::Anova(fi.f.senglob4)

fi.f.senglob5 <- update(fi.f.senglob, ~ . - branched:starting_substrate_new - 
                          size_cat_adjusted_to_flume_exp:starting_substrate_new -
                          CORRECTED_fastest_peak_U:branched)

car::Anova(fi.f.senglob5)

summary(fi.f.senglob5)
write.csv(car::Anova(fi.f.senglob5), file = "car::Anova(fi.f.senglob5).csv")
r.squaredGLMM(fi.f.senglob5)

plot(simulateResiduals(fittedModel = fi.f.senglob5)) # good

cols <- c("4-8 cm" = "#AAADAA", "9-15 cm" = "#64CC7A", "16-23 cm" = "#7DB9EF")
fils <- c("4-8 cm" = "#AAADAA", "9-15 cm" = "#64CC7A", "16-23 cm" = "#7DB9EF")
lins <- c("4-8 cm" = 3, "9-15 cm" = 2, "16-23 cm" = 1)

fi.f.senglob5.grid <- with(over.1, list(CORRECTED_fastest_peak_U = seq(min(CORRECTED_fastest_peak_U), max(CORRECTED_fastest_peak_U), len = 100),
                                        size_cat_adjusted_to_flume_exp = levels(size_cat_adjusted_to_flume_exp),
                                        branched=levels(branched),
                                        starting_substrate_new = levels(starting_substrate_new)))
fi.f.senglob5.grid

fi.f.senglob5.preds <- emmeans(fi.f.senglob5, ~ CORRECTED_fastest_peak_U | size_cat_adjusted_to_flume_exp | branched, 
                               at=fi.f.senglob5.grid, type = 'response') %>% as.data.frame
head(fi.f.senglob5.preds)

pd <- position_dodge(0.5)

fi.f.senglob5.plot <- fi.f.senglob5.preds %>% ggplot(aes(x = CORRECTED_fastest_peak_U, y = prob, 
                                                         colour = size_cat_adjusted_to_flume_exp)) +
  geom_ribbon(aes(ymin=asymp.LCL, ymax=asymp.UCL, 
                  fill=size_cat_adjusted_to_flume_exp), alpha =0.3) +
  geom_line(aes(color = size_cat_adjusted_to_flume_exp, 
                linetype = size_cat_adjusted_to_flume_exp)) +
  theme_classic() + 
  facet_wrap(~branched) +
  scale_color_manual(name = "Rubble length", values= cols) +
  scale_fill_manual(name = "Rubble length", values= fils) +
  scale_linetype_manual(name = "Rubble length", values= lins) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), 
                     limits = c(0,1)) +
  scale_x_continuous(limits = c(0,0.6)) +
  xlab(bquote('Peak bottom orbital velocity ('*m/s*')')) +
  ylab(bquote('Probability of Flipping')) +
  theme(plot.margin = unit(c(9, 9, 9, 9),"points"),
        #panel.grid.major.y = element_line(color = "gray63",
        #                                 size = 0.2, linetype = "dotted"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=12,family = "sans"),
        legend.text = element_text(size=10,family = "sans"),
        legend.position= c(0.66,0.98), #0,0 is for bottom,left
        legend.justification = c(0.66,0.98), #1,1 is for top,right
        axis.text.x.bottom = element_text(size=11,
                                          vjust = -0.05),
        axis.text.y.left = element_text(size=11),
        axis.title.x=element_text(size=14,family = "sans",
                                  color="#000000",vjust = -1),  
        axis.title.y=element_text(size=14,family = "sans",
                                  color="#000000", vjust = 3),
        strip.text.x= element_text(size = 11))

fi.f.senglob5.plot


# Figure 6b -----

fi.f.senglob5.plot <- fi.f.senglob5.preds %>% ggplot(aes(x = CORRECTED_fastest_peak_U, y = prob, 
                                                         colour = size_cat_adjusted_to_flume_exp)) +
  geom_ribbon(aes(ymin=asymp.LCL, ymax=asymp.UCL, 
                  fill=size_cat_adjusted_to_flume_exp), alpha =0.3) +
  geom_line(aes(color = size_cat_adjusted_to_flume_exp, 
                linetype = size_cat_adjusted_to_flume_exp)) +
  theme_classic() + 
  facet_wrap(~branched) +
  scale_color_manual(name = "Rubble length", values= cols) +
  scale_fill_manual(name = "Rubble length", values= fils) +
  scale_linetype_manual(name = "Rubble length", values= lins) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), 
                     limits = c(0,1)) +
  xlab(bquote('Peak bottom orbital velocity ('*m/s*')')) +
  ylab(bquote('Probability of Flipping')) +
  theme_classic() +
  theme(text = element_text(size=25, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position="right") +
  theme(legend.text=element_text(size=20)) +
  theme(axis.title.x=element_text(size=25,face = "bold",color="black")) +
  theme(axis.title.y=element_text(size=25,face = "bold",color="black")) +
  theme(panel.background =element_rect(colour = "black", linewidth=1)) + 
  theme(axis.ticks.length=unit(.2,"cm"))

fi.f.senglob5.plot

ggsave("fi.f.senglob5.plot.pdf", 
       fi.f.senglob5.plot, width = 1170/90, 
       height = 536/90, dpi = 900)

# Table S41 ----

preds <- emmeans(fi.f.senglob5, pairwise ~ branched | size_cat_adjusted_to_flume_exp,
        type = "response")

write.csv(preds$contrasts, file = "emmeans(fi.f.senglob5, pairwise ~ branched | size_cat_adjusted_to_flume_exp.csv")

# Table S42

preds2 <- emmeans(fi.f.senglob5, pairwise ~  size_cat_adjusted_to_flume_exp | branched,
                 type = "response")

write.csv(preds2$contrasts, file = "emmeans(fi.f.senglob5, pairwise ~  size_cat_adjusted_to_flume_exp | branched.csv")

# Table S43 ----

fi.f.senglob5.grid <- with(over.1, list(CORRECTED_fastest_peak_U = c(0.1, 0.2, 0.3, 0.4),
                                        size_cat_adjusted_to_flume_exp = levels(size_cat_adjusted_to_flume_exp)))

fi.f.senglob5.grid

preds <- emmeans(fi.f.senglob5, pairwise ~ size_cat_adjusted_to_flume_exp | CORRECTED_fastest_peak_U,
                 at=fi.f.senglob5.grid, type = "response")

preds

write.csv(preds$contrasts, file = "emmeans(fi.f.senglob5, pairwise ~ size_cat_adjusted_to_flume_exp | CORRECTED_fastest_peak_U.csv")


# Slope

# Table S44 ----

flip.compnewSLONLY <- glmmTMB(flipped_from_previous_day_bin ~ 
                                CORRECTED_fastest_peak_U * 
                                avg_slope_angle_transect.1 +
                                (1|siteUn),
                              data = over.1,
                              family = binomial(link='logit'))

summary(flip.compnewSLONLY) # nothing is significant

car::Anova(flip.compnewSLONLY)
r.squaredGLMM(flip.compnewSLONLY)

# Rubble & substrate characteristics: Distance ----

overall2 <- filter(overall, diag_dist > 1) # picking out those that moved over 1cm in distance only
over.2 <- filter(overall2, day == "1") # now, the data is cut down to 446 rows, If all days, there are 959 rows.

# Table S45 ----

fi.d.add2 <- lme(log(diag_dist) ~
                   CORRECTED_fastest_peak_U + 
                   size_cat_adjusted_to_flume_exp +
                   branched +
                   starting_substrate_new, 
                 random = ~1|siteUn, 
                 data = over.2, method = 'ML') # cannot do an interaction model due to the low replication

plot(fi.d.add2) 
qqnorm(fi.d.add2) 

r.squaredGLMM(fi.d.add2) 

car::Anova(fi.d.add2)
write.csv(car::Anova(fi.d.add2), file = "car::Anova(fi.d.add2).csv")

summary(fi.d.add2)

# Table S46 -----

preds <- emmeans(fi.d.add2, pairwise ~ starting_substrate_new, type = "response")

write.csv(preds$contrasts, file = "emmeans(fi.d.add2, pairwise ~ starting_substrate_new.csv")


# Table S47 -----

fi.d.sl2 <- glmmTMB(log(diag_dist) ~ 
                      CORRECTED_fastest_peak_U *
                      avg_slope_angle_transect.1 +
                      (1|siteUn),
                    data = over.2,
                    family=Gamma (link ='log'))

ggplot() +
  geom_point(data=NULL, aes(y=resid(fi.d.sl2), x = fitted(fi.d.sl2)))

car::Anova(fi.d.sl2)

write.csv(car::Anova(fi.d.sl2), file = "car::Anova(fi.d.sl2).csv")
r.squaredGLMM(fi.d.sl2)

# Table S48 ----

fi.d.sl2.grid <- with(over.2, list(CORRECTED_fastest_peak_U = c(0.1, 0.2, 0.3, 0.4),
                                   avg_slope_angle_transect.1 = c(3, 13, 22)))

fi.d.sl2.preds2 <- emmeans(fi.d.sl2, pairwise ~ avg_slope_angle_transect.1 | CORRECTED_fastest_peak_U, 
                           at=fi.d.sl2.grid, type = 'response')

fi.d.sl2.preds2

write.csv(fi.d.sl2.preds2$emmeans, file = "fi.d.sl2.preds2.csv")
