library(tidyverse)
library(lme4)
library(refund)
library(future)
library(future.apply)
library(reshape)
library(viridis)
library(FunInf)

#--------------------------------------------------------------------------------
# 1. Read in Data
#--------------------------------------------------------------------------------
data = readRDS("Demo Data.rds")

cols = rep(sample(turbo(20)), each = 5)
par(mfrow = c(1,2))
matplot(data$System1, type = 'l', main = "System 1", ylab = "Degrees Flexion", xlab = "% Gait Cycle", lty = 1, col = cols)
matplot(data$System2, type = 'l', main = "System 2", ylab = "Degrees Flexion", xlab = "% Gait Cycle", lty = 1, col = cols)

#--------------------------------------------------------------------------------
# 2. Discrete Limit of Agreement
#--------------------------------------------------------------------------------

# set up data
mod_data = data$Reference %>% mutate(Participant = factor(Participant), Trial = factor(Trial))

# calculate peak knee flexions
mod_data$PKF1 = apply(data$System1, 2, max)
mod_data$PKF2 = apply(data$System2, 2, max)

# compute difference
mod_data$Diff = mod_data$PKF1-mod_data$PKF2

# model expected difference
mod = lmer(Diff ~ 1 + (1|Participant), data = mod_data)

# extract the bias estimate (the fixed intercept)
bias = fixef(mod)

# get the model error variance for construction of the LoA
sigE = sigma(mod)

# add to original data frame for plotting
mod_data$bias = bias

# create the LoA based on a Normal distribution
mod_data$LoA_lower = bias - 1.96*sigE
mod_data$LoA_upper = bias + 1.96*sigE

#plot
ggplot(mod_data)+
  geom_rect(aes(ymin = LoA_lower, ymax = LoA_upper, xmin = -Inf,   xmax = Inf), fill = "lightgrey", alpha = 0.1)+
  geom_jitter(aes(y = Diff, x= 1, colour = Participant))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_hline(aes(yintercept = bias), colour = "red")+
  theme_light()+
  scale_x_continuous(labels = NULL)+
  labs(x = " ", y = "Degrees Difference (System 1 - System 2)",
       title = "Peak Knee Flexion - Limits of Agreement",
       subtitle = paste0("Bias = ", round(bias,2), ", 95% LoA = [", round(bias-1.96*sigE, 2),",",round(bias+1.96*sigE, 2),"]"))+
  theme(legend.position = 'none',
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x =  element_blank(),
        plot.title = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))

#--------------------------------------------------------------------------------
# 1. Functional Limit of Agreement
#--------------------------------------------------------------------------------

# set up functional data for modelling
mod_data$Sys1 = t(data$System1)
mod_data$Sys2 = t(data$System2)

# calculate difference
mod_data$funDif = t(data$System1)-t(data$System2)

# estimate the functional random effects model
t = seq(0,1,0.01)
mod2 = pffr(funDif ~ 1 + s(Participant, k = 20, bs = 're'), data = mod_data,
            bs.int = list(k = 20 , bs = 'ps', m = 2), yind = t, algorithm = 'bam')

# extract functional bias (functional fixed intercept)
bias = coef(mod2)$`smterms`[[1]]$coef

# the way pffr models there is constant intercept that needs to be added in as well
bias$value = bias$value + coef(mod2)$`pterms`[1]

# extract the residuals
resids = matrix(mod2$residuals, 101, 100)

# estimate pointwise error variance via fpca
sigEt = apply(resids, 1, sd)

# get simultaneous confidence bands (Degras 2017)
set.seed(2026)
boot.list = level.boot.list(df = mod_data, id = "Participant", B = 1000)

# note: this step takes time and RAM
fits = fit.boots(model = mod2, boot.list = boot.list, ncores = 10)
gc()
CMA = estimate.CMA(boot.grid = fits[[1]]$boots, est = bias$value)

cor.fac = CMA$adj.fac[1]


# revaluate at same points as the bias
sigEt = spline(x = t, y = sigEt, xout = bias$yindex.vec)$y

bias = bias %>% mutate(LoA_upper = value + cor.fac*sigEt, LoA_lower = value - cor.fac*sigEt,
                       LoApwL = value - 1.96*sigEt, LoApwU = value + 1.96*sigEt, )

# rearrange raw difference data for plotting
diff = data.frame(melt(data$System1-data$System2)) %>% mutate(Participant = factor(rep(mod_data$Participant, each = 101)),
                                                              Trial = factor(rep(1:100, each = 101)),
                                                              t = rep(t, 100)) %>% select(-c(X1, X2))

ggplot()+
  geom_line(data = diff, aes(x = t, y = value, group = Trial, colour = Participant), alpha = 0.3)+
  geom_ribbon(data = bias, aes(x = yindex.vec, ymin = LoA_lower, ymax = LoA_upper), fill = "lightgrey", alpha = 0.8)+
  geom_ribbon(data = bias, aes(x = yindex.vec, ymin = LoApwL, ymax = LoApwU), colour = "darkgrey", alpha = 0, lty = 2)+
  geom_line(data = bias, aes(x = yindex.vec, y = value), colour = "red")+
  geom_hline(yintercept = 0, lty = 2)+
  labs(x = "% Gate Cyle", y = "Degrees Difference (System 1 - System 2)",
       title = "Knee Flexion - Limits of Agreement",
       subtitle = "Bias (red line), Point wise LoA (dashed grey), Simultaneous LoA (shaded grey)")+
  theme_light()+
  theme(legend.position = 'none',
        plot.title = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"))
