############
# Packages
############
library(car)
library(ggplot2)
library(pwr)
library(Rmisc)
library(plyr)
library(visreg)
library(multcomp)
library(lme4)
library(nlme)
library(dplyr)
library(tidyr)
library(EML)
library(plotrix)
library(beanplot)
library(R2admb)
library(glmmADMB)
library(MCMCglmm)
library(ggthemes)
library(RColorBrewer)
library(pegas)
library(rJava)
library(reshape2)
library(xtable)
library(stargazer)
library(ggfortify)
library(AER)
library(glmm)
library(DHARMa)
library(lmerTest)
library(texreg)
library(RVAideMemoire)
library(MASS)
library(emmeans)
library(afex)
library(gridExtra)

library(extrafont)
#font_import()
loadfonts()

############
# Load data
############

S.dat <<-read.csv("~/GitHub/Hookeri-Gardens/Raw Data/seed_2018_Jamal.csv", stringsAsFactors = FALSE, strip.white = TRUE, na.string = c("NA",""))

############
# Data prep
############
S.dat$source <- as.factor(S.dat$source)
S.dat$ms <- as.factor(S.dat$ms)
S.dat$pop <- as.factor(S.dat$pop)
S.dat$garden <- as.factor(S.dat$garden)

S.dat$pop <- factor(S.dat$pop, levels=c("C27", "C88", "C85", "C86", "L17", "L16", "L06", "L12", "L11", "B49", "B46", "B42", "B53")) #sorted by lat, from high to low
S.dat$garden <- factor(S.dat$garden, levels = rev(levels(factor(S.dat$garden))))
S.dat$ms <- factor(S.dat$ms, levels = rev(levels(factor(S.dat$ms)))) #reverse order for figures (south to north)
S.dat$pop <- factor(S.dat$pop, levels = rev(levels(factor(S.dat$pop)))) #reverse order for figures (south to north)

S.dat$seed.mass <- as.numeric(S.dat$seed.mass)

S.dat$good.per.bud <- S.dat$good.seed / S.dat$bud.num

S.dat.garden <- S.dat %>%
  filter(source=="garden")

S.dat.natural <- S.dat %>%
  filter(source=="natural")

S.dat.mass.nat.ms <- S.dat %>%
  filter(source=="natural") %>%
  filter(!is.na(seed.mass)) %>%
  group_by(ms) %>%
  summarize(mean.mass=mean(seed.mass))

S.dat.mass.nat.pop <- S.dat %>%
  filter(source=="natural") %>%
  filter(!is.na(seed.mass)) %>%
  group_by(ms,pop) %>%
  summarize(mean.mass=mean(seed.mass), se.mass=std.error(seed.mass))


############
# Models
############

S.prop.all.glm <- glm(good.ratio~ms, weights= bud.num, data=S.dat, family="binomial")
summary(S.prop.all.glm)

S.prop.garden.glmer <- glm(good.ratio~ms+(1|pop), weights= bud.num, data=S.dat.garden, family="binomial")
summary(S.prop.garden.glm)

S.prop.garden.glmer <- glmer(good.ratio~ms +(1|pop), weights= bud.num, data=S.dat.garden, family="binomial")
summary(S.prop.garden.glmer)

S.prop.natural.glm <- glm(good.ratio~ms, weights= bud.num, data=S.dat.natural, family="binomial")
summary(S.prop.natural.glm)

S.prop.natural.glmer <- glmer(good.ratio~ms+(1|pop), weights= bud.num, data=S.dat.natural, family="binomial")
summary(S.prop.natural.glmer)

S.all.garden.lm <- lm(seed.per.bud~ms, data=S.dat.garden)
summary(S.all.garden.lm)

S.all.natural.lm <- lm(seed.per.bud~ms, data=S.dat.natural)
summary(S.all.natural.lm)

S.good.natural.lmer <- lmer(good.per.bud~ms+(1|pop), data=S.dat.natural)
summary(S.good.natural.lmer)
anova(S.good.natural.lmer)

S.good.garden.lmer <- lmer(good.per.bud~ms+(1|pop), data=S.dat.garden)
summary(S.good.garden.lmer)

S.mass.natural.lm <- lm(seed.mass~ms, data=S.dat.natural)
summary(S.mass.natural.lm)


############
# Plots
############

gg.good.both.pop <- ggplot(S.dat, aes(x=pop, y=good.ratio, fill=pop))+
  geom_boxplot()

gg.good.both.ms <- ggplot(S.dat, aes(x=ms, y=good.ratio))+
  geom_boxplot(width=0.5)

gg.good.natural.pop <- ggplot(S.dat.natural, aes(x=pop, y=good.ratio, fill=ms))+
  geom_boxplot()+
  geom_jitter(alpha=1/3, width=0.1)+
  labs(y="proportion of good seeds produced", x="population")+
  ggtitle("Proportion of good seeds produced \n in natural populations in 2018")

gg.gpb.natural.pop <- ggplot(S.dat.natural, aes(x=pop, y=good.per.bud, fill=ms))+
  geom_boxplot()+
  geom_jitter(alpha=1/3, width=0.1)+
  labs(y="mean number of good seeds per bud", x="population")+
  ggtitle("Mean number of good seeds per bud produced \n in natural populations in 2018")

gg.mass.natural.pop <- ggplot(S.dat.natural, aes(x=pop, y=seed.mass, fill=ms))+
  geom_boxplot()+
  geom_jitter(alpha=1/3, width=0.1)+
  labs(y="Mean mass of seeds (natural pops)",x="population")

gg.mass.garden.pop <- ggplot(S.dat.garden, aes(x=pop, y=seed.mass, fill=ms))+
  geom_boxplot()+
  geom_jitter(alpha=1/3, width=0.1)+
  labs(y="Mean mass of seeds (gardens)", x="population")

gg.good.natural.ms <- ggplot(S.dat.natural, aes(x=ms, y=good.ratio, fill=ms))+
  geom_boxplot(width=0.5)

gg.good.garden.pop <- ggplot(S.dat.garden, aes(x=pop, y=good.ratio, fill=pop))+
  geom_boxplot()+
  geom_jitter(alpha=1/3, width=0.1)+
  labs(y="proportion of good seeds produced", x="population")+
  ggtitle("Proportion of good seeds produced \n in common gardens in 2018")

gg.good.garden.ms <- ggplot(S.dat.garden, aes(x=ms, y=good.ratio, fill=ms))+
  geom_boxplot(width=0.5)

gg.all.garden.ms <- ggplot(S.dat.garden, aes(x=ms, y=seed.per.bud))+
  geom_boxplot()

gg.all.garden.pop <- ggplot(S.dat.garden, aes(x=pop, y=seed.per.bud))+
  geom_boxplot()

gg.all.natural.pop <- ggplot(S.dat.natural, aes(x=pop, y=seed.per.bud))+
  geom_boxplot()

gg.all.natural.ms <- ggplot(S.dat.natural, aes(x=ms, y=seed.per.bud))+
  geom_boxplot()

ggplot(S.dat.natural, aes(x=ms, y=good.per.bud))+
  geom_boxplot()

ggplot(S.dat.garden, aes(x=ms, y=good.per.bud))+
  geom_boxplot()
