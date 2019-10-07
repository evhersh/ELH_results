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


############
# Load data
############

S.dat <-read.csv(file.choose(), stringsAsFactors = FALSE, strip.white = TRUE, na.string = c("NA",""))

############
# Data prep
############

# Sets factors
S.dat$source <- as.factor(S.dat$source)
S.dat$ms <- as.factor(S.dat$ms)
S.dat$pop <- as.factor(S.dat$pop)
S.dat$garden <- as.factor(S.dat$garden)

# Sets levels, and sorts from south to north
S.dat$pop <- factor(S.dat$pop, levels=c("B53", "B42", "B46", "B49", "L11", "L12", "L06", "L16", "L17", "C86", "C85", "C88", "C27"))
S.dat$garden <- factor(S.dat$garden, levels = rev(levels(factor(S.dat$garden)))) # reverse default order to list gardens from South to North
S.dat$ms <- factor(S.dat$ms, levels = rev(levels(factor(S.dat$ms)))) # reverse default order to list sexuals first

# sets seed mass to be a numeric variable
S.dat$seed.mass <- as.numeric(S.dat$seed.mass)

# create new variable for mass per seed
S.dat$mass.per.seed <- S.dat$seed.mass / S.dat$number.weighed

#####################
#### Subset data ####
#####################

# filter data to only keep garden seeds, put in a new data frame
S.dat.garden <- S.dat %>%
  filter(source=="garden")

# same as above for natural seeds
S.dat.natural <- S.dat %>%
  filter(source=="natural")

# calculate means with standard errors for mating system and source (garden and natural)
S.dat.mass.means.ms <- S.dat %>%
  filter(!is.na(seed.mass)) %>%
  group_by(source, ms) %>%
  summarize(mean.mass=mean(seed.mass), se.mass=std.error(seed.mass), n=n())

# same as above but down to the population level
S.dat.mass.means.pop <- S.dat %>%
  filter(!is.na(seed.mass)) %>%
  group_by(source,ms,pop) %>%
  summarize(mean.mass=mean(seed.mass), se.mass=std.error(seed.mass), n=n())

##################
#### Models ######
##################

# Linear model of mating system with no random effects
S.mass.natural.lm <- lm(seed.mass~ms, data=S.dat.natural)
summary(S.mass.natural.lm)
anova(S.mass.natural.lm)

# linear model of population, no random effects
S.mass.natural.pop.lm <- lm(seed.mass~pop, data=S.dat.natural)
summary(S.mass.natural.pop.lm)
anova(S.mass.natural.pop.lm)

# linear model of mating system with population as a random effect
S.mass.natural.lmer <- lmer(seed.mass~ms+(1|pop), data=S.dat.natural)
summary(S.mass.natural.lmer)

# linear model of garden by mating system, no random effects 
S.mass.garden.lm <- lm(seed.mass~ms, data=S.dat.garden)
summary(S.mass.garden.lm)

#################
###### Plots ####
#################

# boxplot of natural populations
gg.mass.natural.pop <- ggplot(S.dat.natural, aes(x=pop, y=seed.mass, fill=ms))+
  geom_boxplot()+
  geom_jitter(alpha=1/3, width=0.1)+
  labs(y="Mean mass of seeds (natural pops)",x="population")

# boxplot of garden seeds by population
gg.mass.garden.pop <- ggplot(S.dat.garden, aes(x=pop, y=seed.mass, fill=ms))+
  geom_boxplot()+
  geom_jitter(alpha=1/3, width=0.1)+
  labs(y="Mean mass of seeds (gardens)", x="population")

# boxplot of garden seeds by garden
gg.mass.garden <- ggplot(S.dat.garden, aes(x=garden, y=seed.mass, fill=ms))+
  geom_boxplot()+
  labs(y="Mean mass of seeds (gardens)", x="garden")