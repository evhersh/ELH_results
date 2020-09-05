############
# Packages
############
library(tidyverse)
library(lme4)
library(ggeffects)
library(ggpmisc)
library(gridExtra)
library(see)
library(plotrix)
library(lmtest)
library(car)
library(ggpubr)
library(lmerTest)

##### DIASPORE #####

### Terminal velocity ###

# ms
Ber.tv <- lmerTest::lmer(terminal.velocity~ms + (1|pop/mom), data=Ber.dat, REML=FALSE)
Ber.tv2 <- lmerTest::lmer(terminal.velocity~  (1|pop/mom), data=Ber.dat, REML=FALSE)
lrtest(Ber.tv, Ber.tv2) # differ
anova(Ber.tv, Ber.tv2)
anova(Ber.tv)
ranova(Ber.tv)

# pop
Ber.tv.pop <- lmer(terminal.velocity ~ pop + (1|mom), data=Ber.dat)
Ber.tv.pop2 <- lmer(terminal.velocity ~ (1|mom), data=Ber.dat)
lrtest(Ber.tv.pop, Ber.tv.pop2) # differ
anova(Ber.tv.pop, Ber.tv.pop2) # differ


### TV regressions ###
tv.aot <- lm(terminal.velocity ~ angle.of.attack.scaled, data=Ber.dat)
summary(tv.aot)
tv.bl <- lm(terminal.velocity ~ bristle.length.scaled, data=Ber.dat)
summary(tv.bl)
tv.nob <- lm(terminal.velocity ~ number.of.bristles.scaled, data=Ber.dat)
summary(tv.nob)
tv.weight <- lm(terminal.velocity ~ weight.scaled, data=Ber.dat)
summary(tv.weight)

### other diaspore traits ###

# angle of attack
Ber.aot <- lmer(angle.of.attack~ms + (1|pop/mom), data=Ber.dat)
Ber.aot2 <- lmer(angle.of.attack~ (1|pop/mom), data=Ber.dat)
lrtest(Ber.aot, Ber.aot2) # differ
anova(Ber.aot, Ber.aot2) # differ
plot(aot.ggfx)

# bristle length
Ber.bl <- lmer(bristle.length~ms + (1|pop/mom), data=Ber.dat)
Ber.bl2 <- lmer(bristle.length~(1|pop/mom), data=Ber.dat)
lrtest(Ber.bl, Ber.bl2) # differ
anova(Ber.bl, Ber.bl2)
plot(bl.ggfx)

# number of bristles
Ber.nob <- glmer(number.of.bristles~ms + (1|pop/mom), family=poisson(link="log"), data=Ber.dat)
Ber.nob2 <- glmer(number.of.bristles~ (1|pop/mom), family=poisson(link="log"), data=Ber.dat)
lrtest(Ber.nob, Ber.nob2) # differ
anova(Ber.nob, Ber.nob2)
plot(nob.ggfx) # large overlap

# weight
Ber.weight <- lmer(weight~ ms + (1|pop/mom), data=Ber.dat)
Ber.weight2 <- lmer(weight~ (1|pop/mom), data=Ber.dat)
lrtest(Ber.weight, Ber.weight2) # differ?
anova(Ber.weight, Ber.weight2)
plot(weight.ggfx) # large overlap

##### LAB STUFF #####

### GERM SUCCESS ###

#ms
L.germ.ms.glmer <- glmer(germ ~ ms + (1|pop/mom)+(1|plate), family=binomial(link="logit"), data=L.dat)
L.germ.ms.glmer2 <- glmer(germ ~ (1|pop/mom)+(1|plate), family=binomial(link="logit"), data=L.dat)
lrtest(L.germ.ms.glmer, L.germ.ms.glmer2)
anova(L.germ.ms.glmer, L.germ.ms.glmer2)
anova(L.germ.ms.glmer)

#pop
L.germ.pop.glmer <- glmer(germ ~ pop + (1|mom)+(1|plate), family=binomial(link="logit"), data=L.dat)
L.germ.pop.glmer2 <- glmer(germ ~ (1|mom)+(1|plate), family=binomial(link="logit"), data=L.dat,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
lrtest(L.germ.pop.glmer, L.germ.pop.glmer2)
anova(L.germ.pop.glmer, L.germ.pop.glmer2)

### GERM SPEED ###

#ms
L.days.ms.glmer <- glmer(days ~ ms + (1|pop/mom)+(1|plate), family=poisson(), data=L.dat.days)
L.days.ms.glmer2 <- glmer(days ~ (1|pop/mom)+(1|plate), family=poisson(), data=L.dat.days)
lrtest(L.days.ms.glmer, L.days.ms.glmer2)
anova(L.days.ms.glmer, L.days.ms.glmer2)


#pop
L.days.pop.glmer <- glmer(days ~ pop + (1|mom)+(1|plate), family=poisson(), data=L.dat.days,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
L.days.pop.glmer2 <- glmer(days ~ (1|mom)+(1|plate), family=poisson(), data=L.dat.days, control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
lrtest(L.days.pop.glmer, L.days.pop.glmer2)
anova(L.days.pop.glmer, L.days.pop.glmer2)


### SEEDLING SURVIVAL ###

#ms
L.surv.ms.glmer <- glmer(surv ~ ms + (1|pop/mom)+(1|rack), family=binomial(link="logit"), data=L.surv)
L.surv.ms.glmer2 <- glmer(surv ~ (1|pop/mom)+(1|rack), family=binomial(link="logit"), data=L.surv)
lrtest(L.surv.ms.glmer, L.surv.ms.glmer.2)
anova(L.surv.ms.glmer, L.surv.ms.glmer.2)


#pop
L.surv.pop.glmer <- glmer(surv ~ pop + (1|mom)+(1|rack), family=binomial(link="logit"), data=L.surv,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
L.surv.pop.glmer2 <- glmer(surv ~ (1|mom)+(1|rack), family=binomial(link="logit"), data=L.surv,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
lrtest(L.surv.pop.glmer, L.surv.pop.glmer2)
anova(L.surv.pop.glmer, L.surv.pop.glmer2)


### LEAF NUM ###

#ms
H.y0.num.ms.glmer <- glmer(leaf.num ~ ms + (1|pop/mom), data=H.y0, family=poisson(link="log"))
H.y0.num.ms.glmer2 <- glmer(leaf.num ~ (1|pop/mom), data=H.y0, family=poisson(link="log"))
lrtest(H.y0.num.ms.glmer, H.y0.num.ms.glmer2)
anova(H.y0.num.ms.glmer, H.y0.num.ms.glmer2)


#pop
H.y0.num.pop.glmer <- glmer(leaf.num ~ pop + (1|mom), data=H.y0, family=poisson(link="log"),control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
H.y0.num.pop.glmer2 <- glmer(leaf.num ~ (1|mom), data=H.y0, family=poisson(link="log"),control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
lrtest(H.y0.num.pop.glmer,H.y0.num.pop.glmer2)
anova(H.y0.num.pop.glmer,H.y0.num.pop.glmer2)


### LEAF LENGTH ###

#ms
H.y0.length.ms.lmer <- lmer(leaf.length ~ ms + (1|pop/mom), data=H.y0)
H.y0.length.ms.lmer2 <- lmer(leaf.length ~ (1|pop/mom), data=H.y0)
lrtest(H.y0.length.ms.lmer,H.y0.length.ms.lmer2)
anova(H.y0.length.ms.lmer,H.y0.length.ms.lmer2)


#pop
H.y0.length.pop.lmer <- lmer(leaf.length ~ pop + (1|mom), data=H.y0)
H.y0.length.pop.lmer2 <- lmer(leaf.length ~ (1|mom), data=H.y0)
lrtest(H.y0.length.pop.lmer, H.y0.length.pop.lmer2)
anova(H.y0.length.pop.lmer, H.y0.length.pop.lmer2)


