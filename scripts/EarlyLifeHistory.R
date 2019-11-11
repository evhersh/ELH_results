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
library(multcompView)

library(extrafont)
#font_import()
#loadfonts()

# check for overdispersion in glmm's
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

############
# Load data
############

L.dat <<-read.csv("~/GitHub/Hookeri-Gardens/RAW DATA/mastergerm_lab.csv", stringsAsFactors = FALSE, strip.white = TRUE, na.string = c("NA",""))
L.surv <<-read.csv("~/GitHub/Hookeri-Gardens/RAW DATA/labseedlingsurvival.csv", stringsAsFactors = FALSE, strip.white = TRUE, na.string = c("NA",""))
Berto.dat <<-read.csv("~/GitHub/Hookeri-Gardens/RAW DATA/alberto.seed.R.csv", stringsAsFactors = FALSE, strip.white = TRUE, na.string = c("NA",""))
H.dat <<-read.csv("~/GitHub/Hookeri-Gardens/RAW DATA/greenhouse.csv", stringsAsFactors = FALSE, strip.white = TRUE, na.string = c("NA",""))
G.dat <<-read.csv("~/GitHub/Hookeri-Gardens/RAW DATA/mastergermv3.csv", stringsAsFactors = FALSE, strip.white = TRUE, na.string = c("NA",""))




##########
# Data Prep
##########
pd <- position_dodge(0.5)

  # GARDEN DATA

# calculate size
H.dat$size <- H.dat$leaf.num * H.dat$leaf.length
H.dat$log.size <- log(H.dat$size)
H.dat$log.length <- log(H.dat$leaf.length)

# Garden data
H.dat$pop <- as.factor(H.dat$pop)
H.dat$ms <- as.factor(H.dat$ms)
H.dat$s.region <- as.factor(H.dat$s.region)
H.dat$mom <- as.factor(H.dat$mom)


H.dat$year <- as.factor(H.dat$year)
H.dat$leaf.num <- as.numeric(H.dat$leaf.num)
H.dat$leaf.length <- as.numeric(H.dat$leaf.length)
H.dat$size <- as.numeric(H.dat$size)
H.dat$log.size <- as.numeric(H.dat$log.size)




#set levels
H.dat$s.region <- factor(H.dat$s.region, levels=c("A.s", "AO.s", "SO.s", "S.s"))
H.dat$pop <- factor(H.dat$pop, levels=c("C27", "C85", "C86", "L17", "L16", "L06", "L12", "L11", "B49", "B46", "B42", "B53")) #sorted by lat, from high to low

H.dat$s.region <- factor(H.dat$s.region, levels = rev(levels(factor(H.dat$s.region)))) #reverse order for figures (south to north)
H.dat$ms <- factor(H.dat$ms, levels = rev(levels(factor(H.dat$ms)))) #reverse order for figures (south to north)
H.dat$pop <- factor(H.dat$pop, levels = rev(levels(factor(H.dat$pop)))) #reverse order for figures (south to north)



#subsetting

H.y0 <- H.dat %>%
  filter(year==0) %>%
  filter(!is.na(log.size)) %>%
  filter(log.size!=0)

H.y0.mom <- H.dat %>%
  filter(year==0) %>%
  filter(!is.na(log.size)) %>%
  filter(log.size!=0) %>%
  group_by(ms, pop, mom) %>%
  summarize(mean.length=mean(leaf.length), mean.num=mean(leaf.num))

H.y0.length.mom <- H.dat %>%
  filter(year==0) %>%
  filter(!is.na(leaf.length)) %>%
  filter(leaf.length!=0) %>%
  group_by(ms, pop, mom) %>%
  summarize(mean.length=mean(leaf.length), se.length=std.error(leaf.length))

H.y0.mean.length.s.region <- H.dat %>%
  filter(year==0) %>%
  filter(!is.na(leaf.length)) %>% #remove NAs from leaf.num
  group_by(ms, s.region) %>%
  summarize(mean.length=mean(leaf.length),
            se.length=std.error(leaf.length))

H.y0.mean.length.pop <- H.dat %>%
  filter(year==0) %>%
  filter(!is.na(leaf.length)) %>% #remove NAs from leaf.num
  group_by(ms, s.region, pop) %>%
  summarize(mean.length=mean(leaf.length),
            se.size=std.error(leaf.length), n())

H.y0.mean.length.ms <- H.dat %>%
  filter(year==0) %>%
  filter(!is.na(leaf.length)) %>% #remove NAs from leaf.num
  group_by(ms) %>%
  summarize(mean.length=mean(leaf.length),
            se.size=std.error(leaf.length), n())

H.y0.mean.num.pop <- H.dat %>%
  filter(year==0) %>%
  filter(!is.na(leaf.num)) %>% #remove NAs from leaf.num
  group_by(ms, s.region, pop) %>%
  summarize(mean.length=mean(leaf.num),
            se.size=std.error(leaf.num), n())

H.y0.mean.num.ms <- H.dat %>%
  filter(year==0) %>%
  filter(!is.na(leaf.num)) %>% #remove NAs from leaf.num
  group_by(ms) %>%
  summarize(mean.length=mean(leaf.num),
            se.size=std.error(leaf.num), n())

H.y1 <- H.dat %>%
  filter(year==1) %>%
  filter(!is.na(leaf.length)) %>%
  filter(leaf.length!=0)

H.y1.mean.length.s.region <- H.dat %>%
  filter(year==1) %>%
  filter(!is.na(leaf.length)) %>% #remove NAs from leaf.num
  group_by(ms, s.region) %>%
  summarize(mean.length=mean(leaf.length),
            se.length=std.error(leaf.length))

H.y1.mean.length.pop <- H.dat %>%
  filter(year==1) %>%
  filter(!is.na(leaf.length)) %>% #remove NAs from leaf.num
  group_by(ms, s.region, pop) %>%
  summarize(mean.length=mean(leaf.length),
            se.size=std.error(leaf.length))

H.y0.length.to.surv.y1 <- H.dat %>%
  filter(year==0) %>%
  filter(!is.na(surv.to.y1)) %>%
  filter(!is.na(leaf.length)) %>%
  group_by(surv.to.y1) %>%
  summarize(mean.length=mean(leaf.length), se.length=std.error(leaf.length))

H.y0.length.to.surv.y1.ms <- H.dat %>%
  filter(year==0) %>%
  filter(!is.na(surv.to.y1)) %>%
  filter(!is.na(leaf.length)) %>%
  group_by(ms, surv.to.y1) %>%
  summarize(mean.length=mean(leaf.length), se.length=std.error(leaf.length))


  # LAB GERMINATION


# factors
L.dat$pop <- as.factor(L.dat$pop)
L.dat$ms <- as.factor(L.dat$ms)
L.dat$mom <- as.factor(L.dat$mom)
L.dat$plate <- as.factor(L.dat$plate)

L.surv$pop <- as.factor(L.surv$pop)
L.surv$ms <- as.factor(L.surv$ms)
L.surv$s.region <- as.factor(L.surv$s.region)
L.surv$rep <- as.factor(L.surv$rep)
L.surv$mom <- as.factor(L.surv$mom)
L.surv$rack <- factor(L.surv$rack)

# levels

L.dat$ms <- factor(L.dat$ms, levels=c("A", "S"))
L.dat$pop <- factor(L.dat$pop, levels=c("C27", "C85", "C86", "L17", "L16", "L06", "L12", "L11", "B49", "B46", "B42", "B53"))
L.dat$ms <- factor(L.dat$ms, levels = rev(levels(factor(L.dat$ms)))) #reverse order for figures (south to north)
L.dat$pop <- factor(L.dat$pop, levels = rev(levels(factor(L.dat$pop)))) #reverse order for figures (south to north)

L.surv$s.region <- factor(L.surv$s.region, levels=c("A.s", "AO.s", "SO.s", "S.s"))
L.surv$ms <- factor(L.surv$ms, levels=c("A", "S"))
L.surv$pop <- factor(L.surv$pop, levels=c("C27", "C85", "C86", "L17", "L16", "L06", "L12", "L11", "B49", "B46", "B42", "B53"))
L.surv$s.region <- factor(L.surv$s.region, levels = rev(levels(factor(L.surv$s.region)))) #reverse order for figures (south to north)
L.surv$ms <- factor(L.surv$ms, levels = rev(levels(factor(L.surv$ms)))) #reverse order for figures (south to north)
L.surv$pop <- factor(L.surv$pop, levels = rev(levels(factor(L.surv$pop)))) #reverse order for figures (south to north)



  # Means, subsets, etc.

L.dat.days <- L.dat %>%
  filter(!is.na(days))

L.germ <- L.dat %>%
  filter(!is.na(germ)) %>%
  group_by(ms) %>%
  summarize(mean.germ=mean(germ), se.germ=std.error(germ), n())

L.germ.pop <- L.dat %>%
  filter(!is.na(germ)) %>%
  group_by(ms, pop) %>%
  summarize(mean.germ=mean(germ), se.germ=std.error(germ), n())

L.germ.mom <- L.dat %>%
  filter(!is.na(germ)) %>%
  group_by(ms, pop, mom) %>%
  summarize(mean.germ=mean(germ), se.germ=std.error(germ))

L.germ.ms <- L.dat %>%
  filter(!is.na(germ)) %>%
  group_by(ms) %>%
  summarize(mean.germ=mean(germ), se.germ=std.error(germ))

L.days <- L.dat %>%
  filter(!is.na(days)) %>%
  group_by(ms) %>%
  summarize(mean.days=mean(days), se.days=std.error(days), n())

L.days1 <- L.dat %>%
  filter(!is.na(days)) %>%
  group_by(ms, pop) %>%
  summarize(mean.days=mean(days), se.days=std.error(days))

L.days.ms <- L.dat %>%
  filter(!is.na(days)) %>%
  group_by(ms) %>%
  summarize(mean.days=mean(days), se.days=std.error(days))

L.days.pop <- L.dat %>%
  filter(!is.na(days)) %>%
  group_by(ms, pop) %>%
  summarize(mean.days=mean(days), se.days=std.error(days), n())

L.days.mom <- L.dat %>%
  filter(!is.na(days)) %>%
  group_by(ms, pop, mom) %>%
  summarize(mean.days=mean(days), se.days=std.error(days))


L.surv.ms <- L.surv %>%
  filter(!is.na(surv)) %>%
  group_by(ms) %>%
  summarize(mean.surv=mean(surv), se.surv=std.error(surv), n())

L.surv.pop <- L.surv %>%
  filter(!is.na(surv)) %>%
  group_by(ms, pop) %>%
  summarize(mean.surv=mean(surv), se.surv=std.error(surv), n())

L.surv.mom <- L.surv %>%
  filter(!is.na(surv)) %>%
  group_by(ms, pop, mom) %>%
  summarize(mean.surv=mean(surv), se.surv=std.error(surv))

L.y0.days.length <- left_join(L.days.mom, H.y0.length.mom)

L.surv.to.seedling.pop <- L.surv %>%
  group_by(ms, pop) %>%
  summarize(mean.surv=mean(surv.to.seedling), se.surv=std.error(surv.to.seedling))

L.surv.to.seedling.mom <- L.surv %>%
  group_by(ms, pop,mom) %>%
  summarize(mean.surv=mean(surv.to.seedling), se.surv=std.error(surv.to.seedling))

L.surv.to.seedling.ms <- L.surv %>%
  group_by(ms) %>%
  summarize(mean.surv=mean(surv.to.seedling), se.surv=std.error(surv.to.seedling))


  # FIELD GERMINATION #

#create factors
G.dat$pop <- as.factor(G.dat$pop)
G.dat$ms <- as.factor(G.dat$ms)
G.dat$s.region <- as.factor(G.dat$s.region)
G.dat$garden <- as.factor(G.dat$garden)
G.dat$g.region <- as.factor(G.dat$g.region)
G.dat$year <- as.numeric(G.dat$year)
G.dat$mom <- as.factor(G.dat$mom)


#set levels
G.dat$s.region <- factor(G.dat$s.region, levels=c("A.s", "AO.s", "SO.s", "S.s"))
G.dat$g.region <- factor(G.dat$g.region, levels=c("A.g", "AO.g", "SO.g", "S.g"))
G.dat$pop <- factor(G.dat$pop, levels=c("C27", "C85", "C86", "L17", "L16", "L06", "L12", "L11", "B49", "B46", "B42", "B53"))

# reverse levels for figures (south to north)
G.dat$s.region <- factor(G.dat$s.region, levels = rev(levels(factor(G.dat$s.region)))) #reverse order for figures (south to north)
G.dat$ms <- factor(G.dat$ms, levels = rev(levels(factor(G.dat$ms)))) #reverse order for figures (south to north)
G.dat$pop <- factor(G.dat$pop, levels = rev(levels(factor(G.dat$pop)))) #reverse order for figures (south to north)
G.dat$g.region <- factor(G.dat$g.region, levels = rev(levels(factor(G.dat$g.region)))) #reverse order for figures (south to north)
G.dat$garden <- factor(G.dat$garden, levels = rev(levels(factor(G.dat$garden)))) #reverse order for figures (south to north)


# Separate sex and asex
G.asex <- G.dat %>%
  filter(ms=="A")

G.sex <- G.dat %>%
  filter(ms=="S")

G.dat.y2 <- G.dat %>%
  filter(!is.na(germ)) %>%
  filter(year==2)

G.dat.y2.mom <- G.dat %>%
  filter(!is.na(germ)) %>%
  group_by(ms,pop,mom) %>%
  summarize(mean.germ=mean(germ))

G.dat.y2.filt <- G.dat %>%
  filter(!is.na(germ)) %>%
  filter(year==2) %>%
  filter(garden!="SS2")

G.dat.S.g <- G.dat %>%
  filter(!is.na(germ)) %>%
  filter(year==2) %>%
  filter(g.region=="S.g")

# subset different means
G.asex.reg <- G.asex %>%
  filter(!is.na(germ)) %>%
  group_by(s.region, g.region) %>%
  summarize(mean.germ=mean(germ), se.germ=std.error(germ))

G.asex.reg2 <- G.asex %>%
  filter(!is.na(germ)) %>%
  group_by(g.region) %>%
  summarize(mean.germ=mean(germ), se.germ=std.error(germ))

G.sex.reg <- G.sex %>%
  filter(!is.na(germ)) %>%
  group_by(s.region, g.region) %>%
  summarize(mean.germ=mean(germ), se.germ=std.error(germ))

G.sex.reg2 <- G.sex %>%
  filter(!is.na(germ)) %>%
  group_by(g.region) %>%
  summarize(mean.germ=mean(germ), se.germ=std.error(germ))

G.all1 <- G.dat %>%
  filter(!is.na(germ)) %>%
  group_by(s.region, g.region) %>%
  summarize(mean.germ=mean(germ), se.germ=std.error(germ))

G.all2 <- G.dat %>%
  filter(!is.na(germ)) %>%
  group_by(g.region) %>%
  summarize(mean.germ=mean(germ), se.germ=std.error(germ))

G.all.ms <- G.dat %>%
  filter(!is.na(germ)) %>%
  group_by(ms, s.region) %>%
  summarize(mean.germ=mean(germ), se.germ=std.error(germ))


G.ms <- G.dat %>%
  filter(!is.na(germ)) %>%
  group_by(ms) %>%
  summarize(mean.germ=mean(germ), se.germ=std.error(germ))

G.all.garden <- G.dat %>%
  filter(!is.na(germ)) %>%
  group_by(ms, garden) %>%
  summarize(mean.germ=mean(germ), se.germ=std.error(germ))


G.all.pop <- G.dat %>%
  filter(!is.na(germ)) %>%
  group_by(ms,s.region,pop) %>%
  summarize(mean.germ=mean(germ), se.germ=std.error(germ))

G.pop.garden <- G.dat %>%
  filter(!is.na(germ)) %>%
  group_by(ms,garden,pop) %>%
  summarize(mean.germ=mean(germ), se.germ=std.error(germ))

G.garden.means <- G.dat %>%
  filter(!is.na(germ)) %>%
  group_by(garden) %>%
  summarize(mean.germ=mean(germ))

G.pop.garden.means <- G.dat %>%
  filter(!is.na(germ)) %>%
  filter(year==2) %>%
  group_by(garden) %>%
  mutate(mean.germ.garden=mean(germ)) %>%
  group_by(garden, pop) %>%
  mutate(mean.germ.pop=mean(germ)) %>%
  mutate(mean.diff.gp=mean.germ.pop-mean.germ.garden)

G.pop.garden.means2 <- G.pop.garden.means %>%
  group_by(garden,pop) %>%
  summarize(mean.germ.garden=mean(mean.germ.garden), mean.germ.pop=mean(mean.germ.pop), diffs=mean(mean.diff.gp))

G.pgm <- G.pop.garden.means %>%
  group_by(garden,pop) %>%
  summarize(diffs=mean(mean.diff.gp))
  
  
  
  # Berto Dispersal Data


# selecting just the main columns that I'm interested in
Ber.dat <- Berto.dat %>%
  filter(region!="WestWY") %>%
  filter(mom!="Bulk")
  

# factors and levels
Ber.dat$pop <- factor(Ber.dat$pop, levels=c("B42", "B42.alt", "B46", "B49", "L12", "C59", "B52", "L16", "L40", "L41", "C21", "C43", "C54"))
Ber.dat$region <- as.factor(Ber.dat$region)
Ber.dat$garden <- as.factor(Ber.dat$garden)
Ber.dat$ms <- factor(Ber.dat$ms)
Ber.dat$lat <- as.factor(Ber.dat$lat)

#Ber.dat$ms <- factor(Ber.dat$ms, levels = rev(levels(factor(Ber.dat$ms)))) #reverse order for figures (south to north)
#Ber.dat$pop <- factor(Ber.dat$pop, levels = rev(levels(factor(Ber.dat$pop)))) #reverse order for figures (south to north)


# recreating table of means for all the pops, sorted by ms and pop
Ber.moms.all <- Ber.dat %>%
  select(pop:number.of.bristles, lat) %>%
  group_by(ms,region,pop,mom,lat) %>%
  summarize_all(funs(mean)) %>%
  arrange(desc(ms), lat)

Ber.moms.garden <- Ber.dat %>% # Only one asex population from here was used in the gardens...
  filter(garden=="yes") %>%
  select(pop:number.of.bristles, lat) %>%
  group_by(pop, ms,region, mom) %>%
  summarize_all(funs(mean))

Ber.ms.all <- Ber.dat %>%
  select(pop:number.of.bristles, lat) %>%
  group_by(ms,region,pop,mom,lat) %>%
  summarize_all(funs(mean)) %>%
  arrange(desc(ms), lat)

# same table but just comparing means of sex and asex
Ber.means.ms <- Ber.dat %>%
  group_by(ms) %>%
  select(-pop) %>%
  summarize_all(funs(mean))

# Terminal velocity means
Ber.means.TV <- Ber.dat %>%
  group_by(ms) %>%
  summarize(mean.terminal.velocity=mean(terminal.velocity), se.terminal.velocity=std.error(terminal.velocity))

Ber.means.TV.pop <- Ber.dat %>%
  group_by(ms, pop) %>%
  summarize(mean.terminal.velocity=mean(terminal.velocity), se.terminal.velocity=std.error(terminal.velocity), n())

# log transform and separate variables
Ber.log <- log(Ber.dat[c(3,4,6,7)])
Ber.pops <- Ber.dat[, 1]

# correlation matrix
Ber.matrix <- cor(Ber.log)
Ber.matrix.4 <- round(Ber.matrix, 4)


##########
# MODELS
##########
  
  # For fixing convergence issues
# control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))

  
  # LAB GERM SUCCESS

# This one has the lowest AIC...
L.germ.ms.glmer <- glmer(germ ~ ms + (1|pop/mom)+(1|plate), family=binomial(link="logit"), data=L.dat)
L.germ.ms.afex <- afex::mixed(germ ~ ms + (1|pop/mom)+(1|plate), family=binomial(link="logit"), data=L.dat, method="LRT")

L.germ.pop.glmer <- glmer(germ ~ pop + (1|mom)+(1|plate), family=binomial(link="logit"), data=L.dat, control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
L.germ.pop.afex <- afex::mixed(germ ~ pop + (1|mom)+(1|plate), family=binomial(link="logit"), data=L.dat, method="LRT", control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))

L.germ.ms.glmer.nob53 <- glmer(germ ~ ms + (1|pop/mom)+(1|plate), family=binomial(link="logit"), data=subset(L.dat, pop!="B53"))
L.germ.pop.glmer.nob53 <- glmer(germ ~ pop + (1|mom)+(1|plate), family=binomial(link="logit"), data=subset(L.dat, pop!="B53"), control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))


summary(L.germ.ms.glmer.nob53)
summary(L.germ.pop.glmer.nob53)

L.germ.ms.glmer2 <- glmer(germ ~ ms + (1|pop/mom), family=binomial(link="logit"), data=L.dat)
L.germ.ms.glmer3 <- glmer(germ ~ ms + (1|pop), family=binomial(link="logit"), data=L.dat)
L.germ.ms.glmer4 <- glmer(germ ~ ms + (1|pop)+(1|plate), family=binomial(link="logit"), data=L.dat)
L.germ.ms.glmer5 <- glmer(germ ~ ms + (1|mom)+(1|plate), family=binomial(link="logit"), data=L.dat)
L.germ.ms.glmer6 <- glmer(germ ~ ms + (1|mom), family=binomial(link="logit"), data=L.dat)
L.germ.ms.glmer7 <- glmer(germ ~ ms + (1|plate), family=binomial(link="logit"), data=L.dat)

L.germ.ms.glmer.AIC <- AIC(L.germ.ms.glmer, L.germ.ms.glmer2, L.germ.ms.glmer3, L.germ.ms.glmer4, L.germ.ms.glmer5, L.germ.ms.glmer6, L.germ.ms.glmer7)


  # glmer ms comparisons
summary(L.germ.ms.glmer)
anova(L.germ.ms.afex) # this for factor level p-values
overdisp.glmer(L.germ.ms.glmer)
L.germ.ms.glmer.texreg <- texreg(L.germ.ms.glmer, float.pos = "h")
L.germ.ms.glmer.emmeans <- emmeans(L.germ.ms.glmer, c("ms"), data=L.dat)
pairs(L.germ.ms.glmer.emmeans)
germ.ms.clds<-cld(L.germ.ms.glmer.emmeans)
pwpp(L.germ.ms.glmer.emmeans)
plot(L.germ.ms.glmer.emmeans, comparisons = TRUE, horizontal=FALSE)+
  theme_classic()

stargazer(anova(L.germ.ms.afex), flip=TRUE)

# glmer pop comparisons
summary(L.germ.pop.glmer)
anova(L.germ.pop.afex) # this for factor level p-values
L.germ.pop.glmer.texreg <- texreg(L.germ.ms.glmer, float.pos = "h")
L.germ.pop.glmer.emmeans <- emmeans(L.germ.pop.glmer, c("pop"), data=L.dat)
pairs(L.germ.pop.glmer.emmeans)
CLD(L.germ.pop.glmer.emmeans)
L.germ.pop.clds <- CLD(L.germ.pop.glmer.emmeans, details=FALSE, sort=FALSE)
plot(L.germ.pop.glmer.emmeans, comparisons = TRUE, horizontal=FALSE)+
  theme_classic()

  # glm regs
L.germ.pop.glm <- glm(germ ~ pop+mom, family=binomial(link="logit"), data=L.dat)
summary(L.germ.pop.glm)
L.germ.pop.glm.texreg <- texreg(L.germ.pop.glm, float.pos="h")
emmeans(L.germ.pop.glm, "pop")  



 # LAB GERM DAYS
L.days.ms.glm <- glm(days~ms, family=quasipoisson(link="log"), data=L.dat.days)  
summary(L.days.ms.glm)

  # GLMM for MS
L.days.ms.glmer <- glmer(days ~ ms + (1|pop/mom)+(1|plate), family=poisson(), data=L.dat.days)
summary(L.days.ms.glmer)
overdisp.glmer(L.days.ms.glmer)
drop1(L.days.ms.glmer)

L.days.ms.glmer.emmeans <- emmeans(L.days.ms.glmer, c("ms"), data=L.dat)
pairs(L.days.ms.glmer.emmeans)
days.ms.clds<-CLD(L.days.ms.glmer.emmeans)
plot(L.days.ms.glmer.emmeans, comparisons = TRUE)

  # afex for ms
L.days.ms.afex <- afex::mixed(days ~ ms + (1|pop/mom)+(1|plate), family=poisson(), data=L.dat.days, method="LRT")
summary(L.days.ms.afex)
anova(L.days.ms.afex)

  # glmer for pop
L.days.pop.glmer <- glmer(days ~ pop + (1|mom) + (1|plate), family=poisson(), data=L.dat.days, control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
L.days.pops.glmer.emmeans<- emmeans(L.days.pop.glmer, c("pop"), data=L.dat.days)
summary(L.days.pops.glmer.emmeans)
pairs(L.days.pops.glmer.emmeans)
CLD(L.days.pops.glmer.emmeans)
L.days.pop.clds <- CLD(L.days.pops.glmer.emmeans, details=FALSE, sort=FALSE)
plot(L.days.pops.glmer.emmeans, comparisons = TRUE, horizontal= FALSE)

  #afex for pop
L.days.pop.afex <- afex::mixed(days ~ pop + (1|mom) + (1|plate), family=poisson(), data=L.dat.days, method="LRT", control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
anova(L.days.pop.afex)
L.days.pops.afex.emmeans<- emmeans(L.days.pop.afex, c("pop"), data=L.dat.days)
pairs(L.days.pops.afex.emmeans)
cld(L.days.pops.afex.emmeans)
plot(L.days.pops.afex.emmeans, comparisons = TRUE, horizontal=FALSE)+
  theme_classic()+
  geom_text(aes(label=L.days.pop.clds$.group), vjust=-2, size=5)


L.days.ms.glmer2 <- glmer(days ~ ms + (1|pop), family=poisson(), data=L.dat.days)
summary(L.days.ms.glmer2)
overdisp.glmer(L.days.ms.glmer2)

#L.days.ms.glmer.nb <- glmer.nb(days~ ms +(1|pop/mom), data=L.dat.days)

L.days.ms.glmer.texreg <- texreg(L.days.ms.glmer, float.pos='h')

#overdisp.glmer(L.days.ms.glmer.nb)

#L.days.s.region.glmer <- glmer(days ~ s.region + (1|pop/mom), family=poisson(), data=na.omit(L.dat))
#summary(L.days.s.region.glmer)

AIC(L.days.ms.glm,L.days.ms.glmer,L.days.ms.glmer2)  

# GREENHOUSE SEEDLING DATA

  #ms glmer
L.surv.ms.glmer <- glmer(surv ~ ms + (1|pop/mom)+(1|rack), family=binomial(link="logit"), data=L.surv, control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
summary(L.surv.ms.glmer)
surv.ms.emms<-emmeans(L.surv.ms.glmer,"ms")
surv.ms.clds<-cld(surv.ms.emms)

L.surv.ms.afex <- afex::mixed(surv ~ ms + (1|pop/mom)+(1|rack), family=binomial(link="logit"), data=L.surv,method="LRT", control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
anova(L.surv.ms.afex)
L.surv.pop.afex <- afex::mixed(surv ~ pop + (1|mom)+(1|rack), family=binomial(link="logit"), data=L.surv,method="LRT", control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
anova(L.surv.pop.afex)


  #pop glmer
L.surv.pop.glmer <- glmer(surv ~ pop + (1|mom)+(1|rack), family=binomial(link="logit"), data=L.surv, control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
summary(L.surv.pop.glmer)

  #pop emmeans
L.surv.pop.glmer.emmeans<- emmeans(L.surv.pop.glmer, c("pop"), data=L.surv)
pairs(L.surv.pop.glmer.emmeans)
L.surv.pop.clds<-cld(L.surv.pop.glmer.emmeans, details=FALSE, sort=FALSE)
plot(L.surv.pop.glmer.emmeans, comparisons = TRUE, horizontal=FALSE)+
  theme_classic()+
  geom_text(aes(label=L.surv.pop.clds$.group), vjust=-2, size=5)

L.surv.ms.glmer.texreg <- texreg(L.surv.ms.glmer, float.pos = 'h')

L.surv.ms.glmer.noL06 <- glmer(surv ~ ms + (1|pop/mom), family=binomial(link="logit"), data=subset(L.surv, pop=!"L06"))
summary(L.surv.ms.glmer.noL06)

L.surv.s.region.glmer <- glmer(surv ~ s.region + (1|pop/mom), family=binomial(link="logit"), data=L.surv)
summary(L.surv.s.region.glmer)

L.surv.to.seedling.glmer <- glmer(surv.to.seedling ~ ms + (1|pop/mom)+(1|rack), family=binomial(link="logit"), data=L.surv)
summary(L.surv.to.seedling.glmer)
surv.to.seedling.emms <- emmeans(L.surv.to.seedling.glmer, 'ms')


L.surv.to.seedling.pop.glmer <- glmer(surv.to.seedling ~ pop + (1|mom)+(1|rack), family=binomial(link="logit"), data=L.surv, control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
summary(L.surv.to.seedling.pop.glmer)

H.y0.length.ms.lmer <- lmer(leaf.length ~ ms + (1|pop/mom)+(1|rack), data=H.y0)
summary(H.y0.length.ms.lmer)
anova(H.y0.length.ms.lmer)
length.ms.emms <- emmeans(H.y0.length.ms.lmer, 'ms')
length.ms.clds <- cld(length.ms.emms)

# lmer and emmeans for leaf length (pop)
H.y0.length.pop.lmer <- lmer(leaf.length ~ pop + (1|mom)+(1|rack), data=H.y0)
anova(H.y0.length.pop.lmer)
length.pop.lmer.emmeans <- emmeans(H.y0.length.pop.lmer, c("pop"), data=H.y0)
pairs(length.pop.lmer.emmeans)
length.pop.clds <- cld(length.pop.lmer.emmeans, details=FALSE, sort=FALSE)

# lmer and emmeans for leaf num (pop)
H.y0.num.ms.afex <- afex::mixed(leaf.num ~ ms + (1|pop/mom)+(1|rack), family=poisson(link="log"), method="LRT", data=H.y0, control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
summary(H.y0.num.ms.afex)
anova(H.y0.num.ms.afex)
num.ms.emms<-emmeans(H.y0.num.ms.afex, "ms")
num.ms.clds<-cld(num.ms.emms)

H.y0.num.pop.afex <- afex::mixed(leaf.num ~ pop + (1|mom)+(1|rack), family=poisson(link="log"), method="LRT", data=H.y0, control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
anova(H.y0.num.pop.afex)
num.pop.emmeans <- emmeans(H.y0.num.pop.afex, c("pop"), data=H.y0)
num.pop.clds <- cld(num.pop.emmeans, details=FALSE, sort=FALSE)

H.y0.num.pop.lmer <- lmer(leaf.num ~ pop + (1|mom), data=H.y0)
summary(H.y0.num.pop.lmer)

H.y0.length.lmer.texreg <- texreg(list(H.y0.length.ms.lmer,H.y0.length.pop.lmer),float.pos = 'h')

H.y0.length.mom.lm <- lm(leaf.length ~ mom, data=H.y0)
anova(H.y0.length.mom.lm)
### visreg(H.y0.length.mom.lm)
summary(H.y0.length.mom.lm)
H.y0.length.mom.lm.texreg <- texreg(H.y0.length.mom.lm, float.pos = "h")


L.length.surv.y1.glmer <- glmer(surv.to.y1 ~ (leaf.length) + (1|pop), family=binomial(link="logit"), data=H.y0)
summary(L.length.surv.y1.glmer)



 # FIELD ESTABLISHMENT #
glm.field.germ.ms <- glm(germ~ms, family = binomial(link="logit"), data=G.dat.y2)
summary(glm.field.germ.ms)
emmeans(glm.field.germ.ms, 'ms')

glm.field.germ.ms2 <- glm(germ~ms, family = binomial(link="logit"), data=G.dat.y2.filt)
summary(glm.field.germ.ms2)

G.germ.ms.glm <- glm(germ ~ ms, family=binomial(link="logit"), data=G.dat.y2)
summary(G.germ.ms.glm)

G.germ.s.region.glm <- glm(germ ~ s.region, family=binomial(link="logit"), data=G.dat.y2)
summary(G.germ.s.region.glm)

G.ms.glmer <- glmer(germ~ms+(1|pop/mom)+(1|garden), family=binomial(link="logit"), data=G.dat)
summary(G.ms.glmer)

glmer.msx <- glmer(germ ~ ms*g.region+(1|pop/mom)+(1|garden), family=binomial(link="logit"), data=G.dat.y2, control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
summary(glmer.msx)
anova(glmer.msx)
emmeans(glmer.msx, 'ms')

glmer.msx2 <- glmer(germ ~ ms*garden+(1|pop/mom), family=binomial(link="logit"), data=G.dat,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
afex.msx2 <- afex::mixed(germ ~ ms*garden+(1|pop/mom), family=binomial(link="logit"), data=G.dat, method="LRT", control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))

summary(glmer.msx2)
anova(afex.msx2)

emmip(glmer.msx2,  ~ ms | garden)

emmeans(glmer.msx2, "ms")
emmeans(glmer.msx2,  ~ ms | garden)
#glmer.msx2 <- glmer(germ ~ ms*garden+(1|pop/census.mom), family=binomial(link="logit"), data=G.dat.y2)
#summary(glmer.msx2)

G.germ.ms.glmer <- glmer(germ~ms+(1|garden)+(1|pop/mom), family=binomial(link="logit"), data=G.dat)
summary(G.germ.ms.glmer)
emmeans(G.germ.ms.glmer, 'ms')


  # BERTO DISPERSAL DATA
  

Ber.pca <- prcomp(Ber.log, center=TRUE, scale=TRUE)
#plot(Ber.pca, type="l")
Ber.pca.sum<- summary(Ber.pca)

# plot terminal velocity against angle of attack

Ber.TVxAA <- lm(terminal.velocity ~ angle.of.attack, data=Ber.dat)
summary(Ber.TVxAA)

Ber.lm.tv <- lm(terminal.velocity ~ ms, data=Ber.dat)
summary(Ber.lm.tv)

Ber.tv.lm.all <- lm(terminal.velocity ~ angle.of.attack + weight + bristle.length + number.of.bristles, data=Ber.dat)
Ber.tv.lm.all.texreg <- texreg(Ber.tv.lm.all, float.pos = "h", include.pvalues=TRUE)
Ber.dat.matrix <- Ber.dat[,3:length(Ber.dat)]
Ber.cor <- cor(Ber.dat.matrix)

Ber.weight <- lmer(weight~ms + (1|pop/mom), data=Ber.dat)
summary(Ber.weight)
Ber.aot <- lmer(angle.of.attack~ms + (1|pop/mom), data=Ber.dat)
summary(Ber.aot)
Ber.bl <- lmer(bristle.length~ms + (1|pop/mom), data=Ber.dat)
summary(Ber.bl)
Ber.nob <- lmer(number.of.bristles~ms + (1|pop/mom), data=Ber.dat)
summary(Ber.nob)

##########
# FIGURES
##########

  # LAB GERMINATION SUCCESS
gg.germ.box.ms <- ggplot(data=L.germ.mom, aes(y=mean.germ, x=ms))+
  geom_boxplot(aes(fill=ms), width=0.3, outlier.shape=NA)+
  scale_fill_manual(values=c("white","darkgrey"))+
  labs(x="mating system", y="mean germination success (moms)")+
  theme_classic()

gg.germ.box.pop <- ggplot(data=L.germ.mom, aes(y=mean.germ, x=pop))+
  geom_boxplot(aes(fill=ms), width=0.3, outlier.shape = NA, position=position_dodge(4))+
  geom_jitter(alpha=1/5, (aes(as.numeric(as.factor(pop)) + 0.3, mean.germ)), width=0.1)+
  scale_fill_manual(values=c("white","darkgrey"))+
  labs(x="population", y="mean germination success (moms)")+
  stat_summary(geom='text', label=L.germ.pop.clds$.group, vjust=-3, size=5)+
  theme_classic()

grid.arrange(gg.germ.box.ms, gg.germ.box.pop)

gg.lab.germ.ms <- ggplot(data=L.germ, aes(x=ms, y=mean.germ, fill=ms))+
  geom_bar(stat="identity", width=0.5, color="black")+
  geom_errorbar(aes(ymin=mean.germ-se.germ, ymax=mean.germ+se.germ),position=position_dodge(.2), width=.1)+
  scale_y_continuous(limits=c(0,1))

gg.lab.germ.pop <- ggplot(data=L.germ.pop, aes(x=pop,y=mean.germ, fill=ms))+
  geom_bar(stat="identity", width=0.5, color="black")+
  geom_errorbar(aes(ymin=mean.germ-se.germ, ymax=mean.germ+se.germ),position=position_dodge(.2), width=.1)+
  scale_y_continuous(limits=c(0,1))+
  geom_text(aes(label=L.germ.pop.clds$.group), vjust=-2, size=5)+
  labs(x="populations", y="mean germination success")+
  theme_classic()


  # LAB GERMINATION DAYS

gg.days.ms.box <- ggplot(data=L.days.mom, aes(y=mean.days, x=ms))+
  geom_boxplot(aes(fill=ms), width=0.3, outlier.shape=NA)+
  scale_fill_manual(values=c("white", "darkgrey"))+
  labs(x="mating system", y="mean days to germination")+
  ylim(0, 15)+
  theme_classic()

gg.days.pop.box <- ggplot(data=L.days.mom, aes(y=mean.days, x=factor(pop)))+
  geom_boxplot(aes(fill=ms), width=0.3, outlier.shape = NA, position=position_dodge(4))+
  geom_jitter(alpha=1/5, (aes(as.numeric(as.factor(pop)) + 0.3, mean.days)), width=0.1)+
  scale_fill_manual(values=c("white","darkgrey"))+
  stat_summary(geom='text', label=L.days.pop.clds$.group, vjust=-5, size=5)+
  labs(x="population", y="mean number of days until germination (mom)")+
  ylim(0,15)+
  theme(legend.position = "none")+
  theme_classic()
  
grid.arrange(gg.days.ms.box, gg.days.pop.box)

#gg.lab.germ.days.s.region <- ggplot(data=L.days, aes(x=ms,y=mean.days, group=s.region))+
#geom_point(aes(color=s.region), position=position_dodge(.2), size=3)+
#geom_errorbar(aes(ymin=mean.days-se.days, ymax=mean.days+se.days),position=position_dodge(.2), width=.1)

gg.days.pop <- ggplot(data=L.days.pop, aes(x=pop,y=mean.days, group=ms))+
  geom_point(aes(color=ms), position=position_dodge(.2), size=3)+
  geom_errorbar(aes(ymin=mean.days-se.days, ymax=mean.days+se.days),position=position_dodge(.2), width=.1)

tbl1 <- with(L.dat, table(ms, days))
tbl2 <- with(L.dat, table(ms,pop, days))

gg.days.freq1 <- ggplot(as.data.frame(tbl1), aes(factor(days), Freq, fill=ms))+
  geom_col(position="dodge")

gg.days.freq2 <- ggplot(as.data.frame(tbl2), aes(factor(days), Freq, fill=pop))+
  geom_col(position="dodge")


  # GREENHOUSE SEEDLING DATA

gg.surv.box.ms <- ggplot(data=L.surv.mom, aes(y=mean.surv, x=ms))+
  geom_boxplot(aes(fill=ms), width=0.3, outlier.shape=NA)+
  scale_fill_manual(values=c("white","darkgrey"))+
  labs(x="mating system", y="mean survival (moms)")+
  ylim(.45,1)+
  theme_classic()

gg.surv.box.pop <- ggplot(data=L.surv.mom, aes(y=mean.surv, x=pop))+
  geom_boxplot(aes(fill=ms), width=0.3, outlier.shape = NA)+
  geom_jitter(alpha=1/5, (aes(as.numeric(as.factor(pop)) + 0.3, mean.surv)), width=0.1)+
  scale_fill_manual(values=c("white","darkgrey"))+
  stat_summary(geom='text', label=L.surv.pop.clds$.group, vjust=4, size=5)+
  labs(x="population", y="mean survival (moms)")+
  ylim(.45,1)+
  theme_classic()

grid.arrange(gg.surv.box.ms,gg.surv.box.pop)

L.surv.to.seedling.mom
L.surv.mom

gg.survtoseed.box.ms <- ggplot(data=L.surv.to.seedling.mom, aes(y=mean.surv, x=ms))+
  geom_boxplot(aes(fill=ms))+
  geom_jitter(alpha=1/5, width=0.15)+
  scale_fill_manual(values=c("white","darkgrey"))+
  labs(x="mating system", y="mean survival (moms)")+
  theme_classic()

gg.survtoseed.box.pop <- ggplot(data=L.surv.to.seedling.mom, aes(y=mean.surv, x=pop))+
  geom_boxplot(aes(fill=ms))+
  geom_jitter(alpha=1/5, width=0.15)+
  scale_fill_manual(values=c("white","darkgrey"))+
  labs(x="population", y="mean survival (moms)")+
  stat_summary(geom='text', label=L.surv.pop.clds$.group, vjust=-3, size=5)+
  theme_classic()





# figures for length and leaf number

gg.box.length.ms.all <- ggplot(data=H.y0.mom, aes(y=mean.length, x=ms))+
  geom_boxplot(aes(fill=ms), width=0.3, outlier.shape=NA)+
  scale_fill_manual(values=c("white", "darkgrey"))+
  labs(x="mating system", y="mean leaf length")+
  theme_classic()

gg.box.length.pop <- ggplot(data=H.y0.mom, aes(y=mean.length, x=pop))+
  geom_boxplot(aes(fill=ms), width=0.3, outlier.shape = NA)+
  geom_jitter(alpha=1/5, (aes(as.numeric(as.factor(pop)) + 0.3, mean.length)), width=0.1)+
  scale_fill_manual(values=c("white","darkgrey"))+
  #stat_summary(geom='text', label=length.pop.clds$.group, vjust=-6, size=5)+
  labs(x="population", y="mean leaf length (moms)")+
  theme_classic()

grid.arrange(gg.box.length.ms.all, gg.box.length.pop)

gg.box.num.ms.all <- ggplot(data=H.y0.mom, aes(y=mean.num, x=ms))+
  geom_boxplot(aes(fill=ms), width=0.3, outlier.shape=NA)+
  scale_fill_manual(values=c("white", "darkgrey"))+
  labs(x="mating system", y="mean leaf number")+
  theme_classic()

gg.box.num.pop <- ggplot(data=H.y0.mom, aes(y=mean.num, x=pop))+
  geom_boxplot(aes(fill=ms), width=0.3, outlier.shape = NA)+
  geom_jitter(alpha=1/5, (aes(as.numeric(as.factor(pop)) + 0.3, mean.num)), width=0.1)+
  scale_fill_manual(values=c("white","darkgrey"))+
  #stat_summary(geom='text', label=num.pop.clds$.group, vjust=-6, size=5)+
  labs(x="population", y="mean number of leaves (moms)")+
  theme_classic()

grid.arrange(gg.box.num.ms.all, gg.box.num.pop)

y0.length.violin <- ggplot(data=H.y0, aes(factor(s.region), leaf.length))+
  geom_violin(data=H.y0, scale="count", aes(fill=factor(pop)), adjust=0.5)+
  geom_errorbar(data=H.y0.mean.length.s.region, aes(y=mean.length, ymin=mean.length-se.length, ymax=mean.length+se.length), colour="black")

y0.length.violin.mom <- ggplot(data=H.y0, aes(factor(mom), leaf.length))+
  geom_violin(data=H.y0, scale="count", aes(fill=factor(mom)))+
  facet_grid(~pop, scales = 'free')+
  theme(legend.position="none")

y1.length.violin <- ggplot(data=H.y1, aes(factor(s.region), leaf.length))+
  geom_violin(data=H.y1, scale="count", aes(fill=factor(pop)), adjust=0.5)+
  geom_errorbar(data=H.y1.mean.length.s.region, aes(y=mean.length, ymin=mean.length-se.length, ymax=mean.length+se.length), colour="black")

y0.y1.size.surv.violin <- ggplot(data=H.y0, aes(factor(surv.to.y1), leaf.length))+
  geom_violin(data=H.y0, scale="count", aes(fill=factor(surv.to.y1)), adjust=0.5)+
  geom_errorbar(data=H.y0.length.to.surv.y1, aes(y=mean.length, ymin=mean.length-se.length, ymax=mean.length+se.length), colour="black")+
  geom_point(data=H.y0.length.to.surv.y1, aes(y=mean.length), position=pd, colour="black", size=.5)

y0.y1.length.surv.boxplot <- ggplot(data=H.y0, aes(factor(surv.to.y1), leaf.length, fill=ms))+
  geom_boxplot()

gg.surv.to.seedling.pop <- ggplot(data=L.surv.to.seedling.pop, aes(x=pop, y=mean.surv, group=ms))+
  geom_point(aes(color=ms), position=position_dodge(.2), size=3)+
  geom_errorbar(aes(ymin=mean.surv-se.surv, ymax=mean.surv+se.surv),position=position_dodge(.2), width=.1)

  # FIELD ESTABLISHMENT

gg.est.box.ms <- ggplot(data=G.dat.y2.mom, aes(y=mean.germ, x=ms))+
  geom_boxplot(aes(fill=ms))+
  geom_jitter(alpha=1/5, width=0.15)+
  scale_fill_manual(values=c("white","darkgrey"))+
  labs(x="mating system", y="mean establishment success (moms)")+
  scale_y_continuous(limits=c(0,1))+
  theme_classic()

# establishment by ms and s.region
gg.germ.ms <- ggplot(data=G.all.ms, aes(x=ms, y=mean.germ, group=s.region))+
  geom_point(aes(shape=s.region, fill=s.region), position=position_dodge(.2), size=6)+
  geom_errorbar(aes(ymin=mean.germ-se.germ, ymax=mean.germ+se.germ),position=position_dodge(.2), width=.1)+
  scale_fill_manual(labels=c("S.s", "SO.s", "AO.s", "A.s"), values=c('black', 'black', 'grey80', 'grey80'), name="Source\nRegion")+
  scale_shape_manual(labels=c("S.s", "SO.s", "AO.s", "A.s"), values=c(22, 23, 25, 24), name="Source\nRegion")+
  theme_tufte()+
  geom_rangeframe()+
  scale_y_continuous(limits=c(0.25,1))+
  theme(text = element_text(size=20, family="CMU Serif"))+
  labs(x="Mating System", y="Proportion of seedlings established")+
  scale_linetype(guide=FALSE)+
  guides(fill=guide_legend(keywidth=0.35,keyheight=0.35,default.unit="inch"))

gg.germ.all <- ggplot(data=G.all1, aes(x=g.region, y=mean.germ, group=s.region))+
  geom_line(aes(linetype=s.region), position=position_dodge(.2))+
  geom_point(aes(shape=s.region, fill=s.region), position=position_dodge(.2), size=6)+
  geom_errorbar(aes(ymin=mean.germ-se.germ, ymax=mean.germ+se.germ),position=position_dodge(.2), width=.1)+
  scale_fill_manual(labels=c("S.s", "SO.s", "AO.s", "A.s"), values=c('black', 'black', 'grey80', 'grey80'), name="Source\nRegion")+
  scale_shape_manual(labels=c("S.s", "SO.s", "AO.s", "A.s"), values=c(22, 23, 25, 24), name="Source\nRegion")+
  theme_tufte()+
  geom_rangeframe()+
  scale_y_continuous(limits=c(0,1))+
  theme(text = element_text(size=20, family="CMU Serif"))+
  labs(x="Garden Region", y="Establishment Success")+
  scale_linetype(guide=FALSE)


gg.germ.all.pop <- ggplot(data=G.all.pop, aes(x=pop, y=mean.germ, group=s.region))+
  geom_point(aes(color=s.region), position=position_dodge(.2), size=3)+
  geom_errorbar(aes(ymin=mean.germ-se.germ, ymax=mean.germ+se.germ),position=position_dodge(.2), width=.1)

gg.germ.ms.garden <- ggplot(data=G.all.garden, aes(x=garden, y=mean.germ, group=ms))+
  geom_point(aes(color=ms), position=position_dodge(.2), size=3)+
  geom_errorbar(aes(ymin=mean.germ-se.germ, ymax=mean.germ+se.germ),position=position_dodge(.2), width=.1)

gg.germ.pop.gardens <- ggplot(data=G.pop.garden, aes(factor(pop), mean.germ))+
  geom_point(aes(color=pop), position=position_dodge(.2), size=3)+
  geom_errorbar(aes(ymin=mean.germ-se.germ, ymax=mean.germ+se.germ),position=position_dodge(.2), width=.1)+
  facet_grid(garden~., scales = 'free')+
  theme(legend.position="none")+
  scale_y_continuous(limits=c(0,1))

gg.germ.pops <- ggplot(data=G.pop.garden, aes(factor(garden), mean.germ))+
  geom_point(aes(color=pop), position=position_dodge(.2), size=3)+
  geom_errorbar(aes(ymin=mean.germ-se.germ, ymax=mean.germ+se.germ),position=position_dodge(.2), width=.1)+
  facet_grid(pop~., scales = 'free')+
  theme(legend.position="none")+
  scale_y_continuous(limits=c(0,1))+
  geom_line(aes(group=pop), linetype=3)

gg.germ.diffs <- ggplot(data=G.pgm, aes(factor(pop), diffs))+
  geom_point(aes(color=pop), position=position_dodge(.2), size=3)+
  geom_hline(yintercept = 0)+
  facet_grid(garden~., scales = 'free')+
  theme(legend.position="none")+
  scale_y_continuous(limits=c(-0.5,0.75))+
  geom_linerange(aes(x=pop, ymax=diffs, ymin=0))


  # BERTO DISPERSAL DATA


#Ber.means.TV$ms <- factor(Ber.means.TV$ms, levels = rev(levels(factor(Ber.means.TV$ms))))

Ber.moms.all

gg.box.ber.pops <- ggplot(data=Ber.moms.all, aes(y=terminal.velocity, x=factor(pop)))+
  geom_boxplot(aes(fill=ms))+
  labs(x="population", y="mean terminal velocity")+
  theme_classic()

#stat_summary(geom='text', label=L.days.pop.clds$.group, vjust=-3, size=5)+

# mean TV
gg.Ber.tv <- ggplot(data=Ber.means.TV, aes(x=ms, y=mean.terminal.velocity, group=ms))+
  geom_point(aes(color=ms), position=position_dodge(.2), size=3)+
  geom_errorbar(aes(ymin=mean.terminal.velocity-se.terminal.velocity, ymax=mean.terminal.velocity+se.terminal.velocity),position=position_dodge(.2), width=.1)+
  scale_y_continuous(limits=c(0,1.5))

# plot of relationship between terminal velocity and angle of attack / bristle length

#function that allows annotation of ggplot with lm equation and r2
lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}

gg.tv.angle<- ggplot(data=Ber.dat, aes(x=angle.of.attack, y=terminal.velocity, colour=ms))+
  geom_point()+
  geom_smooth(aes(group=1), method=lm)+
  annotate("text", x = 150, y = 2.5, label = lm_eqn(lm(terminal.velocity ~ angle.of.attack, Ber.dat)),colour="black",size=5, parse = TRUE)+
  scale_colour_few()+
  theme_few()

gg.tv.length<- ggplot(data=Ber.dat, aes(x=bristle.length, y=terminal.velocity, colour=ms))+
  geom_point()+
  geom_smooth(aes(group=1), method=lm)+
  annotate("text", x = 9, y = 2.5, label = lm_eqn(lm(terminal.velocity ~ bristle.length, Ber.dat)),colour="black",size=5, parse = TRUE)+
  scale_colour_few()+
  theme_few()

# Quick PCA
Ber.plot <- autoplot(Ber.pca, data=Ber.dat, colour='ms', loadings=TRUE, loadings.colour="blue", loadings.label=TRUE, loadings.label.colour="black", loadings.label.size=5, frame=TRUE, frame.type='norm')


# leaf length x leaf num
gg.num.x.length.pop<- ggplot(data=H.y0, aes(x=leaf.num, y=leaf.length, colour=pop))+
  geom_point()

gg.num.x.length.ms<- ggplot(data=H.y0, aes(x=leaf.num, y=leaf.length, colour=ms))+
  geom_point()

gg.num.x.length.s.region<- ggplot(data=H.y0, aes(x=leaf.num, y=leaf.length, colour=s.region))+
  geom_point()

gg.days.x.length <- ggplot(data=L.y0.days.length, aes(x=mean.days, y=mean.length, colour=ms))+
  geom_point()+
  geom_smooth(aes(group=1), method=lm)+
  annotate("text", x = 18, y = 18, label = lm_eqn(lm(mean.days ~ mean.length, L.y0.days.length)),colour="black",size=5, parse = TRUE)

# seedling surv and field est comparison plot
# comparisons between survival to seedling stage (y0) and field establishment 
LG.est <- left_join(L.surv.to.seedling.pop, G.all.pop)

gg.LG.est <- ggplot(LG.est, aes(y=mean.surv, x=mean.germ))+
  geom_point(aes(color=pop), size=4)+
  geom_text(label=LG.est$pop)+
  geom_smooth(method=lm, se=FALSE)+
  labs(x="mean field establishment", y="mean survival to seedling")