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

############
# Load data
############

L.dat <<-read.csv("./data/mastergerm_lab.csv", stringsAsFactors = FALSE, strip.white = TRUE, na.string = c("NA",""))
L.surv <<-read.csv("./data/labseedlingsurvival.csv", stringsAsFactors = FALSE, strip.white = TRUE, na.string = c("NA",""))
Berto.dat <<-read.csv("./data/alberto.seed.R.csv", stringsAsFactors = FALSE, strip.white = TRUE, na.string = c("NA",""))
H.dat <<-read.csv("./data/greenhouse.csv", stringsAsFactors = FALSE, strip.white = TRUE, na.string = c("NA",""))
G.dat <<-read.csv("./data/mastergermv3.csv", stringsAsFactors = FALSE, strip.white = TRUE, na.string = c("NA",""))

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

# H.y0.length.to.surv.y1 <- H.dat %>%
#   filter(year==0) %>%
#   filter(!is.na(surv.to.y1)) %>%
#   filter(!is.na(leaf.length)) %>%
#   group_by(surv.to.y1) %>%
#   summarize(mean.length=mean(leaf.length), se.length=std.error(leaf.length))

# H.y0.length.to.surv.y1.ms <- H.dat %>%
#   filter(year==0) %>%
#   filter(!is.na(surv.to.y1)) %>%
#   filter(!is.na(leaf.length)) %>%
#   group_by(ms, surv.to.y1) %>%
#   summarize(mean.length=mean(leaf.length), se.length=std.error(leaf.length))


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

# selecting just the main columns that I'm interested in
Ber.dat <- Berto.dat %>%
  filter(region!="WestWY") %>%
  filter(mom!="Bulk")


# factors and levels
Ber.dat$pop <- factor(Ber.dat$pop, levels=c("B42", "B42.alt", "B46", "B49", "L12", "C59", "B52", "L16", "L40", "L41", "C21", "C43", "C54"))
Ber.dat$region <- as.factor(Ber.dat$region)
Ber.dat$garden <- as.factor(Ber.dat$garden)
Ber.dat$ms <- factor(Ber.dat$ms, levels=c("S", "A"))
Ber.dat$lat <- as.factor(Ber.dat$lat)

Ber.dat$terminal.velocity.log <- log(Ber.dat$terminal.velocity)
Ber.dat$weight.log <- log(Ber.dat$weight)
Ber.dat$angle.of.attack.log <- log(Ber.dat$angle.of.attack)
Ber.dat$bristle.length.log <- log(Ber.dat$bristle.length)
Ber.dat$number.of.bristles.log <- log(Ber.dat$number.of.bristles)
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
#Ber.log <- log(Ber.dat[c(3,4,6,7)])
#Ber.pops <- Ber.dat[, 1]

# correlation matrix
#Ber.matrix <- cor(Ber.log)
#Ber.matrix.4 <- round(Ber.matrix, 4)

##########
# Models
##########

# lab germination

L.germ.ms.glmer <- glmer(germ ~ ms + (1|pop/mom)+(1|plate), family=binomial(link="logit"), data=L.dat)
summary(L.germ.ms.glmer)
germ.ggfx <- ggpredict(L.germ.ms.glmer, terms = c("ms"))
plot(germ.ggfx)
germ.ggfx.ms.df <- as.data.frame(germ.ggfx)

L.germ.ms.glmer2 <- glmer(germ ~ pop + (1|mom)+(1|plate), family=binomial(link="logit"), data=L.dat, control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
summary(L.germ.ms.glmer2)
germ.ggfx2 <- ggpredict(L.germ.ms.glmer2, terms = c("pop"))
plot(germ.ggfx2)
germ.ggfx.pop.df <- as.data.frame(germ.ggfx2)
germ.ggfx.pop.df$group <- c("1","1","1","1","1","1","2","2","2","2","2","2")

lrtest(L.germ.ms.glmer,L.germ.ms.glmer2)

var.test(germ ~ ms, data=L.dat, alternative = "two.sided")

bartlett.test(germ ~ ms, data=L.dat)

leveneTest(germ ~ ms, data=L.dat)
# germ days

L.days.ms.glmer <- glmer(days ~ ms + (1|pop/mom)+(1|plate), family=poisson(), data=L.dat.days)
summary(L.days.ms.glmer)
days.ggfx <- ggpredict(L.days.ms.glmer, terms = c("ms"), type="re")
days.ggfx.randpop <- ggpredict(L.days.ms.glmer, terms = c("pop"), type="re")
plot(days.ggfx)
days.ggfx.ms.df <- as.data.frame(days.ggfx)

L.days.pop.glmer <- glmer(days ~ pop + (1|mom)+(1|plate), family=poisson(), data=L.dat.days, control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
summary(L.days.pop.glmer)
days.ggfx2 <- ggpredict(L.days.pop.glmer, terms = c("pop"))
plot(days.ggfx2)
days.ggfx.pop.df <- as.data.frame(days.ggfx2)
days.ggfx.pop.df$group <- c("1","1","1","1","1","1","2","2","2","2","2","2")

# greenhouse surv

L.surv.ms.glmer <- glmer(surv ~ ms + (1|pop/mom)+(1|rack), family=binomial(link="logit"), data=L.surv)
summary(L.surv.ms.glmer)
gsurv.ggfx <- ggpredict(L.surv.ms.glmer, terms = c("ms"))
plot(gsurv.ggfx)
gsurv.ggfx.ms.df <- as.data.frame(gsurv.ggfx)

L.surv.ms.glmer.2 <- glmer(surv ~ (1|pop/mom)+(1|rack), family=binomial(link="logit"), data=L.surv)
lrtest(L.surv.ms.glmer, L.surv.ms.glmer.2)
drop1(L.surv.ms.glmer, test="Chisq")
anova(L.surv.ms.glmer, L.surv.ms.glmer.2)

L.surv.pop.glmer <- glmer(surv ~ pop + (1|mom)+(1|rack), family=binomial(link="logit"), data=L.surv, control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
summary(L.surv.pop.glmer)
gsurv.ggfx2 <- ggpredict(L.surv.pop.glmer, terms = c("pop"))
plot(gsurv.ggfx2)
gsurv.ggfx.pop.df <- as.data.frame(gsurv.ggfx2)
gsurv.ggfx.pop.df$group <- c("1","1","1","1","1","1","2","2","2","2","2","2")

lrtest(L.surv.ms.glmer, L.surv.pop.glmer)

# length 

H.y0.length.ms.lmer <- lmer(leaf.length ~ ms + (1|pop/mom), data=H.y0)
summary(H.y0.length.ms.lmer)
length.ggfx <- ggpredict(H.y0.length.ms.lmer, terms="ms")
plot(length.ggfx)
length.ggfx.ms.df <- as.data.frame(length.ggfx)

H.y0.length.pop.lmer <- lmer(leaf.length ~ pop + (1|mom), data=H.y0)
length.pop.ggfx <- ggpredict(H.y0.length.pop.lmer, terms="pop")
plot(length.pop.ggfx)
length.ggfx.pop.df <- as.data.frame(length.pop.ggfx)
length.ggfx.pop.df$group <- c("1","1","1","1","1","1","2","2","2","2","2","2")

# num

H.y0.num.ms.glmer <- glmer(leaf.num ~ ms + (1|pop/mom), data=H.y0, family=poisson(link="log"))
summary(H.y0.num.ms.glmer)
num.ggfx <- ggpredict(H.y0.num.ms.glmer, terms="ms")
plot(num.ggfx)
num.ggfx.ms.df <- as.data.frame(num.ggfx) 

H.y0.num.pop.glmer <- glmer(leaf.num ~ pop + (1|mom), data=H.y0, family=poisson(link="log"), control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
num.pop.ggfx <- ggpredict(H.y0.num.pop.glmer, terms = "pop")
plot(num.pop.ggfx)
num.ggfx.pop.df <- as.data.frame(num.pop.ggfx)
num.ggfx.pop.df$group <- c("1","1","1","1","1","1","2","2","2","2","2","2")

# berto dat

Ber.weight <- lmer(weight~ms + (1|pop/mom), data=Ber.dat)
summary(Ber.weight)
weight.ggfx <- ggpredict(Ber.weight, terms="ms")
plot(weight.ggfx) # weight does not differ between MS

Ber.aot <- lmer(angle.of.attack~ms + (1|pop/mom), data=Ber.dat)
summary(Ber.aot)
aot.ggfx <- ggpredict(Ber.aot, terms="ms")
plot(aot.ggfx) # angle of attack differs significantly by ms
ggplot(aes(ms, angle.of.attack), data=Ber.dat)+geom_jitter()

Ber.bl <- lmer(bristle.length~ms + (1|pop/mom), data=Ber.dat)
summary(Ber.bl)
bl.ggfx <- ggpredict(Ber.bl, terms="ms")
plot(bl.ggfx) # very slight overlap
ggplot(aes(ms, bristle.length), data=Ber.dat)+geom_jitter()

Ber.bl2 <- lmer(bristle.length~(1|pop/mom), data=Ber.dat)
lrtest(Ber.bl, Ber.bl2) # ms is significant by LRT

Ber.nob <- lmer(number.of.bristles~ms + (1|pop/mom), data=Ber.dat)
summary(Ber.nob)
nob.ggfx <- ggpredict(Ber.nob, terms="ms")
plot(nob.ggfx) # no difference in number of bristles

Ber.tv <- lmer(terminal.velocity~ms + (1|pop/mom), data=Ber.dat)
summary(Ber.tv)
tv.ggfx <- ggpredict(Ber.tv, terms="ms")
tv.ggfx.ms.df <- as.data.frame(tv.ggfx)
plot(tv.ggfx)
ggplot(aes(ms, terminal.velocity, fill=ms), data=Ber.dat)+geom_jitter(shape=21)
var.test(terminal.velocity ~ ms, data=Ber.dat)

Ber.tv.pop <- lmer(terminal.velocity~pop+(1|mom), data=Ber.dat)
tv.ggfx.pop <- ggpredict(Ber.tv.pop, terms="pop")
tv.ggfx.pop.df <- as.data.frame(tv.ggfx.pop)
tv.ggfx.pop.df$group <- c("1","1","1","1","1","2","2","2","2","2","2","2")


Ber.pop.tv <- lmer(terminal.velocity~ pop + (1|mom), data=Ber.dat)
tv.pop.ggfx <- ggpredict(Ber.pop.tv, terms="pop")
plot(tv.pop.ggfx)


Ber.tv.lm.all <- lm(terminal.velocity.log ~ angle.of.attack.log + weight.log + bristle.length.log + number.of.bristles.log, data=Ber.dat)
summary(Ber.tv.lm.all)


########
# figures
########

# germ

gg.germ.ms <-ggplot()+
  geom_boxplot(data=L.germ.mom, aes(y=mean.germ, x=ms), alpha=0, width=0.3, outlier.shape=NA)+
  geom_jitter(data=L.germ.mom, aes(x=ms, y=mean.germ, fill=L.germ.mom$ms), alpha=0.7, shape=21, width=0.1)+
  geom_point(data=germ.ggfx.ms.df, aes(x=as.numeric(as.factor(x))+0.4, y=predicted, fill=x), shape=21, size=3)+
  geom_errorbar(data=germ.ggfx.ms.df, aes(ymax=conf.high, ymin=conf.low, x=as.numeric(as.factor(x))+0.4), width=0.1)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="mating system", labels=c("sexual", "apomict"))+
  labs(x="mating system", y="mean germination success")+
  theme_classic()

# plot(germ.ggfx)
# 
# gg.germ.box.ms <- ggplot(data=L.germ.mom, aes(y=mean.germ, x=ms))+
#   geom_boxplot(aes(fill=ms), width=0.3, outlier.shape=NA)+
#   scale_fill_manual(values=c("white","darkgrey"))+
#   labs(x="mating system", y="mean germination success (moms)")+
#   theme_classic()

#tiff("ELH_newplot.tiff", height=6, width=8, res=300, units="in")

#dev.off()

gg.germ.box.pop <- ggplot(data=L.germ.mom, aes(y=mean.germ, x=pop))+
  geom_boxplot(aes(fill=ms), width=0.3, outlier.shape = NA, position=position_dodge(4))+
  geom_jitter(alpha=1/5, (aes(as.numeric(as.factor(pop)) + 0.3, mean.germ)), width=0.1)+
  scale_fill_manual(values=c("white","darkgrey"))+
  labs(x="population", y="mean germination success (moms)")+
  #stat_summary(geom='text', label=L.germ.pop.clds$.group, vjust=-3, size=5)+
  theme_classic()

#grid.arrange(gg.germ.box.ms, gg.germ.box.pop)

# days

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
  #stat_summary(geom='text', label=L.days.pop.clds$.group, vjust=-5, size=5)+
  labs(x="population", y="mean number of days until germination (mom)")+
  ylim(0,15)+
  theme(legend.position = "none")+
  theme_classic()

#grid.arrange(gg.days.ms.box, gg.days.pop.box)

# surv

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
  #stat_summary(geom='text', label=L.surv.pop.clds$.group, vjust=4, size=5)+
  labs(x="population", y="mean survival (moms)")+
  ylim(.45,1)+
  theme_classic()

#grid.arrange(gg.surv.box.ms,gg.surv.box.pop)

# tv

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}


gg.tv.weight<-ggplot(aes(y=terminal.velocity, x=weight, fill=ms), data=Ber.dat)+
  geom_point(shape=21)+
  geom_smooth(method='lm', se=FALSE, aes(color=ms))+
  theme_classic()+
  labs(x="seed mass (g)", y="terminal velocity (cm/s)")+
  scale_colour_manual(values=c("coral3","cornflowerblue"), name="mating system", labels=c("sexual", "apomict"))+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="mating system", labels=c("sexual", "apomict"))+
  theme(legend.position = "none")

gg.tv.aoa<-ggplot(aes(y=terminal.velocity, x=angle.of.attack, fill=ms), data=Ber.dat)+
  geom_point(shape=21)+
  geom_smooth(method='lm', se=FALSE, aes(color=ms))+
  theme_classic()+
  labs(x="angle of attack", y="terminal velocity (cm/s)")+
  scale_colour_manual(values=c("coral3","cornflowerblue"), name="mating system", labels=c("sexual", "apomict"))+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="mating system", labels=c("sexual", "apomict"))+
  theme(legend.position = "none")

gg.tv.bl<-ggplot(aes(y=terminal.velocity, x=bristle.length, fill=ms), data=Ber.dat)+
  geom_point(shape=21)+
  geom_smooth(method='lm', se=FALSE, aes(color=ms))+
  theme_classic()+
  labs(x="bristle length (mm)", y="terminal velocity (cm/s)")+
  scale_colour_manual(values=c("coral3","cornflowerblue"), name="mating system", labels=c("sexual", "apomict"))+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="mating system", labels=c("sexual", "apomict"))+
  theme(legend.position = "none")


gg.tv.nb<-ggplot(aes(y=terminal.velocity, x=number.of.bristles, fill=ms), data=Ber.dat)+
  geom_point(shape=21)+
  geom_smooth(method='lm', se=FALSE, aes(color=ms))+
  theme_classic()+
  labs(x="number of bristles", y="terminal velocity (cm/s)")+
  scale_colour_manual(values=c("coral3","cornflowerblue"), name="mating system", labels=c("sexual", "apomict"))+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="mating system", labels=c("sexual", "apomict"))+
  theme(legend.position = "bottom")

ggarrange(gg.tv.aoa,gg.tv.bl,gg.tv.nb,gg.tv.weight, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

# ms
gg.box.ber.ms.all <- ggplot(data=Ber.dat, aes(y=terminal.velocity, x=ms))+
  geom_boxplot(aes(fill=ms), width=0.3, outlier.shape=NA)+
  scale_fill_manual(values=c("white", "darkgrey"))+
  labs(x="mating system", y="mean terminal velocity (cm/s)")+
  theme_classic()

ggplot()+
  geom_boxplot(data=Ber.dat, aes(y=terminal.velocity, x=ms), alpha=0, width=0.3, outlier.shape=NA)+
  geom_jitter(data=Ber.dat, aes(x=ms, y=terminal.velocity, fill=ms), alpha=0.7, shape=21, width=0.1)+
  geom_point(data=tv.ggfx.ms.df, aes(x=as.numeric(as.factor(x))+0.4, y=predicted, fill=x), shape=21, size=3)+
  geom_errorbar(data=tv.ggfx.ms.df, aes(ymax=conf.high, ymin=conf.low, x=as.numeric(as.factor(x))+0.4), width=0.1)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="mating system", labels=c("sexual", "apomict"))+
  labs(x="mating system", y="mean terminal velocity")+
  theme_classic()

# pops
gg.box.ber.allpops <- ggplot(data=Ber.dat, aes(y=terminal.velocity, x=pop))+
  geom_boxplot(aes(fill=ms), width=0.3, outlier.shape = NA)+
  geom_jitter(alpha=1/5, (aes(as.numeric(as.factor(pop)) + 0.3, terminal.velocity)), width=0.1)+
  scale_fill_manual(values=c("white","darkgrey"))+
  labs(x="population", y="mean terminal velocity (cm/s)")+
  theme_classic()

#grid.arrange(gg.box.ber.ms.all,gg.box.ber.allpops)
