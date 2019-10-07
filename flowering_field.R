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
library(xtable)
library(stargazer)
library(ggfortify)

library(extrafont)
#font_import()
#loadfonts()

############
# Load data
############

H.dat <<-read.csv("~/GitHub/Hookeri-Gardens/RAW DATA/mastergarden2018.csv", stringsAsFactors = FALSE, strip.white = TRUE, na.string = c("NA",""))

############
# Data Prep
############

# calculate size
H.dat$size <- H.dat$leaf.num * H.dat$leaf.length
H.dat$log.size <- log(H.dat$size+1)

# Garden data
H.dat$pop <- as.factor(H.dat$pop)
H.dat$ms <- as.factor(H.dat$ms)
H.dat$s.region <- as.factor(H.dat$s.region)
H.dat$garden <- as.factor(H.dat$garden)
H.dat$g.region <- as.factor(H.dat$g.region)
H.dat$year <- as.numeric(H.dat$year)
H.dat$leaf.num <- as.numeric(H.dat$leaf.num)
H.dat$leaf.length <- as.numeric(H.dat$leaf.length)
H.dat$size <- as.numeric(H.dat$size)
H.dat$log.size <- as.numeric(H.dat$log.size)
H.dat$surv.to.y1 <- as.factor(H.dat$surv.to.y1)
H.dat$surv.to.y2 <- as.factor(H.dat$surv.to.y2)
H.dat$surv.years <- as.numeric(H.dat$surv.years)
H.dat$flower.years <- as.numeric(H.dat$flower.years)

# set levels
H.dat$s.region <- factor(H.dat$s.region, levels=c("A.s", "AO.s", "SO.s", "S.s"))
H.dat$g.region <- factor(H.dat$g.region, levels=c("A.g", "AO.g", "SO.g", "S.g"))
H.dat$pop <- factor(H.dat$pop, levels=c("C27", "C85", "C86", "L17", "L16", "L06", "L12", "L11", "B49", "B46", "B42", "B53")) #sorted by lat, from high to low
H.dat$garden <- factor(H.dat$garden, levels = rev(levels(factor(H.dat$garden))))
H.dat$ms <- factor(H.dat$ms, levels = rev(levels(factor(H.dat$ms))))

############
# Means, subsets, etc.
############

# count flowers y2
H.flwr.inds.y2 <- H.dat %>%
  filter(year==2) %>%
  filter(!is.na(flower.y2.10)) %>%
  group_by(g.region, ms) %>%
  summarize(mean=mean(flower.y2.10))

# count flowers y3
H.flwr.inds.y3 <- H.dat %>%
  filter(year==3) %>%
  filter(!is.na(flower.y3.10)) %>%
  group_by(g.region, ms) %>%
  summarize(mean=mean(flower.y3.10))

# count flowers y4
H.flwr.inds.y4 <- H.dat %>%
  filter(year==4) %>%
  filter(flower.y4.10==1) %>%
  group_by(garden) %>%
  summarize(n=n())

# mean size in y0 of inds flowering in y2
H.y0.size.flwr <- H.dat %>%
  filter(year==0) %>%
  filter(!is.na(flower.y2)) %>%
  filter(!is.na(log.size)) %>%
  group_by(flower.y2) %>%
  summarize(mean.size=mean(log.size), se.size=std.error(log.size))

# mean size y1 of inds flowering in y2
H.y1.size.flwr <- H.dat %>%
  filter(year==1) %>%
  filter(!is.na(flower.y2)) %>%
  filter(!is.na(log.size)) %>%
  group_by(flower.y2) %>%
  summarize(mean.size=mean(log.size), se.size=std.error(log.size))

# size y2 of inds flowering y2
H.y2.size.flwr <- H.dat %>%
  filter(year==2) %>%
  filter(!is.na(flower.y2)) %>%
  filter(!is.na(log.size)) %>%
  group_by(flower.y2) %>%
  summarize(mean.size=mean(log.size), se.size=std.error(log.size))

# mean number of inds flowering by year 3
H.flwr.10 <- H.dat %>%
  filter(year==3) %>%
  filter(!is.na(flower.10)) %>%
  group_by(g.region, ms) %>%
  summarize(mean=mean(flower.10), se=std.error(flower.10))

# mean number of inds flowering by year 3 and survived in y2?
H.flwr.10.surv.y3 <- H.dat %>%
  filter(year==3) %>%
  filter(surv==1) %>%
  filter(!is.na(flower.y3.10)) %>%
  group_by(g.region, ms) %>%
  summarize(mean=mean(flower.y3.10), se=std.error(flower.y3.10), count=sum(flower.y3.10))

H.flwr.10.surv.y2 <- H.dat %>%
  filter(year==2) %>%
  filter(surv==1) %>%
  filter(!is.na(flower.y2.10)) %>%
  group_by(g.region, ms) %>%
  summarize(mean=mean(flower.y2.10), se=std.error(flower.y2.10), count=sum(flower.y2.10))



# number of inds that flowered twice
H.flwr.yrs <- H.dat %>%
  filter(year==3) %>%
  filter(!is.na(flower.years)) %>%
  filter(flower.years==2) %>%
  group_by(g.region, ms) %>%
  summarize(count=n_distinct(rep))

# mean number of flowers from y3
H.flwr.num <- H.dat %>%
  filter(year==3) %>%
  filter(flower.y3.10==1) %>%
  group_by(g.region, ms) %>%
  summarize(mean=mean(bud.num))

H.flwer.count.y2 <- H.dat %>%
  filter(year==2) %>%
  filter(!is.na(flower.y2.10)) %>%
  filter(surv==1) %>%
  group_by(ms) %>%
  summarize(flower.count=sum(flower.y2.10), surv.count=sum(surv))

H.flwer.count.y3 <- H.dat %>%
  filter(year==3) %>%
  filter(!is.na(flower.y3.10)) %>%
  filter(surv==1) %>%
  group_by(ms) %>%
  summarize(flower.count=sum(flower.y3.10), surv.count=sum(surv))

H.flwer.count.y3.2 <- H.dat %>%
  filter(year==3) %>%
  filter(!is.na(flower.y3.10)) %>%
  filter(surv==1) %>%
  group_by(garden, ms) %>%
  summarize(flower.count=sum(flower.y3.10), surv.count=sum(surv))

H.flwr.mom <-H.dat %>%
  filter(!is.na(flower)) %>%
  filter(year>=2) %>%
  group_by(year, garden, ms) %>%
  summarize(flower.mom=sum(flower))

H.flwr.num <-H.dat %>%
  filter(!is.na(bud.num)) %>%
  filter(year>=2) %>%
  group_by(year, garden, ms) %>%
  summarize(flower.num=sum(bud.num))

H.mean.flwr.num <-H.dat %>%
  filter(!is.na(bud.num)) %>%
  filter(year>=2) %>%
  filter(bud.num>=1) %>%
  group_by(year, garden, ms) %>%
  summarize(mean.flower.num=mean(bud.num))

############
# Models
############

# #this one has the lowest AIC, and tons of sig.
# H.mod.flwr.1 <- glmer(flower.y3.10 ~ ms*g.region+(1|garden)+(1|pop/s.region), family=binomial(link="logit"), data=subset(H.dat, year==3))
# summary(H.mod.flwr.1)
# 
# H.mod.flwr.10.1 <- glmer(flower.y3.10 ~ ms*g.region+(1|garden)+(1|pop/s.region), family=binomial(link="logit"), data=subset(H.dat, year==3 & surv==1))
# summary(H.mod.flwr.10.1)
# 
# H.mod.flwr.10.2 <- glmer(flower.y3.10 ~ ms*g.region+(1|garden)+(1|s.region), family=binomial(link="logit"), data=subset(H.dat, year==3 & surv==1))
# summary(H.mod.flwr.10.2)
# 
# H.mod.flwr.10.3 <- glmer(flower.y3.10 ~ ms*g.region+(1|garden), family=binomial(link="logit"), data=subset(H.dat, year==3 & surv==1))
# summary(H.mod.flwr.10.3)
# 
# H.mod.flwr.10.4 <- glmer(flower.y3.10 ~ s.region*g.region+(1|garden)+(1|pop), family=binomial(link="logit"), data=subset(H.dat, year==3 & surv==1))
# summary(H.mod.flwr.10.4)
############
# Figures
############

# reverse levels for figures (south to north)
H.flwr.inds.y2$ms <- factor(H.flwr.inds.y2$ms, levels = rev(levels(factor(H.flwr.inds.y2$ms))))
H.flwr.inds.y2$g.region <- factor(H.flwr.inds.y2$g.region, levels = rev(levels(factor(H.flwr.inds.y2$g.region))))

H.flwr.inds.y3$ms <- factor(H.flwr.inds.y3$ms, levels = rev(levels(factor(H.flwr.inds.y3$ms))))
H.flwr.inds.y3$g.region <- factor(H.flwr.inds.y3$g.region, levels = rev(levels(factor(H.flwr.inds.y3$g.region))))

H.flwr.10$ms <- factor(H.flwr.10$ms, levels = rev(levels(factor(H.flwr.10$ms))))
H.flwr.10$g.region <- factor(H.flwr.10$g.region, levels = rev(levels(factor(H.flwr.10$g.region))))

# H.flwr.10.surv$ms <- factor(H.flwr.10.surv$ms, levels = rev(levels(factor(H.flwr.10.surv$ms))))
# H.flwr.10.surv$g.region <- factor(H.flwr.10.surv$g.region, levels = rev(levels(factor(H.flwr.10.surv$g.region))))

H.flwr.10.surv.y2$ms <- factor(H.flwr.10.surv.y2$ms, levels = rev(levels(factor(H.flwr.10.surv.y2$ms))))
H.flwr.10.surv.y2$g.region <- factor(H.flwr.10.surv.y2$g.region, levels = rev(levels(factor(H.flwr.10.surv.y2$g.region))))

H.flwr.10.surv.y3$ms <- factor(H.flwr.10.surv.y3$ms, levels = rev(levels(factor(H.flwr.10.surv.y3$ms))))
H.flwr.10.surv.y3$g.region <- factor(H.flwr.10.surv.y3$g.region, levels = rev(levels(factor(H.flwr.10.surv.y3$g.region))))

H.flwr.yrs$ms <- factor(H.flwr.yrs$ms, levels = rev(levels(factor(H.flwr.yrs$ms))))
H.flwr.yrs$g.region <- factor(H.flwr.yrs$g.region, levels = rev(levels(factor(H.flwr.yrs$g.region))))


# barplot of # of flowering inds per region
gg.flwr.y2 <- ggplot(H.flwr.inds.y2, aes(x=g.region, y=mean, fill=ms))+
  geom_bar(stat="identity", colour="black", width=0.5)+
  theme_tufte()+
  geom_rangeframe()+
  theme(text = element_text(size=20, family="CMU Serif"))+
  scale_fill_manual(values=c("black", "grey80"), name="Mating\nSystem")+
  scale_y_continuous(limits=c(0,.15))+
  labs(x="Garden Region", y="Proportion of individuals flowering in y2")

gg.flwr.y3 <- ggplot(H.flwr.inds.y3, aes(x=g.region, y=mean, fill=ms))+
  geom_bar(stat="identity", colour="black", width=0.5)+
  theme_tufte()+
  geom_rangeframe()+
  theme(text = element_text(size=20, family="CMU Serif"))+
  scale_fill_manual(values=c("black", "grey80"), name="Mating\nSystem")+
  scale_y_continuous(limits=c(0,.15))+
  labs(x="Garden Region", y="Proportion of individuals flowering in y3")

gg.flwr.10.surv.y2 <- ggplot(H.flwr.10.surv.y2, aes(x=g.region, y=mean, fill=ms))+
  geom_bar(stat="identity",position="dodge", colour="black", width=0.5)+
  theme_tufte()+
  geom_rangeframe()+
  theme(text = element_text(size=20, family="CMU Serif"))+
  scale_fill_manual(values=c("black", "grey80"), name="Mating\nSystem")+
  scale_y_continuous(limits=c(0,1))+
  labs(x="Garden Region", y="Proportion of individuals flowering by y2")

gg.flwr.10.surv.y3 <- ggplot(H.flwr.10.surv.y3, aes(x=g.region, y=mean, fill=ms))+
  geom_bar(stat="identity",position="dodge", colour="black", width=0.5)+
  theme_tufte()+
  geom_rangeframe()+
  theme(text = element_text(size=20, family="CMU Serif"))+
  scale_fill_manual(values=c("black", "grey80"), name="Mating\nSystem")+
  scale_y_continuous(limits=c(0,1))+
  labs(x="Garden Region", y="Proportion of individuals flowering by y3")

gg.flwr.mom.yr <- ggplot(H.flwr.mom, aes(x=garden, y=flower.mom, fill=ms))+
  geom_bar(stat="identity", position="dodge", colour="black", width=0.5)+
  facet_grid(year~., scales = 'free')+
  scale_y_continuous(limits=c(0,15))+
  theme_classic()+
  

gg.flwr.num.yr <- ggplot(H.flwr.num, aes(x=garden, y=flower.num, fill=ms))+
  geom_bar(stat="identity", position="dodge", colour="black", width=0.5)+
  facet_grid(year~., scales = 'free')+
  scale_y_continuous(limits=c(0,30))

gg.mean.flwr.num.yr <- ggplot(H.mean.flwr.num, aes(x=garden, y=mean.flower.num, fill=ms))+
  geom_bar(stat="identity", position="dodge", colour="black", width=0.5)+
  facet_grid(year~., scales = 'free')+
  scale_y_continuous(limits=c(0,3.5))

# gg.flwr.10.surv <- ggplot(data=H.flwr.10.surv, aes(x=g.region, y=mean, group=ms))+
#   geom_point(aes(shape=ms, fill=ms), position=position_dodge(.2), size=6)+
#   geom_line(aes(linetype=ms), position=position_dodge(.2))+
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(.2), width=.1)+
#   scale_fill_manual(labels=c("S", "A"), values=c('black', 'grey80'), name="Mating\nSystem")+
#   scale_shape_manual(labels=c("S", "A"), values=c(22, 24), name="Mating\nSystem")+
#   theme_tufte()+
#   geom_rangeframe()+
#   scale_y_continuous(limits=c(0,1))+
#   theme(text = element_text(size=20, family="CMU Serif"))+
#   labs(title="Flowering by year 3", x="Garden Region", y="Proportion of inds flowering by year 3")+
#   scale_linetype(guide=FALSE)


# gg.flwr.years <- ggplot(H.flwr.yrs, aes(x=g.region, y=count, fill=ms))+
#   geom_bar(stat="identity", colour="black", width=0.5)+
#   theme_tufte()+
#   geom_rangeframe()+
#   theme(text = element_text(size=20, family="CMU Serif"))+
#   scale_fill_manual(values=c("black", "grey80"), name="Mating\nSystem")+
#   scale_y_continuous(limits=c(0,40))+
#   labs(x="Garden Region", y="Number of individuals flowering twice")


