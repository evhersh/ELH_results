############
# Packages
############
library(car)
library(ggplot2)
library(pwr)
library(Rmisc)
#library(plyr)
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
library(ggpmisc)

library(extrafont)
#font_import()
#loadfonts()








############
# Data Prep
############
Berto.dat <<-read.csv("~/GitHub/Hookeri-Gardens/RAW DATA/alberto.seed.R.csv", stringsAsFactors = FALSE, strip.white = TRUE, na.string = c("NA",""))


# selecting just the main columns that I'm interested in
df<-data.frame("")

Ber.dat <- Berto.dat %>%
  filter(mom!="Bulk")

Ber.sub <- Berto.dat%>%
  filter(mom!="Bulk")%>%
  filter(pop %in% c("B42", "B46", "B49", "L12", "B52", "L16","L18","L20", "L40", "L41"))


# factors and levels
Ber.dat$pop <- factor(Ber.dat$pop, levels=c("B42", "B46", "B49", "L12", "C59", "B52", "L16","L18","L20", "L40", "L41", "C21", "C43", "C54"))
Ber.dat$region <- as.factor(Ber.dat$region)
Ber.dat$garden <- as.factor(Ber.dat$garden)
Ber.dat$ms <- factor(Ber.dat$ms)
Ber.dat$ms <- factor(Ber.dat$ms, levels = c("S", "A"))
Ber.dat$mom <- as.factor(Ber.dat$mom)


#log transform
Ber.dat$terminal.velocity.log <- log(Ber.dat$terminal.velocity)
Ber.dat$weight.log <- log(Ber.dat$weight)
Ber.dat$angle.of.attack.log <- log(Ber.dat$angle.of.attack)
Ber.dat$bristle.length.log <- log(Ber.dat$bristle.length)
Ber.dat$number.of.bristles.log <- log(Ber.dat$number.of.bristles)

#Ber.dat$pop <- factor(Ber.dat$pop, levels = rev(levels(factor(Ber.dat$pop)))) #reverse order for figures (south to north)

# remove pops
Ber.dat.sub <- Ber.dat %>%
  group_by(ms,region,pop,mom,lat) %>%
  filter(pop %in% c("B42", "B46", "B49", "L12", "B52", "L16","L18","L20", "L40", "L41")) %>%
  select(pop:number.of.bristles, lat) %>%
  summarize_all(funs(mean))

Ber.dat.sub.pop <- Ber.dat %>%
  group_by(ms,region,pop,lat) %>%
  filter(pop %in% c("B42", "B46", "B49", "L12", "B52", "L16","L18","L20", "L40", "L41")) %>%
  select(pop:number.of.bristles, lat) %>%
  summarize_all(funs(mean))

# Means for moms
Ber.dat.all.pop <- Ber.dat %>%
  group_by(ms,region,pop,lat) %>%
  select(pop:number.of.bristles, lat) %>%
  summarize_all(funs(mean))

Ber.dat.all.mom <- Ber.dat %>%
  group_by(ms,region,pop,mom,lat) %>%
  select(pop:number.of.bristles, lat) %>%
  summarize_all(funs(mean))

Ber.pop.sub <- Ber.dat.sub %>%
  group_by(ms,region,pop,lat) %>%
  summarize_all(funs(mean)) %>%
  arrange(desc(ms), lat)

Ber.dat.ms <- Ber.dat %>%
  group_by(ms) %>%
  filter(!is.na(terminal.velocity))%>%
  filter(terminal.velocity!=0)%>%
  summarize(terminal.velocity=mean(terminal.velocity), se=std.error(terminal.velocity, na.rm=TRUE))


############
# Models
############
# multiple regression
ber.MR <- lmer(terminal.velocity.log ~ angle.of.attack.log + bristle.length.log + number.of.bristles.log + weight.log + (1|pop/mom), data=Ber.dat)
summary(ber.MR)
anova(ber.MR)
coefficients(ber.MR)
resid(ber.MR)
coef(summary(ber.MR))[,"Estimate"]

visreg(ber.MR)

layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page 
plot(ber.MR)
VarCorr(ber.MR)
confint(ber.MR)
fitted(ber.MR)

# ms all
ber.ms.lmer <- lmer(terminal.velocity ~ ms + (1|pop/mom), data=Ber.dat)
summary(ber.ms.lmer)
anova(ber.ms.lmer)
emmeans(ber.ms.lmer, 'ms')

summary(lmer(angle.of.attack~ms+(1|pop/mom), data=Ber.dat))
summary(lmer(bristle.length~ms+(1|pop/mom), data=Ber.dat))
summary(lmer(number.of.bristles~ms+(1|pop/mom), data=Ber.dat))
summary(lmer(weight~ms+(1|pop/mom), data=Ber.dat))


ber.lat.lmer <- lmer(terminal.velocity.log ~ lat + (1|pop/mom), data=subset(Ber.dat, region!="Yukon"))
summary(ber.lat.lmer)
plot(ber.lat.lmer)

ber.lat.lmer.sex <- lmer(terminal.velocity.log ~ lat + (1|pop/mom), data=subset(Ber.dat, region!="Yukon" & ms=="S"))
summary(ber.lat.lmer.sex)
plot(ber.lat.lmer.sex)

ber.lat.lmer.asex <- lmer(terminal.velocity.log ~ lat + (1|pop/mom), data=subset(Ber.dat, region!="Yukon" & ms=="A"))
summary(ber.lat.lmer.asex)
plot(ber.lat.lmer.asex)

# ms subset
ber.ms.sub.lmer <- lmer(terminal.velocity ~ ms + (1|pop/mom), data=Ber.sub)
summary(ber.ms.sub.lmer)
ber.ms.emms<-emmeans(ber.ms.sub.lmer, 'ms')
ber.ms.clds<-cld(ber.ms.emms)

ber.pop.lmer <- lmer(terminal.velocity.log~pop+(1|mom), data=Ber.dat)
anova(ber.pop.lmer)
ber.pop.lmer.emms <- emmeans(ber.pop.lmer, 'pop')
pairs(ber.pop.lmer.emms)
ber.pop.lmer.clds <-cld(ber.pop.lmer.emms, sort=FALSE)
plot(ber.pop.lmer.emms, comparisons = TRUE, horizontal=FALSE)+
  theme_classic()+
  geom_text(aes(label=ber.pop.lmer.clds$.group), vjust=-2, size=5)

# pop subset
ber.pop.sub.lmer <- lmer(terminal.velocity~pop+(1|mom), data=Ber.dat.sub)
anova(ber.pop.sub.lmer)
ber.pop.sub.lmer.emms <- emmeans(ber.pop.sub.lmer, 'pop')
pairs(ber.pop.sub.lmer.emms)
ber.pop.sub.lmer.clds <-cld(ber.pop.sub.lmer.emms, sort=FALSE)
plot(ber.pop.sub.lmer.emms, comparisons = TRUE, horizontal=FALSE)+
  theme_classic()+
  geom_text(aes(label=ber.pop.sub.lmer.clds$.group), vjust=-2, size=5)


#############
# Figures
#############
#ms subset


gg.box.ber.ms.all <- ggplot(data=Ber.dat, aes(y=terminal.velocity, x=ms))+
  geom_boxplot(aes(fill=ms), width=0.3, outlier.shape=NA)+
  scale_fill_manual(values=c("white", "darkgrey"))+
  labs(x="mating system", y="mean terminal velocity (cm/s)")+
  theme_classic()

#all pops
gg.box.ber.allpops <- ggplot(data=Ber.dat, aes(y=terminal.velocity, x=pop))+
  geom_boxplot(aes(fill=ms), width=0.3, outlier.shape = NA)+
  geom_jitter(alpha=1/5, (aes(as.numeric(as.factor(pop)) + 0.3, terminal.velocity)), width=0.1)+
  scale_fill_manual(values=c("white","darkgrey"))+
  labs(x="population", y="mean terminal velocity (cm/s)")+
  theme_classic()+
  stat_summary(geom='text', label=ber.pop.lmer.clds$.group, vjust=-9, size=5)
  

grid.arrange(gg.box.ber.ms.all, gg.box.ber.allpops)



#try to make space between ms
gg.box.ber.allpops2 <- ggplot(data=Ber.dat, aes(y=terminal.velocity, x=pop))+
  geom_boxplot(aes(fill=ms), width=0.3, outlier.shape = NA, position=position_dodge(4))+
  geom_jitter(alpha=1/5, (aes(as.numeric(as.factor(pop)) + 0.3, terminal.velocity)), width=0.1)+
  scale_fill_manual(values=c("white","darkgrey"))+
  labs(x="population", y="mean terminal velocity (cm/s)")+
  stat_summary(geom='text', label=ber.pop.lmer.clds$.group, vjust=-9, size=5)+
  scale_x_discrete(breaks=Ber.dat$pop[nchar(as.character(Ber.dat$pop))!=1])+
  theme_classic()

gg.box.ber.ms.sub <- ggplot(data=Ber.dat.sub.pop, aes(y=terminal.velocity, x=ms))+
  geom_boxplot(aes(fill=ms), width=0.5)+
  geom_jitter(alpha=1/5, width=.25)+
  scale_fill_manual(values=c("white", "darkgrey"))+
  labs(x="mating system", y="mean terminal velocity")+
  scale_y_continuous(limits=c(0,2))+
  theme_classic()

#pop subset
gg.box.ber.sub.pops <- ggplot(data=Ber.dat.sub, aes(y=terminal.velocity, x=pop))+
  geom_boxplot(aes(fill=ms, color=region))+
  geom_jitter(alpha=1/5, width=.15)+
  scale_fill_manual(values=c("white", "darkgrey"))+
  labs(x="population", y="mean terminal velocity")+
  stat_summary(geom='text', label=ber.pop.sub.lmer.clds$.group, vjust=-3, size=5)+
  scale_y_continuous(limits=c(0,3))+
  theme_classic()

#####facet for regressions#####
#function that allows annotation of ggplot with lm equation and r2
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}


gg.tv.weight<-ggplot(aes(y=terminal.velocity.log, x=weight.log, shape=ms, linetype=ms), data=Ber.dat)+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)+
  scale_shape_manual(values=c(1,16))+
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  theme_classic()+
  labs(x="log seed mass (g)", y="log terminal velocity (cm/s)")+
  scale_linetype_manual(values=c(2,1))

gg.tv.aoa<-ggplot(aes(y=terminal.velocity.log, x=angle.of.attack.log, shape=ms, linetype=ms), data=Ber.dat)+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)+
  scale_shape_manual(values=c(1,16))+
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  theme_classic()+
  labs(x="log angle of attack", y="log terminal velocity (cm/s)")+
  scale_linetype_manual(values=c(2,1))

gg.tv.bl<-ggplot(aes(y=terminal.velocity.log, x=bristle.length.log, shape=ms, linetype=ms), data=Ber.dat)+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)+
  scale_shape_manual(values=c(1,16))+
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  theme_classic()+
  labs(x="log bristle length (mm)", y="log terminal velocity (cm/s)")+
  scale_linetype_manual(values=c(2,1))

gg.tv.nb<-ggplot(aes(y=terminal.velocity.log, x=number.of.bristles.log, shape=ms, linetype=ms), data=Ber.dat)+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)+
  scale_shape_manual(values=c(1,16))+
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  theme_classic()+
  labs(x="log number of bristles", y="log terminal velocity (cm/s)")+
  scale_linetype_manual(values=c(2,1))

gg.tv.lat <- ggplot(aes(y=terminal.velocity, x=lat), data=subset(Ber.dat, region!="Yukon"))+
  geom_point(aes(shape=ms))+
  geom_smooth(method='lm', se=FALSE)+
  scale_shape_manual(values=c(1,16))+
  theme_classic()

gg.tv.lat.sex <- ggplot(aes(y=terminal.velocity, x=lat), data=subset(Ber.dat, ms=="S"))+
  geom_point(aes(shape=ms))+
  geom_smooth(method='lm', se=FALSE)+
  scale_shape_manual(values=c(1,16))+
  theme_classic()

gg.tv.lat.asex <- ggplot(aes(y=terminal.velocity, x=lat), data=subset(Ber.dat, ms=="A"))+
  geom_point(aes(shape=ms))+
  geom_smooth(method='lm', se=FALSE)+
  scale_shape_manual(values=c(1,16))+
  theme_classic()

########

gg.tv.regressions <- grid.arrange(gg.tv.aoa,gg.tv.bl,gg.tv.nb,gg.tv.weight)

# # Quick PCA
Ber.log <- log(Ber.dat[, 8:12])
Ber.pca <- prcomp(Ber.log, center=TRUE, scale=TRUE)
Ber.plot <- autoplot(Ber.pca, data=Ber.dat, colour='ms', loadings=TRUE, loadings.colour="blue", loadings.label=TRUE, loadings.label.colour="black", loadings.label.size=5, frame=TRUE, frame.type='norm')

##### OLD SCRIPT STUFF TO MINE #######

# ############
# # Packages
# ############
# library(car)
# library(ggplot2)
# library(pwr)
# library(Rmisc)
# library(plyr)
# library(visreg)
# library(multcomp)
# library(lme4)
# library(nlme)
# library(dplyr)
# library(tidyr)
# library(EML)
# library(plotrix)
# library(beanplot)
# library(R2admb)
# library(glmmADMB)
# library(MCMCglmm)
# library(ggthemes)
# library(RColorBrewer)
# library(pegas)
# library(rJava)
# library(xtable)
# library(stargazer)
# library(ggfortify)
# library(cluster)
#
# library(extrafont)
# #font_import()
# #loadfonts()
#
# ############
# # Load Berto's data
# ############
# Berto.dat <<-read.csv("~/GitHub/Hookeri-Gardens/RAW DATA/alberto.seed.R.csv", stringsAsFactors = FALSE, strip.white = TRUE, na.string = c("NA",""))
#
#
# ############
# # Organizing data
# ############
#
# Ber.moms.garden <- Ber.dat %>% # Only one asex population from here was used in the gardens...
#   filter(garden=="yes") %>%
#   select(pop:number.of.bristles, lat) %>%
#   group_by(pop, ms,region, mom) %>%
#   summarize_all(funs(mean))

# # selecting just the main columns that I'm interested in
# Ber.dat <- Berto.dat %>%
#   select(pop:number.of.bristles) %>%
#   select(-ploidy, -rep)
#
# # factors and levels
# Ber.dat$pop <- as.factor(Ber.dat$pop)
# Ber.dat$ms <- factor(Ber.dat$ms, levels=c("S", "A"))
#
#
# # recreating table of means for all the pops, sorted by ms and pop
# Ber.means <- Ber.dat %>%
#   group_by(pop, ms) %>%
#   summarize_all(funs(mean)) %>%
#   arrange(ms, pop)
#
# # same table but just comparing means of sex and asex
# Ber.means.ms <- Ber.dat %>%
#   group_by(ms) %>%
#   select(-pop) %>%
#   summarize_all(funs(mean)) %>%
#   arrange(desc(ms))
#
# # Terminal velocity means
# Ber.means.TV <- Ber.dat %>%
#   group_by(ms) %>%
#   summarize(mean.terminal.velocity=mean(terminal.velocity), se.terminal.velocity=std.error(terminal.velocity))
#
# # log transform and separate variables
# Ber.log <- log(Ber.dat[, 3:7])
# Ber.pops <- Ber.dat[, 1]
#
# # correlation matrix
# Ber.matrix <- cor(Ber.log)
# Ber.matrix.4 <- round(Ber.matrix, 4)
#
# ###########
# Models (PCA, etc)
# ###########

# Ber.pca <- prcomp(Ber.log, center=TRUE, scale=TRUE)
# #plot(Ber.pca, type="l")
# Ber.pca.sum<- summary(Ber.pca)
#
# # plot terminal velocity against angle of attack
#
# Ber.tv.lm.all <- lm(terminal.velocity ~ angle.of.attack + weight + bristle.length + number.of.bristles, data=Ber.dat)
# summary(Ber.tv.lm.all, correlation=TRUE)
# anova(Ber.tv.lm.all)
# Ber.dat.matrix <- Ber.dat[,3:length(Ber.dat)]
# Ber.cor <- round(cor(Ber.dat.matrix), 3)
#
# Ber.tv.lm.ms <- lm(terminal.velocity ~ ms, data=Ber.dat)
# summary(Ber.lm.tv)
#
# ############
# # Figures
############

# Ber.means.TV$ms <- factor(Ber.means.TV$ms, levels = rev(levels(factor(Ber.means.TV$ms))))
#
# # mean TV
# gg.Ber.tv <- ggplot(data=Ber.means.TV, aes(x=ms, y=mean.terminal.velocity, group=ms))+
#   geom_point(aes(shape=ms, fill=ms), position=position_dodge(.2), size=6)+
#   geom_errorbar(aes(ymin=mean.terminal.velocity-se.terminal.velocity, ymax=mean.terminal.velocity+se.terminal.velocity),position=position_dodge(.2), width=.1)+
#   scale_fill_manual(labels=c("Sex", "Asex"), values=c('black', 'grey80'), name="Mating\nSystem")+
#   scale_shape_manual(labels=c("Sex", "Asex"), values=c(22, 25), name="Mating\nSystem")+
#   theme_tufte()+
#   geom_rangeframe()+
#   scale_y_continuous(limits=c(0,1.5))+
#   theme(text = element_text(size=20, family="CMU Serif"))+
#   labs(x="Mating System", y="Mean Terminal Velocity")+
#   scale_linetype(guide=FALSE)
#
# # plot of relationship between terminal velocity and angle of attack / bristle length
#
# #function that allows annotation of ggplot with lm equation and r2
# lm_eqn = function(m) {
#
#   l <- list(a = format(coef(m)[1], digits = 2),
#             b = format(abs(coef(m)[2]), digits = 2),
#             r2 = format(summary(m)$r.squared, digits = 3));
#
#   if (coef(m)[2] >= 0)  {
#     eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
#   } else {
#     eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
#   }
#
#   as.character(as.expression(eq));
# }
#
#
# gg.tv.angle<- ggplot(data=Ber.dat, aes(x=angle.of.attack, y=terminal.velocity, colour=ms))+
#   geom_point()+
#   geom_smooth(aes(group=1), method=lm)+
#   annotate("text", x = 150, y = 2.5, label = lm_eqn(lm(terminal.velocity ~ angle.of.attack, Ber.dat)),colour="black",size=5, parse = TRUE)+
#   scale_colour_few()+
#   theme_few()
#
# gg.tv.length<- ggplot(data=Ber.dat, aes(x=bristle.length, y=terminal.velocity, colour=ms))+
#   geom_point()+
#   geom_smooth(aes(group=1), method=lm)+
#   annotate("text", x = 9, y = 2.5, label = lm_eqn(lm(terminal.velocity ~ bristle.length, Ber.dat)),colour="black",size=5, parse = TRUE)+
#   scale_colour_few()+
#   theme_few()
#
#
#
#
#
# # Quick PCA
 Ber.plot <- autoplot(Ber.pca, data=Ber.dat, colour='ms', loadings=TRUE, loadings.colour="blue", loadings.label=TRUE, loadings.label.colour="black", loadings.label.size=5, frame=TRUE, frame.type='norm')
#
#
#
# # # same table but just comparing means of sex and asex
# Ber.means.ms <- Ber.dat %>%
#   group_by(ms) %>%
#   select(-pop) %>%
#   summarize_all(funs(mean)) %>%
#   arrange(desc(ms))
#
# # Terminal velocity means
# Ber.means.TV <- Ber.dat %>%
#   group_by(ms) %>%
#   summarize(mean.terminal.velocity=mean(terminal.velocity), se.terminal.velocity=std.error(terminal.velocity))
#
# # log transform and separate variables
# Ber.log <- log(Ber.dat[c(3,4,6,7)])
# Ber.pops <- Ber.dat[, 1]
#
# # correlation matrix
# Ber.matrix <- cor(Ber.log)
# Ber.matrix.4 <- round(Ber.matrix, 4)
# #
