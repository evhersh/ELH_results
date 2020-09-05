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

##### GERM SUCCESS #####

# plot ms
gg.germ.ms <- ggplot()+
  geom_boxplot(data=L.germ.mom, aes(y=mean.germ, x=ms), alpha=0, width=0.3, outlier.shape=NA)+
  geom_jitter(data=L.germ.mom, aes(x=ms, y=mean.germ, fill=L.germ.mom$ms), alpha=0.7, shape=21, width=0.1, size=2)+
  geom_errorbar(data=germ.ggfx.ms.df, aes(ymax=conf.high, ymin=conf.low, x=as.numeric(as.factor(x))+0.4), width=0.1)+
  geom_point(data=germ.ggfx.ms.df, aes(x=as.numeric(as.factor(x))+0.4, y=predicted, fill=x), shape=21, size=3)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("Sexual", "Apomictic"))+
  labs(x="Mating system", y="Mean germination success")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("S" = "Sexual", "A" = "Apomictic"))

# plot pop
gg.germ.pop <- ggplot()+
  geom_errorbar(data=germ.ggfx.pop.df, aes(ymax=conf.high, ymin=conf.low, x=x, width=0.1))+
  geom_point(data=germ.ggfx.pop.df, aes(x=x, y=predicted, fill=group), shape=21, size=3)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  labs(x="", y="Predicted germination success")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##### GERM DAYS #####

# plot ms
gg.days.ms<-ggplot()+
  geom_boxplot(data=L.days.mom, aes(y=mean.days, x=ms), alpha=0, width=0.3, outlier.shape=NA)+
  geom_jitter(data=L.days.mom, aes(x=ms, y=mean.days, fill=L.days.mom$ms), alpha=0.7, shape=21, width=0.1)+
  geom_errorbar(data=days.ggfx.ms.df, aes(ymax=conf.high, ymin=conf.low, x=as.numeric(as.factor(x))+0.4), width=0.1)+
  geom_point(data=days.ggfx.ms.df, aes(x=as.numeric(as.factor(x))+0.4, y=predicted, fill=x), shape=21, size=3)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  labs(x="", y="Mean germination speed (days)")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("S" = "Sexual", "A" = "Apomictic"))

# plot pop
gg.days.pop <- ggplot()+
  geom_errorbar(data=days.ggfx.pop.df, aes(ymax=conf.high, ymin=conf.low, x=x, width=0.1))+
  geom_point(data=days.ggfx.pop.df, aes(x=x, y=predicted, fill=group), shape=21, size=3)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  labs(x="", y="Predicted germination speed (days)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##### GSURV #####

# plot ms
gg.gsurv.ms<-ggplot()+
  geom_boxplot(data=L.surv.mom, aes(y=mean.surv, x=ms), alpha=0, width=0.3, outlier.shape=NA)+
  geom_jitter(data=L.surv.mom, aes(x=ms, y=mean.surv, fill=L.surv.mom$ms), alpha=0.7, shape=21, width=0.1)+
  geom_errorbar(data=gsurv.ggfx.ms.df, aes(ymax=conf.high, ymin=conf.low, x=as.numeric(as.factor(x))+0.4), width=0.1)+
  geom_point(data=gsurv.ggfx.ms.df, aes(x=as.numeric(as.factor(x))+0.4, y=predicted, fill=x), shape=21, size=3)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  labs(x="", y="Mean seedling survival")+
  theme_bw()+
  ylim(0.65,1)+
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("S" = "Sexual", "A" = "Apomictic"))

# plot pop
gg.gsurv.pop <- ggplot()+
  geom_errorbar(data=gsurv.ggfx.pop.df, aes(ymax=conf.high, ymin=conf.low, x=x, width=0.1))+
  geom_point(data=gsurv.ggfx.pop.df, aes(x=x, y=predicted, fill=group), shape=21, size=3)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  labs(x="Population", y="Predicted greenhouse survival")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##### LENGTH #####

# plot ms
gg.length.ms<-ggplot()+
  geom_boxplot(data=H.y0.mom, aes(y=mean.length, x=ms), alpha=0, width=0.3, outlier.shape=NA)+
  geom_jitter(data=H.y0.mom, aes(x=ms, y=mean.length, fill=H.y0.mom$ms), alpha=0.7, shape=21, width=0.1)+
  geom_errorbar(data=length.ggfx.ms.df, aes(ymax=conf.high, ymin=conf.low, x=as.numeric(as.factor(x))+0.4), width=0.1)+
  geom_point(data=length.ggfx.ms.df, aes(x=as.numeric(as.factor(x))+0.4, y=predicted, fill=x), shape=21, size=3)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  labs(x="Mating system", y="Mean seedling leaf length (mm)")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("S" = "Sexual", "A" = "Apomictic"))

# plot pop
gg.length.pop <- ggplot()+
  geom_errorbar(data=length.ggfx.pop.df, aes(ymax=conf.high, ymin=conf.low, x=x, width=0.1))+
  geom_point(data=length.ggfx.pop.df, aes(x=x, y=predicted, fill=group), shape=21, size=3)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  labs(x="Population", y="Predicted seedling leaf length")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##### NUM #####

# plot ms
gg.num.ms<-ggplot()+
  geom_boxplot(data=H.y0.mom, aes(y=mean.num, x=ms), alpha=0, width=0.3, outlier.shape=NA)+
  geom_jitter(data=H.y0.mom, aes(x=ms, y=mean.num, fill=H.y0.mom$ms), alpha=0.7, shape=21, width=0.1)+
  geom_errorbar(data=num.ggfx.ms.df, aes(ymax=conf.high, ymin=conf.low, x=as.numeric(as.factor(x))+0.4), width=0.1)+
  geom_point(data=num.ggfx.ms.df, aes(x=as.numeric(as.factor(x))+0.4, y=predicted, fill=x), shape=21, size=3)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  labs(x="Mating system", y="Mean seedling leaf number")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("S" = "Sexual", "A" = "Apomictic"))

# plot pop
gg.num.pop <- ggplot()+
  geom_errorbar(data=num.ggfx.pop.df, aes(ymax=conf.high, ymin=conf.low, x=x, width=0.1))+
  geom_point(data=num.ggfx.pop.df, aes(x=x, y=predicted, fill=group), shape=21, size=3)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  labs(x="Population", y="Predicted leaf number")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##### TERMINAL VELOCITY #####

# plot ms
gg.tv.ms<-ggplot()+
  geom_boxplot(data=Ber.dat, aes(y=terminal.velocity, x=ms), alpha=0, width=0.3, outlier.shape=NA)+
  geom_jitter(data=Ber.dat, aes(x=ms, y=terminal.velocity, fill=Ber.dat$ms), alpha=0.7, shape=21, width=0.1, size=2)+
  geom_errorbar(data=tv.ggfx.ms.df, aes(ymax=conf.high, ymin=conf.low, x=as.numeric(as.factor(x))+0.4), width=0.1)+
  geom_point(data=tv.ggfx.ms.df, aes(x=as.numeric(as.factor(x))+0.4, y=predicted, fill=x), shape=21, size=3)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("Sexual", "Apomictic"))+
  labs(x="Mating system", y="Terminal velocity (cm/s)")+
  scale_x_discrete(labels=c("S" = "Sexual", "A" = "Apomictic"))+
  theme_bw()

# plot pop
gg.tv.pop <- ggplot()+
  geom_errorbar(data=tv.ggfx.pop.df, aes(ymax=conf.high, ymin=conf.low, x=x, width=0.1))+
  geom_point(data=tv.ggfx.pop.df, aes(x=x, y=predicted, fill=group), shape=21, size=3)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  labs(x="", y="Predicted terminal velocity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##### Diaspore traits #####

gg.tv.weight<-ggplot(aes(y=terminal.velocity, x=weight.scaled, fill=ms), data=Ber.dat)+
  geom_point(shape=21)+
  geom_smooth(method='lm', se=FALSE, aes(color=ms))+
  theme_classic()+
  labs(x="Seed mass (scaled)", y="Terminal velocity (cm/s)")+
  scale_colour_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  theme(legend.position = "none")

gg.tv.aoa<-ggplot(aes(y=terminal.velocity, x=angle.of.attack.scaled, fill=ms), data=Ber.dat)+
  geom_point(shape=21)+
  geom_smooth(method='lm', se=FALSE, aes(color=ms))+
  theme_classic()+
  labs(x="Angle of attack (scaled)", y="Terminal velocity (cm/s)")+
  scale_colour_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  theme(legend.position = "none")

gg.tv.bl<-ggplot(aes(y=terminal.velocity, x=bristle.length.scaled, fill=ms), data=Ber.dat)+
  geom_point(shape=21)+
  geom_smooth(method='lm', se=FALSE, aes(color=ms))+
  theme_classic()+
  labs(x="Bristle length (scaled)", y="Terminal velocity (cm/s)")+
  scale_colour_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  theme(legend.position = "none")


gg.tv.nb<-ggplot(aes(y=terminal.velocity, x=number.of.bristles.scaled, fill=ms), data=Ber.dat)+
  geom_point(shape=21)+
  geom_smooth(method='lm', se=FALSE, aes(color=ms))+
  theme_classic()+
  labs(x="Number of bristles (scaled)", y="Terminal velocity (cm/s)")+
  scale_colour_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("sexual", "apomict"))+
  theme(legend.position = "none")


#### diaspore by ms plots ####

gg.aoa.ms <- ggplot()+
  geom_boxplot(data=Ber.dat, aes(y=angle.of.attack, x=ms), alpha=0, width=0.3, outlier.shape=NA)+
  geom_jitter(data=Ber.dat, aes(x=ms, y=angle.of.attack, fill=Ber.dat$ms), alpha=0.7, shape=21, width=0.1, size=2)+
  geom_errorbar(data=aot.ggfx, aes(ymax=conf.high, ymin=conf.low, x=as.numeric(as.factor(x))+0.4), width=0.1)+
  geom_point(data=aot.ggfx, aes(x=as.numeric(as.factor(x))+0.4, y=predicted, fill=x), shape=21, size=3)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("Sexual", "Apomictic"))+
  labs(x="", y="Angle of attack")+
  scale_x_discrete(labels=c("S" = "Sexual", "A" = "Apomictic"))+
  theme_bw()

gg.bl.ms <- ggplot()+
  geom_boxplot(data=Ber.dat, aes(y=bristle.length, x=ms), alpha=0, width=0.3, outlier.shape=NA)+
  geom_jitter(data=Ber.dat, aes(x=ms, y=bristle.length, fill=Ber.dat$ms), alpha=0.7, shape=21, width=0.1, size=2)+
  geom_errorbar(data=bl.ggfx, aes(ymax=conf.high, ymin=conf.low, x=as.numeric(as.factor(x))+0.4), width=0.1)+
  geom_point(data=bl.ggfx, aes(x=as.numeric(as.factor(x))+0.4, y=predicted, fill=x), shape=21, size=3)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("Sexual", "Apomictic"))+
  labs(x="", y="Bristle length (mm)")+
  scale_x_discrete(labels=c("S" = "Sexual", "A" = "Apomictic"))+
  theme_bw()

gg.nob.ms <- ggplot()+
  geom_boxplot(data=Ber.dat, aes(y=number.of.bristles, x=ms), alpha=0, width=0.3, outlier.shape=NA)+
  geom_jitter(data=Ber.dat, aes(x=ms, y=number.of.bristles, fill=Ber.dat$ms), alpha=0.7, shape=21, width=0.1, size=2)+
  geom_errorbar(data=nob.ggfx, aes(ymax=conf.high, ymin=conf.low, x=as.numeric(as.factor(x))+0.4), width=0.1)+
  geom_point(data=nob.ggfx, aes(x=as.numeric(as.factor(x))+0.4, y=predicted, fill=x), shape=21, size=3)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("Sexual", "Apomictic"))+
  labs(x="Mating system", y="Number of bristles")+
  scale_x_discrete(labels=c("S" = "Sexual", "A" = "Apomictic"))+
  theme_bw()

gg.weight.ms <- ggplot()+
  geom_boxplot(data=Ber.dat, aes(y=weight, x=ms), alpha=0, width=0.3, outlier.shape=NA)+
  geom_jitter(data=Ber.dat, aes(x=ms, y=weight, fill=Ber.dat$ms), alpha=0.7, shape=21, width=0.1, size=2)+
  geom_errorbar(data=weight.ggfx, aes(ymax=conf.high, ymin=conf.low, x=as.numeric(as.factor(x))+0.4), width=0.1)+
  geom_point(data=weight.ggfx, aes(x=as.numeric(as.factor(x))+0.4, y=predicted, fill=x), shape=21, size=3)+
  scale_fill_manual(values=c("coral3","cornflowerblue"), name="Mating system", labels=c("Sexual", "Apomictic"))+
  labs(x="Mating system", y="Seed mass (g)")+
  scale_x_discrete(labels=c("S" = "Sexual", "A" = "Apomictic"))+
  theme_bw()

##### arrange #####

png("diaspore.png", height=7, width=9, res=300, units="in")
ggarrange(gg.tv.ms, nrow=1, labels="A", ggarrange(gg.tv.aoa,gg.tv.bl,gg.tv.nb,gg.tv.weight, ncol=2, nrow=2, labels=c("B", "C", "D", "E")),common.legend = TRUE, legend="bottom")
dev.off()

# ms big plot
png("lab.png", height=7, width=8, res=300, units="in")
#ggarrange(gg.germ.ms, nrow=1, labels="A", ggarrange(gg.days.ms,gg.gsurv.ms,gg.length.ms,gg.num.ms, ncol=2, nrow=2, labels=c("B", "C", "D", "E")),common.legend = TRUE, legend="bottom")
ggarrange(gg.germ.ms,gg.days.ms,gg.gsurv.ms,gg.length.ms,gg.num.ms, ncol=3, nrow=2, labels=c("A", "B", "C", "D", "E"),common.legend = TRUE, legend="bottom")
dev.off()

# ms big plot 2
png("lab2.png", height=5, width=8, res=300, units="in")
#ggarrange(gg.germ.ms, nrow=1, labels="A", ggarrange(gg.days.ms,gg.gsurv.ms,gg.length.ms,gg.num.ms, ncol=2, nrow=2, labels=c("B", "C", "D", "E")),common.legend = TRUE, legend="bottom")
ggarrange(gg.germ.ms,gg.days.ms,gg.gsurv.ms, ncol=3, nrow=1, labels=c("A", "B", "C"),common.legend = TRUE, legend="bottom")
dev.off()

# pop big plot
png("pops.png", height=7, width=8, res=300, units="in")
ggarrange(gg.tv.pop, gg.germ.pop, gg.days.pop, gg.gsurv.pop, gg.length.pop, gg.num.pop, common.legend = TRUE, legend="bottom", labels=c("A", "B", "C", "D", "E", "F"))
dev.off()

# diaspore by ms
png("diaspore_ms.png", height=7, width=8, res=300, units="in")
ggarrange(gg.aoa.ms, gg.bl.ms, gg.nob.ms, gg.weight.ms, common.legend = TRUE, legend="bottom", labels=c("A", "B", "C", "D"))
dev.off()

# arrange for defence
png("TV_germ.png", height=7, width=9, res=300, units="in")
ggarrange(gg.tv.ms, gg.germ.ms, common.legend = TRUE, legend="bottom", labels=c("A", "B"))
dev.off()