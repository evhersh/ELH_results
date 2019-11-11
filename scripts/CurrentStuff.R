##### Current stuff #####
###to show Dolph#####

# diaspore traits
  # multiple regression
summary(ber.MR)
anova(ber.MR)
  #regressions
grid.arrange(gg.tv.aoa,gg.tv.bl,gg.tv.nb,gg.tv.weight)
  #PCA
Ber.plot
  # lmer for mating system
summary(ber.ms.lmer)
anova(ber.ms.lmer)
anova(ber.pop.lmer)
emmeans(ber.ms.lmer, 'ms')
emmeans(ber.pop.lmer, 'pop')
  # boxplots for ms and pops
gg.box.ber.ms.all
gg.box.ber.allpops
grid.arrange(gg.box.ber.ms.all, gg.box.ber.allpops)
  # lat (Yukon removed)
summary(ber.lat.lmer)
summary(ber.lat.lmer.sex)
summary(ber.lat.lmer.asex)

# Lab germ
  # glmer for ms
summary(L.germ.ms.glmer)
anova(L.germ.ms.afex)
L.germ.ms.glmer.emmeans
  # boxplots for ms and pops
gg.germ.box.ms
gg.germ.box.pop
grid.arrange(gg.germ.box.ms, gg.germ.box.pop)

  # glmer for days, ms and pop
summary(L.days.ms.afex)
anova(L.days.ms.afex)
anova(L.days.pop.afex)
  # boxplots for days (pops)
gg.days.pop.box
grid.arrange(gg.days.ms.box, gg.days.pop.box)
# Greenhouse surv
  #surv (post germ)
summary(L.surv.ms.glmer)
  #surv (including germ)
summary(L.surv.to.seedling.glmer)
  #boxplots for surv (ms and pop)
gg.surv.box.ms
gg.surv.box.pop
grid.arrange(gg.surv.box.ms,gg.surv.box.pop)
# Greenhouse leaf.length
  # models
anova(H.y0.length.ms.lmer)
anova(H.y0.length.pop.lmer)
  # boxplot for pops
grid.arrange(gg.box.length.ms.all, gg.box.length.pop)

# Greenhouse leaf.num
  # models
anova(H.y0.num.ms.afex)
anova(H.y0.num.pop.afex)
  # figure
grid.arrange(gg.box.num.ms.all, gg.box.num.pop)




# field establishment
  # model for ms*g.region
summary(glmer.msx)
emmeans(glmer.msx, 'ms')
  # model for ms*garden
summary(glmer.msx2)
anova(afex.msx2)
  # figures
gg.est.box.ms
gg.germ.all
emmip(glmer.msx2,  ~ ms | garden)






#### alberto seed ######


gg.box.ber.pops <- ggplot(data=Ber.moms.all, aes(y=terminal.velocity, x=pop))+
  geom_boxplot(aes(fill=ms))+
  labs(x="population", y="mean terminal velocity")+
  theme_classic()

gg.box.ber.ms <- ggplot(data=Ber.moms.all, aes(y=terminal.velocity, x=ms))+
  geom_boxplot(aes(fill=ms))+
  labs(x="mating system", y="mean terminal velocity")+
  theme_classic()

#ms model
ber.ms.lmer <- lmer(terminal.velocity ~ ms + (1|pop/mom), data=Ber.dat)
summary(ber.ms.lmer)
emmeans(ber.ms.lmer, 'ms')

ber.pop.lmer <- lmer(terminal.velocity~pop+(1|mom), data=Ber.dat)
anova(ber.pop.lmer)
ber.pop.lmer.emms <- emmeans(ber.pop.lmer, 'pop')
pairs(ber.pop.lmer.emms)
ber.pop.lmer.clds <-cld(ber.pop.lmer.emms)
plot(ber.pop.lmer.emms, comparisons = TRUE, horizontal=FALSE)+
  theme_classic()+
  geom_text(aes(label=ber.pop.lmer.clds$.group), vjust=-2, size=5)

ggplot(aes(y=terminal.velocity, x=weight), data=Ber.dat)+
  geom_point()

ggplot(aes(y=terminal.velocity, x=angle.of.attack), data=Ber.dat)+
  geom_point()

ggplot(aes(y=terminal.velocity, x=bristle.length), data=Ber.dat)+
  geom_point()
