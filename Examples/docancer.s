#
# Test out all of the routines on a more complex data set
#
attach(cancer)
stime <- Surv(time, status)
temp <- survfit(stime ~ ph.ecog)
summary(temp, times=c(30*1:11, 365*1:3))
print(temp[2:3])

temp <- survfit(stime, type='fleming',
		   conf.int=.9, conf.type='log-log', error='tsiatis')
summary(temp, times=30 *1:5)

temp <- survdiff(stime ~ inst, rho=.5)
print(temp, digits=6)

temp <- coxph(stime ~ ph.ecog + ph.karno + pat.karno + wt.loss + sex +
		       age + meal.cal + strata(inst))
summary(temp)
cox.zph(temp)
cox.zph(temp, transform='identity')

coxph(Surv(rep(0,length(time)), time, status) ~ ph.ecog + ph.karno + pat.karno
		+ wt.loss + sex + age + meal.cal + strata(inst))
detach(w=2)
