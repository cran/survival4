#
# Run a test that can be verified using SAS's LIFEREG
#
fit1w <- survreg(Surv(time, status) ~x, test1, link='log', dist='extreme')
fit1w
summary(fit1w)

fit1e <- survreg(Surv(time, status) ~x, test1, link='log', dist='exp')
fit1e
summary(fit1e)

fit1l <- survreg(Surv(time, status) ~x, test1, link='log', dist='logistic')
fit1l
summary(fit1l)

fit1g <- survreg(Surv(time, status) ~x, test1, link='log', dist='gaussian')
summary(fit1g)
#
#  Do a test with the ovarian data
#
fitfw <- survreg(Surv(futime, fustat) ~ age + ecog.ps, ovarian,
	link='log', dist='extreme')
fitfw

fitfl <- survreg(Surv(futime, fustat) ~ age + ecog.ps, ovarian,
	link='log', dist='logistic')
fitfl

flem2 <- scan("fleming.data2", what=list(ltime=0, rtime=0))
flsurv<- Surv(flem2$ltime, flem2$rtime, type='interval2')

fitfw2 <- survreg(flsurv ~ age + ecog.ps, ovarian,
	 link='log', dist='extreme')
summary(fitfw2)

fitfl2 <- survreg(flsurv ~ age + ecog.ps, ovarian,
	 link='log', dist='logistic')
summary(fitfl2)

fitfg2 <- survreg(flsurv ~ age + ecog.ps, ovarian,
	 link='log', dist='gaussian')
summary(fitfg2)
