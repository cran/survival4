#
#   Do the test on the simple agreg data set
#
temp <- Surv(test2$start, test2$stop, test2$event)
fit <-coxph(temp~ x, test2, method='breslow')
fit
fit0 <-coxph(temp~ x, test2, iter=0)
fit0$coef
coxph(temp~ x, test2, iter=1, method='breslow')$coef
coxph(temp~ x, test2, iter=2, method='breslow')$coef
coxph(temp~ x, test2, iter=3, method='breslow')$coef

coxph(temp ~ test2$x, method='efron')
coxph(temp ~ test2$x, method='exact')

resid(fit0)
resid(fit0, 'scor')
resid(fit0, 'scho')

resid(fit)
resid(fit, 'scor')
resid(fit, 'scho')

resid(coxph(temp~ x, test2, iter=0, init=log(2)), 'score')

sfit <-survfit(fit)
sfit
summary(sfit)

temp <- Surv(rep(test2$start,2), rep(test2$stop,2), rep(test2$event,2))
ss <- rep(0:1, rep(length(test2$x),2))
options(contrasts=c("contr.treatment", "contr.treatment"))
fitx <- coxph(temp ~ c(test2$x, test2$x^2)*strata(ss), method='breslow')
sfit <- survfit(fitx, c(fitx$means[1], 0) )
sfit
summary(sfit)
#
# Even though everyone in strata 1 has x2==0, I won't get the same survival
#  curve above if survfit is called without forcing predicted x2 to be
#  zero-- otherwise I am asking for a different baseline than the
#  simple model did.  In this particular case coef[2] is nearly zero, so
#  the curves are the same, but the variances differ.
#

# This mimics 'byhand3' and the documentation
fit <- coxph(Surv(start, stop, event) ~x, test2, method='breslow')
tdata <- data.frame( start=c(0,20, 6,0,5),
		     stop =c(5,23,10,5,15),
		     event=rep(0,5),
		     x=c(0,0,1,0,2) )

temp <- survfit(fit, tdata, individual=T)
temp2 <- byhand3(fit$coef)
all.equal(temp$surv, temp2$surv)
all.equal(temp2$std, temp$std.err)

temp2
