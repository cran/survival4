#
#   Do the test on the simple data set
#
temp <- Surv(test1$time, test1$status)
fit <-coxph(temp~ test1$x, method='breslow')
fit
fit0 <-coxph(temp~ test1$x, iter=0)
fit0$coef
coxph(temp~ test1$x, iter=1, method='breslow')$coef
coxph(temp~ test1$x, iter=2, method='breslow')$coef
coxph(temp~ test1$x, iter=3, method='breslow')$coef

coxph(Surv(time, status) ~ x, test1, method='efron')
coxph(Surv(time, status) ~ x, test1, method='exact')

resid(fit0)
resid(coxph(temp~1))
resid(fit0, 'scor')
resid(fit0, 'scho')

resid(fit)
resid(fit, 'scor')
resid(fit, 'scho')

predict(fit, type='lp')
predict(fit, type='risk')
predict(fit, type='expected')
predict(fit, type='terms')
predict(fit, type='lp', se.fit=T)
predict(fit, type='risk', se.fit=T)
predict(fit, type='expected', se.fit=T)
predict(fit, type='terms', se.fit=T)

summary(survfit(fit))
summary(survfit(fit, list(test1=list(x=2))))
