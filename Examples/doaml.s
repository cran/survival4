#
# These results can be found in Miller
#
fit <- coxph(Surv(aml$time, aml$status) ~ aml$x, method='breslow')
fit
resid(fit, type='mart')
resid(fit, type='score')
resid(fit, type='scho')

fit <- survfit(Surv(aml$time, aml$status) ~ aml$x)
fit
summary(fit)
survdiff(Surv(aml$time, aml$status)~ aml$x)

#
# Test out the weighted K-M
#
#  First, equal case weights- shouldn't change the survival, but will
#    halve the variance
temp2 <-survfit(Surv(aml$time, aml$status), type='kaplan', weight=rep(2,23))
temp  <-survfit(Surv(time, status)~1, aml)
temp$surv/temp2$surv
(temp$std.err/temp2$std.err)^2

# Risk weights-- use a null Cox model
summary(survfit(coxph(Surv(aml$time, aml$status) ~ offset(log(1:23)))))

# Lots of ties, so its a good test case
x1 <- coxph(Surv(time, status)~x, aml, method='efron')
x1
x2 <- coxph(Surv(rep(0,23),time, status) ~x, aml, method='efron')
x1$coef - x2$coef

