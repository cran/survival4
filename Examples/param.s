times <- c(1,3,3,6,8,9,10)
e <- c(1,1,1,0,0,1,0)
f <- survreg(Surv(times,e)~1)
print(f)
