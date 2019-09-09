# problem1
data.bioassay = data.frame(dose = c(0, 1, 2, 3, 4), 
                  n = c(2, 8, 15, 23, 27))
resp=cbind(data.bioassay$n, 30 - data.bioassay$n)
pred=data.bioassay$dose

## logit
fit=glm(resp~pred,family=binomial(link='logit'))

# CI for beta
vcov(fit) # covariance of beta MLE (fisher information inverse)
beta=fit$coefficients[2]
se=sqrt(vcov(fit)[2,2]) # (same as in summary(glm_logit))
beta+c(qnorm(0.025),0,-qnorm(0.025))*se

# diviance
dev=sum(residuals(fit,type='deviance')^2);dev

# predict
predi = predict(fit, data.frame(pred=0.01), se.fit=TRUE,type='response')$fit
exp(predi)/(1+exp(predi))

# LD50 est. and CI
beta0=fit$coefficients[1]
beta1=fit$coefficients[2]
betacov=vcov(fit) # inverse fisher information
x0fit=-beta0/beta1
exp(x0fit) # point estimate of LD50
varx0=betacov[1,1]/(beta1^2)+betacov[2,2]*(beta0^2)/(beta1^4)-2*betacov[1,2]*beta0/(beta1^3)
c(x0fit,sqrt(varx0)) # point est and se
exp(x0fit+c(qnorm(0.05),-qnorm(0.05))*sqrt(varx0)) # 90% CI for LD50

## probit
fit=glm(resp~pred,family=binomial(link='probit'))

# CI for beta
vcov(fit) # covariance of beta MLE (fisher information inverse)
beta=fit$coefficients[2]
se=sqrt(vcov(fit)[2,2]) # (same as in summary(glm_logit))
beta+c(qnorm(0.025),0,-qnorm(0.025))*se

# diviance
dev=sum(residuals(fit,type='deviance')^2);dev

# predict
predict(fit, data.frame(pred=0.01), se.fit=TRUE,type='response')$fit

# LD50 est. and CI
beta0=fit$coefficients[1]
beta1=fit$coefficients[2]
betacov=vcov(fit) # inverse fisher information
x0fit=-beta0/beta1
exp(x0fit) # point estimate of LD50
varx0=betacov[1,1]/(beta1^2)+betacov[2,2]*(beta0^2)/(beta1^4)-2*betacov[1,2]*beta0/(beta1^3)
c(x0fit,sqrt(varx0)) # point est and se
exp(x0fit+c(qnorm(0.05),-qnorm(0.05))*sqrt(varx0)) # 90% CI for LD50

## cloglog
fit=glm(resp~pred,family=binomial(link='cloglog'))

# CI for beta
vcov(fit) # covariance of beta MLE (fisher information inverse)
beta=fit$coefficients[2]
se=sqrt(vcov(fit)[2,2]) # (same as in summary(glm_logit))
beta+c(qnorm(0.025),0,-qnorm(0.025))*se

# diviance
dev=sum(residuals(fit,type='deviance')^2);dev

# predict
predict(fit, data.frame(pred=0.01), se.fit=TRUE,type='response')$fit

# LD50 est. and CI
beta0=fit$coefficients[1]
beta1=fit$coefficients[2]
betacov=vcov(fit) # inverse fisher information
x0fit=(log(-log(1-0.5))-beta0)/beta1
exp(x0fit) # point estimate of LD50
varx0=betacov[1,1]/(beta1^2)+betacov[2,2]*((beta0 - log(-log(1-0.5)))^2)/(beta1^4)-2*betacov[1,2]*(beta0 - log(-log(1-0.5)))/(beta1^3)
c(x0fit,sqrt(varx0)) # point est and se
exp(x0fit+c(qnorm(0.05),-qnorm(0.05))*sqrt(varx0)) # 90% CI for LD50

##############

data.mph = data.frame(amount = seq(from = 10000, to = 90000, by = 5000), 
                      offers = c(4, 6, 10, 12, 39, 36, 22, 14, 10, 12, 8, 9, 3, 1, 5, 2, 1), 
                      enrolls = c(0, 2, 4, 2, 12, 14, 10, 7, 5, 5, 3, 5, 2, 0, 4, 2, 1))
resp = cbind(data.mph$enrolls, data.mph$offers - data.mph$enrolls)
pred = data.mph$amount
fit_mph = glm(resp~pred,family=binomial(link='logit'))
summary(fit_mph)

pea = sum(residuals(fit_mph, type = 'pearson')^2) # pearson chisq follows chi-square(df = 8-2)
dev = sum(residuals(fit_mph, type = 'deviance')^2) # deviance (or obtain from summary(glm_logit)) 
pval = 1 - pchisq(pea,15);pval # compare with chisq(8-2), fit is ok, fails to reject

beta=fit_mph$coefficients[2]
se=sqrt(vcov(fit_mph)[2,2]) # (same as in summary(glm_logit))
beta+c(qnorm(0.025),0,-qnorm(0.025))*se

# 0.5 est. and CI
beta0=fit_mph$coefficients[1]
beta1=fit_mph$coefficients[2]
betacov=vcov(fit_mph) # inverse fisher information
x0fit=(log(0.4/0.6)-beta0)/beta1
varx0=betacov[1,1]/(beta1^2)+betacov[2,2]*((beta0 - log(0.4/0.6))^2)/(beta1^4)-2*betacov[1,2]*(beta0 - log(0.4/0.6))/(beta1^3)
c(x0fit,sqrt(varx0)) # point est and se
x0fit+c(qnorm(0.025),-qnorm(0.025))*sqrt(varx0) # 95% CI for 0.4
