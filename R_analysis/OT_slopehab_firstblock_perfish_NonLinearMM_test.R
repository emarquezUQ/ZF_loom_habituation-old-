
##### testing how to fit a nonlinear model 


###### basing how to fit a esponential model from https://rpubs.com/mengxu/exponential-model

######################## Prepare fitting a nonlinear model ####################

#Select an approximate $\theta$, since theta must be lower than min(y), and greater than zero
theta.0 <- min(Firstblock$resp) * 0.5 

# Estimate the rest parameters using a linear model
model.0 <- lm(log(resp - theta) ~ loom, data=Firstblock)  
alpha.0 <- exp(coef(model)[1])
beta.0 <- coef(model)[2]

# Starting parameters
start1 <- list(alpha=alpha.0,beta=beta.0,theta=theta.0)
start1

######################Fit the model (with estimated starting parameters) #################

model <- nls(resp ~ alpha*exp(beta*loom) + theta , data = Firstblock, start = start1) ### not working...

# Plot fitted curve
plot(Firstblock$loom, Firstblock$resp)
lines(Firstblock$loom, predict(model, list(loom = Firstblock$loom)), col = 'skyblue', lwd = 3)


summary(model)

##### i dont manage to make it work... 
fit<-nlmer(resp~alpha*exp(beta*loom) + theta~~loom*cond
           +  (1|fishID),
           family="inverse.gaussian", 
           control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
           data=Firstblock)



#### apparently there is not much documentation on how to do it per groups... it seems hard to do