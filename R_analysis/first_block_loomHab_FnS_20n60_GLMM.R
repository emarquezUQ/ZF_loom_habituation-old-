
#### this is to try a GLMM again on the FnS_20n60 dataset but taking into account the fish ID

getwd()


##setwd("C:/Users/uqemarqu/Documents/R/behaviour analysis")
##setwd("R:/EMLPHD-Q0556/R backup 20190118/behaviour analysis")

## to import the data
Firstblock<-read.csv("FnS_20n60_forGLMM.csv",header=T,sep=",") 


### trying the looms as ordinal. 
#Firstblock$loom<-ordered(Firstblock$loom)

install.packages("MuMIn")  ### package to calculate the Rsquare

library(lme4)
require(MuMIn)

fit <- glmer(Firstblock$resp~Firstblock$loom+Firstblock$speed+Firstblock$ISI
              +  (1|fishID),
             family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
             data=Firstblock)


summary(fit)

#### calculating the Rsquare. i get marginal which epresents the variance explained by the fixed effects
#### and conditional wich is interpreted as a variance explained by the entire model,  including bothfixed and random effects
r.squaredGLMM(fit) 


#### the interaction models fails to converge... not sure why. i manage to make it converge by adding the control parameters
fit2 <- glmer(Firstblock$resp~Firstblock$loom*Firstblock$speed*Firstblock$ISI
             +  (1|fishID),
             family="binomial", 
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
             data=Firstblock)


summary(fit2)
r.squaredGLMM(fit2) 

#### comparing models. it seems that as independent variables i get a slightly higher R2m (0.1 more). but that is significant
#### I think this means that the independent (not interacting) variables model is better? 
anova(fit,fit2)




#### taking away some variables. its always better with the 3 variables. and most variance is explained by loom. as expected
fit3 <- glmer(Firstblock$resp~Firstblock$ISI+Firstblock$speed
              +  (1|fishID),
              family="binomial", 
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
              data=Firstblock)


summary(fit3)
r.squaredGLMM(fit3) 

