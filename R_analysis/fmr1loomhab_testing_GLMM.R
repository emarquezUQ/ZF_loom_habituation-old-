

##### this script is to try a GLMM analysis at the behavioural responses of the 
#### of the fmr1 loom habituation experiment. 


### I based the formula on the script of my behavioural experiments from Lizard Island. 

### I tested with all the looms and with just the first 5 and I dont find that the effect of the genotye is significant
### but I am not sure I am doing it the right way... for instance, I dont get convergence and I think that might be relevant. 
### update: I used -optimizer="bobyqa",- to make it converge. 

## for all looms i gives another warning saying that i should rescale... but I tried it and it didnt work. I think that the problem is that all my variables are factors
## for looms 1-5 it worked but still not significant.
##it worked for looms 2-4 but still not significant. 

### I am still not sure the approach is the best... I think I might not be uploading the data properly


### another option could be the GEE that looks at population effects instead of individual ones. 



getwd()


setwd("C:/Users/uqemarqu/Documents/R/behaviour_Fmr1_loomhab")


## to import the data
fmr1loomhab<-read.csv("fmr1loomhab_for_GLMM.csv",header=T,sep=",") 

summary(fmr1loomhab);

fmr1loomhab$loom <- factor(fmr1loomhab$loom)
fmr1loomhab$response <- factor(fmr1loomhab$response)
fmr1loomhab$fishID <- factor(fmr1loomhab$fishID)

install.packages("lme4")



library(lme4)


fit <- glmer(response~loom*genotype 
               + (1|fishID) ,
             family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=100000)),
             data=fmr1loomhab)

fit2 <- glmer(response~loom+genotype
                 + (1|fishID) ,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
              data=fmr1loomhab)



anova(fit2, fit)

### No significant interaction

fit3 <- glmer(response~loom
              + (1|fishID) ,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
              data=fmr1loomhab)


anova(fit2, fit3)

### no effect of genotype... but this is with the 20 looms. maybe I need to test with the first 5

fit4 <- glmer(response~genotype
              + (1|fishID) ,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
              data=fmr1loomhab)


anova(fit2, fit4)


### sigificant effect of the loom

install.packages("multcomp")
library(multcomp)
tests <- glht(fit2, linfct=mcp(loom="Tukey"))

summary(tests)




summary(fit2)



############ trying now with only the first 5 looms

## to import the data
fmr1loomhab<-read.csv("fmr1loomhab_for_GLMM_short.csv",header=T,sep=",") 

summary(fmr1loomhab);

fmr1loomhab$loom <- factor(fmr1loomhab$loom)
fmr1loomhab$response <- factor(fmr1loomhab$response)
fmr1loomhab$fishID <- factor(fmr1loomhab$fishID)

install.packages("lme4")



library(lme4)


fit <- glmer(response~loom*genotype 
             + (1|fishID) ,
             family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
             data=fmr1loomhab)

fit2 <- glmer(response~loom+genotype
              + (1|fishID) ,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
              data=fmr1loomhab)



anova(fit2, fit)

### No significant interaction

fit3 <- glmer(response~loom
              + (1|fishID) ,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
              data=fmr1loomhab)


anova(fit2, fit3)

### no effect of genotype... 

fit4 <- glmer(response~genotype
              + (1|fishID) ,
              family="binomial",
              control=glmerControl(optCtrl = list(maxfun=5000)),
              data=fmr1loomhab)


anova(fit2, fit4)


### sigificant effect of the loom

install.packages("multcomp")
library(multcomp)
tests <- glht(fit2, linfct=mcp(loom="Tukey"))

summary(tests)




summary(fit2)

##### another attempt with just looms 2-4

### no effect of genotype... but got closer... 

## to import the data
fmr1loomhab<-read.csv("fmr1loomhab_for_GLMM_short2.csv",header=T,sep=",") 

summary(fmr1loomhab);

fmr1loomhab$loom <- factor(fmr1loomhab$loom)
fmr1loomhab$response <- factor(fmr1loomhab$response)
fmr1loomhab$fishID <- factor(fmr1loomhab$fishID)

install.packages("lme4")



library(lme4)


fit <- glmer(response~loom*genotype 
             + (1|fishID) ,
             family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
             data=fmr1loomhab)

fit2 <- glmer(response~loom+genotype
              + (1|fishID) ,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
              data=fmr1loomhab)



anova(fit2, fit)

### No significant interaction

fit3 <- glmer(response~loom
              + (1|fishID) ,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
              data=fmr1loomhab)


anova(fit2, fit3)

### no effect of genotype... but got closer

fit4 <- glmer(response~genotype
              + (1|fishID) ,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
              data=fmr1loomhab)


anova(fit2, fit4)


### sigificant effect of the loom

install.packages("multcomp")
library(multcomp)
tests <- glht(fit2, linfct=mcp(loom="Tukey"))

summary(tests)




summary(fit2)



##### another attempt... this time changing a bit the upload data and the formula...



##### with just looms 2-4

### no effect of genotype... and same result as before

## to import the data
fmr1loomhab<-read.csv("fmr1loomhab_for_GLMM_short2_alt.csv",header=T,sep=",") 



fmr1loomhab$loom <- factor(fmr1loomhab$loom)

fmr1loomhab$fishID <- factor(fmr1loomhab$fishID)

summary(fmr1loomhab);

install.packages("lme4")



library(lme4)


fit <- glmer(cbind(response,response2)~loom*genotype 
             + (1|fishID) ,
             family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
             data=fmr1loomhab)

fit2 <- glmer(cbind(response,response2)~loom+genotype
              + (1|fishID) ,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
              data=fmr1loomhab)



anova(fit2, fit)

### No significant interaction

fit3 <- glmer(cbind(response,response2)~loom
              + (1|fishID) ,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
              data=fmr1loomhab)


anova(fit2, fit3)

### no effect of genotype... and same result as before

fit4 <- glmer(cbind(response,response2)~genotype
              + (1|fishID) ,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
              data=fmr1loomhab)


anova(fit2, fit4)


### sigificant effect of the loom

install.packages("multcomp")
library(multcomp)
tests <- glht(fit2, linfct=mcp(loom="Tukey"))

summary(tests)




summary(fit2)



#### yet another attempt... this one with the average of the responses instead of the individual fish
### it doesnt work... I think the mean in proportions is not accepted in a binominal distribution

## to import the data
fmr1loomhab<-read.csv("fmr1loomhab_for_GLMM_means.csv",header=T,sep=",") 

summary(fmr1loomhab);

fmr1loomhab$loom <- factor(fmr1loomhab$loom)



install.packages("lme4")



library(lme4)


fit <- glmer(response~loom*genotype 
             + (1|loom) ,
             family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
             data=fmr1loomhab)

fit2 <- glmer(response~loom+genotype
              + (1|fishID) ,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
              data=fmr1loomhab)



anova(fit2, fit)

### No significant interaction

fit3 <- glmer(response~loom
              + (1|fishID) ,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
              data=fmr1loomhab)


anova(fit2, fit3)

### no effect of genotype... but got closer

fit4 <- glmer(response~genotype
              + (1|fishID) ,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=5000)),
              data=fmr1loomhab)


anova(fit2, fit4)


### sigificant effect of the loom

install.packages("multcomp")
library(multcomp)
tests <- glht(fit2, linfct=mcp(loom="Tukey"))

summary(tests)




summary(fit2)


