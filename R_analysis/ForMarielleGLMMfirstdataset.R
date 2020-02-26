


#### this is to try a GLMM again on the FnS_20n60 dataset but taking into account the fish ID
getwd()

##setwd("C:/Users/uqemarqu/Documents/R/behaviour analysis")
##setwd("R:/EMLPHD-Q0556/R backup 20190118/behaviour analysis")
## to import the data
Firstblock<-read.csv("Behaviour_larct_startlefrequ_forGLMM.csv",header=T,sep=",") 

### trying the stimuluss as ordinal. 
#Firstblock$stimulus<-ordered(Firstblock$stimulus)
install.packages("MuMIn")  ### package to calculate the Rsquare
install.packages("lme4") 
install.packages("Matrix") 
library(lme4)
require(MuMIn)




fit <- glmer(Firstblock$resp~Firstblock$stimulus*Firstblock$cond
             +  (1|fishID),
             family="binomial",
             control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=5e5)),
             data=Firstblock)

summary(fit)
#### calculating the Rsquare. i get marginal which epresents the variance explained by the fixed effects
#### and conditional wich is interpreted as a variance explained by the entire model,  including bothfixed and random effects
r.squaredGLMM(fit) 

#### the interaction models fails to converge... not sure why. i manage to make it converge by adding the control parameters
fit2 <- glmer(Firstblock$resp~Firstblock$stimulus+Firstblock$cond
              +  (1|fishID),
              family="binomial", 
              control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=5e5)),
              data=Firstblock)

summary(fit2)
r.squaredGLMM(fit2) 
#### comparing models. it seems that as independent variables i get a slightly higher R2m (0.1 more). but that is significant
#### I think this means that the independent (not interacting) variables model is better? 
anova(fit,fit2)


#### taking away some variables. its always better with the 3 variables. and most variance is explained by stimulus. as expected
fit3 <- glmer(Firstblock$resp~Firstblock$ISI+Firstblock$speed
              +  (1|fishID),
              family="binomial", 
              control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)),
              data=Firstblock)

summary(fit3)
r.squaredGLMM(fit3) 


##################################
###Chganging the reference relevel()

fit <- glmer(Firstblock$resp~Firstblock$stimulus*Firstblock$cond
             +  (1|fishID),
             family="binomial",
             control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)),
             data=Firstblock)

summary(fit)
#### calculating the Rsquare. i get marginal which epresents the variance explained by the fixed effects
#### and conditional wich is interpreted as a variance explained by the entire model,  including bothfixed and random effects
r.squaredGLMM(fit) 

Firstblock$cond<-relevel(Firstblock$cond,ref="Loom")


fit <- glmer(Firstblock$resp~Firstblock$stimulus*Firstblock$cond
             +  (1|fishID),
             family="binomial",
             control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)),
             data=Firstblock)

summary(fit)

################################ to do multiple comparisons with in a model

#### we need to make the looms ordinals if we want to compare between looms
Firstblock$stimulus<-ordered(Firstblock$stimulus)


install.packages("multcomp")
library(multcomp)
tests <- glht(fit, linfct=mcp(stimulus="Tukey"))
summary(tests)


############## testing fiting the formula in a different way cause the glht doenst find it. 



fit <- glmer(resp ~ stimulus*cond
             +  (1|fishID),data = Firstblock,
             family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))


#### comparing between stimulus presentations and adjusting multiple comparisons. however this doesnt do it per group
tests <- glht(fit, linfct=mcp(stimulus="Tukey"))

summary(tests, test = adjusted("bonferroni"))


### now with condition. so it shows that they are significantly different in general.
tests <- glht(fit, linfct=mcp(cond="Tukey"))

summary(tests, test = adjusted("bonferroni"))

#### but I want to be able to ask in particular looms, specially the first one. 


#### I will first try to fit a model with just the 1st loom. it doesnt work. i need 2 or more levels

Firstblock2<-subset(Firstblock,Firstblock$stimulus %in% c(1))

fit3 <- glmer(resp ~ stimulus*cond
             +  (1|fishID),data = Firstblock2,
             family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))



### trying with 3 looms (first, 5th and 10th). ### but still getting everthing mixed

Firstblock2<-subset(Firstblock,Firstblock$stimulus %in% c(1,5,10))

fit3 <- glmer(resp ~ stimulus*cond
              +  (1|fishID),data = Firstblock2,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

summary(fit3)

tests <- glht(fit3, linfct=mcp(stimulus="Tukey"))

summary(tests, test = adjusted("bonferroni"))



#### aparantly there is many ways to do this. 

## one is using the lsmeans package. see:
## https://www.researchgate.net/post/How_to_test_Tukey_for_interaction_between_categorical_and_continuos_variables_in_R


## another one is creating a contrast matrix. see: 
## https://thebiobucket.blogspot.com/2011/06/glmm-with-custom-multiple-comparisons.html#more

## finally, there are a couple of ways using some R functions see page 10 of (although they also need a contrast matrix):
## https://cran.r-project.org/web/packages/multcomp/vignettes/multcomp-examples.pdf


#### I will first try the last one

### making the contrast matrix
Tukey<-contrMat(table(Firstblock2$stimulus),"Tukey")
K1<-cbind(Tukey,matrix(0,nrow=nrow(Tukey),ncol=ncol(Tukey)))
rownames(K1)<-paste(levels(Firstblock2$cond)[1],rownames(K1),sep=":")
K2<-cbind(matrix(0,nrow=nrow(Tukey),ncol=ncol(Tukey)),Tukey)
rownames(K2)<-paste(levels(Firstblock2$cond)[2],rownames(K2),sep=":")
K<-rbind(K1,K2)
colnames(K)<-c(colnames(Tukey),colnames(Tukey))


## making an interaction variable
Firstblock$SnC <- with(Firstblock, interaction(stimulus, cond))

### then fiting a model
fit3 <- glmer(resp ~ SnC -1 
              +  (1|fishID),data = Firstblock,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

summary(fit3)

#### the glht didnt work with Firstblock2 because the matrix is not the same size. by some reason K is made with all the stimulus and not only the subset that I selected...
#### so I tried again but with the whole Firstblock to see if having the same size make it work.
#### but it still didnt work... 
tests <- glht(fit3, linfct=K)

summary(tests, test = adjusted("bonferroni"))


#### testing the lsmeans package. It seems that it works!!!

install.packages("lsmeans")
library(lsmeans)



fit3 <- glmer(resp ~ stimulus*cond 
              +  (1|fishID),data = Firstblock2,
              family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

summary(fit3)

leastsquare = lsmeans(fit3,
                      pairwise ~ stimulus:cond,
                      adjust = "bon")
leastsquare$contrasts


