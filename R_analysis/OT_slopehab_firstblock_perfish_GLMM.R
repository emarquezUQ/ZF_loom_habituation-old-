


#### this is to try a GLMM again on the FnS_20n60 max responses IN THE TECTUM dataset but taking into account the fish ID

getwd()

##setwd("C:/Users/uqemarqu/Documents/E Marquez/graph_testing_post_nature_review/loom-hab_paper-progress_nov-dec_2019/Max response SPIM analysis")
##setwd("C:/Users/Siebeck_lab.GCI-PF0EZR63/Documents/ethan lab/Max response SPIM analysis")
##setwd("R:/EMLPHD-Q0556/R backup 20190118/behaviour analysis")

## to import the data. at the moment testing for the slopehab cluster
Firstblock<-read.csv("OT_slopehab_firstblock_perfish.csv",header=T,sep=",") 


### trying the looms as ordinal. ##### I need to do this if I want to do comparisons within looms
#Firstblock$loom<-ordered(Firstblock$loom)


##### to plot it and check the structure is correct. 
#install.packages("ggplot2")
library(ggplot2)
#library(lattice)

meanresp1<-aggregate(Firstblock[, 1], list(Firstblock$loom,Firstblock$cond), mean)

#qplot(meanresp1, aes(x = Group.1,y = x),xlim=c(0,10),ylim=c(0,1)) #### this one is not working yet...

##plot(meanresp1,xlim=c(0,10),ylim=c(0,1)) #### this one is not working yet...

#xyplot(meanresp1$x ~ meanresp1$Group.1, groups=meanresp1$Group.2,xlim=c(0,10),ylim=c(0,1))

### OR better

ggplot(meanresp1, aes(meanresp1$Group.1, meanresp1$x),xlim=c(0,10),ylim=c(0,1)) +
  geom_line(aes(linetype = meanresp1$Group.2, group = meanresp1$Group.2))

#####

#install.packages("MuMIn")  ### package to calculate the Rsquare
#install.packages("lme4")  ### package to do the gl

library(lme4)
require(MuMIn)

##### checking with distribution to use

hist(Firstblock$resp)
shapiro.test(Firstblock$resp)  #### the result is sig... this means it is not normal. 

#### is skewed but there are a lot of 1 values too... 
#### i dont think doing the logaritmic corraction would be a good idea
#### i think i can use an inverse gaussian. 
#### however the r.squaredGLMM doesnt work with inverse gaussians... so I will try to log correction too

Firstblock$resp_trans<-log(Firstblock$resp)
hist(Firstblock$resp_trans)
shapiro.test(Firstblock$resp_trans) #### still not normal, but looks much better.but still not convinced on using it 


#####

fit <- glmer(Firstblock$resp~Firstblock$loom+Firstblock$cond
             +  (1|fishID),
             family="inverse.gaussian",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
             data=Firstblock)


summary(fit)

#### calculating the Rsquare. it doesnt work with inverse.gaussian

#r.squaredGLMM(fit) 


#### the interaction models fails to converge... not sure why. i manage to make it converge by adding the control parameters
fit2 <- glmer(Firstblock$resp~Firstblock$loom*Firstblock$cond
              +  (1|fishID),
              family="inverse.gaussian", 
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
              data=Firstblock)


summary(fit2)

#### interesting. it seems that as interactive model f20 is sig different from the 60ISIs and and with s20
#### however with s20 is just sig. it might not stay after multiple comparisons corrections. 

#r.squaredGLMM(fit2) 

#### comparing models. they are significantly different. however I cannot see how much they explain...
#### so i am not sure which one explains more. 
anova(fit,fit2)

#######

##### checking that the condition is significant for the model

fit3 <- glmer(Firstblock$resp~Firstblock$loom
              +  (1|fishID),
              family="inverse.gaussian", 
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
              data=Firstblock)


summary(fit3)

#### comparing models. fit vs fi 3 are not significantly different... but fit2 and fit3 yes. this means that the interaction between looms and condition is significant
anova(fit2,fit3)


####### testing for other groups to see different comparisons. 

Firstblock$cond <- relevel(Firstblock$cond,"f60") # changes baseline factor level


#####

fit <- glmer(Firstblock$resp~Firstblock$loom+Firstblock$cond
             +  (1|fishID),
             family="inverse.gaussian",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
             data=Firstblock)


summary(fit)


####### or a muliple comparison... not working yet... not sure why. it seems that is is something related to how the formula is writen
library("multcomp")
mc_group<-glht(fit, linfct=mcp(loom="Tukey"))

summary(mc_group)


############## testing fiting the formula in a different way cause the glht doenst find it. 
####### now it works
fit <- glmer(resp ~ loom*cond  #### apparantly if i do with in interaction I can calculate the differences between slopes
             +  (1|fishID),data = Firstblock,
             family="inverse.gaussian",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

summary(fit)

mc_group<-glht(fit, linfct=mcp(cond="Tukey"))

summary(mc_group, test = adjusted("bonferroni"))

###### this is comparing the slopes! so it is relevant 

####################


##### but this is mixing the results of all the looms. If I need to do it for only some
##### it seems that I found a solution in the ForMarielleGLMMfirstdataset script



install.packages("lsmeans")
library(lsmeans)

##### but this needs that I make the looms categorical 
### trying the looms as ordinal. ##### I need to do this if I want to do comparisons within looms
#Firstblock$loom_cat<-ordered(Firstblock$loom)


Firstblock2<-subset(Firstblock,Firstblock$loom_cat %in% c(1,5,10))


fit4 <- glmer(resp ~ loom_cat*cond 
              +  (1|fishID),data = Firstblock2,
              family="inverse.gaussian",
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

summary(fit4)

leastsquare = lsmeans(fit4,
                      pairwise ~ loom_cat:cond,
                      adjust = "bon")
leastsquare$contrasts

#### the problem with this is that it also tries to compare across looms so I could be over correcting... 
#### still loom 5 seems to have significant differences between f60-s20


#####



##### for a kruskal-wallis test between groups at the 5th loom
#### this is based on: http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r

Firstblock2<-subset(Firstblock,Firstblock$loom %in% c(5)) ###### selecting a loom

kruskal.test(resp ~ cond, data = Firstblock2)

pairwise.wilcox.test(Firstblock2$resp, Firstblock2$cond,
                     p.adjust.method = "bonferroni")


### loom 5 seems to have significant differences between f20-f60 and f60-s20





######## testing with a gaussian distribution to try compare variance explained


##### it doesnt seem to work... probably is not the right approach. 

fit <- glmer(Firstblock$resp~Firstblock$loom+Firstblock$cond
             +  (1|fishID),
             family="gaussian",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
             data=Firstblock)


summary(fit)

#### calculating the Rsquare. it doesnt work with inverse.gaussian

#r.squaredGLMM(fit) 




########


####### what if i do it as for behaviour with separeting the stimulus features
####### although this is not the question in this case, the question is if they differn between groups...


#####

fit <- glmer(Firstblock$resp~Firstblock$loom+Firstblock$speed+Firstblock$ISI
             +  (1|fishID),
             family="inverse.gaussian",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
             data=Firstblock)


summary(fit)



####### then maybe what would be interesting would be also to look at the differences at the 11th loom


##### for that I need to load the data of the 2nd block

## to import the data. at the moment testing for the slopehab cluster
Secondblock<-read.csv("OT_slopehab_secondblock_perfish.csv",header=T,sep=",") 


meanresp2<-aggregate(Secondblock[, 1], list(Secondblock$loom,Secondblock$cond), mean)


ggplot(meanresp2, aes(meanresp2$Group.1, meanresp2$x),xlim=c(0,10),ylim=c(0,1)) +
  geom_line(aes(linetype = meanresp2$Group.2, group = meanresp2$Group.2))


##### but this needs that I make the looms categorical 
### trying the looms as ordinal. ##### I need to do this if I want to do comparisons within looms
Secondblock$loom_cat<-ordered(Secondblock$loom)


##### for a kruskal-wallis test between groups at the 5th loom
#### this is based on: http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r

Secondblock2<-subset(Secondblock,Secondblock$loom_cat %in% c(1)) ###### selecting a loom

kruskal.test(resp ~ cond, data = Secondblock2) #### sig different

pairwise.wilcox.test(Secondblock2$resp, Secondblock2$cond,
                     p.adjust.method = "bonferroni")


### loom 11 seems to have significant differences between f20-f60 and f60-s20

