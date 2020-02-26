##### this script is to run a linear regression analysis at the max responses of the Optic Tectum 
####(SPIM) of the of the first block of the loom habituation FnS_20n60 experiment. 



getwd()


setwd("C:/Users/uqemarqu/Documents/R/Max response SPIM analysis")


#### for OT slopehab data

## to import the data
Firstblock<-read.csv("OT_slopehab_firstblock.csv",header=T,sep=",") 


N<-colnames(Firstblock)


## looking at the structure. 
str(Firstblock)
summary(Firstblock)

#######

### to visualize the distribution of the data
hist(Firstblock$response) ### it is right skewed... 

### to test if the distribution is normal. if it is significant is means it is not normal
shapiro.test(Firstblock$response) 

### because is right skewed i do a log transformation
Firstblock$logresponse<-log(Firstblock$response+1) ###  ### the +1 is to avoid 0s 

hist(Firstblock$logresponse)
shapiro.test(Firstblock$logresponse) ### is still not normal...  the distribution looks similar... 




#######################

### so we test with a correlation first to see if there is a trend

### to look at the correlation

cor.test(Firstblock$loom,Firstblock$logresponse,method=c("spearman"))
cor(Firstblock$loom,Firstblock$response,method=c("spearman")) ### it gives the same results as it should. 

plot(Firstblock$loom,Firstblock$logresponse) ### (X,Y)
plot(Firstblock$loom~Firstblock$logresponse)  ### when I use a ~ it changes the order to (Y,X)



##### there is a negative correlation. so i can try my multiple linear regression. but first i will do a simple one.


#########################


##### linear regression in categorical data ###############

#### is response related to speed?



### now to visualize, because we are using cathegorical data is better to use a boxplot.
boxplot(Firstblock$logresponse~Firstblock$speed)

### linear regression
result2<-lm(Firstblock$logresponse~Firstblock$speed)


summary(result2)

coef(result2)

#### slow b=-0.01, and it is not significant..., but in only explains 0.2% of the variance...

### in categorical data the first variable is taken as a reference, that why we only get one

N2<-c("loom" , "speed" ,"isi" )

#### doing a loop

for (i in 1:length(N2)){
  print(N2[i])
  resulttest<-lm(Firstblock$logresponse~Firstblock[,N2[i]])
  print(summary(resulttest))
  print(anova(resulttest))
  print(coef(resulttest))
}


### it seems that the only significant one is the looms presentations... not the other two. 


######## now trying with multiple linear regression

### we first do a linear regression 

result3<-lm(Firstblock$logresponse~Firstblock$loom+Firstblock$speed+Firstblock$isi)

step<-step(result3,direction="both")  #### be mindful that if there are NA's the stepwise can fail. 

summary(step)

anova(step)



#### to check for the assumptions. 
plot(step)

#### looking at the assumptions it seems that is not a linear system so this is not an accurate test...


##### to check combinations

result4<-lm(Firstblock$logresponse~Firstblock$loom+Firstblock$isi)

step2<-step(result4,direction="both")  #### be mindful that if there are NA's the stepwise can fail. 

summary(step2)

#####  to check individual variables again...
result4<-lm(Firstblock$logresponse~Firstblock$isi)


summary(result4)



##### non linear regression#####

### i am now going to try with a nonlinear model. the idea is to make different models and compare them. 

##### its tricky cause I need to put my own formula...

### sorting by loom

y<-Firstblock$response

x<-Firstblock$loom

z<-Firstblock$speed

##### i will use the formula and initial values from the one phase decay in prism. i take the average of K and Y0
#####  Y=(Y0 - Plateau)*exp(-K*X) + Plateau


result5<-nls(y~(a-c)*exp(-b*x)+c,start=list(a=2.25,b=1.28,c=0)) ### this is my simple model with only the looms


#result6<-nls(y~(a-c)*exp(-b*x)+(a-c)*exp(-b*z)+c,start=list(a=2.25,b=1.28,c=0)) ### second model... is not working, i think because of the categorical data

#get some estimation of goodness of fit
cor(y,predict(result5))


#plot
plot(x,y)
lines(x,predict(result5),lty=1,col="red",lwd=1)

summary(result5)


#anova(result5,result6)

######## what about gam models?

library("mgcv") 

gam1<-gam(logresponse~s(loom)+0, data=Firstblock,method = "REML")

plot(gam1)

gam.check

#### the assumptions look terrible... not sure what to do now...
