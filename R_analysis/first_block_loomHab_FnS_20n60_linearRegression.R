##### this script is to run a linear regression analysis at the behavioural responses of the 
#### of the first block of the loom habituation FnS_20n60 experiment. 



getwd()


##setwd("C:/Users/uqemarqu/Documents/R/behaviour analysis")
##setwd("R:/EMLPHD-Q0556/R backup 20190118/behaviour analysis")

## to import the data
Firstblock<-read.csv("first_block_loomHab_FnS_20n60.csv",header=T,sep=",") 


N<-colnames(Firstblock)


## looking at the structure. 
str(Firstblock)
summary(Firstblock)

#######

### to visualize the distribution of the data
hist(Firstblock$response) ### it is left skewed... 

### to test if the distribution is normal. if it is significant is means it is not normal
shapiro.test(Firstblock$response) 

### because is left skewed i do a square transformation
Firstblock$sqresponse<-Firstblock$response^2 ### 

hist(Firstblock$sqresponse)
shapiro.test(Firstblock$sqresponse) ### it worked 

#### now it is better.


#######################

### so we test with a correlation first to see if there is a trend

### to look at the correlation

cor.test(Firstblock$loom,Firstblock$sqresponse,method=c("spearman"))
cor(Firstblock$loom,Firstblock$response,method=c("spearman")) ### it gives the same results as it should. 

plot(Firstblock$loom,Firstblock$sqresponse) ### (X,Y)
plot(Firstblock$loom~Firstblock$sqresponse)  ### when I use a ~ it changes the order to (Y,X)



##### there is a negative correlation. so i can try my multiple linear regression. but first i will do a simple one.


#########################


##### linear regression in categorical data ###############

#### is response related to speed?



### now to visualize, because we are using cathegorical data is better to use a boxplot.
boxplot(Firstblock$sqresponse~Firstblock$speed)

### linear regression
result2<-lm(Firstblock$sqresponse~Firstblock$speed)


summary(result2)

coef(result2)

#### slow b=-0.24, and it is significant, but in only explains 31% of the variance

### in categorical data the first variable is taken as a reference, that why we only get one

N2<-c("loom" , "speed" ,"isi" )

#### doing a loop

for (i in 1:length(N2)){
  print(N2[i])
  resulttest<-lm(Firstblock$sqresponse~Firstblock[,N2[i]])
  print(summary(resulttest))
  print(anova(resulttest))
}


### it seems that the significant ones are the looms presentations and the speed!! not the ISI. 


######## now trying with multiple linear regression

### we first do a linear regression with the variables that were sig or borderline

result3<-lm(Firstblock$sqresponse~Firstblock$loom+Firstblock$speed+Firstblock$isi)

step<-step(result3,direction="both")  #### be mindful that if there are NA's the stepwise can fail. 

summary(step)

### taken together ISI does appear significant...

#### to check for the assumptions. 
plot(step)



##### to check combinations

result4<-lm(Firstblock$sqresponse~Firstblock$isi+Firstblock$loom)

step2<-step(result4,direction="both")  #### be mindful that if there are NA's the stepwise can fail. 

summary(step2)

#####  to check individual variables again...
result4<-lm(Firstblock$sqresponse~Firstblock$loom)


summary(result4)


#####  to check for interactions... not sure how to interpret them yet


result4<-lm(Firstblock$sqresponse~Firstblock$speed*Firstblock$loom+Firstblock$isi*Firstblock$loom) ### this seemt o be the model that explains more of the data.Multiple R-squared:  0.8756,	Adjusted R-squared:  0.8573  

#result5<-lm(Firstblock$sqresponse~Firstblock$speed+Firstblock$isi+Firstblock$loom)
result5<-lm(Firstblock$sqresponse~Firstblock$speed*Firstblock$loom*Firstblock$isi)  ### this one is R2=0.8795. the difference between them is not significant though. 


summary(result4)
summary(result5)

anova(result4,result5)

