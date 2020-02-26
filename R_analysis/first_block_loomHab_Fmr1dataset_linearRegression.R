
##### this script is to try to analyze the fmr1 loom hab behavioura data as I did with the wiltype dataset




getwd()


##setwd("C:/Users/uqemarqu/Documents/R/behaviour analysis")
##setwd("R:/EMLPHD-Q0556/R backup 20190118/behaviour analysis")

## to import the data
Firstblock<-read.csv("first_block_loomHab_Fmr1dataset.csv",header=T,sep=",") 


N<-colnames(Firstblock)


## looking at the structure. 
str(Firstblock)
summary(Firstblock)

#######

### to visualize the distribution of the data
hist(Firstblock$response) ### it is right skewed... 

### to test if the distribution is normal. if it is significant is means it is not normal
shapiro.test(Firstblock$response) 

### because is rigth skewed i do a log transformation
Firstblock$logresponse<-log(Firstblock$response) ### 

hist(Firstblock$logresponse)
shapiro.test(Firstblock$logresponse) ### it worked... although by just a bit 



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

#### is response related to genotype?



### now to visualize, because we are using cathegorical data is better to use a boxplot.
boxplot(Firstblock$logresponse~Firstblock$genotype)

### linear regression
result2<-lm(Firstblock$logresponse~Firstblock$genotype)


summary(result2)

coef(result2)

#### not significant... and even got a negative R-squared

### in categorical data the first variable is taken as a reference, that why we only get one

N2<-c("loom" , "genotype")

#### doing a loop

for (i in 1:length(N2)){
  print(N2[i])
  resulttest<-lm(Firstblock$logresponse~Firstblock[,N2[i]])
  print(summary(resulttest))
  print(anova(resulttest))
}


### it seems that the significant ones are the looms, not the genotype. 


######## now trying with multiple linear regression

### we first do a linear regression with the variables that were sig or borderline

result3<-lm(Firstblock$logresponse~Firstblock$loom+Firstblock$genotype)

step<-step(result3,direction="both")  #### be mindful that if there are NA's the stepwise can fail. 

summary(step)

### taken together genotype does appear significant... at least wt vs fmr1

#### to check for the assumptions. 
plot(step)

#### they dont look very good... 



#####  to check individual variables again...
result6<-lm(Firstblock$logresponse~Firstblock$loom)


summary(result6)


#####  to check for interactions... not sure how to interpret them yet


result4<-lm(Firstblock$logresponse~Firstblock$genotype*Firstblock$loom) ### 
#### only loom alone appears sig... not even in the interactions. 


summary(result3)
summary(result6)

anova(result3,result6)
 ### I am not getting p values. I guess i dont have enough data. 


#### it seems that the genotype is not sigificant seen this way. 

