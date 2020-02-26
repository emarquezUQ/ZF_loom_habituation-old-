#### testing the multiple linear regression but taking into account 
#### the binomial distribution

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



##### now if I change the loom as factors the correlations dont work of course... so I will go on with the model

#########################

#### i made a new dataframe to take the escape responses
escapes<-read.csv("escapes.csv",header=T,sep=",") 

y<-cbind(escapes$escapes,escapes$noescapes)

##### glm model for binomial data ###############

fit1<-glm(y~Firstblock$loom+Firstblock$speed+Firstblock$isi, family = "binomial",data = Firstblock)

summary(fit1)


fit2<-glm(y~Firstblock$loom+Firstblock$speed, family = "binomial",data = Firstblock)

summary(fit2)

anova(fit1,fit2)

#### i dont get sig values... apparantly this can happen when you saturate the models. 



#### now I will try a multiple comparisons for the looms

library(multcomp)

tests <- glht(fit1, linfct=mcp(loom="Tukey"))

summary(tests)


################

##### glm model for binomial data ###############

fit1<-glm(y~loom+speed+isi, family = "binomial",data = Firstblock)

summary(fit1)


fit2<-glm(y~loom+speed, family = "binomial",data = Firstblock)

summary(fit2)

anova(fit1,fit2)

#### i dont get sig values... apparantly this can happen when you saturate the models. 



#### now I will try a multiple comparisons for the looms. to do this I do need to make the looms a factor

### i think I also need to change the loom variable as a factor... not sure

Firstblock$loom<-as.factor(Firstblock$loom)

fit1<-glm(y~loom+speed+isi, family = "binomial",data = Firstblock)

summary(fit1)

library(multcomp)

tests <- glht(fit1, linfct=mcp(loom="Tukey"))

summary(tests)

#### but this gets me the comparisns between looms but not between groups per loom. 

#### not sure what to do now... 