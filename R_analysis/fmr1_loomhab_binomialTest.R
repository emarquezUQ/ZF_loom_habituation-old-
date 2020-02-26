
#### this script is to do a binomial test on the fmr1 responses to the looms based on the wt probability of response



getwd()


##setwd("C:/Users/uqemarqu/Documents/R/behaviour analysis")
##setwd("R:/EMLPHD-Q0556/R backup 20190118/behaviour analysis")

## to import the data
Firstblock<-read.csv("loomHab_Fmr1dataset_binomialTest_forR.csv",header=T,sep=",") 


N<-colnames(Firstblock)


## looking at the structure. 
str(Firstblock)
summary(Firstblock)


### to do a binomial test
test<-binom.test(26, 42, 0.307692307692308, alternative="greater")

summary(test)

test$p.value



#### doing a loop for fmr1

for (i in 1:length(Firstblock$loom)){
  
  resulttest<-binom.test(Firstblock$Fmr1Escap[i], Firstblock$Fmr1N[i], Firstblock$WTp[i], alternative="greater")
  #resulttestAll[i]<-resulttest$p.value
  print(c(Firstblock$loom[i],resulttest$p.value))
}

### it seems that with this test only the 2nd (p=3.056486e-05), 3rd (p=0.03403846) 
### and 16th (p=0.03969664) loom are sig. the 11th is marginal (p=0.05599378).



#### doing a loop for hets

for (i in 1:length(Firstblock$loom)){
  
  resulttest<-binom.test(Firstblock$HetsEscap[i], Firstblock$HetsN[i], Firstblock$WTp[i], alternative="greater")
  #resulttestAll[i]<-resulttest$p.value
  print(c(Firstblock$loom[i],resulttest$p.value))
}

### it seems that with this test only the 2nd (p=0.0005916218) and 9th (p=0.0395049) loom are sig.


########## getting the powers now ######

#### based on:Power of a Binomial Test T.L. Scofield 9/23/2015
# https://sites.calvin.edu/scofield/courses/m343/F15/handouts/binomialTestPower.pdf
##### I actually found a function for it!!

install.packages("MESS")

library(MESS)


power_binom_test(42,0.307692307692308,26/42,sig.level = 0.00125,alternative = "greater")

## power = 0.866487


#### another options based on :https://rpubs.com/alex-lev/46271



install.packages("binom")
library(binom)

binom.power(p.alt=26/42,n=42,p=0.307692307692308,alpha = 0.00125,alternative="greater")

##### slightly different power = 0.8869561

##### finally

install.packages("pwr")
library(pwr)

pwr.p.test(ES.h(26/42,0.307692307692308),n = 42,alternative = "greater",sig.level = 0.00125)

## similar to the first one power = 0.8628459

####### now i need to make a loop. 


#### doing a loop for fmr1

for (i in 1:length(Firstblock$loom)){
  
  resulttest<-binom.test(Firstblock$Fmr1Escap[i], Firstblock$Fmr1N[i], Firstblock$WTp[i], alternative="greater")
  
  resulttestPower<-binom.power(p.alt=Firstblock$Fmr1Escap[i]/Firstblock$Fmr1N[i],n=Firstblock$Fmr1N[i],p=Firstblock$WTp[i],alpha = 0.00125,alternative="greater")
  
  print(c(Firstblock$loom[i],resulttest$p.value,resulttestPower))
}

### with a bonferroni correction it seems that only the 2nd (p=3.056486e-05, pwr=8.869561e-01) is sig. other "close" ones are, 3rd (p=0.03403846, pwr=0.17783262) 
### and 16th (p=0.03969664,pwr=0.07943119) loom are sig. the 11th is marginal (p=0.05599378,pwr=0.12894134).



#### doing a loop for hets

for (i in 1:length(Firstblock$loom)){
  
  resulttest<-binom.test(Firstblock$HetsEscap[i], Firstblock$HetsN[i], Firstblock$WTp[i], alternative="greater")
  
  resulttestPower<-binom.power(p.alt=Firstblock$HetsEscap[i]/Firstblock$HetsN[i],n=Firstblock$HetsN[i],p=Firstblock$WTp[i],alpha = 0.00125,alternative="greater")
  
  
  print(c(Firstblock$loom[i],resulttest$p.value,resulttestPower))
}

### it seems that with this test only the 2nd (p=0.0005916218, pwr=0.6689282146) and 9th (p=0.0395049, pwr=0.09904768) loom are sig.


##### to check the n we would need to find a sig difference at alpha; 0.00125 at the 3rd and 11th loom for fmr1

i=3
pwr.p.test(ES.h(Firstblock$Fmr1Escap[i]/Firstblock$Fmr1N[i],Firstblock$WTp[i]),,alternative = "greater",power = 0.8,sig.level = 0.00125)

### we would need 167 fmr1 mutant fish with such a distribution. 

i=11
pwr.p.test(ES.h(Firstblock$Fmr1Escap[i]/Firstblock$Fmr1N[i],Firstblock$WTp[i]),,alternative = "greater",power = 0.8,sig.level = 0.00125)

### we would need 213 fmr1 mutant fish with such a distribution. 
