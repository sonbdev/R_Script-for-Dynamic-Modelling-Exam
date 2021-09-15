
################################################################## ######
#####################  MODELING POPULATION DYNAMIC         ##############
###### TOPIC:  Role of semelparity, delayed reproduction and ############
#### density dependence on population dynamics               ############

##Student: Bienvenu sonon|sonbien92@gmail.com
##Prof: Orou GAOUE
##Date: 14 September 2021
## Based on the code by Orou G. Gaoue, 2016
  
     ## Clear all
rm(list=ls(all=TRUE))
##Creating projection matrix function called "projection" using popbio package
projection <- function(sigma1,sigma2,gama,phi){
  A <- matrix(c(sigma1*(1-gama),phi,sigma1*gama, sigma2), 
              byrow = T, ncol = 2)
  library(popbio)
  e.a <- eigen.analysis(A) 
  as <- e.a$lambda1
  tr <- -log(as)
  Dynamic <- list(projection.matrix = A, asymtotic.dynamic=as,
                  transient.dynamic=tr, Elasticity=e.a$elasticities)
  return(Dynamic)
}
    ##############################
     ##      TASK 1            ##
   ###############################
#  simulate the effects of delayed reproduction on the asymptotic
#  and transient dynamics of the study species.
Task1 <- function(Gamma){
  n <- length(as.vector(Gamma))
  g <- list()
  for (i in 1:n){g[[i]]<-projection(sigma1 = 0.5,sigma2 =  0.9,
                                    gama = Gamma[i],phi = 1.5)}
  return(g)
} 
Gama <- seq(from=0, to=1, by= 0.2) # this step for variation of the gama

resG <- Task1(Gama)
resG
asymptotic_lambda0 <- sapply(resG, "[[", "asymtotic.dynamic") 
## compute the variation of lambda0 asymtotic.dynamic 
vari.lambda0 <- abs(asymptotic_lambda0[-1] - asymptotic_lambda0[-length(asymptotic_lambda0)])
vari.lambda0
#compute the meaneffect of the delayed reproduction on the asymptotic dynamic
summary(vari.lambda0)
##compute lambda transient
transient_lambda1 <-  sapply(resG, "[[", "transient.dynamic")
# compute the variation of lambda1 transient.dynamic
vari.lambda1 <- abs(transient_lambda1[-1] - transient_lambda1[-length(transient_lambda1)])
vari.lambda1
#compute the mean effect of the delayed reproduction on the transient dynamic
summary(vari.lambda1)

Table1.data <- data.frame(asymptotic_lambda0, transient_lambda1)
Table1.data
library(nlme)
t1 <- lm(asymptotic_lambda0~Gama)
summary(t1)

par(mfrow=c(1, 2))
matplot(Gama, asymptotic_lambda0,xlab = "Gama", ylab = "asymptotic_lambda0", type = "l", lty = 1, col=rainbow(10),
        main = 'delayed reproduction:asymptotic')
matplot(Gama, transient_lambda1,xlab = "Gama", ylab = "transient_lambda1", type = "l", lty = 1, col=rainbow(10),
        main = 'delayed reproduction:Transient')

#########################################
##                                    ##
##      TASK 2                        ##
########################################

#simulate the effects of iteroparity on the asymptotic and transient dynamics of the study species
Task2 <- function(Sigma2){
  n <- length(Sigma2)
  g <- list()
  for (i in 1:n){g[[i]]<-projection(sigma1 = 0.5,sigma2=Sigma2[i],
                                    gama=0.1,phi = 1.5)}
  return(g)
} 
Sgma <- seq(from=0, to=1, by=0.2)
resS <- Task2(Sgma)
resS
asymptotic.lambda0.2 <- sapply(resS, "[[", "asymtotic.dynamic") # values the asymptotic lambda
transient.lambda1.2 <-  sapply(resS, "[[", "transient.dynamic")# values the transient lambda
Table2.data <- data.frame(asymptotic.lambda0.2, transient.lambda1.2)
Table2.data
# compute the variation of lambda02 asymptotic dynamic
vari.lambda02 <- abs(asymptotic.lambda0.2[-1] - asymptotic.lambda0.2[-length(asymptotic.lambda0.2)])
vari.lambda02
#compute the mean effect of the iteroparity on the asymptotic dynamic
summary(vari.lambda02)
# compute the variation of lambda1.2 transient dynamic
vari.lambda1.2 <- abs(transient.lambda1.2[-1] - transient.lambda1.2[-length(transient.lambda1.2)])
vari.lambda1.2
#compute the mean effect of the iteroparity on the asymptotic dynamic
summary(vari.lambda1.2)
## test the significant of the sigmma2.
library(nlme)
t2 <- lm(asymptotic.lambda0.2~Sgma)
summary(t2)

par(mfrow=c(1, 2))
matplot(Sgma, asymptotic.lambda0.2,xlab = "Sgma", ylab = "asymptotic.lambda0.2", 
        type = "l", lty = 1, col=rainbow(10),main = 'iteroparity:asymptotic')
matplot(Sgma, transient.lambda1.2,xlab = "Sgma", ylab = "transient.lambda1.2", type = "l", lty = 1, col=rainbow(10),verbose = getOption("verbose"),
        main = 'iteroparity:Transient')

  
#########################################
##                                    ##
##      TASK 3                       ##
########################################

#let's simulate the effect of over-compensatory density dependence on the 
# asymptotic and transient population dynamics for semelparous 
#(sigma2 = 0.1) and iteroparous (sigma2 = 0.9) species.

# Semelparous
Task3.semelparous <- function(N){
  n <- length(N)
  g <- list()
  for (i in 1:n){g[[i]]<-projection(sigma1 = 0.5, sigma2 = 0.1,
                                    gama= 0.1, phi = 1.5*exp(-N[i]))}
  return(g)
} 
Numbers <- c(seq(from=1,to=20,by=1))
resS.N <- Task3.semelparous(Numbers)
resS.N
asymptotic3.lambda03 <- sapply(resS.N, "[[", "asymtotic.dynamic") 
transient3.lambda1.3 <-  sapply(resS.N, "[[", "transient.dynamic")

# compute the variation of lambda03 asymptotic dynamic
vari.lambda03 <- abs(asymptotic3.lambda03[-1] - asymptotic3.lambda03[-length(asymptotic3.lambda03)])
vari.lambda03
#compute the mean effect of the density dependent  on the asymptotic dynamic
summary(vari.lambda03)
# compute the variation of lambda1.3 transient dynamic
vari.lambda1.3 <- abs(transient3.lambda1.3[-1] - transient3.lambda1.3[-length(transient3.lambda1.3)])
vari.lambda1.3
#compute the mean effect of the density dependent on the asymptotic dynamic
summary(vari.lambda1.3)
## test the significant of the numbers.
library(nlme)
t3 <- lm(asymptotic3.lambda03~Numbers)
summary(t3)
Table3.data.s <- data.frame(asymptotic3.lambda03, transient3.lambda1.3)
# plot
par(mfrow=c(1,2))
matplot(Numbers, asymptotic3.lambda03,xlab = "Numbers", ylab = "asymptotic3.lambda03", 
        type = "l", lty = 1, col=rainbow(10),main = 'density dependence:asymptotic')
matplot(Numbers, transient3.lambda1.3,xlab = "Numbers", ylab = "transient3.lambda1.3", type = "l", lty = 1, col=rainbow(10),verbose = getOption("verbose"),
        main = 'density dependence:Transient')

# Iteroparous
Task3.iteroparous <- function(N){
  n <- length(N)
  g <- list()
  for (i in 1:n){g[[i]]<-projection(sigma1 = 0.5, sigma2 = 0.1,
                                    gama= 0.1, phi = 1.5*exp(-N[i]))}
  return(g)
}
Numbers <- c(seq(from=1,to=20,by=1))
resI.N <- Task3.iteroparous(Numbers)
resI.N
asymptotic3.lambda03.2 <- sapply(resI.N, "[[", "asymtotic.dynamic") 
transient3.lambda3.2 <-  sapply(resI.N, "[[", "transient.dynamic")
# compute the variation of lambda03 asymptotic dynamic
vari.lambda04 <- abs(asymptotic3.lambda03.2[-1] - asymptotic3.lambda03.2[-length(asymptotic3.lambda03.2)])
vari.lambda04
#compute the mean effect of the density dependent  on the asymptotic dynamic
summary(vari.lambda04)
# compute the variation of lambda1.3 transient dynamic
vari.lambda1.4 <- abs(transient3.lambda3.2[-1] - transient3.lambda3.2[-length(transient3.lambda3.2)])
vari.lambda1.4
#compute the mean effect of the density dependent on the asymptotic dynamic
summary(vari.lambda1.4)
## test the significance of numbers on rate growth (lambda) using nlme 
library(nlme)
t4 <- lm(asymptotic3.lambda03.2~Numbers)
summary(t4)

Table3.1.data.i <- data.frame(asymptotic3.lambda03.2, transient3.lambda3.2)

par(mfrow=c(1,2))
matplot(Numbers, asymptotic3.lambda03.2,xlab = "Numbers", ylab = "asymptotic3.lambda03.2", 
        type = "l", lty = 1, col=rainbow(10),main = 'density dependence:asymptotic')
matplot(Numbers, transient3.lambda3.2,xlab = "Numbers", ylab = "transient3.lambda3.2", type = "l", lty = 1, col=rainbow(10),verbose = getOption("verbose"),
        main = 'density dependence:Transient')



