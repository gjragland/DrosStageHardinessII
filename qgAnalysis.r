#R
#GJR 6/14/2018
#QG analysis of Phil's factorial experiment rearing DGRP lines at 18 and 25 C and measuring both cold and heat hardiness using survival assays after acute temperature stress
#Two identical sets of analyses, one for cold, one for heat, FOCUSING HERE ON PLASTICITY (h2 and rg for cold and heat hardiness at a particular temperature are elsewhere) 


#######################################################################################
############################ COLD HARDINESS ANALYSIS ##################################
#######################################################################################




setwd('/media/raglandlab/ExtraDrive1/DrosDgrp/ConstraintIIMs')

library(nlme)
require(lme4)

countToBinom<-function(x) {

    library(dplyr)
    library(tidyr)
    a<- x %>% mutate(fdead=s_female-f_female) %>% select(-n) %>% 
        gather(status, count, c(f_female,fdead)) %>% 
        slice(rep(1:n(), .$count)) %>% select(-count) %>% 
        transform(surv=ifelse(status=="f_female",1,0), status=NULL) %>%
        arrange(line, rep, surv) %>% mutate(sex='female') %>% select(line, rep, sex, surv) 
    b<- x %>% mutate(mdead=s_male-f_male) %>% select(-n) %>% 
        gather(status, count, c(f_male,mdead)) %>% 
        slice(rep(1:n(), .$count)) %>% select(-count) %>% 
        transform(surv=ifelse(status=="f_male",1,0), status=NULL) %>%
        arrange(line, rep, surv) %>% mutate(sex='male') %>% select(line, rep, sex, surv)
    out<-rbind(a,b)
    out$line<-factor(out$line)
    out$rep<-factor(out$rep)
    return(out)
    
}


A18C<-read.csv('Adult18C.csv')
A18C<-countToBinom(A18C)
A25C<-read.csv('Adult25C.csv')
A25C<-countToBinom(A25C)

L18C<-read.csv('Larvae18C.csv')
L18C <- L18C %>% mutate(dead=n-alive) %>% select(-n) %>% 
  gather(status, count, c(alive,dead)) %>% 
  slice(rep(1:n(), .$count)) %>% select(-count) %>% 
  transform(surv=ifelse(status=="alive",1,0), status=NULL) %>%
    arrange(line, rep, surv)  %>% select(line, rep, surv)
L18C$line<-factor(L18C$line)
L18C$rep<-factor(L18C$rep)
L25C<-read.csv('Larvae25C.csv')
L25C <- L25C %>% mutate(dead=n-alive) %>% select(-n) %>% 
  gather(status, count, c(alive,dead)) %>% 
  slice(rep(1:n(), .$count)) %>% select(-count) %>% 
  transform(surv=ifelse(status=="alive",1,0), status=NULL) %>%
  arrange(line, rep, surv)  %>% select(line, rep, surv)
L25C$line<-factor(L25C$line)
L25C$rep<-factor(L25C$rep)


A18C$temp<-18
A25C$temp<-25
adat<-rbind(A18C,A25C)
adat$temp<-factor(adat$temp)

L18C$temp<-18
L25C$temp<-25
ldat<-rbind(L18C,L25C)
ldat$temp<-factor(ldat$temp)

adat$stage='A'
ldat$stage='L'
ldat$sex=NA
alldat<-rbind(adat,ldat)


#significant GxE (genetic variance for plasticity) for both adults and larvae
fma<-glmer(surv ~ factor(sex) + factor(temp) + (1| line/rep) + (1| line:temp) , family=binomial,data=adat) #full model AIC = 5607.822
fm<-glmer(surv ~ factor(sex) + factor(temp) + (1| line/rep), family=binomial,data=adat) #reduced model AIC = 5694.337
fm<-glmer(surv ~ factor(sex)*factor(temp) + (1| line/rep) +  (1| line:temp) +  (1| line:sex) + (1| line:sex:temp), family=binomial,data=adat) # 5497.199
fm<-glmer(surv ~ factor(sex) + factor(temp) + (1| line/rep) +  (1| line:temp) +  (1| line:sex) + (1| line:sex:temp), family=binomial,data=adat) # 5495.467
fm<-glmer(surv ~ factor(temp) + (1| line/rep) +  (1| line:temp) +  (1| line:sex) + (1| line:sex:temp), family=binomial,data=adat) # 5503.220
fm<-glmer(surv ~ factor(sex)  + (1| line/rep) +  (1| line:temp) +  (1| line:sex) + (1| line:sex:temp), family=binomial,data=adat) #non - convergent
fm<-glmer(surv ~ factor(sex) + factor(temp) + (1| line/rep) +  (1| line:temp) +  (1| line:sex) , family=binomial,data=adat) # 5510.907
fm<-glmer(surv ~ factor(sex) + factor(temp) + (1| line/rep) +  (1| line:temp)  , family=binomial,data=adat)
fm<-glmer(surv ~ factor(sex) + factor(temp) + (1| line/rep) +  (1| line:temp) + (1| line:sex:temp), family=binomial,data=adat) #5495.278 
fm<-glmer(surv ~ factor(sex) + factor(temp) + (1| line/rep) + (1| line:sex) + (1| line:sex:temp), family=binomial,data=adat) #5498.213
fm<-glmer(surv ~ factor(sex) + factor(temp) + (1| line/rep) +  (1| line:sex:temp), family=binomial,data=adat) #5496.221

#however, because we can't take sex into account for the estimation of rg, need to use a model without sex
fma<-glmer(surv ~ factor(temp) + (1| line/rep) + (1| line:temp) , family=binomial,data=adat) #reduced model AIC = 5671.501

fml<-glmer(surv ~ factor(temp) + (1| line/rep) + (1| line:temp) , family=binomial,data=ldat) #full model AIC = 8381.639
fm<-glmer(surv ~ factor(temp) + (1| line/rep) , family=binomial,data=ldat) #reduced model AIC = 8469.425


#significant GxDxE using the entire data set, though the two-way lineXtemp could be dropped
fm3<-glmer(surv ~ stage*temp + (1| line:stage:temp) + (1| line:stage) + (1| line:temp) + (1| line/rep) , family=binomial,data=alldat) #full model AIC = 14143.088
fm<-glmer(surv ~ stage*temp + (1| line:stage) + (1| line:temp) + (1| line/rep) , family=binomial,data=alldat) #reduced model AIC = 14203.415
fm<-glmer(surv ~ stage*temp + (1| line:stage:temp) + (1| line:temp) + (1| line/rep) , family=binomial,data=alldat) #reduced model AIC = 14159.664
fmb<-glmer(surv ~ stage*temp + (1| line:stage:temp) + (1| line:stage) + (1| line/rep) , family=binomial,data=alldat) #reduced model AIC = 14141.52 **BEST MODEL**
fm<-glmer(surv ~ stage*temp + (1| line:stage:temp) + (1| line/rep) , family=binomial,data=alldat) #reduced model AIC = 14157.664
fm<-glmer(surv ~ stage*temp + (1| line/rep) , family=binomial,data=alldat) #reduced model AIC = 14288.325
fm<-glmer(surv ~ stage + temp + (1| line:stage:temp) + (1| line:stage) + (1| line/rep) , family=binomial,data=alldat) #reduced model AIC = 14163 


lineXtempVar<-as.data.frame(VarCorr(fma))[2,4]
lineXtempVar2<-as.data.frame(VarCorr(fml))[2,4]
#covar<-as.data.frame(VarCorr(fm3))[5,4]
covar<-as.data.frame(VarCorr(fm3))[3,4]
genCor<-covar/sqrt(lineXtempVar*lineXtempVar2)



library(reshape2)



# Have to make "rep" unique
s<-alldat[alldat$temp==18 & alldat$stage=='A',]
rep<-as.numeric(s$line)*10 + as.numeric(s$rep)
s$rep<-rep
dat<-s
s<-alldat[alldat$temp==25 & alldat$stage=='A',]
rep<-as.numeric(s$line)*100 + as.numeric(s$rep)
s$rep<-rep
dat<-rbind(dat,s)
s<-alldat[alldat$temp==18 & alldat$stage=='L',]
rep<-as.numeric(s$line)*1000 + as.numeric(s$rep)
s$rep<-rep
dat<-rbind(dat,s)
s<-alldat[alldat$temp==25 & alldat$stage=='L',]
rep<-as.numeric(s$line)*10000 + as.numeric(s$rep)
s$rep<-rep
dat<-rbind(dat,s)


#permutation test of correlation b/t betaA and betaL

testCor<-function(data) {
    shuffle<-function(x) {
        reps<-unique(x$rep)
        reps<-sample(reps,length(reps))
        bootDat<-x[1,]
        bootDat[,]<-NA
        i=1
        for (line in unique(x$line)) {
            nreps<-length(unique(x$rep[x$line==line]))
            chunk<-x[x$rep %in% reps[i:(i+nreps-1)],]
            chunk$line<-line
            bootDat<-rbind( bootDat,chunk )
            i=i+nreps
        }
        bootDat<-bootDat[-1,]
        return(bootDat)
    }
    l18<-shuffle(data[data$temp==18 & data$stage=='L',])
    l25<-shuffle(data[data$temp==25 & data$stage=='L',])
    a18<-shuffle(data[data$temp==18 & data$stage=='A',])
    a25<-shuffle(data[data$temp==25 & data$stage=='A',])

    dat<-rbind(l18,l25,a18,a25)
    a<-aggregate(surv ~ line:temp:stage,dat,mean)
    b<-a %>% spread(stage,surv) %>% select(-L) %>% spread(temp,A)
    A<-data.frame(b$line,betaA=(b[,2] - b[,3]))

    b<-a %>% spread(stage,surv) %>% select(-A) %>% spread(temp,L)
    L<-data.frame(b$line,betaL=(b[,2] - b[,3]))

    dat<-merge(A,L)
    b<-a %>% spread(stage,surv) %>% select(-L) %>% spread(temp,A)
    Am<-data.frame(b$line,Amean=rowMeans(b[,2:3]))

    b<-a %>% spread(stage,surv) %>% select(-A) %>% spread(temp,L)
    Lm<-data.frame(b$line,Lmean=rowMeans(b[,2:3]))

    dat<-merge(dat,Am)
    dat<-merge(dat,Lm)
    pcor<-cor(dat$betaA,dat$betaL)
    pcor<-c(pcor,cor(dat$betaA,dat$Amean))
    pcor<-c(pcor,cor(dat$betaL,dat$Lmean))
    pcor
}

library(parallel)
nCores=10
nIter=10000 
#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(nCores)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("testCor","dat","%>%","spread","select"),envir=environment())
hvec<-parLapply(cl, 1:nIter, function(x,y) testCor(y) ,y=dat )
hmat<-matrix(unlist(hvec), ncol = 3, byrow = TRUE)
stopCluster(cl)


a<-aggregate(surv ~ line:temp:stage,alldat,mean)
b<-a %>% spread(stage,surv) %>% select(-L) %>% spread(temp,A)
A<-data.frame(b$line,betaA=(b[,2] - b[,3]))
b<-a %>% spread(stage,surv) %>% select(-A) %>% spread(temp,L)
L<-data.frame(b$line,betaL=(b[,2] - b[,3]))
dat<-merge(A,L)
b<-a %>% spread(stage,surv) %>% select(-L) %>% spread(temp,A)
Am<-data.frame(b$line,Amean=rowMeans(b[,2:3]))
b<-a %>% spread(stage,surv) %>% select(-A) %>% spread(temp,L)
Lm<-data.frame(b$line,Lmean=rowMeans(b[,2:3]))
dat<-merge(dat,Am)
dat<-merge(dat,Lm)

pcor<-cor(dat$betaA,dat$betaL)
# 0.3433245
#two-tailed p
(1-ecdf(hmat[,1])(pcor))*2
# p < 0.0338

pcor<-cor(dat$betaA,dat$Amean)
# -0.8560066
#two-tailed p
(ecdf(hmat[,2])(pcor))*2
# p < 0.0242
#however, the median of the null distribution is -0.7401147, so correlation seems to be driven mainly by the magnitude of the temperature effect


pcor<-cor(dat$betaL,dat$Lmean)
# 0.1132921
#two-tailed p
(1-ecdf(hmat[,3])(pcor))*2
# p < 0.1144

#from above, a reasonably high and significant phenotypic correlation for plasticity of adult and larval cold hardiness
# tried pretty hard, but haven't been able to find a good method for genetic correlation estimate, may just have to let that one go


bootCor<-function(data) {

    sampledLines<-sample(unique(data$line),length(unique(data$line)),replace=T)
    bootDat<-data.frame()
    for (j in 1:length(sampledLines)) {
        lineDat<-data[data$line==sampledLines[j],]
        lineDat$line <- j
        bootDat<-rbind(bootDat,lineDat)
    }
    a<-aggregate(surv ~ line:temp:stage,bootDat,mean)
    b<-a %>% spread(stage,surv) %>% select(-L) %>% spread(temp,A)
    A<-data.frame(b$line,betaA=(b[,2] - b[,3]))
    b<-a %>% spread(stage,surv) %>% select(-A) %>% spread(temp,L)
    L<-data.frame(b$line,betaL=(b[,2] - b[,3]))
    dat<-merge(A,L)
    pcor<-cor(dat$betaA,dat$betaL)
    return(pcor)
}

library(parallel)
nCores=10
nIter=1000 
#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(nCores)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("bootCor","alldat","%>%","spread","select"),envir=environment())
hvec<-unlist(parLapply(cl, 1:nIter, function(x,y) bootCor(y) ,y=alldat ))
stopCluster(cl)

quantile(hvec,c(0.025,0.975))
#     2.5%     97.5% 
#0.0906036 0.5650792 

library(multcomp)


#get estimate of betaA and betaL from model, and their difference
#model.matrix(surv ~ stage*temp + (1| line:stage:temp) + (1| line:stage) + (1| line/rep),data=alldat)
# Fixed Effects:
#  (Intercept)         stageL         temp25  stageL:temp25

#betaA
#Adult 18 - 25
# 1 0 0 0
# 1 0 1 0
cont<-c(0,0,-1,0)
#so, just the negative of the parameter estimate for 25C, which is significant
#betaL
#Larva 18 - 25
# 1 1 0 0
# 1 1 1 1
cont<-c(0,0,-1,-1)
summary(glht(fmc, linfct=t(as.matrix(cont))))
#       Estimate Std. Error z value Pr(>|z|)    
#1 == 0   1.1773     0.1566   7.516 5.66e-14 ***
#also significant

#betaA - betaL
# 0,0,-1,0
# 0,0,-1,-1
cont<-c(0,0,0,-1)
#so, just the interaction term, which is sig

#mean hardiness adult
# 1 0 0.5 0
#plasticity adult
# 0 0 1 0
cont<-c(1,0,-0.5,0)
summary(glht(fmc, linfct=t(as.matrix(cont))))

fma<-lme()

#having weird issues with lme4, this is the only way I can get it to load
library(nlme)
require(lme4)
library(parallel)



################################################################################################
########################### BOOTSTRAP LINES FOR C.I. ESTIMATE OF H^2 ########################### 
################################################################################################

####NOTE tecnically, all of these h2 estimates are on the latent (logit) scale, see https://www.biorxiv.org/content/biorxiv/early/2016/02/08/026377.full.pdf

###############################
####### Larvae  Cold #######
###############################




bootH2<-function(data) {
    ##define residual variance
    ##residual in logit model with no overdispersion parameter is fixed, see:
    ## http://stats.stackexchange.com/questions/134020/interpretation-of-variance-in-multilevel-logistic-regression
    ## and, Hox, J. (2010). Multilevel Analysis: Techniques and Applications. New York: Routledge.
    ##Especially, see: 
    ##Nakagawa, Shinichi, et Holger Schielzeth. « Repeatability for Gaussian and Non Gaussian Data: a Practical Guide for Biologists
    ##Biological Reviews 85, (novembre 1, 2010): 935-956. 
    ##doi:10.1111/j.1469-185X.2010.00141.x.

    residVar<-pi^2/3
    #set warnings to errors so that they can be caught by 'try'; continure bootstrapping until nIter error-free iterations are reached
    options( warn = 2 )
    
    sampledLines<-sample(unique(data$line),length(unique(data$line)),replace=T)
    bootDat<-data.frame()
    for (j in 1:length(sampledLines)) {
        lineDat<-data[data$line==sampledLines[j],]
        lineDat$line <- j
        bootDat<-rbind(bootDat,lineDat)
    }
    tt<-try(fm<-glmer(surv ~ temp + (1| line/rep) + (1| line:temp) , family=binomial,data=bootDat),silent=T )
    if ( is(tt,"try-error") ) {
        h2<-NA
    } else {
        lineByTempVar<-as.data.frame(VarCorr(fm))[2,4]
        h2<-lineByTempVar/(lineByTempVar+residVar + as.data.frame(VarCorr(fm))[3,4])
    }
    return(h2)  
}


library(parallel)
nCores=10
nIter=600 #really 500, but run 600 to account for some iterations not converging
#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(nCores)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("bootH2","ldat","glmer","VarCorr"),envir=environment())
hvec<-unlist(parLapply(cl, 1:nIter, function(x,y) bootH2(y) ,y=ldat ) )
stopCluster(cl)
hvec<-hvec[!is.na(hvec)]
hvec<-hvec[1:500]

quantile(hvec,c(0.025,0.975))
#      2.5%      97.5% 
#0.03685444 0.10895572  


#point estimate of larval Cold h2
residVar<-pi^2/3
fm<-glmer(surv ~ temp + (1| line/rep) + (1| line:temp) , family=binomial,data=ldat)
lineByTempVar<-as.data.frame(VarCorr(fm))[2,4]
h2<-lineByTempVar/(lineByTempVar+residVar+as.data.frame(VarCorr(fm))[3,4])
h2
# 0.06855459


###############################
####### Adult  Cold #######
###############################

library(parallel)
nCores=10
nIter=600 #really 500, but run 600 to account for some iterations not converging
#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(nCores)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("bootH2","ldat","glmer","VarCorr"),envir=environment())
hvec<-unlist(parLapply(cl, 1:nIter, function(x,y) bootH2(y) ,y=adat ) )
stopCluster(cl)
hvec<-hvec[!is.na(hvec)]
hvec<-hvec[1:500]

quantile(hvec,c(0.025,0.975))
#      2.5%      97.5% 
#0.05829152 0.20355290 


#point estimate of adult Cold h2
residVar<-pi^2/3
fm<-glmer(surv ~ temp + (1| line/rep) + (1| line:temp) , family=binomial,data=adat)
lineByTempVar<-as.data.frame(VarCorr(fm))[2,4]
h2<-lineByTempVar/(lineByTempVar+residVar+as.data.frame(VarCorr(fm))[3,4])
h2
# 0.1250851





######################################################################################
######################### PERMUTATION TESTS FOR H2 ESTIMATES #########################
######################################################################################

###############################
####### Larvae Cold #######
###############################


#rep has to be unique within each lineXtemp combination
# shuffles replicates among lines WITHIN temperatures, to preserve the temperature fixed effect
testH2<-function(data) {
    options( warn = 2 )
    #define residual variance
    residVar<-pi^2/3
    d18<-data[data$temp==18,]
    d25<-data[data$temp==25,]

    data<-d18
    reps<-unique(data$rep)
    reps<-sample(reps,length(reps))
    bootDat<-data[1,]
    bootDat[,]<-NA
    i=1
    for (line in unique(data$line)) {
        nreps<-length(unique(data$rep[data$line==line]))
        chunk<-data[data$rep %in% reps[i:(i+nreps-1)],]
        chunk$line<-line
        bootDat<-rbind( bootDat,chunk )
        i=i+nreps
    }
    bootDat18<-bootDat[-1,]

    data<-d25
    reps<-unique(data$rep)
    reps<-sample(reps,length(reps))
    bootDat<-data[1,]
    bootDat[,]<-NA
    i=1
    for (line in unique(data$line)) {
        nreps<-length(unique(data$rep[data$line==line]))
        chunk<-data[data$rep %in% reps[i:(i+nreps-1)],]
        chunk$line<-line
        bootDat<-rbind( bootDat,chunk )
        i=i+nreps
    }
    bootDat25<-bootDat[-1,]
    bootDat<-rbind(bootDat18,bootDat25)
    
    tt<-try(fm<-glmer(surv ~ temp + (1| line/rep) + (1| line:temp) , family=binomial,data=bootDat),silent=T )
    if ( is(tt,"try-error") ) {
        h2<-NA
    } else {
        lineByTempVar<-as.data.frame(VarCorr(fm))[2,4]
        h2<-lineByTempVar/(lineByTempVar+residVar+as.data.frame(VarCorr(fm))[3,4])
    }
    return(h2)  
}


# Have to make "rep" unique
dat18<-ldat[ldat$temp==18,]
rep<-as.numeric(dat18$line)*10 + as.numeric(dat18$rep)
dat18$rep<-rep
dat25<-ldat[ldat$temp==25,]
rep<-as.numeric(dat25$line)*100 + as.numeric(dat25$rep)
dat25$rep<-rep
dat<-rbind(dat18,dat25)

library(parallel)
nCores=10
nIter=1200 #really 1000, but run 1200 to account for some iterations not converging
#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(nCores)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("testH2","dat","glmer","VarCorr"),envir=environment())
hvec<-unlist(parLapply(cl, 1:nIter, function(x,y) testH2(y) ,y=dat ) )
stopCluster(cl)
hvec<-hvec[!is.na(hvec)]
hvec<-hvec[1:1000]


#point estimate of larval cold h2:
residVar<-pi^2/3
fm<-glmer(surv ~ temp + (1| line/rep) + (1| line:temp) , family=binomial,data=ldat)
lineByTempVar<-as.data.frame(VarCorr(fm))[2,4]
#h2<-lineByTempVar/(lineByTempVar+residVar+as.data.frame(VarCorr(fm))[3,4])
h2<-lineByTempVar/(lineByTempVar+residVar+as.data.frame(VarCorr(fm))[3,4])
h2
# 0.06855459
#two-tailed p
(1-ecdf(hvec)(h2))*2
# p < 0.001


###############################
####### Adult Cold #######
###############################

# Have to make "rep" unique
dat<-adat
rep<-as.numeric(dat$line)*10 + as.numeric(dat$rep)
dat$rep<-rep

library(parallel)
nCores=10
nIter=1200 #really 1000, but run 1200 to account for some iterations not converging
#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(nCores)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("testH2","dat","glmer","VarCorr"),envir=environment())
hvec<-unlist(parLapply(cl, 1:nIter, function(x,y) testH2(y) ,y=dat ) )
stopCluster(cl)
hvec<-hvec[!is.na(hvec)]
hvec<-hvec[1:1000]


#point estimate of larval cold h2:
residVar<-pi^2/3
fm<-glmer(surv ~ temp + (1| line/rep) + (1| line:temp) , family=binomial,data=adat)
lineByTempVar<-as.data.frame(VarCorr(fm))[2,4]
h2<-lineByTempVar/(lineByTempVar+residVar+as.data.frame(VarCorr(fm))[3,4])
h2
# 0.1250851
#two-tailed p
(1-ecdf(hvec)(h2))*2
# p < 0.001





#R
#GJR 6/14/2018
#QG analysis of Phil's factorial experiment rearing DGRP lines at 18 and 25 C and measuring both cold and heat hardiness using survival assays after acute temperature stress
#Two identical sets of analyses, one for cold, one for heat, FOCUSING HERE ON PLASTICITY (h2 and rg for cold and heat hardiness at a particular temperature are elsewhere) 


#######################################################################################
############################ HEAT HARDINESS ANALYSIS ##################################
#######################################################################################




setwd('/media/raglandlab/ExtraDrive1/DrosDgrp/ConstraintIIMs')

countToBinom<-function(x) {

    library(dplyr)
    library(tidyr)
    a<- x %>% mutate(fdead=s_female-f_female) %>% select(-n) %>% 
        gather(status, count, c(f_female,fdead)) %>% 
        slice(rep(1:n(), .$count)) %>% select(-count) %>% 
        transform(surv=ifelse(status=="f_female",1,0), status=NULL) %>%
        arrange(line, rep, surv) %>% mutate(sex='female') %>% select(line, rep, sex, surv) 
    b<- x %>% mutate(mdead=s_male-f_male) %>% select(-n) %>% 
        gather(status, count, c(f_male,mdead)) %>% 
        slice(rep(1:n(), .$count)) %>% select(-count) %>% 
        transform(surv=ifelse(status=="f_male",1,0), status=NULL) %>%
        arrange(line, rep, surv) %>% mutate(sex='male') %>% select(line, rep, sex, surv)
    out<-rbind(a,b)
    out$line<-factor(out$line)
    out$rep<-factor(out$rep)
    return(out)
    
}


A18H<-read.csv('Adult18H.csv')
A18H<-countToBinom(A18H)
A25H<-read.csv('Adult25H.csv')
A25H<-countToBinom(A25H)

L18H<-read.csv('Larvae18H.csv')
L18H <- L18H %>% mutate(dead=n-alive) %>% select(-n) %>% 
  gather(status, count, c(alive,dead)) %>% 
  slice(rep(1:n(), .$count)) %>% select(-count) %>% 
  transform(surv=ifelse(status=="alive",1,0), status=NULL) %>%
    arrange(line, rep, surv)  %>% select(line, rep, surv)
L18H$line<-factor(L18H$line)
L18H$rep<-factor(L18H$rep)
L25H<-read.csv('Larvae25H.csv')
L25H <- L25H %>% mutate(dead=n-alive) %>% select(-n) %>% 
  gather(status, count, c(alive,dead)) %>% 
  slice(rep(1:n(), .$count)) %>% select(-count) %>% 
  transform(surv=ifelse(status=="alive",1,0), status=NULL) %>%
  arrange(line, rep, surv)  %>% select(line, rep, surv)
L25H$line<-factor(L25H$line)
L25H$rep<-factor(L25H$rep)


A18H$temp<-18
A25H$temp<-25
adat<-rbind(A18H,A25H)
adat$temp<-factor(adat$temp)

L18H$temp<-18
L25H$temp<-25
ldat<-rbind(L18H,L25H)
ldat$temp<-factor(ldat$temp)

adat$stage='A'
ldat$stage='L'
ldat$sex=NA

alldat$trait<-'cold'
a<-alldat
alldat<-rbind(adat,ldat)
alldat$trait<-'heat'
alldatTrait<-rbind(a,alldat)







#significant GxE (genetic variance for plasticity) for both adults and larvae
fma<-glmer(surv ~ factor(sex) + factor(temp) + (1| line/rep) + (1| line:temp) , family=binomial,data=adat) #full model AIC = 7352.412
fm<-glmer(surv ~ factor(sex) + factor(temp) + (1| line/rep), family=binomial,data=adat) #reduced model AIC = 7737.931


fml<-glmer(surv ~ factor(temp) + (1| line/rep) + (1| line:temp) , family=binomial,data=ldat) #full model AIC = 8381.639
fm<-glmer(surv ~ factor(temp) + (1| line/rep) , family=binomial,data=ldat) #reduced model AIC = 9163.448


#significant GxDxE using the entire data set, though the two-way lineXtemp could be dropped
fm3<-glmer(surv ~ stage*temp + (1| line:stage:temp) + (1| line:stage) + (1| line:temp) + (1| line/rep) , family=binomial,data=alldat) #full model AIC = 17178.523
fm<-glmer(surv ~ stage*temp + (1| line:stage) + (1| line:temp) + (1| line/rep) , family=binomial,data=alldat) #reduced model AIC = 17420.462
fm<-glmer(surv ~ stage*temp + (1| line:stage:temp) + (1| line:temp) + (1| line/rep) , family=binomial,data=alldat) #reduced model AIC = 17193.358
fm<-glmer(surv ~ stage*temp + (1| line:stage:temp) + (1| line:stage) + (1| line/rep) , family=binomial,data=alldat) #reduced model AIC = 17177.499 
fm<-glmer(surv ~ stage*temp + (1| line:stage:temp) + (1| line/rep) , family=binomial,data=alldat) #reduced model AIC = 17191.358
fm<-glmer(surv ~ stage*temp + (1| line/rep) , family=binomial,data=alldat) #reduced model AIC = 18938.885
fm<-glmer(surv ~ stage + temp + (1| line:stage:temp) + (1| line:stage) + (1| line/rep) , family=binomial,data=alldat) #reduced AIC = 17176 **BEST MODEL**


lineXtempVar<-as.data.frame(VarCorr(fma))[2,4]
lineXtempVar2<-as.data.frame(VarCorr(fml))[2,4]
#covar<-as.data.frame(VarCorr(fm3))[5,4]
covar<-as.data.frame(VarCorr(fm3))[3,4]
genCor<-covar/sqrt(lineXtempVar*lineXtempVar2)


fmAll<-glmer(surv ~ stage*temp*trait + (1| line:stage:temp) + (1| line:stage) + (1| line:temp) + (1| line/rep) , family=binomial,data=alldatTrait)


a<-aggregate(surv ~ line:temp, adat, mean)
b<-aggregate(surv ~ line:temp, ldat, mean)

betaA<-a$surv[a$temp==25] - a$surv[a$temp==18]
betaL<-b$surv[b$temp==25] - b$surv[b$temp==18]
phenCor<-cor(betaA,betaL)
# 0.2444697


#from above, a reasonably high phenotypic correlation for plasticity of adult and larval cold hardiness
# tried pretty hard, but haven't been able to find a good method for genetic correlation estimate, may just have to let that one go

library(dplyr)
library(tidyr)
library(reshape2)



# Have to make "rep" unique
s<-alldat[alldat$temp==18 & alldat$stage=='A',]
rep<-as.numeric(s$line)*10 + as.numeric(s$rep)
s$rep<-rep
dat<-s
s<-alldat[alldat$temp==25 & alldat$stage=='A',]
rep<-as.numeric(s$line)*100 + as.numeric(s$rep)
s$rep<-rep
dat<-rbind(dat,s)
s<-alldat[alldat$temp==18 & alldat$stage=='L',]
rep<-as.numeric(s$line)*1000 + as.numeric(s$rep)
s$rep<-rep
dat<-rbind(dat,s)
s<-alldat[alldat$temp==25 & alldat$stage=='L',]
rep<-as.numeric(s$line)*10000 + as.numeric(s$rep)
s$rep<-rep
dat<-rbind(dat,s)


#permutation test of correlation b/t betaA and betaL

testCor<-function(data) {
    shuffle<-function(x) {
        reps<-unique(x$rep)
        reps<-sample(reps,length(reps))
        bootDat<-x[1,]
        bootDat[,]<-NA
        i=1
        for (line in unique(x$line)) {
            nreps<-length(unique(x$rep[x$line==line]))
            chunk<-x[x$rep %in% reps[i:(i+nreps-1)],]
            chunk$line<-line
            bootDat<-rbind( bootDat,chunk )
            i=i+nreps
        }
        bootDat<-bootDat[-1,]
        return(bootDat)
    }
    l18<-shuffle(data[data$temp==18 & data$stage=='L',])
    l25<-shuffle(data[data$temp==25 & data$stage=='L',])
    a18<-shuffle(data[data$temp==18 & data$stage=='A',])
    a25<-shuffle(data[data$temp==25 & data$stage=='A',])

    dat<-rbind(l18,l25,a18,a25)
    a<-aggregate(surv ~ line:temp:stage,dat,mean)
    b<-a %>% spread(stage,surv) %>% select(-L) %>% spread(temp,A)
    A<-data.frame(b$line,betaA=(b[,2] - b[,3]))

    b<-a %>% spread(stage,surv) %>% select(-A) %>% spread(temp,L)
    L<-data.frame(b$line,betaL=(b[,2] - b[,3]))

    dat<-merge(A,L)
    b<-a %>% spread(stage,surv) %>% select(-L) %>% spread(temp,A)
    Am<-data.frame(b$line,Amean=rowMeans(b[,2:3]))

    b<-a %>% spread(stage,surv) %>% select(-A) %>% spread(temp,L)
    Lm<-data.frame(b$line,Lmean=rowMeans(b[,2:3]))

    dat<-merge(dat,Am)
    dat<-merge(dat,Lm)
    pcor<-cor(dat$betaA,dat$betaL)
    pcor<-c(pcor,cor(dat$betaA,dat$Amean))
    pcor<-c(pcor,cor(dat$betaL,dat$Lmean))
    pcor
}

library(parallel)
nCores=10
nIter=10000 
#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(nCores)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("testCor","dat","%>%","spread","select"),envir=environment())
hvec<-parLapply(cl, 1:nIter, function(x,y) testCor(y) ,y=dat )
hmat<-matrix(unlist(hvec), ncol = 3, byrow = TRUE)
stopCluster(cl)


a<-aggregate(surv ~ line:temp:stage,alldat,mean)
b<-a %>% spread(stage,surv) %>% select(-L) %>% spread(temp,A)
A<-data.frame(b$line,betaA=(b[,2] - b[,3]))
b<-a %>% spread(stage,surv) %>% select(-A) %>% spread(temp,L)
L<-data.frame(b$line,betaL=(b[,2] - b[,3]))
dat<-merge(A,L)
b<-a %>% spread(stage,surv) %>% select(-L) %>% spread(temp,A)
Am<-data.frame(b$line,Amean=rowMeans(b[,2:3]))
b<-a %>% spread(stage,surv) %>% select(-A) %>% spread(temp,L)
Lm<-data.frame(b$line,Lmean=rowMeans(b[,2:3]))
dat<-merge(dat,Am)
dat<-merge(dat,Lm)

pcor<-cor(dat$betaA,dat$betaL)
# 0.2444697
#two-tailed p
(1-ecdf(hmat[,1])(pcor))*2
# p < 0.1354

pcor<-cor(dat$betaA,dat$Amean)
# -0.1070893
#two-tailed p
(ecdf(hmat[,2])(pcor))*2
# p < 0.3532


pcor<-cor(dat$betaL,dat$Lmean)
# -0.203455
#two-tailed p
(ecdf(hmat[,3])(pcor))*2
# p < 0.418



bootCor<-function(data) {

    sampledLines<-sample(unique(data$line),length(unique(data$line)),replace=T)
    bootDat<-data.frame()
    for (j in 1:length(sampledLines)) {
        lineDat<-data[data$line==sampledLines[j],]
        lineDat$line <- j
        bootDat<-rbind(bootDat,lineDat)
    }
    a<-aggregate(surv ~ line:temp:stage,bootDat,mean)
    b<-a %>% spread(stage,surv) %>% select(-L) %>% spread(temp,A)
    A<-data.frame(b$line,betaA=(b[,2] - b[,3]))
    b<-a %>% spread(stage,surv) %>% select(-A) %>% spread(temp,L)
    L<-data.frame(b$line,betaL=(b[,2] - b[,3]))
    dat<-merge(A,L)
    pcor<-cor(dat$betaA,dat$betaL)
    return(pcor)
}

library(parallel)
nCores=10
nIter=1000 
#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(nCores)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("bootCor","alldat","%>%","spread","select"),envir=environment())
hvec<-unlist(parLapply(cl, 1:nIter, function(x,y) bootCor(y) ,y=alldat ))
stopCluster(cl)

quantile(hvec,c(0.025,0.975))
#      2.5%      97.5% 
#0.01398381 0.48972900 




fma<-lme()

#having weird issues with lme4, this is the only way I can get it to load
library(nlme)
require(lme4)
library(parallel)



################################################################################################
########################### BOOTSTRAP LINES FOR C.I. ESTIMATE OF H^2 ########################### 
################################################################################################

####NOTE tecnically, all of these h2 estimates are on the latent (logit) scale, see https://www.biorxiv.org/content/biorxiv/early/2016/02/08/026377.full.pdf

###############################
####### Larvae  Heat #######
###############################




bootH2<-function(data) {
    ##define residual variance
    ##residual in logit model with no overdispersion parameter is fixed, see:
    ## http://stats.stackexchange.com/questions/134020/interpretation-of-variance-in-multilevel-logistic-regression
    ## and, Hox, J. (2010). Multilevel Analysis: Techniques and Applications. New York: Routledge.
    ##Especially, see: 
    ##Nakagawa, Shinichi, et Holger Schielzeth. « Repeatability for Gaussian and Non Gaussian Data: a Practical Guide for Biologists
    ##Biological Reviews 85, (novembre 1, 2010): 935-956. 
    ##doi:10.1111/j.1469-185X.2010.00141.x.

    residVar<-pi^2/3
    #set warnings to errors so that they can be caught by 'try'; continure bootstrapping until nIter error-free iterations are reached
    options( warn = 2 )
    
    sampledLines<-sample(unique(data$line),length(unique(data$line)),replace=T)
    bootDat<-data.frame()
    for (j in 1:length(sampledLines)) {
        lineDat<-data[data$line==sampledLines[j],]
        lineDat$line <- j
        bootDat<-rbind(bootDat,lineDat)
    }
    tt<-try(fm<-glmer(surv ~ temp + (1| line/rep) + (1| line:temp) , family=binomial,data=bootDat),silent=T )
    if ( is(tt,"try-error") ) {
        h2<-NA
    } else {
        lineByTempVar<-as.data.frame(VarCorr(fm))[2,4]
        h2<-lineByTempVar/(lineByTempVar+residVar + as.data.frame(VarCorr(fm))[3,4])
    }
    return(h2)  
}


library(parallel)
nCores=10
nIter=600 #really 500, but run 600 to account for some iterations not converging
#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(nCores)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("bootH2","ldat","glmer","VarCorr"),envir=environment())
hvec<-unlist(parLapply(cl, 1:nIter, function(x,y) bootH2(y) ,y=ldat ) )
stopCluster(cl)
hvec<-hvec[!is.na(hvec)]
hvec<-hvec[1:500]

quantile(hvec,c(0.025,0.975))
#      2.5%      97.5% 
#0.06543181 0.20790028   


#point estimate of larval Cold h2
residVar<-pi^2/3
fm<-glmer(surv ~ temp + (1| line/rep) + (1| line:temp) , family=binomial,data=ldat)
lineByTempVar<-as.data.frame(VarCorr(fm))[2,4]
h2<-lineByTempVar/(lineByTempVar+residVar+as.data.frame(VarCorr(fm))[3,4])
h2
# 0.1333803


###############################
####### Adult  Heat #######
###############################

library(parallel)
nCores=10
nIter=600 #really 500, but run 600 to account for some iterations not converging
#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(nCores)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("bootH2","ldat","glmer","VarCorr"),envir=environment())
hvec<-unlist(parLapply(cl, 1:nIter, function(x,y) bootH2(y) ,y=adat ) )
stopCluster(cl)
hvec<-hvec[!is.na(hvec)]
hvec<-hvec[1:500]

quantile(hvec,c(0.025,0.975))
#      2.5%      97.5% 
#0.08568465 0.26870317 


#point estimate of adult Cold h2
residVar<-pi^2/3
fm<-glmer(surv ~ temp + (1| line/rep) + (1| line:temp) , family=binomial,data=adat)
lineByTempVar<-as.data.frame(VarCorr(fm))[2,4]
h2<-lineByTempVar/(lineByTempVar+residVar+as.data.frame(VarCorr(fm))[3,4])
h2
#  0.1693904




######################################################################################
######################### PERMUTATION TESTS FOR H2 ESTIMATES #########################
######################################################################################

###############################
####### Larvae Heat #######
###############################


#rep has to be unique within each lineXtemp combination
# shuffles replicates among lines WITHIN temperatures, to preserve the temperature fixed effect
testH2<-function(data) {
    options( warn = 2 )
    #define residual variance
    residVar<-pi^2/3
    d18<-data[data$temp==18,]
    d25<-data[data$temp==25,]

    data<-d18
    reps<-unique(data$rep)
    reps<-sample(reps,length(reps))
    bootDat<-data[1,]
    bootDat[,]<-NA
    i=1
    for (line in unique(data$line)) {
        nreps<-length(unique(data$rep[data$line==line]))
        chunk<-data[data$rep %in% reps[i:(i+nreps-1)],]
        chunk$line<-line
        bootDat<-rbind( bootDat,chunk )
        i=i+nreps
    }
    bootDat18<-bootDat[-1,]

    data<-d25
    reps<-unique(data$rep)
    reps<-sample(reps,length(reps))
    bootDat<-data[1,]
    bootDat[,]<-NA
    i=1
    for (line in unique(data$line)) {
        nreps<-length(unique(data$rep[data$line==line]))
        chunk<-data[data$rep %in% reps[i:(i+nreps-1)],]
        chunk$line<-line
        bootDat<-rbind( bootDat,chunk )
        i=i+nreps
    }
    bootDat25<-bootDat[-1,]
    bootDat<-rbind(bootDat18,bootDat25)
    
    tt<-try(fm<-glmer(surv ~ temp + (1| line/rep) + (1| line:temp) , family=binomial,data=bootDat),silent=T )
    if ( is(tt,"try-error") ) {
        h2<-NA
    } else {
        lineByTempVar<-as.data.frame(VarCorr(fm))[2,4]
        h2<-lineByTempVar/(lineByTempVar+residVar+as.data.frame(VarCorr(fm))[3,4])
    }
    return(h2)  
}


# Have to make "rep" unique
dat18<-ldat[ldat$temp==18,]
rep<-as.numeric(dat18$line)*10 + as.numeric(dat18$rep)
dat18$rep<-rep
dat25<-ldat[ldat$temp==25,]
rep<-as.numeric(dat25$line)*100 + as.numeric(dat25$rep)
dat25$rep<-rep
dat<-rbind(dat18,dat25)

library(parallel)
nCores=10
nIter=1200 #really 1000, but run 1200 to account for some iterations not converging
#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(nCores)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("testH2","dat","glmer","VarCorr"),envir=environment())
hvec<-unlist(parLapply(cl, 1:nIter, function(x,y) testH2(y) ,y=dat ) )
stopCluster(cl)
hvec<-hvec[!is.na(hvec)]
hvec<-hvec[1:1000]


#point estimate of larval cold h2:
residVar<-pi^2/3
fm<-glmer(surv ~ temp + (1| line/rep) + (1| line:temp) , family=binomial,data=dat)
lineByTempVar<-as.data.frame(VarCorr(fm))[2,4]
#h2<-lineByTempVar/(lineByTempVar+residVar+as.data.frame(VarCorr(fm))[3,4])
h2<-lineByTempVar/(lineByTempVar+residVar+as.data.frame(VarCorr(fm))[3,4])
h2
# 0.1333803
#two-tailed p
(1-ecdf(hvec)(h2))*2
# p < 0.001


###############################
####### Adult Cold #######
###############################

# Have to make "rep" unique
dat18<-adat[adat$temp==18,]
rep<-as.numeric(dat18$line)*10 + as.numeric(dat18$rep)
dat18$rep<-rep
dat25<-adat[adat$temp==25,]
rep<-as.numeric(dat25$line)*100 + as.numeric(dat25$rep)
dat25$rep<-rep
dat<-rbind(dat18,dat25)

library(parallel)
nCores=10
nIter=1200 #really 1000, but run 1200 to account for some iterations not converging
#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(nCores)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("testH2","dat","glmer","VarCorr"),envir=environment())
hvec<-unlist(parLapply(cl, 1:nIter, function(x,y) testH2(y) ,y=dat ) )
stopCluster(cl)
hvec<-hvec[!is.na(hvec)]
hvec<-hvec[1:1000]


#point estimate of larval cold h2:
residVar<-pi^2/3
fm<-glmer(surv ~ temp + (1| line/rep) + (1| line:temp) , family=binomial,data=adat)
lineByTempVar<-as.data.frame(VarCorr(fm))[2,4]
h2<-lineByTempVar/(lineByTempVar+residVar+as.data.frame(VarCorr(fm))[3,4])
h2
# 0.1693904
#two-tailed p
(1-ecdf(hvec)(h2))*2
# p < 0.001






######################################################################################
#TESTS FOR CORRELATION BETWEEN PLASTICITY OF COLD AND HEAT STRESS WITHIN LIFE STAGES #
######################################################################################





library(dplyr)
library(tidyr)
library(reshape2)



a<-aggregate(surv ~ line:temp:stage:trait,alldatTrait,mean)

Ah<-a[a$stage=='A' & a$trait=='heat',]
b<- Ah %>% spread(temp,surv)
Ah<-data.frame( line=b$line,betaAh=(b[,4] - b[,5]) )

Ac<-a[a$stage=='A' & a$trait=='cold',]
b<- Ac %>% spread(temp,surv)
Ac<-data.frame( line=b$line,betaAc=(b[,4] - b[,5]) )

Lh<-a[a$stage=='L' & a$trait=='heat',]
b<- Lh %>% spread(temp,surv)
Lh<-data.frame( line=b$line,betaLh=(b[,4] - b[,5]) )

Lc<-a[a$stage=='L' & a$trait=='cold',]
b<- Lc %>% spread(temp,surv)
Lc<-data.frame( line=b$line,betaLc=(b[,4] - b[,5]) )

adat<-merge(Ah,Ac)
pcorA<-cor(adat$betaAh,adat$betaAc)
# 0.01624054

ldat<-merge(Lh,Lc)
pcorL<-cor(ldat$betaLh,ldat$betaLc)
# 0.170839


# Have to make "rep" unique
hdat<-alldatTrait[alldatTrait$trait=='heat',]
s<-hdat[hdat$temp==18 & hdat$stage=='A',]
rep<-as.numeric(s$line)*10 + as.numeric(s$rep)
s$rep<-rep
dat<-s
s<-hdat[hdat$temp==25 & hdat$stage=='A',]
rep<-as.numeric(s$line)*100 + as.numeric(s$rep)
s$rep<-rep
dat<-rbind(dat,s)
s<-hdat[hdat$temp==18 & hdat$stage=='L',]
rep<-as.numeric(s$line)*1000 + as.numeric(s$rep)
s$rep<-rep
dat<-rbind(dat,s)
s<-hdat[hdat$temp==25 & hdat$stage=='L',]
rep<-as.numeric(s$line)*10000 + as.numeric(s$rep)
s$rep<-rep
dat<-rbind(dat,s)
hdat<-dat

cdat<-alldatTrait[alldatTrait$trait=='cold',]
s<-cdat[cdat$temp==18 & cdat$stage=='A',]
rep<-as.numeric(s$line)*10 + as.numeric(s$rep)
s$rep<-rep
dat<-s
s<-cdat[cdat$temp==25 & cdat$stage=='A',]
rep<-as.numeric(s$line)*100 + as.numeric(s$rep)
s$rep<-rep
dat<-rbind(dat,s)
s<-cdat[cdat$temp==18 & cdat$stage=='L',]
rep<-as.numeric(s$line)*1000 + as.numeric(s$rep)
s$rep<-rep
dat<-rbind(dat,s)
s<-cdat[cdat$temp==25 & cdat$stage=='L',]
rep<-as.numeric(s$line)*10000 + as.numeric(s$rep)
s$rep<-rep
dat<-rbind(dat,s)
cdat<-dat

dat<-rbind(hdat,cdat)

testCor<-function(data) {
    shuffle<-function(x) {
        reps<-unique(x$rep)
        reps<-sample(reps,length(reps))
        bootDat<-x[1,]
        bootDat[,]<-NA
        i=1
        for (line in unique(x$line)) {
            nreps<-length(unique(x$rep[x$line==line]))
            chunk<-x[x$rep %in% reps[i:(i+nreps-1)],]
            chunk$line<-line
            bootDat<-rbind( bootDat,chunk )
            i=i+nreps
        }
        bootDat<-bootDat[-1,]
        return(bootDat)
    }
    l18<-shuffle(data[data$temp==18 & data$stage=='L',])
    l25<-shuffle(data[data$temp==25 & data$stage=='L',])
    a18<-shuffle(data[data$temp==18 & data$stage=='A',])
    a25<-shuffle(data[data$temp==25 & data$stage=='A',])

    dat<-rbind(l18,l25,a18,a25)
    a<-aggregate(surv ~ line:temp:stage:trait,dat,mean)

    Ah<-a[a$stage=='A' & a$trait=='heat',]
    b<- Ah %>% spread(temp,surv)
    Ah<-data.frame( line=b$line,betaAh=(b[,4] - b[,5]) )

    Ac<-a[a$stage=='A' & a$trait=='cold',]
    b<- Ac %>% spread(temp,surv)
    Ac<-data.frame( line=b$line,betaAc=(b[,4] - b[,5]) )

    Lh<-a[a$stage=='L' & a$trait=='heat',]
    b<- Lh %>% spread(temp,surv)
    Lh<-data.frame( line=b$line,betaLh=(b[,4] - b[,5]) )

    Lc<-a[a$stage=='L' & a$trait=='cold',]
    b<- Lc %>% spread(temp,surv)
    Lc<-data.frame( line=b$line,betaLc=(b[,4] - b[,5]) )

    adat<-merge(Ah,Ac)
    pcor<-cor(adat$betaAh,adat$betaAc)
                                        # 0.01624054

    ldat<-merge(Lh,Lc)
    pcor<-c(pcor,cor(ldat$betaLh,ldat$betaLc))
    return(pcor)
}

library(parallel)
nCores=10
nIter=10000 
#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(nCores)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("testCor","dat","%>%","spread","select"),envir=environment())
hvec<-parLapply(cl, 1:nIter, function(x,y) testCor(y) ,y=dat )
stopCluster(cl)
hmat<-matrix(unlist(hvec), ncol = 2, byrow = TRUE)



#two-tailed p
(ecdf(hmat[,1])(pcorA))*2
# p < 0.731

#two-tailed p
(1-ecdf(hmat[,2])(pcorL))*2
# p < 0.8678




bootCor<-function(data) {

    sampledLines<-sample(unique(data$line),length(unique(data$line)),replace=T)
    bootDat<-data.frame()
    for (j in 1:length(sampledLines)) {
        lineDat<-data[data$line==sampledLines[j],]
        lineDat$line <- j
        bootDat<-rbind(bootDat,lineDat)
    }
    a<-aggregate(surv ~ line:temp:stage:trait,bootDat,mean)

    Ah<-a[a$stage=='A' & a$trait=='heat',]
    b<- Ah %>% spread(temp,surv)
    Ah<-data.frame( line=b$line,betaAh=(b[,4] - b[,5]) )

    Ac<-a[a$stage=='A' & a$trait=='cold',]
    b<- Ac %>% spread(temp,surv)
    Ac<-data.frame( line=b$line,betaAc=(b[,4] - b[,5]) )

    Lh<-a[a$stage=='L' & a$trait=='heat',]
    b<- Lh %>% spread(temp,surv)
    Lh<-data.frame( line=b$line,betaLh=(b[,4] - b[,5]) )

    Lc<-a[a$stage=='L' & a$trait=='cold',]
    b<- Lc %>% spread(temp,surv)
    Lc<-data.frame( line=b$line,betaLc=(b[,4] - b[,5]) )

    adat<-merge(Ah,Ac)
    pcor<-cor(adat$betaAh,adat$betaAc)
                                        # 0.01624054

    ldat<-merge(Lh,Lc)
    pcor<-c(pcor,cor(ldat$betaLh,ldat$betaLc))
    pcor
    return(pcor)
}

library(parallel)
nCores=10
nIter=10000 
#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(nCores)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("bootCor","alldatTrait","%>%","spread","select"),envir=environment())
hvec<-unlist(parLapply(cl, 1:nIter, function(x,y) bootCor(y) ,y=alldatTrait ))
stopCluster(cl)
hmat<-matrix(unlist(hvec), ncol = 2, byrow = TRUE)

#adult
quantile(hmat[,1],c(0.025,0.975))
#      2.5%      97.5% 
#-0.2267015  0.2532972 

#larva
quantile(hmat[,2],c(0.025,0.975))
#      2.5%      97.5% 
#-0.1303040  0.4003278 
















######################## OTHER STUFF, MOSTLY ATTEMPTS TO ESTIMATE rg ##################################

















library(parallel)
totalCores<-detectCores()

#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(4)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("testH2","lar18C","glmer","VarCorr"),envir=environment())
#run more than 1000 because some error out
lar18CPermVec<-unlist(parLapply(cl, 1:1500, function(x,y) testH2(y) ,y=lar18C) )
stopCluster(cl)
lar18CPermVec<-lar18CPermVec[!is.na(lar18CPermVec)]


#point estimate of larval 18C COld Shock h2:
residVar<-pi^2/3
fm<-glmer(ld ~ (1| line/rep) , family=binomial,data=lar18C)
lineVar<-as.data.frame(VarCorr(fm))[2,4]
h2<-lineVar/(lineVar+residVar)
h2
# 0.1266239
ecdf(lar18CPermVec)(h2)
# 1
#two-tailed p
(1-ecdf(genCorPermVec)(genCor))*2
# 0
#one-tailed p
(1-ecdf(genCorPermVec)(genCor))
# 0
# Actual p-value if ecdf = 1:
1/1000
# 0.001

save.image('lar18C_permh2.R')






ecdf(larvH2PermVec)(h2)
# 1
#one or two-tailed p is essentially zero
save.image(file='lar18C.Rdat')




























##### We'll run the bootstrap replicates in parallel usin the parallel package

data<-adu25C
data2<-lar25C

#define residual variance

#set warnings to errors so that they can be caught by 'try'; continure bootstrapping until nIter error-free iterations are reached
options( warn = 2 )
library(nlme)
require(lme4)

estGenCor<-function(data,data2) {
  sampledLines<-sample(unique(data$line),length(unique(data$line)),replace=T)
  bootDat<-data.frame()
  bootDat2<-data.frame()
  for (j in 1:length(sampledLines)) {
    lineDat<-data[data$line==sampledLines[j],]
    lineDat$line <- j
    bootDat<-rbind(bootDat,lineDat)
  }
  for (j in 1:length(sampledLines)) {
    lineDat<-data2[data2$line==sampledLines[j],]
    lineDat$line <- j
    bootDat2<-rbind(bootDat2,lineDat)
  }
  
  tt<-try(fm<-glmer(surv ~ factor(temp) + (1| line/rep) + (1| line:temp) , family=binomial,data=bootDat),silent=T )
  tt2<-try(fm2<-glmer(surv ~ factor(temp) + (1| line/rep) + (1| line:temp), family=binomial,data=bootDat2),silent=T )
  alldat<-rbind(bootDat,ldat)
  
  mns<-aggregate(surv ~ line:stage,alldat,mean)
  cor(mns$surv[mns$stage=='A'],mns$surv[mns$stage=='L'])
  
  corIsPos=T
  ##### detect if correlation is negative, and if so flip and store state (cor has to be positive for this method to work)
  if ( cor(mns$surv[mns$stage=='A'],mns$surv[mns$stage=='L']) < 0 ) {
    corIsPos=F
    i1<-alldat$stage=='L' & alldat$surv==0
    i2<-alldat$stage=='L' & alldat$surv==1
    alldat$surv[i1]<-1
    alldat$surv[i2]<-0
  }
  
  mns<-aggregate(surv ~ line:stage,alldat,mean)
  cor(mns$surv[mns$stage=='A'],mns$surv[mns$stage=='L'])
  tt3<-try(fm3<-glmer(surv ~ factor(stage) + (1| line:stage) + (1| line/rep) , family=binomial,data=alldat),silent=T )
  if ( (is(tt,"try-error")) | (is(tt2,"try-error")) |  (is(tt3,"try-error")) ) {
    genCor<-NA
  } else {
    lineVar<-as.data.frame(VarCorr(fm))[2,4]
    lineVar2<-as.data.frame(VarCorr(fm2))[2,4]
    covar<-as.data.frame(VarCorr(fm3))[3,4]
    genCor<-covar/sqrt(lineVar*lineVar2)
    #make correlation negative if the values were flipped above
    if (corIsPos==F) {genCor<--genCor}
  }
  return(genCor)
}

library(parallel)
totalCores<-detectCores()

#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(4)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("estGenCor","data","data2","glmer","VarCorr"),envir=environment())
genCorVec<-unlist(parLapply(cl, 1:500, function(x,y,z) estGenCor(y,z) ,y=data, z=data2 ) )
stopCluster(cl)
quantile(genCorVec,c(0.025,0.5,0.975))
#      2.5%          50%        97.5% 
# -0.517519435  0.006706088  0.594440260

save.image('adu25Cvlar25C_CI_GenCor.R')

########## Permutation test for significance of genetic correlation ############

testGenCor<-function(data,data2) {
  #keep names from bootstrap function to avoid changing throughout
  shuffle<-function(dat) {
    lines<-unique(dat$line)
    newOrder<-sample(lines,length(lines),replace=F)
    datLines<-dat$line
    for (i in 1:length(newOrder) ) {
      dat$line[datLines==lines[i]]<-newOrder[i]
    }
    return(dat)
  }
  
  bootDat<-shuffle(data)
  bootDat2<-shuffle(data2)
  
  tt<-try(fm<-glmer(ld ~ factor(sex) + (1| line/rep) , family=binomial,data=bootDat),silent=T )
  tt2<-try(fm2<-glmer(ld ~ (1| line/rep) , family=binomial,data=bootDat2),silent=T )
  adat<-bootDat
  adat$sex<-rep('A',nrow(adat))
  names(adat)[1]<-'stage'
  ldat<-bootDat2
  #makes nesting simpler
  ldat$rep<-ldat$rep*10000
  ldat<-data.frame(stage=rep('L',nrow(ldat)),ldat)
  alldat<-rbind(adat,ldat)
  
  mns<-aggregate(ld ~ line:stage,alldat,mean)
  cor(mns$ld[mns$stage=='A'],mns$ld[mns$stage=='L'])
  
  corIsPos=T
  ##### detect if correlation is negative, and if so flip and store state (cor has to be positive for this method to work)
  if ( cor(mns$ld[mns$stage=='A'],mns$ld[mns$stage=='L']) < 0 ) {
    corIsPos=F
    i1<-alldat$stage=='L' & alldat$ld==0
    i2<-alldat$stage=='L' & alldat$ld==1
    alldat$ld[i1]<-1
    alldat$ld[i2]<-0
  }
  
  mns<-aggregate(ld ~ line:stage,alldat,mean)
  cor(mns$ld[mns$stage=='A'],mns$ld[mns$stage=='L'])
  tt3<-try(fm3<-glmer(ld ~ factor(stage) + (1| line:stage) + (1| line/rep) , family=binomial,data=alldat),silent=T )
  if ( (is(tt,"try-error")) | (is(tt2,"try-error")) |  (is(tt3,"try-error")) ) {
    genCor<-NA
  } else {
    lineVar<-as.data.frame(VarCorr(fm))[2,4]
    lineVar2<-as.data.frame(VarCorr(fm2))[2,4]
    covar<-as.data.frame(VarCorr(fm3))[3,4]
    genCor<-covar/sqrt(lineVar*lineVar2)
    #make correlation negative if the values were flipped above
    if (corIsPos==F) {genCor<--genCor}
  }
  return(genCor)
}

library(parallel)
totalCores<-detectCores()

#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(4)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("testGenCor","data","data2","glmer","VarCorr"),envir=environment())
genCorPermVec<-unlist(parLapply(cl, 1:1000, function(x,y,z) testGenCor(y,z) ,y=data, z=data2 ) )
stopCluster(cl)
quantile(genCorPermVec,c(0.025,0.5,0.975))

## point estimate of genetic correlation
adat<-data
adat$sex<-rep('A',nrow(adat))
names(adat)[1]<-'stage'
ldat<-data2
#makes nesting simpler
ldat$rep<-ldat$rep*10000
ldat<-data.frame(stage=rep('L',nrow(ldat)),ldat)
alldat<-rbind(adat,ldat)
fm<-glmer(ld ~ factor(stage) + (1| line:stage) + (1| line/rep) , family=binomial,data=alldat) 

lineVar<-as.data.frame(VarCorr(fmadu25C))[2,4]
lineVar2<-as.data.frame(VarCorr(fmlarv25C))[2,4]
covar<-as.data.frame(VarCorr(fm))[3,4]
genCor<-covar/sqrt(lineVar*lineVar2)
genCor
# 0.0437575
ecdf(genCorPermVec)(genCor)
# 0.604
# two-tailed p
(1-ecdf(genCorPermVec)(genCor))*2
# 0.792
save.image('adu25Cvslar25C_GenCor.R')


############################
########## NEW #############
############################

################################################################################################
######################## GENETIC CORRELATIONS FOR PLASTICITYS (BETAS) ##########################
################################################################################################

# Import data sets:
larC_Beta <- read.csv("C:/Users/phili/Google Drive/K-State/Research/Dmel Temp Tolerance/Chapters/Chapter 2 - 18v25 HvC/Input files/heritability/larC_Beta.csv")
aduC_Beta <- read.csv("C:/Users/phili/Google Drive/K-State/Research/Dmel Temp Tolerance/Chapters/Chapter 2 - 18v25 HvC/Input files/heritability/aduC_Beta.csv")

################################################################
####### Larval Cold Plasticity vs. Adult Cold Plasticity #######
################################################################

##### We'll run the bootstrap replicates in parallel usin the parallel package

data<-larC_Beta
data2<-aduC_Beta

#define residual variance

#set warnings to errors so that they can be caught by 'try'; continure bootstrapping until nIter error-free iterations are reached
options( warn = 2 )
library(nlme)
require(lme4)

estGenCor<-function(data,data2) {
  sampledLines<-sample(unique(data$line),length(unique(data$line)),replace=T)
  bootDat<-data.frame()
  bootDat2<-data.frame()
  for (j in 1:length(sampledLines)) {
    lineDat<-data[data$line==sampledLines[j],]
    lineDat$line <- j
    bootDat<-rbind(bootDat,lineDat)
  }
  for (j in 1:length(sampledLines)) {
    lineDat<-data2[data2$line==sampledLines[j],]
    lineDat$line <- j
    bootDat2<-rbind(bootDat2,lineDat)
  }
  
  tt<-try(fm<-lm(ld ~ (1| line/rep) ,data=bootDat),silent=T )
  tt2<-try(fm2<-lm(ld ~ (1| line/rep) ,data=bootDat2),silent=T )
  adat<-bootDat
  adat<-data.frame(stage=rep('L',nrow(adat)),adat)
  #adat$sex<-rep('A',nrow(adat))
  #adat$rep<-adat$rep*10000
  #names(adat)[1]<-'stage'
  ldat<-bootDat2
  ldat$rep<-ldat$rep*10000
  ldat<-data.frame(stage=rep('A',nrow(ldat)),ldat)
  alldat<-rbind(adat,ldat)
  alldat<-alldat[c(1,2,4,3)]
  
  mns<-aggregate(ld ~ line:stage,alldat,mean)
  cor(mns$ld[mns$stage=='L'],mns$ld[mns$stage=='A'])
  
  corIsPos=T
  ##### detect if correlation is negative, and if so flip and store state (cor has to be positive for this method to work)
  #  if ( cor(mns$ld[mns$stage=='L'],mns$ld[mns$stage=='A']) < 0 ) {
  #    corIsPos=F
  #    i1<-alldat$stage=='A' & alldat$ld==0
  #    i2<-alldat$stage=='A' & alldat$ld==1
  #    alldat$ld[i1]<-1
  #    alldat$ld[i2]<-0
  #  }
  
  mns<-aggregate(ld ~ line:stage,alldat,mean)
  cor(mns$ld[mns$stage=='L'],mns$ld[mns$stage=='A'])
  tt3<-try(fm3<-lm(ld ~ factor(stage) + (1| line:stage) + (1| line/rep) ,data=alldat),silent=T )
  if ( (is(tt,"try-error")) | (is(tt2,"try-error")) |  (is(tt3,"try-error")) ) {
    genCor<-NA
  } else {
    lineVar<-as.data.frame(VarCorr(fm))[2,4]
    lineVar2<-as.data.frame(VarCorr(fm2))[2,4]
    covar<-as.data.frame(VarCorr(fm3))[3,4]
    genCor<-covar/sqrt(lineVar*lineVar2)
    #make correlation negative if the values were flipped above
    if (corIsPos==F) {genCor<--genCor}
  }
  return(genCor)
}

library(parallel)
totalCores<-detectCores()

#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(4)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("estGenCor","data","data2","glmer","VarCorr"),envir=environment())
genCorVec<-unlist(parLapply(cl, 1:500, function(x,y,z) estGenCor(y,z) ,y=data, z=data2 ) )
stopCluster(cl)
quantile(genCorVec,c(0.025,0.5,0.975))
#     2.5%        50%      97.5% 
# -0.4067212  0.1318012  0.7146185 

save.image('larC_BetavaduC_Beta_CI_GenCor.R')

########## Permutation test for significance of genetic correlation ############

testGenCor<-function(data,data2) {
  #keep names from bootstrap function to avoid changing throughout
  shuffle<-function(dat) {
    lines<-unique(dat$line)
    newOrder<-sample(lines,length(lines),replace=F)
    datLines<-dat$line
    for (i in 1:length(newOrder) ) {
      dat$line[datLines==lines[i]]<-newOrder[i]
    }
    return(dat)
  }
  
  bootDat<-shuffle(data)
  bootDat2<-shuffle(data2)
  
  tt<-try(fm<-lm(ld ~ (1| line/rep) ,data=bootDat),silent=T )
  tt2<-try(fm2<-lm(ld ~ (1| line/rep) ,data=bootDat2),silent=T )
  adat<-bootDat
  #adat$sex<-rep('A',nrow(adat))
  #names(adat)[1]<-'stage'
  #adat$rep<-adat$rep*10000
  adat<-data.frame(stage=rep('L',nrow(adat)),adat)
  ldat<-bootDat2
  #makes nesting simpler
  ldat$rep<-ldat$rep*10000
  ldat<-data.frame(stage=rep('A',nrow(ldat)),ldat)
  alldat<-rbind(adat,ldat)
  alldat<-alldat[c(1,2,4,3)]
  
  mns<-aggregate(ld ~ line:stage,alldat,mean)
  cor(mns$ld[mns$stage=='L'],mns$ld[mns$stage=='A'])
  
  corIsPos=T
  ##### detect if correlation is negative, and if so flip and store state (cor has to be positive for this method to work)
  if ( cor(mns$ld[mns$stage=='L'],mns$ld[mns$stage=='A']) < 0 ) {
    corIsPos=F
    i1<-alldat$stage=='A' & alldat$ld==0
    i2<-alldat$stage=='A' & alldat$ld==1
    alldat$ld[i1]<-1
    alldat$ld[i2]<-0
  }
  
  mns<-aggregate(ld ~ line:stage,alldat,mean)
  cor(mns$ld[mns$stage=='L'],mns$ld[mns$stage=='A'])
  tt3<-try(fm3<-lm(ld ~ factor(stage) + (1| line:stage) + (1| line/rep) ,data=alldat),silent=T )
  if ( (is(tt,"try-error")) | (is(tt2,"try-error")) |  (is(tt3,"try-error")) ) {
    genCor<-NA
  } else {
    lineVar<-as.data.frame(VarCorr(fm))[2,4]
    lineVar2<-as.data.frame(VarCorr(fm2))[2,4]
    covar<-as.data.frame(VarCorr(fm3))[3,4]
    genCor<-covar/sqrt(lineVar*lineVar2)
    #make correlation negative if the values were flipped above
    if (corIsPos==F) {genCor<--genCor}
  }
  return(genCor)
}

library(parallel)
totalCores<-detectCores()

#create a virtual cluster, reserving computational space for a parallel loop
#stopCluster(cl)
#cl <- makeCluster( (totalCores/2) )
cl <- makeCluster(4)
##need to export variable/functions so that they are visible in the virtual machine
clusterExport(cl=cl, varlist=c("testGenCor","data","data2","glmer","VarCorr"),envir=environment())
genCorPermVec<-unlist(parLapply(cl, 1:1000, function(x,y,z) testGenCor(y,z) ,y=data, z=data2 ) )
stopCluster(cl)
quantile(genCorPermVec,c(0.025,0.5,0.975))

## point estimate of genetic correlation
adat<-data
#adat$sex<-rep('A',nrow(adat))
#names(adat)[1]<-'stage'
adat<-data.frame(stage=rep('L',nrow(adat)),adat)
ldat<-data2
#makes nesting simpler
ldat$rep<-ldat$rep*10000
ldat<-data.frame(stage=rep('A',nrow(ldat)),ldat)
alldat<-rbind(adat,ldat)
alldat<-alldat[c(1,2,4,3)]
fm<-lm(ld ~ factor(stage) + (1| line:stage) + (1| line/rep) ,data=alldat) 

lineVar<-as.data.frame(VarCorr(fmlar25C))[2,4]
lineVar2<-as.data.frame(VarCorr(fmlar25H))[2,4]
covar<-as.data.frame(VarCorr(fm))[3,4]
genCor<-covar/sqrt(lineVar*lineVar2)
genCor
# 0.1458375
ecdf(genCorPermVec)(genCor)
# 0.751
# two-tailed p
(1-ecdf(genCorPermVec)(genCor))*2
# 0.498

save.image('larC_BetavaduC_Beta_GenCor.R')
