
##########

#Simulation to produce Tables 1 and 2 of manuscript
#Author: Carrie Caswell
#This code runs a model with one covariate (subject to measurement error)
#NOTE: Due to parallel packages used, the code may not run properly on a Mac computer.
#This code is intended for PC use.


#The output of this code is one 5 x 2 section of Table 1 (or 2).
#User: change the following inputs to reproduce different sections of Tables 1 and 2
b.11<-log(2) #true log hazard ratio for event of interest (beta_0)
sigma2<-0.2 #measurement error variance (sigma^2)
#set cov.dist to "normal" for Table 1
#set cov.dist to "uniform" for Table 2
cov.dist<-"normal" 

#NSR is sigma^2/k, and k = 4. Set sigma2 to 0.2 for NSR = 0.05, set to 1 for NSR = 0.25,
#and set to 4 for NSR = 1.


##############################################################################
if(!require(doParallel)){install.packages("doParallel")}
if(!require(doRNG)){install.packages("doRNG")}
if(!require(foreach)){install.packages("foreach")}
if(!require(parallel)){install.packages("parallel")}
if(!require(survival)){install.packages("survival")}
if(!require(cmprsk)){install.packages("cmprsk")}
if(!require(rootSolve)){install.packages("rootSolve")}
library(doParallel)
library(doRNG)
library(foreach)
library(parallel)


#The following parameters need not be changed to reproduce Tables 1 and 2
n<-200 #number of subjects
reps<-1000 #number of simulations
maxk<-4 #number of replicates per subject (k)
b.21<-b.11 #true log hazard ratio for competing event (equal to beta_0)


if(cov.dist=="normal"){
  p<-0.3 #parameter fixed to generate data according to subdistribution (see Fine & Gray 1999)
a<-0.5 #ensures approx. 30% censoring
b<-2
}

if(cov.dist=="uniform"){
  p<-0.1
  zlimita<-0
  zlimitb<-sqrt(12)
  a<-0.2
  b<-1
}
sseed<-4567
set.seed(sseed)

print(paste("covariate distribution: ",cov.dist))
print(paste("reps: ",reps))
print(paste("n: ",n))
print(paste("censor limit a: ",a))
print(paste("censor limit b: ",b))
print(paste("seed: ",sseed))
print(paste("true beta of interest: ",b.11))
print(paste("true other beta: ",b.21))
print(paste("ME var: ",sigma2))
print(paste("p: ",p,sep=""))

#function to be parallelized across as many cores as machine has
fixthatME<-function(k)
    {
  
  library(survival)
  library(cmprsk)
  library(rootSolve)
  
  
 
   
     dat<-data.frame(subjectnum=1:n)
  
      if(cov.dist=="normal"){
        dat$z1<-rnorm(n,0,1)
      }
     
      if(cov.dist=="uniform"){
        dat$z1<-runif(n,zlimita,zlimitb)
      }
      
      
      #Data generation process:
      dat$zB1<-(dat$z1*b.11)
      dat$zB2<-(dat$z1*b.21)
      
      #Calculate P(eps=1|Z) based on equation given in Fine & Gray 1999, for each value of Z
      #P(eps=1|Z)=lim(t->Inf)P(T<t,eps=1|Z) = 1-(1-p)^exp(zB1)
      dat$p.eps1z<-1-(1-p)^exp(dat$zB1)
      
      #Generate P(T<t|eps=1,Z) from Unif(0,1)
      dat$u.1<-runif(n,0,dat$p.eps1z)
      
      #Now P(T<t,eps=1|Z)=P(T<t|eps=1,Z)*P(eps=1|Z)
      #F.1<-u.1*p.eps1z
      #Now use equation given above to solve for t
      dat$t.1<-(-log(1-((1-((1-dat$u.1)^exp(-dat$zB1)))/p)))
     
      
      #Repeat the process for t.2, event time of type 2
      dat$u.2<-runif(n,0,1)
      dat$p.eps2z<-(1-dat$p.eps1z)
      
      
      #P(T<t,eps=2|Z) = P(T<t|eps=2,Z)*P(eps=2|Z) 
      # u.2 = (1-exp(-lamba.2*t.2))
      #Solve for t.2 using this equation
      dat$t.2<-log(1-dat$u.2)/(-exp(dat$zB2))
      
      #Generate censoring times
      dat$cens<-runif(n,a,b)
      
      
      #Now assess whether each subject experienced event 1, event 2, or was censored
   
      dat$delta<-vector(length=n)
      dat$times<-vector(length=n)
      dat$bern<-vector(length=n)
      for(i in 1:n){
        dat[i,"bern"]<-rbinom(1,size=1,prob=dat[i,"p.eps1z"])
        if(dat[i,"bern"]==1){dat[i,"times"]<-dat[i,"t.1"]
        dat[i,"delta"]<-1}
        if(dat[i,"bern"]==0){dat[i,"times"]<-dat[i,"t.2"]
        dat[i,"delta"]<-2}
        
        if(dat[i,"cens"]<dat[i,"times"]){dat[i,"times"]<-dat[i,"cens"]
        dat[i,"delta"]<-0}
        
      }
      
   
    
      #Now implement estimating equation instead of using crr package

      
      dat$delta.1<-ifelse(dat$delta==1,1,0)
      
      tots<-as.vector(table(dat$delta))
      
      #output percent of subjects dead(competing event), censored, and had event of interest
      percents<-tots/n
      ratio<-percents[2]/percents[3]
    
      
      #Now add measurement error replicates and implement RRC
  
      
      ###############################################################################################################
      #For JASA submission:
      #Make numreps equal for all subjects
      numreps<-rep(maxk,n)
      
      #generate measurement error
      eps.l<-lapply(1:length(numreps), function(i) rnorm(numreps[i],mean=0,sd=sqrt(sigma2)))
    
      #generate true covariates for each subject, replicate according to number of reps
      z1.l<-lapply(1:length(numreps), function(i) rep(dat$z1[i],numreps[i]))
 
      #zstar is covariate replicates with measurement error 
      zstar.l<-lapply(1:length(z1.l), function(i) z1.l[[i]]+eps.l[[i]])
      
      
      #dat$zstarbar is the average observed for each subject (W.i.bar in manuscript)
      dat$zstarbar<-lapply(zstar.l,mean)
  
      #calculate Delta (quantity to estimate measurement error - see appendix)
      dev<-unlist(zstar.l)-rep(unlist(dat$zstarbar),numreps)

      Deltahat<-sum(dev**2)/(sum(numreps)-n)
      
      #Check Deltahat
      #devcheck<-list()
      #for(i in 1:n){
      #  devcheck[[i]]<-list()
      #  for(j in 1:numreps[i]){
      #    devcheck[[i]][j]<-(zstar.l[[i]][j]-dat$zstarbar[[i]])**2
      #  }
      #}
      #sum(unlist(devcheck))/(sum(numreps)-n)
      
      
      #Now identify the at-risk sets for each type 1 failure timepoint
      #vector of type 1 failure times and corresponding true covariates
      f1time<-dat$times[which(dat$delta.1==1)]
      zi.vec<-dat$z1[which(dat$delta==1)]
      
      #zstarbar.vec is the average observed for each subject who had a type 1 failure
      zstarbar.vec<-unlist(dat$zstarbar[which(dat$delta==1)])
    
      
      #censoring distribution for weights
      censind<-ifelse(dat$delta==0,1,0) #This counts only true censors as censors
      G<-survfit(Surv(dat$times,censind)~1,type="kaplan-meier")
      
      #identify risk set at each f1time
      atrisk<-function(i)
      {
        dat$Y.t.i<-ifelse((dat$times>=f1time[i] & dat$delta %in% c(0,1)) | dat$delta==2,1,0)
        return(dat$subjectnum[which(dat$Y.t.i==1)])
      }
      
      
      #identify which subjects were at risk at each type 1 failure time
      subj.list<-lapply(1:length(f1time),atrisk)
      #corresponding covariates of the at-risk subjects
      zj.list<-lapply(1:length(f1time),function(i) (dat$z1[which(dat$subjectnum %in% subj.list[[i]])]))
      zstarbar.list<-lapply(1:length(f1time),function(i) (unlist(dat$zstarbar[which(dat$subjectnum %in% subj.list[[i]])])))
      
      #calculate censoring weight w.i for score equation
      censweight<-function(i)
      {
          w.j<-sapply(subj.list[[i]],function(j)
          ifelse((dat[j,"times"]>=f1time[i]),1,
                 ifelse((dat[j,"times"] < f1time[i] & dat[j,"delta"]!=0), G$surv[which(round(G$time,8)==round(f1time[i],8))]/
                          G$surv[which(round(G$time,8)==round(dat[j,"times"],8))[1]],0)))
        return(w.j)
      }
      
      wj.list<-lapply(1:length(f1time),censweight)
      

      #Now we have the necessary covariates. Assemble into an estimating function and optimize.
 
    
    
    #Find mu hat for each t
    mu.hat.t<-sapply(1:length(subj.list), function(i) (sum(zstarbar.list[[i]])/length(zstarbar.list[[i]])))
    
    #number of replicates for each subject - this quantity is useful when k.i is not the same for each subject
    k.i<-lapply(1:length(subj.list), function(i) (numreps[subj.list[[i]]]))
   
    #Find Sigma for each t
    Sigma.hat.t<-sapply(1:length(subj.list), function(i) ((sum((zstarbar.list[[i]]-mu.hat.t[i])**2)/(length(zstarbar.list[[i]])-1))-
                                              (sum(Deltahat/k.i[[i]])/length(zstarbar.list[[i]]))))

    
    #Find z.hat for each subject at risk at each t.i
    z.hat.list<-lapply(1:length(k.i), function(i) ((((Deltahat/k.i[[i]])*mu.hat.t[i])/(Sigma.hat.t[i]+Deltahat/k.i[[i]]))+
                                              (Sigma.hat.t[i]*zstarbar.list[[i]]/(Sigma.hat.t[i]+Deltahat/k.i[[i]]))))
    #Now take z.hat for each subject who had a type 1 event
    k.i.vec<-numreps[which(dat$delta==1)]
    z.hat.vec<-sapply(1:length(k.i.vec), function(i) ((Deltahat/k.i.vec[i])*mu.hat.t[i]/(Sigma.hat.t[i]+Deltahat/k.i.vec[i])+
                                                        (Sigma.hat.t[i]*zstarbar.vec[i]/(Sigma.hat.t[i]+Deltahat/k.i.vec[i]))))
    
    
    
  
    

    
    
    #Now calculate betahat based on RRC estimates using u3b function (u3b is the score equation)
   
    
    u3b.ME<-function(beta,zj=z.hat.list,zi=z.hat.vec)
    {
      wj<-wj.list
      sum(zi-unlist(lapply(1:length(zj),
                           function(i) sum(zj[[i]]*wj[[i]]*exp(zj[[i]]*beta))/sum(wj[[i]]*exp(zj[[i]]*beta)))))
    }
    

    
    #This should be closest to the truth
    new.beta<-uniroot(u3b.ME,interval=c(-100,100))$root
    

   ############################################################################################
   #Asymptotic variance
   ############################################################################################

 
    #quantities in this section should be consistent with notation used in appendix of manuscript
    
    #Estimate S.hat values
    #Calculate S.hat for each failure timepoint
    
    S0.hat<-sapply(1:length(f1time), function(i) sum(wj.list[[i]]*exp(z.hat.list[[i]]*new.beta))/n)
    S1.hat<-sapply(1:length(f1time), function(i) sum(wj.list[[i]]*z.hat.list[[i]]*exp(z.hat.list[[i]]*new.beta))/n)
    S2.hat<-sapply(1:length(f1time), function(i) sum(wj.list[[i]]*(z.hat.list[[i]]**2)*exp(z.hat.list[[i]]*new.beta))/n)
    
    A.mat<-sum((S2.hat/S0.hat)-((S1.hat/S0.hat)**2))/n
    
    
    ################################
    #Sum of iid variables part

    #Assemble D.ik
    #Split the population into those with the same number of reps together
    Y.i.t<-lapply(1:length(dat$subjectnum), function(i) ifelse((dat$times[i]>=f1time & dat$delta[i] %in% c(0,1)) | dat$delta[i]==2,1,0))

    E.Yt<-sapply(1:length(f1time),function(i) length(wj.list[[i]])/n)
   
    
    
    k.to.l<-sort(unique(numreps))
    a.m<-sapply(1:length(k.to.l), function(m) length(numreps[which(numreps==k.to.l[m])])/n)
    Gamma1.t<-lapply(1:length(k.to.l),function(m) (Sigma.hat.t/(Sigma.hat.t+(Deltahat/k.to.l[[m]])))/E.Yt)

    Gamma2.t<-lapply(1:length(k.to.l),function(m) 1/(Sigma.hat.t+(Deltahat/k.to.l[[m]])))

   
    Delta.i<-sapply(1:n, function(i) (sum((zstar.l[[i]]-unlist(dat$zstarbar)[[i]])**2)-((numreps[i]-1)*Deltahat))/sum(a.m*(k.to.l-1)))
    

    D.ik.t<-lapply(1:length(dat$subjectnum), function(i) Gamma1.t[[which(k.to.l==numreps[i])]]*
      ((Deltahat/(k.to.l[which(k.to.l==numreps[i])]*Sigma.hat.t))*((Y.i.t[[i]]*(unlist(dat$zstarbar)[i]-mu.hat.t)**2)-
      (Delta.i[i]*E.Yt*sum(a.m/k.to.l))-(Deltahat*Y.i.t[[i]]/numreps[i])-(Y.i.t[[i]]*Sigma.hat.t))
      -(Delta.i[i]*E.Yt/numreps[i]))*Gamma2.t[[which(k.to.l==numreps[i])]])

    w.tilde<-lapply(1:length(f1time),function(t) sapply(subj.list[[t]],function(j)
      ifelse((dat[j,"times"]>=f1time[t]),1,
             ifelse((dat[j,"times"] < f1time[t] & dat[j,"delta"]!=0), G$surv[which(round(G$time,8)==round(f1time[t],8))]/
                      G$surv[which(round(G$time,8)==round(dat[j,"times"],8))[1]],0))))
    
  
    #dNjc.t should just be an indicator for whether the subject was censored at time t
    dNjc.t<-lapply(1:length(dat$subjectnum), function(i) sapply(1:length(dat$times), 
                    function(t) ifelse(dat$times[i]==dat$times[t] & dat$delta[i]==1,1,0)))
  

    
    dat$lambda.c.t<-sapply(1:length(dat$times), function(t) length(dat$subjectnum[which(dat$times==dat$times[t] & dat$delta==0)])/
                             length(dat$subjectnum[which(dat$times>=dat$times[t])]))
    
     
    dAjc.t<-lapply(1:length(dat$subjectnum),function(i) {subT<-dat$times[which(dat$subjectnum==dat$subjectnum[i])]
    sapply(1:length(dat$times), function(t) 
      ifelse(subT>=dat$times[t],
             dat$lambda.c.t[which(dat$times==dat$times[t])],0))})
    
    
    
    
    dMjc.t<-lapply(1:length(dat$subjectnum),function(i) dNjc.t[[i]]-dAjc.t[[i]])
    
    
    dMjc.t.rev<-lapply(1:length(dat$times), function(t) sapply(1:length(dat$subjectnum), 
                                                           function(i) dMjc.t[[i]][t]))
    
    EdMc.u<-sapply(1:length(dMjc.t.rev), function(t) sum(dMjc.t.rev[[t]])/length(dat$subjectnum))
    
    
    Y.i.t.rev<-lapply(1:length(f1time), function(t) sapply(1:length(dat$subjectnum),
                                                           function(i) Y.i.t[[i]][t]))
    
    Y.t.i.all<-lapply(1:length(dat$times), function(t) sapply(1:length(dat$subjectnum), function(i)
      ifelse((dat$times[i]<=t & dat$delta[i] %in% c(0,1)) | (dat$delta[i]==2),1,0)))
    
    pi.hat.u<-sapply(1:length(dat$times), function(t) sum(Y.t.i.all[[t]])/length(dat$subjectnum))
    
    
    

    
    #########################################
    #Part 2 of variance
    
    #################################################
    #Parts D and E
    numreps.at.risk<-lapply(1:length(f1time), function(t) numreps[which(dat$subjectnum %in% subj.list[[t]])])
    failtimes.at.risk<-lapply(1:length(f1time), function(t) dat$times[which(dat$subjectnum %in% subj.list[[t]])])
    failsubs<-dat$subjectnum[which(dat$delta==1)]
    
    dR.tilde.t<-1/n
    dat$numreps<-numreps
    E.WwdN.t<-lapply(1:length(k.to.l), function(m) sapply(1:length(f1time), 
        function(t) {num<-dat$numreps[which(dat$subjectnum==failsubs[t])]
    ifelse(num==k.to.l[[m]],unlist(dat$zstarbar)[which(dat$subjectnum==failsubs[t])]/length(dat$numreps[which(dat$numreps==num)]),0)}))
   


    
    T.ik.t.1<-sapply(1:length(dat$subjectnum), function(i) sum(D.ik.t[[i]]*E.WwdN.t[[which(k.to.l==dat$numreps[i])]]*
              a.m[which(k.to.l==dat$numreps[i])]))
   
  
 

    #For parts 2 and 3 of Tik, need to split into groups by number of replicates
    rev.Y.i.t<-lapply(1:length(f1time), function(t) ifelse((dat$times>=f1time[t] & dat$delta %in% c(0,1)) | dat$delta==2,1,0))
 
    num.list<-lapply(1:length(subj.list), function(t) sapply(1:length(subj.list[[t]]), 
            function(i) dat$numreps[which(dat$subjectnum==subj.list[[t]][i])]))
  
    
    
    E.2<-lapply(1:length(f1time), function(t) sapply(1:length(k.to.l), 
          function(m) sum(((z.hat.list[[t]][which(num.list[[t]]==k.to.l[m])]-(S1.hat[t]/S0.hat[t]))*
          unlist(dat$zstarbar)[which(dat$numreps==k.to.l[m] & rev.Y.i.t[[t]]==1)]*new.beta*
          exp(z.hat.list[[t]][which(num.list[[t]]==k.to.l[m])]*new.beta)*
          w.tilde[[t]][which(num.list[[t]]==k.to.l[m])])/(S0.hat[t]))/
            length(dat$numreps[which(dat$numreps==k.to.l[m])])))
    
    #Only keep the expected values whose number of replicates match the failed subjects at the same timepoint
    E.2.rev<-lapply(1:length(k.to.l), function(m) sapply(1:length(f1time), function(t) E.2[[t]][[m]]*
          ifelse(dat$numreps[which(dat$subjectnum==failsubs[t])]==k.to.l[m],1,0)))

    

    #For T.ik.t.2, multiple E.2.rev by D.ik.t, but also divide by n to account for the R.tilde.t component
    T.ik.t.2<-sapply(1:length(dat$subjectnum), function(i) sum(D.ik.t[[i]]*E.2.rev[[which(k.to.l==dat$numreps[i])]]*dR.tilde.t*
              a.m[which(k.to.l==dat$numreps[i])]))
    
    
    E.3<-lapply(1:length(f1time), function(t) sapply(1:length(k.to.l),
          function(m) (sum((unlist(dat$zstarbar)[which(dat$numreps==k.to.l[m] & rev.Y.i.t[[t]]==1)]*
          exp(z.hat.list[[t]][which(num.list[[t]]==k.to.l[m])]*new.beta)*
          w.tilde[[t]][which(num.list[[t]]==k.to.l[m])])/(S0.hat[t]))/
          length(dat$numreps[which(dat$numreps==k.to.l[m])]))))
    

    
    E.3.rev<-lapply(1:length(k.to.l), function(m) sapply(1:length(f1time), function(t) E.3[[t]][[m]]*
           ifelse(dat$numreps[which(dat$subjectnum==failsubs[t])]==k.to.l[m],1,0)))
    
    T.ik.t.3<-sapply(1:length(dat$subjectnum), function(i) sum(D.ik.t[[i]]*E.3.rev[[which(k.to.l==dat$numreps[i])]])*dR.tilde.t*
                       a.m[which(k.to.l==dat$numreps[i])])
    
    
    
    ######################################################################
    #Parts that were already iid without manipulation:
    
    #part3 is only subjects that failed of type 1
    part3.a<-sapply(1:length(dat$subjectnum), function(i) ifelse(dat$delta[i]==1,
            z.hat.vec[which(f1time==dat$times[i])]-(S1.hat[which(f1time==dat$times[i])]/S0.hat[which(f1time==dat$times[i])]),0))

    
    #partA includes only timepoints with a type 1 failure, but add up all subjects at risk at those timepoints
    partA.a<-sapply(1:length(dat$subjectnum), function(i) sum(sapply(1:length(f1time),
          function(t) ifelse(Y.i.t[[i]][t]==1,w.tilde[[t]][which(subj.list[[t]]==dat$subjectnum[i])]*
            exp(z.hat.list[[t]][which(subj.list[[t]]==dat$subjectnum[i])]*new.beta)*
            (z.hat.list[[t]][which(subj.list[[t]]==dat$subjectnum[i])]-(S1.hat[t]/S0.hat[t]))*dR.tilde.t/
            (S0.hat[t]),0))))
    
    
   
    ###################Part A
    
    
    #Start by fixing j
    #Then fix u
    #Then fix subject i with k replicates
    #Then add up all t with t >= u > X_i
    
    #Then take expected value across i
    #Then add up over u
    #Then we have one value for each subject j
    
    
     f.t.a<-lapply(1:length(k.to.l), function(m) sapply(1:length(dat$times), function(u) {
       tlist<-which(f1time>=dat$times[u])
       ifelse(length(tlist>0), sum(sapply(tlist, function(t) 
       {sublist<-which(failtimes.at.risk[[t]]<dat$times[u] & numreps.at.risk[[t]]==k.to.l[m])
       sum(w.tilde[[t]][sublist]*exp(z.hat.list[[t]][sublist]*new.beta)*
             (z.hat.list[[t]][sublist]-(S1.hat[t]/S0.hat[t]))*dR.tilde.t/
             (S0.hat[t]))/
         length(dat$subjectnum==k.to.l[m])})),0)}))
     
     partA.b.2<-sapply(1:length(dat$subjectnum), function(i)
      (-1)*a.m[which(k.to.l==dat$numreps[i])]*sum(f.t.a[[which(k.to.l==dat$numreps[i])]]*dMjc.t[[i]]/pi.hat.u))
    
  
    
    

   
    
   ###################################
    #Part 3
    
    #Same adding up scheme as part A
 

    
     f.t.u<-lapply(1:length(k.to.l), function(m) sapply(1:length(f1time), function(t) {ulist<-which(dat$times<=f1time[t])
       sum(sapply(ulist, 
     function(u){ sublist<-which(failtimes.at.risk[[t]]<dat$times[u] & numreps.at.risk[[t]]==k.to.l[m])
     sum(w.tilde[[t]][sublist]*exp(z.hat.list[[t]][sublist]*new.beta)*
       (((z.hat.list[[t]][sublist]-(S1.hat[t]/S0.hat[t]))*new.beta*unlist(dat$zstarbar)[sublist])+
       unlist(dat$zstarbar)[sublist]))/length(dat$subjectnum[which(numreps==k.to.l[m])])})*EdMc.u[ulist]/
         (pi.hat.u[ulist]*S0.hat[t]))})*dR.tilde.t)
     
     
     partD.E.b.3<-sapply(1:length(dat$subjectnum), function(i) a.m[which(k.to.l==dat$numreps[i])]*
                           sum((-1)*f.t.u[[which(k.to.l==dat$numreps[i])]]*D.ik.t[[i]]))

   

    finalvar<-(sum((T.ik.t.1-T.ik.t.2-T.ik.t.3
                     -partA.b.2-partD.E.b.3
                     +part3.a-partA.a)**2)/n)/(n*A.mat**2)
   sqrt(finalvar)
   
   
    #Fine & Gray estimates using crr function - for simulation Table 1
    crr.cols<-crr(ftime=dat$times, fstatus=dat$delta, cov1=unlist(dat$zstarbar), failcode=1, cencode=0)
    crr.beta<-summary(crr.cols)$coef[,1]
    crr.se<-summary(crr.cols)$coef[,3]
    
    ci<-c(new.beta-(qnorm(0.975)*sqrt(finalvar)),new.beta+(qnorm(0.975)*sqrt(finalvar)))
    coverage<-ifelse(ci[1]<=b.11 & b.11<=ci[2],1,0)
    
    crr.ci<-c(crr.beta-(qnorm(0.975)*crr.se),crr.beta+(qnorm(0.975)*crr.se))
    crr.coverage<-ifelse(crr.ci[1]<=b.11 & b.11<=crr.ci[2],1,0)
    
    mse<-((new.beta-b.11)**2)+finalvar
    
    crr.mse<-((crr.beta-b.11)**2)+((crr.se)**2)
    
  
    return(c(new.beta,sqrt(finalvar),coverage,crr.beta,crr.se,crr.coverage,percents[1],percents[2],percents[3],ratio))
    

}

  registerDoParallel(detectCores()) 
  results.mat<-foreach(k=1:reps, .options.RNG=1234, .combine=rbind) %dorng% {fixthatME(k)}
  stopImplicitCluster()
  
  bias.new<-mean(as.numeric(results.mat[,1])-b.11)
  sd.new<-sd(as.numeric(results.mat[,1]))
  se.new<-mean(as.numeric(results.mat[,2]))
  mse.new<-bias.new**2 + se.new**2
  coverage.new<-mean(as.numeric(results.mat[,3]))
  
  bias.crr<-mean(as.numeric(results.mat[,4])-b.11)
  sd.crr<-sd(as.numeric(results.mat[,4]))
  se.crr<-mean(as.numeric(results.mat[,5]))
  mse.crr<-bias.crr**2 + se.crr**2
  coverage.crr<-mean(as.numeric(results.mat[,6]))
  
  meanvec1<-c(round(bias.new,3),
             round(sd.new,3),
             round(se.new,3),
             round(mse.new,3),
             round(coverage.new,3))
  meanvec2<-c(round(bias.crr,3),
              round(sd.crr,3),
              round(se.crr,3),
              round(mse.crr,3),
              round(coverage.crr,3))

 
  
  labels<-c("Bias","SD","SE","MSE","Coverage")
  print(cbind(labels,meanvec1,meanvec2))
 
  print(paste("Censored: ",round(mean(as.numeric(results.mat[,7])),3),sep=""))
  print(paste("AD: ",round(mean(as.numeric(results.mat[,8])),3),sep=""))
  print(paste("Death: ",round(mean(as.numeric(results.mat[,9])),3),sep=""))
 print(paste("Ratio: ",round(mean(as.numeric(results.mat[,10])),3),sep=""))
