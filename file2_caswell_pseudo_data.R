#Data Example JASA Manuscript Submission
#Author: Carrie Caswell
#September 2017
#Reviewed July 2018


if(!("cmprsk" %in% rownames(installed.packages()))){install.packages("cmprsk")}
if(!("rootSolve" %in% rownames(installed.packages()))){install.packages("rootSolve")}
if(!("survival" %in% rownames(installed.packages()))){install.packages("survival")}
if(!("optimx" %in% rownames(installed.packages()))){install.packages("optimx")}
library(cmprsk)
library(rootSolve)
library(survival)
library(optimx)

#########################################################################
#Generate pseudo dataset
#This dataset should have 70% censoring, which emulates ADNI data
set.seed(5677)

n<-200 #number of subjects
maxk<-4 #number of replicates per subject (k)
b.11<-log(2) #true log hazard ratio for event of interest
b.12<-log(0.7)
b.21<-b.11 #true log hazard ratio for competing event
b.22<-b.12
numvars<-2
sigma2<-0.2


p<-0.3 #parameter fixed to generate data according to subdistribution (see Fine & Gray 1999)
a<-0.5 #ensures approx. 30% censoring
b<-2



dat<-data.frame(subjectnum=1:n)


dat$z1<-rnorm(n,0,1)
dat$z2<-rbinom(n,size=1,prob=0.5)



#Data generation process:
dat$zB1<-(dat$z1*b.11)+(dat$z2*b.12)
dat$zB2<-(dat$z1*b.21)+(dat$z2*b.22)

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

dat$ttime<-dat$times
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
z2.l<-lapply(1:length(numreps), function(i) rep(dat$z2[i],numreps[i]))

#zstar is covariate replicates with measurement error 
zstar.l<-list(lapply(1:length(z1.l), function(i) z1.l[[i]]+eps.l[[i]]),
              z2.l)

#########################################################################
#C-RRC Method Application


#dat$zstarbar is the average observed for each subject (W.i.bar in manuscript)
zstarbar.l<-lapply(1:numvars,function(p) sapply(zstar.l[[p]],mean))
zstarmed.l<-lapply(1:numvars, function(i) sapply(zstar.l[[i]],median))



sublist<-dat$subjectnum

####################################################################################################
####################################################################################################


dev<-lapply(1:numvars, function(p) lapply(1:n, function(i) zstar.l[[p]][[i]]-rep(zstarbar.l[[p]][i],numreps[i])))
dev.prelim<-lapply(1:n, function(i) lapply(1:numreps[i], function(j) sapply(1:numvars, function(h) dev[[h]][[i]][j])))
dev.mat<-lapply(1:n, function(i) lapply(1:numreps[i], 
                                        function(j) dev.prelim[[i]][[j]]%*%t(dev.prelim[[i]][[j]])))

Deltahat<-(Reduce("+",lapply(1:n, function(i) Reduce("+",dev.mat[[i]]))))/(sum(numreps-1))


# #Check Deltahat
# devcheck<-list()
# for(i in 1:n){
# devcheck[[i]]<-list()
# for(j in 1:numreps[i]){
#     devcheck[[i]][j]<-(zstar.l[[1]][[i]][j]-zstarbar.l[[1]][i])*(zstar.l[[1]][[i]][j]-zstarbar.l[[1]][i])
#   }
# }
# sum(unlist(devcheck))/(sum(numreps)-n)



#Now identify the at-risk sets for each type 1 failure timepoint
#For each subject


#vector of type 1 failure times
f1time<-dat$ttime[which(dat$delta.1==1)]

#zstar.vec is the average observed for each subject who had a type 1 failure
zstarbar.fail<-lapply(1:numvars, function(i) zstarbar.l[[i]][which(dat$delta==1)])
zstarbar.fail.forzhat<-lapply(1:length(f1time), function(t) sapply(1:numvars, function(h) zstarbar.fail[[h]][t]))


#censoring distribution for weights
censind<-ifelse(dat$delta==0,1,0) #This counts only true censors as censors
G<-survfit(Surv(dat$ttime,censind)~1,type="kaplan-meier")


atrisk<-function(i)
{
  dat$Y.t.i<-ifelse((dat$ttime>=f1time[i] & dat$delta %in% c(0,1)) | dat$delta==2,1,0)
  return(dat$subjectnum[which(dat$Y.t.i==1)])
}



#identify risk set at each f1time
subj.list<-lapply(1:length(f1time),atrisk)

zstarbar.atrisk.list<-lapply(1:numvars,function(i) lapply(1:length(f1time),
                                                          function(j) (zstarbar.l[[i]][which(dat$subjectnum %in% subj.list[[j]])])))

censweight<-function(i)
{
  #see picture saved in this folder with censoring weight chart
  w.j<-sapply(subj.list[[i]],function(j)
    ifelse((dat[j,"ttime"]>=f1time[i]),1,
           ifelse((dat[j,"ttime"] < f1time[i] & dat[j,"delta"]!=0), G$surv[round(G$time,8)==round(f1time[i],8)]/
                    G$surv[round(G$time,8)==round(dat[j,"ttime"],8)],0)))
  return(w.j)
}

wj.list<-lapply(1:length(f1time),censweight)

#rearrange zstarbar.atrisk.list for zhat calculation
zstarbar.atrisk.forzhat<-lapply(1:length(f1time), 
                                function(j) lapply(1:length(subj.list[[j]]), 
                                                   function(i) sapply(1:numvars, function(h) zstarbar.atrisk.list[[h]][[j]][i])))


#Find mu hat for each t
#mu.hat.t<-lapply(1:length(subj.list), function(j) sapply(1:numvars, 
#       function(i) (sum(zstarbar.atrisk.list[[i]][[j]])/length(zstarbar.atrisk.list[[i]][[j]]))))
if(numvars>1){
mu.hat.t<-lapply(1:length(f1time), function(t) apply(sapply(1:length(subj.list[[t]]),
      function(i) (zstarbar.atrisk.forzhat[[t]][[i]])),1,sum)/(length(subj.list[[t]])))
}
if(numvars==1){
mu.hat.t<-lapply(1:length(f1time), function(t) sum(sapply(1:length(subj.list[[t]]),
          function(i) zstarbar.atrisk.forzhat[[t]][[i]]))/length(subj.list[[t]]))
}
#mu.hat.t.old<-lapply(1:numvars, 
#                 function(i) sapply(1:length(subj.list), function(j) (sum(zstarbar.atrisk.list[[i]][[j]])/length(zstarbar.atrisk.list[[i]][[j]]))))

#Find Sigma for each t
k.i<-lapply(1:length(subj.list), function(j) (numreps[subj.list[[j]]]))
#k.i<-lapply(1:length(subj.list), function(i) (rep(samereps,length(subj.list))))


#05JUN2017 Sigma hat should be a p xp matrix
Sigma.mat<-lapply(1:length(f1time), function(t) lapply(1:length(subj.list[[t]]), 
                                                       function(i) ((zstarbar.atrisk.forzhat[[t]][[i]])-
                                                                      mu.hat.t[[t]])%*%t(zstarbar.atrisk.forzhat[[t]][[i]]-mu.hat.t[[t]])))
Sigma.sum<-lapply(1:length(f1time), function(t) Reduce("+",Sigma.mat[[t]]))



Sigma.hat.t<-lapply(1:length(f1time), 
                    function(t) (Sigma.sum[[t]]/(length(subj.list[[t]])-1))-((Deltahat*sum(1/k.i[[t]]))/length(subj.list[[t]])))



#Find z.hat for each subject at risk at each t.i
z.hat.list<-lapply(1:length(f1time), function(t) lapply(1:length(subj.list[[t]]), 
                                                        function(i) ((((Deltahat)/k.i[[t]][i])%*%solve(Sigma.hat.t[[t]]+(Deltahat/k.i[[t]][i])))%*%mu.hat.t[[t]])+
                                                          (Sigma.hat.t[[t]]%*%solve(Sigma.hat.t[[t]]+(Deltahat/k.i[[t]][i]))%*%unlist(zstarbar.atrisk.forzhat[[t]][i]))))


#Now take z.hat for each subject who had a type 1 event
k.i.vec<-numreps[which(dat$delta==1)]
#z.hat.vec<-lapply(1:numvars, function(i) sapply(1:length(k.i.vec[[i]]), 
#                             function(j) (((Deltahat[i]/k.i.vec[[i]][j])*mu.hat.t[[i]][j])/(Sigma.hat.t[[i]][j]+Deltahat[i]/k.i.vec[[i]][j])+
#                             (Sigma.hat.t[[i]][j]*zstarbar.fail[[i]][j]/(Sigma.hat.t[[i]][j]+Deltahat[i]/k.i.vec[[i]][j])))))

z.hat.vec<-lapply(1:length(f1time), function(t) ((Deltahat/k.i.vec[t])%*%solve(Sigma.hat.t[[t]]+
                                                                                 (Deltahat/k.i.vec[t]))%*%mu.hat.t[[t]])+
                    (Sigma.hat.t[[t]]%*%solve(Sigma.hat.t[[t]]+(Deltahat/k.i.vec[t]))%*%zstarbar.fail.forzhat[[t]]))


#05JUN2017: Rearrange z.hat.list and z.hat.vec for use in the optimizing functions
z.hat.list.2<-lapply(1:numvars, function(h) 
  lapply(1:length(f1time), function(t) sapply(1:length(subj.list[[t]]), function(i) z.hat.list[[t]][[i]][h])))
z.hat.vec.2<-lapply(1:numvars, function(h)
  sapply(1:length(f1time), function(t) z.hat.vec[[t]][h]))

loglik<-function(param)
{
  (-1)*sum(sapply(1:length(f1time), function(t) param%*%z.hat.vec[[t]]-
                    log(sum(sapply(1:length(wj.list[[t]]), function(i) wj.list[[t]][[i]]*exp(param%*%z.hat.list[[t]][[i]]))))))
}

new.zhat<-optimx(par=rep(1,numvars),fn=loglik,method="BFGS")

new.beta<-as.numeric(new.zhat[1:numvars])

new.beta
###############################################################################################
#Confidence interval using asymptotic variance
###############################################################################################

#Dimension of S0 should be 1 x 1 at each timepoint
#Dimension of S1 should be p x 1 at each timepoint
#Dimension of S2 should be p x p at each timepoint


S0.hat<-sapply(1:length(f1time), function(t) 
  sum(wj.list[[t]]*exp(sapply(1:length(z.hat.list[[t]]), function(i) new.beta%*%z.hat.list[[t]][[i]])))/n)
# wj.list[[1]]*z.hat.list[[1]]*exp(sapply(1:length(z.hat.list[[1]]), function(i) new.beta%*%z.hat.list[[1]][[i]]))

#S1.hat.new<-lapply(1:length(f1time), function(t) Reduce("+",lapply(1:length(wj.list[[t]]), 
#       function(i) {exppart<-sapply(1:length(z.hat.list[[t]]), function(i) exp(new.beta%*%z.hat.list[[t]][[i]]))
# wj.list[[t]][i]*z.hat.list[[t]][[i]]*exppart[i]}))/n)

#S1.hat.what<-lapply(1:length(z.hat.list), function(t) Reduce("+",lapply(1:length(z.hat.list[[t]]), 
#       function(i) wj.list[[t]][[i]]*z.hat.list[[t]][[i]]*
#       as.numeric(exp(new.beta%*%z.hat.list[[t]][[i]]))))/n)


S1.hat<-lapply(1:length(f1time), function(t) sapply(1:numvars, 
                                                    function(p) sum(z.hat.list.2[[p]][[t]]*wj.list[[t]]*exp(apply(sapply(1:numvars, 
                                                                                                                         function(g) z.hat.list.2[[g]][[t]]*new.beta[g]),1,sum)))/n))

S2.hat<-lapply(1:length(z.hat.list), function(t) Reduce("+",lapply(1:length(z.hat.list[[t]]), 
                                                                   function(i) wj.list[[t]][[i]]*z.hat.list[[t]][[i]]%*%t(z.hat.list[[t]][[i]])*
                                                                     as.numeric(exp(new.beta%*%z.hat.list[[t]][[i]]))))/n)

A.mat<-Reduce("+",lapply(1:length(f1time), 
  function(t) (S2.hat[[t]]/S0.hat[[t]])-(S1.hat[[t]]/S0.hat[[t]])%*%
    t(S1.hat[[t]]/S0.hat[[t]])))/n


################################333
#Sum of iid variables part

#Assemble D.ik
#Split the population into those with the same number of reps together
Y.i.t<-lapply(1:length(dat$subjectnum), function(i) ifelse((dat$ttime[i]>=f1time & dat$delta[i] %in% c(0,1)) | dat$delta[i]==2,1,0))

E.Yt<-sapply(1:length(f1time),function(i) length(wj.list[[i]])/n)


k.to.l<-sort(unique(numreps))
a.m<-sapply(1:length(k.to.l), function(m) length(numreps[which(numreps==k.to.l[m])])/n)

Gamma1.t<-lapply(1:length(k.to.l), function(m) lapply(1:length(f1time), 
                                                      function(t) Sigma.hat.t[[t]]%*%solve(Sigma.hat.t[[t]]+(Deltahat/k.to.l[[m]]))/E.Yt[t]))

Gamma2.t<-lapply(1:length(k.to.l), function(m) lapply(1:length(f1time),
                                                      function(t) solve(Sigma.hat.t[[t]]+(Deltahat/k.to.l[[m]]))))





Delta.i<-lapply(1:n, function(i) (Reduce("+",dev.mat[[i]]) - ((numreps[i]-1)*Deltahat))/sum(a.m*(k.to.l-1)))


zstarbar.l.2<-lapply(1:n, function(i) sapply(1:numvars, function(p) zstarbar.l[[p]][[i]]))

D.ik.t<-lapply(1:n, function(i) lapply(1:length(f1time), function(t) (Gamma1.t[[which(k.to.l==numreps[i])]][[t]]%*%(((Deltahat/(k.to.l[which(k.to.l==numreps[i])])))%*%solve(Sigma.hat.t[[t]])%*%
  (Y.i.t[[i]][t]*((zstarbar.l.2[[i]]-mu.hat.t[[t]])%*%t(zstarbar.l.2[[i]]-mu.hat.t[[t]])) - 
     (Delta.i[[i]]*E.Yt[t]*sum(a.m/k.to.l)) - (Y.i.t[[i]][t]*Deltahat/numreps[i]) - Y.i.t[[i]][t]*Sigma.hat.t[[t]]) -
  ((Delta.i[[i]]*E.Yt[t])/(k.to.l[which(k.to.l==numreps[i])])))%*%Gamma2.t[[which(k.to.l==numreps[i])]][[t]])))



########################################
#Part 2 of variance

dR.tilde.t<-1/n
dat$numreps<-numreps
failsubs<-dat$subjectnum[which(dat$delta==1)]


E.WwdN.t<-lapply(1:length(k.to.l), function(m) lapply(1:length(f1time),
    function(t) {num<-dat$numreps[which(dat$subjectnum==failsubs[t])]
    ifelse(rep(num==k.to.l[[m]],numvars),zstarbar.l.2[[which(dat$subjectnum==failsubs[t])]]/
             length(dat$numreps[which(dat$numreps==num)]),
    rep(0,numvars))}))
    



T.ik.t.1<-lapply(1:n, function(i) Reduce("+", lapply(1:length(f1time), 
     function(t) D.ik.t[[i]][[t]]%*%E.WwdN.t[[which(k.to.l==dat$numreps[i])]][[t]]*
       a.m[which(k.to.l==dat$numreps[i])])))



#For parts 2 and 3 of Tik, need to split into groups by number of replicates

num.list<-lapply(1:length(subj.list), function(t) sapply(1:length(subj.list[[t]]), 
       function(i) dat$numreps[which(dat$subjectnum==subj.list[[t]][i])]))


w.tilde<-lapply(1:length(f1time),function(t) sapply(subj.list[[t]],function(j)
  ifelse((dat[j,"ttime"]>=f1time[t]),1,
         ifelse((dat[j,"ttime"] < f1time[t] & dat[j,"delta"]!=0), G$surv[round(G$time,8)==round(f1time[t],8)]/
                  G$surv[round(G$time,8)==round(dat[j,"ttime"],8)],0))))




E.2<-lapply(1:length(f1time), function(t) lapply(1:length(k.to.l), 
     function(m) {
    worklist<-z.hat.list[[t]][which(num.list[[t]]==k.to.l[m])]
    if(length(worklist)>0){
     z.minus.s<-lapply(1:length(worklist), function(i) worklist[[i]]-(S1.hat[[t]]/S0.hat[[t]]))
     worklist2<-zstarbar.l.2[which(num.list[[t]]==k.to.l[m])]
     beta.w<-sapply(1:length(worklist2), function(i) worklist2[[i]]%*%new.beta)
     nextpiece<-sapply(1:length(worklist), function(i) exp(t(worklist[[i]])%*%new.beta))
     therest<-w.tilde[[t]][which(num.list[[t]]==k.to.l[m])]/(S0.hat[t])
     all<-lapply(1:length(z.minus.s), function(i) z.minus.s[[i]]*beta.w[i]*nextpiece[i]*therest[i])
     ret<-Reduce("+",all)/length(dat$numreps==k.to.l[m])}
    if(length(worklist)==0){
      ret<-rep(0,numvars)
    }
    ret
     }))

#Only keep the expected values whose number of replicates match the failed subjects at the same timepoint
E.2.rev<-lapply(1:length(k.to.l), function(m) lapply(1:length(f1time), 
   function(t) E.2[[t]][[m]]*ifelse(dat$numreps[which(dat$subjectnum==failsubs[t])]==k.to.l[m],1,0)))


#For T.ik.t.2, multiple E.2.rev by D.ik.t, but also divide by n to account for the R.tilde.t component
T.ik.t.2<-lapply(1:n, 
  function(i) Reduce("+",lapply(1:length(f1time), 
  function(t) D.ik.t[[i]][[t]]%*%E.2.rev[[which(k.to.l==dat$numreps[i])]][[t]]*
    dR.tilde.t*a.m[which(k.to.l==dat$numreps[i])])))



E.3<-lapply(1:length(f1time), function(t) lapply(1:length(k.to.l),
   function(m) {
     w.piece<-zstarbar.l.2[which(num.list[[t]]==k.to.l[m])]
     if(length(w.piece)>0){
       therest<-w.tilde[[t]][which(num.list[[t]]==k.to.l[m])]/(S0.hat[t])
   worklist<-z.hat.list[[t]][which(num.list[[t]]==k.to.l[m])]
   nextpiece<-sapply(1:length(worklist), function(i) exp(t(worklist[[i]])%*%new.beta))
   all<-lapply(1:length(worklist), function(j) therest[j]*nextpiece[j]*w.piece[[j]])
   ret<-Reduce("+",all)/length(dat$numreps==k.to.l[m])}
    if(length(w.piece)==0){
      ret<-rep(0,numvars)
    }
     ret
   }))



E.3.rev<-lapply(1:length(k.to.l), function(m) lapply(1:length(f1time), function(t) E.3[[t]][[m]]*
   ifelse(dat$numreps[which(dat$subjectnum==failsubs[t])]==k.to.l[m],1,0)))

T.ik.t.3<-lapply(1:n, 
  function(i) Reduce("+",lapply(1:length(f1time),
  function(t) D.ik.t[[i]][[t]]%*%E.3.rev[[which(k.to.l==dat$numreps[i])]][[t]]*
    dR.tilde.t*a.m[which(k.to.l==dat$numreps[i])]))) 



######################################################################
#Parts that were already iid without manipulation:

#part3 is only subjects that failed of type 1
part3.a<-lapply(1:n, function(i) ifelse(rep(dat$delta[i]==1,numvars),
      z.hat.vec[[which(failsubs==dat$subjectnum[i])]]-
    (S1.hat[[which(failsubs==dat$subjectnum[i])]]/S0.hat[which(failsubs==dat$subjectnum[i])]),
      rep(0,numvars)))


#partA includes only timepoints with a type 1 failure, but add up all subjects at risk at those timepoints
partA.a<-lapply(1:n, function(i) Reduce("+",lapply(1:length(f1time), function(t) ifelse(rep(Y.i.t[[i]][t]==1,numvars),
  w.tilde[[t]][which(subj.list[[t]]==dat$subjectnum[i])]*
    as.numeric(exp(t(z.hat.list[[t]][[which(subj.list[[t]]==dat$subjectnum[i])]])%*%new.beta))*
    (z.hat.list[[t]][[which(subj.list[[t]]==dat$subjectnum[i])]]-(S1.hat[[t]]/S0.hat[t]))*dR.tilde.t/
    (S0.hat[t]),
  rep(0,numvars)))))

final.expected<-Reduce("+",lapply(1:length(T.ik.t.1), function(i) (T.ik.t.1[[i]]-T.ik.t.2[[i]]-T.ik.t.3[[i]]+
   part3.a[[i]]-partA.a[[i]])%*%t(T.ik.t.1[[i]]-T.ik.t.2[[i]]-T.ik.t.3[[i]]+
                                                                                                      part3.a[[i]]-partA.a[[i]])))/n
finalvar<-(solve(A.mat)%*%final.expected%*%solve(A.mat))/n




se.new<-sqrt(diag(finalvar))
finalvar.old<-finalvar


crr.cols<-crr(ftime=dat$ttime,fstatus=dat$delta,
              cov1=cbind(zstarmed.l[[1]],zstarmed.l[[2]]),
              failcode=1,cencode=0)
# crr.cols<-crr(ftime=dat$ttime,fstatus=dat$delta,cov1=dat$med.n100,
#              failcode=1,cencode=0)
crr.beta<-as.numeric(crr.cols$coef)
crr.se<-as.numeric(summary(crr.cols)$coef[,3])
crr.beta
crr.se

cilower<-new.beta-(qnorm(0.975)*se.new)
ciupper<-new.beta+(qnorm(0.975)*se.new)
crr.lower<-crr.beta-(qnorm(0.975)*crr.se)
crr.upper<-crr.beta+(qnorm(0.975)*crr.se)



labels<-c("C-RRC Estimate","C-RRC SE","C-RRC CI","SHR Estimate","SHR SE","SHR CI")

z1vec<-c(round(new.beta[1],3),round(se.new[1],3),paste("(",round(cilower[1],3),", ",round(ciupper[1],3),")",sep=""),
          round(crr.beta[1],3),round(crr.se[1],3),paste("(",round(crr.lower[1],3),", ",round(crr.upper[1],3),")",sep=""))
z2vec<-c(round(new.beta[2],3),round(se.new[2],3),paste("(",round(cilower[2],3),", ",round(ciupper[2],3),")",sep=""),
             round(crr.beta[2],3),round(crr.se[2],3),paste("(",round(crr.lower[2],3),", ",round(crr.upper[2],3),")",sep=""))
print(rbind(labels,z1vec,z2vec))
