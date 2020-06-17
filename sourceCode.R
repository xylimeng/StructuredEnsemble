# source code for double-spike Dirichlet priors

### Utility function for calculating Metropolis-Hastings acceptance rate
loglike.decouple=function(x,y,A,invsigma2){
  beta_tp=A/sum(A)
  return((sum(x%*%beta_tp*y)-0.5*sum((x%*%beta_tp)^2))*invsigma2)
}

# logACC.decouple calculates the log value of Metropolis-Hastings acceptance rate
logAcc.decouple=function(invsigma2,x,y,gam_new,A_new,gam_old,A_old,theta,K){
  return(loglike.decouple(x,y,A_new,invsigma2)-loglike.decouple(x,y,A_old,invsigma2)+sum(gam_new-gam_old)*(log(theta)-log(1-theta)))
}

### Obtain the posterior with the double-spike Dirichlet prior 
# alpha1, alpha2 are the prior hyperparameters of the variance paramter \sigma^2
gprior.full.decouple=function(x,y,theta,rho1,rho2,alpha1=0.01,alpha2=0.01,iter,prob_add=0.25,prob_delete=0.25,prob_swap=0.25){
  #initialization
  K=ncol(x)
  gam_list=matrix(0,nrow=iter,ncol = K)
  A_list=matrix(nrow=iter,ncol=K)
  #use lm as a good initial value, if possible
  beta.lm=coef(lm(y~x+0))
  if(sum(is.na(beta.lm))==0){
    large.bound=sort(beta.lm,decreasing = TRUE)[min(K/2,10)]
    gam_list[1,which(beta.lm>=large.bound)]=1 
  }
  else{
    gam_list[1,]=rbinom(K,size = 1,prob = theta)
  }
  #avoid assigning all 0 to the vecotr gamma as an initial value
  if(sum(gam_list[1,])==0){
    gam_list[1,sample(1:K,1)]=1
  }
  
  A_list[1,]=rgamma(K,shape = rho1*gam_list[1,]+rho2*(1-gam_list[1,]),rate = 1)
  #avoid assigning all 0 to the vector A_list
  if(sum(A_list[1,])==0){
    A_list[1,]=A_list[1,]+1e-50
  }
  beta_list=matrix(nrow=iter,ncol=K)
  beta_list[1,]=A_list[1,]/sum(A_list[1,])
  invsigma2_list=rep(NA,iter)
  invsigma2_list[1]=1
  acc_list=rep(0,iter)
  
  #start
  for (i in 2:iter) {
    large_ind=which(gam_list[i-1,]>0)
    small_ind=which(gam_list[i-1,]==0)
    #propose a candidate gamma; default is the previous one.
    gam_prop=gam_list[i-1,]
    #random search sampling for the candidate gam
    unif_gam=runif(1)
    #with prob_add, randomly add a gam=1
    if(length(large_ind)==0 | unif_gam<prob_add){
      chosen=sample(small_ind)[1]
      gam_prop[chosen]=1
    }
    #with prob_delete,randomly delete a gam=1
    else if(length(small_ind)==0 |unif_gam<prob_add+prob_delete){
      chosen=sample(large_ind)[1]
      gam_prop[chosen]=0
    }
    #with prob_swap,randomly swap a gam=1 and a gam=0
    else if(unif_gam<prob_add+prob_delete+prob_swap){
      chosen1=sample(large_ind)[1]
      chosen2=sample(small_ind)[1]
      gam_prop[chosen1]=0
      gam_prop[chosen2]=1
    }
    #propose A's based on the gam candidate
    A_prop=rgamma(K,shape = rho1*gam_prop+rho2*(1-gam_prop),rate = 1)
    #avoid assigning all 0 to the vector A_prop
    if(sum(A_prop)==0){
      A_prop=A_prop+1e-50
    }
    #calculate the acceptance rate
    logacc=logAcc.decouple(invsigma2_list[i-1],x,y,gam_prop,A_prop,gam_list[i-1,],A_list[i-1,],theta,K)
    if(runif(1)<exp(logacc)){
      A_list[i,]=A_prop
      gam_list[i,]=gam_prop
      acc_list[i]=1
    }
    else
    {
      A_list[i,]=A_list[i-1,]
      gam_list[i,]=gam_list[i-1,]
    }
    
    beta_list[i,]=A_list[i,]/sum(A_list[i,])
    invsigma2_list[i]=rgamma(1,shape = length(y)/2+alpha1, rate = alpha2+sum((y-x%*%beta_list[i,])^2)/2)
  }
  return(list("betas"=beta_list,"gams"=gam_list,"accs"=acc_list,'invsigma2s'=invsigma2_list))
}
#-------------------end of gprior model-------------------------

