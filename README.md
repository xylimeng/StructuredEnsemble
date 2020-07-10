# Double spike Dirichlet priors
Bayesian procedure using double-spike Dirichlet priors for inferring structured probabiliy simplexes. 



## Examples using Concrete dataset from the UCI repository  

The followling demonstration is also in the "example.R" script.

### Improve random forests predictions using Bayesian procedure with double spike Dirichlet priors

```R
rm(list = ls())
source("sourceCode.R")

### the concrete dataset from the UCI repository
library(readxl)
Concrete_Data <- read_excel("Concrete_Data.xls")
colnames(Concrete_Data)=c('component1','component2','component3','component4','component5','component6','component7','age','ccs')

# train-test split and other parameters
set.seed(123)
theta=0.2
nsample=nrow(Concrete_Data)
train_id=sample(1:nsample, round(nsample *0.5))
test_id=setdiff(1:nsample,train_id)
burning=20000
iteration=30000

# fit a random forest model
library(randomForest)
rf.mdl=randomForest(ccs~.,data=Concrete_Data[train_id,])#use default number of trees, 500.

rf.allpred=predict(rf.mdl,data.frame(Concrete_Data[,-9]),predict.all=TRUE)
#aggregate prediction of random forest
rf.agg=rf.allpred$aggregate
#individual tree prediction of random forest
rf.indiv=rf.allpred$individual


K=ncol(rf.indiv)
rho1=K^1.5
rho2=1/K
# apply the Bayesian method using double-spike Dirichlet prior
gp.results=gprior.full.decouple(x=rf.indiv[train_id,],y=Concrete_Data$ccs[train_id],theta=theta,rho1=rho1,rho2=rho2,
                                   alpha1 = 0.01,alpha2 = 0.01,iter=iteration)
# posterior mean \beta
gp_beta_mean=colMeans(gp.results$betas[burning:iteration,])
# predictions of the Bayesian procedure using double-spike Dirichlet priors
gp.pred.tp=rf.indiv%*%gp_beta_mean
```


#### In and out-of sample predictions  

```R
# Compare root MSE
# out-of-sample rmse
rmse.out.rf=sqrt(mean((rf.agg[test_id]-Concrete_Data$ccs[test_id])^2))
rmse.out.gp=sqrt(mean((gp.pred.tp[test_id]-Concrete_Data$ccs[test_id])^2))
# in-sample rmse
rmse.in.rf=sqrt(mean((rf.agg[train_id]-Concrete_Data$ccs[train_id])^2))
rmse.in.gp=sqrt(mean((gp.pred.tp[train_id,]-Concrete_Data$ccs[train_id])^2))

c(rmse.out.proposed=rmse.out.gp, rmse.out.randomforest=rmse.out.rf, rmse.in.proposed=rmse.in.gp, rmse.in.randomforest=rmse.in.rf)

# output: 
#   rmse.out.proposed rmse.out.randomforest      rmse.in.proposed  rmse.in.randomforest
#            5.875607              6.522016              2.736968              3.414659 
```



## Reference: 

Lin, H. and Li, M. (2020). Double spike Dirichlet priors for structured weighting. Under review. arXiv:2007.04387v1

