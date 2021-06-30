######## create sample data ########

require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,doParallel,tictoc,patchwork,matrixcalc)

n.cores <- parallel::detectCores() - 1

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
clusterCall(my.cluster, function() require(pacman))
clusterCall(my.cluster, function() p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,doParallel))

# initialize some parameters
n <- 10000
m <- 500

truebeta <- c(4,12,-8,0,0,0)
#beta1, beta2, beta3, z, z2, z3
#betas as params in utility function
#z is unobserved factor in selection
#z2, z3 are observed factors in selection

vcov <- matrix(c( 1,  0, 0,.5,.3,.3,
                  NA, 1, 0,.5,.4,.3,
                  NA,NA, 1,0,0,0,
                  NA,NA,NA, 1,.3,.3,
                  NA,NA,NA,NA, 1,.2,
                  NA,NA,NA,NA,NA, 1),nrow=6,ncol=6,byrow=T)
vcov[lower.tri(vcov)] <- t(vcov)[lower.tri(vcov)]
if (!is.positive.definite(vcov)) {cat("\n\n\n\t\t\t!!!! VCOV NOT PD !!!!\n\n\n")}

# utility function for RUM
u <- function (b1,b2,b3,x1,x2,x3) {b1*x1 + b2*x2 + b3*x3}

###### run sim #####
set.seed(666)
tic()
{
  results <- foreach (i = 1:m,
                      .combine = "rbind") %dopar% 
    {
      
      
      subjects <- mvrnorm(n=n,mu=truebeta,Sigma=vcov) %>% as.data.frame()
      names(subjects) <- c("beta1","beta2","beta3","z","z2","z3")
      subjects$id <- 1:n
      
      choices <- data.frame(id=1:n,att1=rnorm(n),att2=rnorm(n),att3=rnorm(n)+1,eps=rlogis(n))
      
      dat <- left_join(choices,subjects,by="id")
      
      dat <- dat %>% mutate(utility = u(beta1,beta2,beta3,att1,att2,att3) + eps)
      
      dat$best <- ifelse(dat$utility > 0, 1, 0)
      dat$select <- ifelse(dat$z + dat$z2 + dat$z3 > 0,1,0)
      
      stage1 <- glm(select ~ z2 + z3,data=dat, family = binomial(link= "probit"))
      
      
      dat$RP <- predict(stage1,dat)
      dat$RP <- dat$RP - mean(dat$RP)
      
      uncor <- glm(best ~ att1 + att2 + att3, family = "binomial", 
                   data = dat[which(dat$select==1),])
      
      
      fullsamp <- glm(best ~ att1 + att2 + att3, family = "binomial", 
                      data = dat[1:sum(dat$select),])
      
      
      corrected <- glm(best ~ att1 + att2 + att3 + att1:RP + att2:RP + att3:RP, family = "binomial", 
                       data = dat[which(dat$select==1),])
      
      # results$unbiased1[i] <- -fullsamp$coefficients[["att1"]]/fullsamp$coefficients[["att3"]]
      # results$unbiased2[i] <- -fullsamp$coefficients[["att2"]]/fullsamp$coefficients[["att3"]]
      # 
      # results$uncorrected1[i] <- -uncor$coefficients[["att1"]]/uncor$coefficients[["att3"]]
      # results$uncorrected2[i] <- -uncor$coefficients[["att2"]]/uncor$coefficients[["att3"]]
      # 
      # results$corrected1[i] <- -corrected$coefficients[["att1"]]/corrected$coefficients[["att3"]]
      # results$corrected2[i] <- -corrected$coefficients[["att2"]]/corrected$coefficients[["att3"]]
      
      c(-fullsamp$coefficients[["att1"]]/fullsamp$coefficients[["att3"]],
        -fullsamp$coefficients[["att2"]]/fullsamp$coefficients[["att3"]],
        
        -uncor$coefficients[["att1"]]/uncor$coefficients[["att3"]],
        -uncor$coefficients[["att2"]]/uncor$coefficients[["att3"]],
        
        -corrected$coefficients[["att1"]]/corrected$coefficients[["att3"]],
        -corrected$coefficients[["att2"]]/corrected$coefficients[["att3"]])
      
      
      
    }
}
toc()

###### plot results #####
{
results <- results %>% as.data.frame()
names(results) <- c("unbiased1","unbiased2","uncorrected1","uncorrected2","corrected1","corrected2")

ggplot(results) +
  geom_density(aes(x=unbiased1),fill="forestgreen",alpha=.5) +
  geom_vline(xintercept=mean(results$unbiased1),color="forestgreen") +
  geom_density(aes(x=corrected1),fill="blue",alpha=.5) +
  geom_vline(xintercept=mean(results$corrected1),color="blue") +
  geom_density(aes(x=uncorrected1),fill="red",alpha=.5) +
  geom_vline(xintercept=mean(results$uncorrected1),color="red") +
  labs(x="Estimated mWTP for attribute 1") +
ggplot(results) +
  geom_density(aes(x=unbiased2),fill="forestgreen",alpha=.5) +
  geom_vline(xintercept=mean(results$unbiased2),color="forestgreen") +
  geom_density(aes(x=corrected2),fill="blue",alpha=.5) +
  geom_vline(xintercept=mean(results$corrected2),color="blue") +
  geom_density(aes(x=uncorrected2),fill="red",alpha=.5) +
  geom_vline(xintercept=mean(results$uncorrected2),color="red") +
  labs(x="Estimated mWTP for attribute 2")
}


############## heckit 2-step doesn't work ##############
############ (using LPM for outcome equation) ##########

stage1 <- glm(select ~ att1 + att2 + att3 + z2, data=dat,family = binomial(link = "probit"))
summary(stage1)
dat$yhat <- predict(stage1)
dat$imr <- dnorm(dat$yhat)/pnorm(dat$yhat)

stage2 <- lm(best ~ att1 + att2 + att3 + imr,data=dat[which(dat$select==1),])
summary(stage2)

naivelpm <- lm(best ~ att1 + att2 + att3,data=dat[which(dat$select==1),])
summary(naivelpm)

#### obsolete code for slow for loop ####
results <- data.frame(unbiased1=rep(NA,m),unbiased2=NA,uncorrected1=NA,uncorrected2=NA, corrected1=NA,corrected2=NA)
set.seed(123)
tic()
{
  results <- for (i in 1:m) {
    
    subjects <- mvrnorm(n=n,mu=truebeta,Sigma=vcov) %>% as.data.frame()
    names(subjects) <- c("beta1","beta2","beta3","z","z2","z3")
    subjects$id <- 1:n
    
    choices <- data.frame(id=1:n,att1=rnorm(n),att2=rnorm(n),att3=rnorm(n),eps=rlogis(n))
    
    dat <- left_join(choices,subjects,by="id")
    
    dat <- dat %>% mutate(utility = u(beta1,beta2,beta3,att1,att2,att3) + eps)
    
    dat$best <- ifelse(dat$utility > 0, 1, 0)
    dat$select <- ifelse(dat$z + dat$z2 + dat$z3 > 0,1,0)
    
    stage1 <- glm(select ~ z2 + z3 + 0,data=dat, family = binomial(link= "probit"))
    summary(stage1)
    
    dat$RP <- predict(stage1,dat)
    dat$RP <- dat$RP - mean(dat$RP)
    
    uncor <- glm(best ~ att1 + att2 + att3, family = "binomial", 
                 data = dat[which(dat$select==1),])
    
    
    fullsamp <- glm(best ~ att1 + att2 + att3, family = "binomial", 
                    data = dat[1:sum(dat$select),])
    
    
    corrected <- glm(best ~ att1 + att2 + att3 + att1:RP + att2:RP + att3:RP + 0, family = "binomial", 
                     data = dat[which(dat$select==1),])
    
    # results$unbiased1[i] <- -fullsamp$coefficients[["att1"]]/fullsamp$coefficients[["att3"]]
    # results$unbiased2[i] <- -fullsamp$coefficients[["att2"]]/fullsamp$coefficients[["att3"]]
    # 
    # results$uncorrected1[i] <- -uncor$coefficients[["att1"]]/uncor$coefficients[["att3"]]
    # results$uncorrected2[i] <- -uncor$coefficients[["att2"]]/uncor$coefficients[["att3"]]
    # 
    # results$corrected1[i] <- -corrected$coefficients[["att1"]]/corrected$coefficients[["att3"]]
    # results$corrected2[i] <- -corrected$coefficients[["att2"]]/corrected$coefficients[["att3"]]
    
    c(-fullsamp$coefficients[["att1"]]/fullsamp$coefficients[["att3"]],
      -fullsamp$coefficients[["att2"]]/fullsamp$coefficients[["att3"]],
      
      -uncor$coefficients[["att1"]]/uncor$coefficients[["att3"]],
      -uncor$coefficients[["att2"]]/uncor$coefficients[["att3"]],
      
      -corrected$coefficients[["att1"]]/corrected$coefficients[["att3"]],
      -corrected$coefficients[["att2"]]/corrected$coefficients[["att3"]])
    
    
    
  }
}
toc()