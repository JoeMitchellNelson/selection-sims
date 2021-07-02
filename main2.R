######## create sample data ########

require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,doParallel,tictoc,patchwork,matrixcalc,survival)

n.cores <- parallel::detectCores() - 1

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
clusterCall(my.cluster, function() require(pacman))
clusterCall(my.cluster, function() p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,doParallel,survival))

# initialize some parameters
n <- 10000
m <- 100

truebeta <- c(4,12,-8,0,0,0)
#beta1, beta2, beta3, z, z2, z3
#betas as params in utility function
#z is unobserved factor in selection
#z2, z3 are observed factors in selection

vcovgood <- matrix(c(   1, 0, 0,.3,.3,.3,
                       NA, 1, 0,.3,.3,.3,
                       NA,NA, 1,.2,.2,.2,
                       NA,NA,NA, 1,.5,.5,
                       NA,NA,NA,NA, 1,.0,
                       NA,NA,NA,NA,NA, 1),nrow=6,ncol=6,byrow=T)

vcovbad <- matrix(c(   1, 0, 0,.3,.0,.0,
                      NA, 1, 0,.3,.0,.0,
                      NA,NA, 1,.2,.0,.0,
                      NA,NA,NA, 1,.0,.0,
                      NA,NA,NA,NA, 1,.0,
                      NA,NA,NA,NA,NA, 1),nrow=6,ncol=6,byrow=T)

vcov <- vcovgood

vcov[lower.tri(vcov)] <- t(vcov)[lower.tri(vcov)]
if (!is.positive.definite(vcov)) {cat("\n\n\n\t\t\t!!!! VCOV NOT PD !!!!\n\n\n")}

# utility function for RUM
u <- function (b1,b2,b3,x1,x2,x3) {b1*x1 + b2*x2 + b3*x3}

###### run sim #####
set.seed(666)
tic()
{
  results <- foreach (i = 1:m,
                     .combine = "rbind"
                      ) %dopar% 
    {

      subjects <- mvrnorm(n=n,mu=truebeta,Sigma=vcov) %>% as.data.frame()
      names(subjects) <- c("beta1","beta2","beta3","z","z2","z3")
      subjects$id <- 1:n
      
      choices <- data.frame(id=1:n,att1=rnorm(n),att2=rnorm(n),att3=rnorm(n),eps=rgumbel(n))
      
      dat <- left_join(choices,subjects,by="id")
      
      dat <- dat %>% mutate(utility = u(beta1,beta2,beta3,att1,att2,att3) + eps)
      
      dat$best <- ifelse(dat$utility > 0, 1, 0)
      prs <- 2*dat$z + 1*dat$z2 + 1*dat$z3 # sum of selection factors
      dat$pr <- pnorm(prs,mean=mean(prs),sd=sd(prs)) # convert sum to [0,1] probability
      #dat$select <- ifelse(dat$z + dat$z2 + dat$z3 > 0,1,0)
      dat$select <- rbinom(nrow(dat),1,prob=dat$pr) # respond with calculated probability
      
      stage1 <- glm(select ~  z2,data=dat, family = "binomial")
      
      dat$RP <- predict(stage1,dat)
      dat$RP <- (dat$RP - mean(dat$RP))/sd(dat$RP)
      #ggplot(dat) + geom_point(aes(x=RP,y=pr))
      
      d2 <- dat %>% mutate(att1=0,att2=0,att3=0,eps=rgumbel(n))
      d2$best <- (dat$best - 1)^2
      
      dat <- rbind(dat,d2)
      
      dat <- dat[order(dat$id),]
      dat$utility <- ifelse(dat$att1==0,dat$eps,dat$utility)
      dat <- dat %>% group_by(id) %>% mutate(best = ifelse(utility==max(utility),1,0)) %>% ungroup()
      head(dat)
      
      uncor <- clogit(best ~ att1 + att2 + att3 + strata(id),  
                   data = dat[which(dat$select==1),])
      
      
      fullsamp <- clogit(best ~ att1 + att2 + att3 + strata(id),
                      data = dat[1:sum(dat$select),])
      
      
      corrected <- clogit(best ~ att1 + att2 + att3 + att1:RP + att2:RP + att3:RP + strata(id), 
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


