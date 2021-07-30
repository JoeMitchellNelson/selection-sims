######## create sample data ########

require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,doParallel,tictoc,patchwork,matrixcalc,survival,plotly)

n.cores <- parallel::detectCores() - 4

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
clusterCall(my.cluster, function() require(pacman))
clusterCall(my.cluster, function() p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,doParallel,survival))

# parameters initialized in sim2

sim2 <- function(covs) {
  
  cov1 <- covs[1]
  cov2 <- covs[2]
  covz <- covs[3]
  
  n <- 2000
  m <- 100
  truebeta <- c(4,-8,0,0)
  vcov <- matrix(c(      1, 0,.25,cov1,
                        NA, 1,-.5,cov2,
                        NA,NA,  1,covz,
                        NA,NA,NA, 1),nrow=4,ncol=4,byrow=T)
  
  
  
  vcov[lower.tri(vcov)] <- t(vcov)[lower.tri(vcov)]
  if (!is.positive.definite(vcov)) {stop("\n\n\n\t\t\t!!!! VCOV NOT PD !!!!\n\n\n")}
  
  # utility function for RUM
  u <- function (b1,b2,x1,x2) {b1*x1 + b2*x2}
  
  ###### run sim #####
  {
    results <- foreach (i = 1:m,
                        .combine = "rbind"
    ) %dopar% 
      {
        subjects <- mvrnorm(n=n,mu=truebeta,Sigma=vcov) %>% as.data.frame()
        names(subjects) <- c("beta1","beta2","z","z2")
        subjects$id <- 1:n
        
        choices <- data.frame(id=1:n,att1=rnorm(n),att2=rnorm(n),eps=rgumbel(n))
        
        dat <- left_join(choices,subjects,by="id")
        
        
        prs <- 1*dat$z # sum of selection factors
        dat$pr <- pnorm(prs,mean=mean(prs),sd=sd(prs)) # convert sum to [0,1] probability
        dat$select <- rbinom(nrow(dat),1,prob=dat$pr) # respond with calculated probability
        
        stage1 <- glm(select ~  z2,data=dat, family = "binomial")
        
        dat$RP <- predict(stage1,dat)
        dat$RP <- (dat$RP - mean(dat$RP))/sd(dat$RP)
        #ggplot(dat) + geom_point(aes(x=RP,y=pr))
        
        d2 <- dat %>% mutate(att1=0,att2=0,eps=rgumbel(n))
        #d2$best <- (dat$best - 1)^2
        
        dat <- bind_rows(dat,d2)
        
        dat <- dat[order(dat$id),]
        dat <- dat %>% mutate(utility = u(beta1,beta2,att1,att2) + eps)
        dat <- dat %>% group_by(id) %>% mutate(best = ifelse(utility==max(utility),1,0)) %>% ungroup()
        # head(dat)
        
        uncor <- clogit(best ~ att1 + att2 +  strata(id),  
                        data = dat[which(dat$select==1),])
        
        
        # fullsamp <- clogit(best ~ att1 + att2 +  strata(id),
        #                    data = dat[which(dat$id %in% sample(n,n/2)),])
        
        
        corrected <- clogit(best ~ att1 + att2 +  att1:RP + att2:RP  + strata(id), 
                            data = dat[which(dat$select==1),])
        
        
        c(-uncor$coefficients[["att1"]]/uncor$coefficients[["att2"]],
          
          -corrected$coefficients[["att1"]]/corrected$coefficients[["att2"]])
        
        
        
      }
  }
  
  c(mean(results[,1]),mean(results[,2]))
  
}

sim3 <- function (covs) {tryCatch({sim2(covs)},error=function(cond) {return(c(NA,NA))})}

gridsearch <- expand.grid(cov1=seq(-1,1,by=.1),cov2=seq(-1,1,by=.1),covz=seq(-1,1,by=.1))

# remove vcovs that are not PD
{
  to_delete <- c()
  
  for (i in 1:nrow(gridsearch)) {
    cov1 <- gridsearch$cov1[i]
    cov2 <- gridsearch$cov2[i]
    covz <- gridsearch$covz[i]
    
    vcov <- matrix(c(     1, 0,.25,cov1,
                          NA, 1,-.5,cov2,
                          NA,NA, 1,covz,
                          NA,NA,NA, 1),nrow=4,ncol=4,byrow=T)
    
    
    
    vcov[lower.tri(vcov)] <- t(vcov)[lower.tri(vcov)]
    if (!is.positive.definite(vcov)) {to_delete <- c(to_delete,i)}
  }
  
  gridsearch <- gridsearch[-to_delete,]
  
}

system.time({
  bigresults <- apply(gridsearch,1,sim3) %>% t() %>% as.data.frame()
})



gridsearch <- cbind(gridsearch,bigresults)

gridsearch$bias_remain <- abs(.5-gridsearch$V2)/abs(.5-gridsearch$V1)

gridsearch$bias_remain_truncated <- ifelse(gridsearch$bias_remain >= 3,3,gridsearch$bias_remain)

gridsearch$improve <- ifelse(gridsearch$bias_remain>=1,"Worse","Better")
saveRDS(gridsearch, file = "sim_results_no_causality_.25_neg.5.rds")

###### plot results #####

#gridsearch <- readRDS("sim_results_no_causality.rds")
#gridsearch <- gridsearch %>% dplyr::filter(bias_remain <= 1)

fig <- plot_ly(gridsearch[which(gridsearch$bias_remain<=1),], x = ~cov1, y = ~cov2, z = ~covz,color=~bias_remain)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Cov(beta1,instrument)'),
                                   yaxis = list(title = 'Cov(beta2,instrument)'),
                                   zaxis = list(title = 'Cov(Z,instrument)')))

fig

