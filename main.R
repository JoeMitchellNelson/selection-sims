######## create sample data ########

require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection)

# initialize some parameters
n <- 10000

# utility function for RUM

u <- function (b1,b2,b3,x1,x2,x3) {
  b1*x1 + b2*x2 + b3*x3
}

{

truebeta <- c(3,9,-6)

vcov <- matrix(c(1,.5,.4,.4,
                 .5,1,.7,.3,
                 .4,.7,1,.6,
                 .4,.3,.6,1),nrow=4,ncol=4,byrow=T)


subjects <- mvrnorm(n=n,mu=c(truebeta,0),Sigma=vcov) %>% as.data.frame()
names(subjects) <- c("beta1","beta2","beta3","z")
subjects$id <- 1:n

choices <- data.frame(id=1:n,att1=rnorm(n),att2=rnorm(n),att3=rnorm(n),eps=rlogis(n))

dat <- left_join(choices,subjects)

dat <- dat %>% mutate(utility = u(beta1,beta2,beta3,att1,att2,att3) + eps)
head(dat)

dat$best <- ifelse(dat$utility > 0, 1, 0)
dat$select <- ifelse(dat$z > 0,1,0)
dat$z2 <- dat$z + rnorm(n,0,.5)

stage1 <- glm(select ~ z2,data=dat, binomial(link = "probit"))

summary(stage1)

dat$RP <- predict(stage1,dat)
dat$RP <- dat$RP - mean(dat$RP)

uncor <- glm(best ~ att1 + att2 + att3, family = "binomial", 
         data = dat[which(dat$select==1),])


fullsamp <- glm(best ~ att1 + att2 + att3, family = "binomial", 
    data = dat)


corrected <- glm(best ~ att1 + att2 + att3 + att1:RP + att2:RP + att3:RP, family = "binomial", 
                 data = dat[which(dat$select==1),])

cat("Average mWTP in full sample for attribute 1:", 
    -fullsamp$coefficients[["att1"]]/fullsamp$coefficients[["att3"]],
    "(should be ",-mean(dat$beta1/dat$beta3),")\n")

cat("Average mWTP in full sample for attribute 2:", 
    -fullsamp$coefficients[["att2"]]/fullsamp$coefficients[["att3"]],
    "(should be ",-mean(dat$beta2/dat$beta3),")\n")

cat("Average mWTP in observed sample for attribute 1 (uncorrected):", 
    -uncor$coefficients[["att1"]]/uncor$coefficients[["att3"]],
    "(should be ",-mean(dat$beta1/dat$beta3),")\n")

cat("Average mWTP in observed sample for attribute 2 (uncorrected):", 
    -uncor$coefficients[["att2"]]/uncor$coefficients[["att3"]],
    "(should be ",-mean(dat$beta2/dat$beta3),")\n")


cat("Average mWTP in observed sample for attribute 1 (corrected):", 
    -corrected$coefficients[["att1"]]/corrected$coefficients[["att3"]],
    "(should be ",-mean(dat$beta1/dat$beta3),")\n")

cat("Average mWTP in observed sample for attribute 2 (corrected):", 
    -corrected$coefficients[["att2"]]/corrected$coefficients[["att3"]],
    "(should be ",-mean(dat$beta2/dat$beta3),")\n")



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


