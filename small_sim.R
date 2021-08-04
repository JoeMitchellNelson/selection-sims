require(pacman)
p_load(tidyverse,evd,MASS)

n = 2000
m = 200

means <- c()
medians <- c()
trues <- c()



for (i in 1:m) {



df <- data.frame(ID = rep(1:n,each=3),                     
                       best=NA,                                
                       beta_cost = -2,
                       beta_x = rep(rnorm(n,mean=6,sd=1)),
                       x=rep(runif(n,min=0,max=3)),       
                       cost=rep(runif(n,min=0,max=3)),   
                       nu = rep(rgumbel(n*3,loc=0,scale=15)))

truewtp <- mean(-df$beta_x/df$beta_cost)

df$u <- df$beta_cost*df$cost + df$beta_x*df$x + df$nu

df <- df %>% group_by(ID) %>% mutate(maxu = max(u)) %>% ungroup()
df$best <- ifelse(df$maxu == df$u,1,0)
summary(res <- clogit(best ~ x + cost + strata(ID),data=df))

getwtp <- mvrnorm(n=10000,mu=res$coefficients,Sigma=res$var) %>% as.data.frame()
getwtp$wtp <- -getwtp$x/getwtp$cost

means <- c(means,mean(getwtp$wtp))
medians <- c(medians,median(getwtp$wtp))
trues <- c(trues,truewtp)

}

summary(means)
summary(medians)
summary(trues)

simres <- data.frame(
  stat=rep(c("mean_WTP","median_WTP","true_WTP"),each=m),
  value=c(means,medians,trues)
  )

ggplot(simres[simres$stat %in% c("median_WTP","mean_WTP"),]) +
  geom_density(aes(x=value,group=stat,fill=stat),alpha=.5) +
  lims(x=c(-10,10))
