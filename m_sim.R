require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,doParallel,tictoc,patchwork,matrixcalc,survival,plotly)

# choices (1 choice per person)
n <- 100000

# alternatives per choice
J <- 3

# parameters
beta10 <- -2
beta20 <- 4
beta11 <- 1
beta21 <- 1

df <- data.frame(alt = rep(1:J,n),
                 id = rep(1:n,each=J),
                 m = rep(rnorm(n=n),each=J),
                 nu = rgumbel(n*J),
                 cost = runif(n=n*J,min=0,max=3),
                 x = runif(n=n*J,min=0,max=4))

# make alternative J the statquo
df$x <- ifelse(df$alt==J,0,df$x)
df$cost <- ifelse(df$alt==J,0,df$cost)

# calculate utility, choose best
df <- df %>% mutate(u = (beta10)*cost + (beta20 + beta21*m)*x + 2*nu)
df <- df %>% group_by(id) %>%  mutate(best = ifelse(u==max(u),1,0)) %>% ungroup()
head(df)
# naive, uncorrected WTP
res1 <- clogit(best ~ cost + x + strata(id),data=df)
res2 <- clogit(best ~ cost + x*m + strata(id),data=df)

WTP_dumb <- -res1$coefficients[["x"]]/res1$coefficients[["cost"]]
WTP_dumb2 <- -res2$coefficients[["x"]]/res2$coefficients[["cost"]]

WTP <- mvrnorm(n=1000,mu=res2$coefficients,Sigma=res2$var) %>% as.data.frame()
WTP$wtp <- -WTP$x/WTP$cost
mean(WTP$wtp)
quantile(WTP$wtp,c(0.025,0.975))
plot(density(WTP$wtp))

mean(-(beta20 + beta21*df$m)/(beta10+beta11*0))

