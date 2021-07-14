require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,doParallel,tictoc,patchwork,matrixcalc,survival,plotly)

set.seed(301)

n <- 100000 # choices (1 choice per person)

J <- 3 # alternatives per choice

# set up data
df <- data.frame(alt = rep(1:J,n),                  # alternative index
                 id = rep(1:n,each=J),              # choice/person index (one choice per person)
                 m = rep(rnorm(n=n),each=J),        # the source of all the trouble
                 nu = rgumbel(n*J),                 # Type I extreme value error for utility
                 cost = runif(n=n*J,min=0,max=3),   # cost of each alternative (bad)
                 x = runif(n=n*J,min=0,max=4),      # some attribute of each choice alternative (good)
                 omegastar= rep(rnorm(n=n),each=J), # error term in selection model
                 w = rep(rnorm(n=n),each=J)         # observable characteristic that affects selection
)

# parameters
beta1 <- -2 # effect for cost (choice model)
beta20 <- 4  # common component of effect for x
beta21 <- .5  # m-component of effect for x

alpha0 <- 0 # constant term, increase to boost the simulated response rate
alpha1 <- 0 # unused coef for w in selection equation
alpha2 <- 3 # coef for m in selection equation



# constructed variables for selection model
df$Z <- alpha0 + alpha1*df$w + alpha2*df$m + df$omegastar
df$respond <- ifelse(df$Z > 0,1,0)

# make alternative J the statquo, (cost = 0, x = 0)
df$x <- ifelse(df$alt==J,0,df$x)
df$cost <- ifelse(df$alt==J,0,df$cost)

# calculate utility, choose best
df <- df %>% mutate(u = (beta1)*cost + (beta20 + beta21*m)*x + nu)
df <- df %>% group_by(id) %>%  mutate(best = ifelse(u==max(u),1,0)) %>% ungroup()

mean(-(beta20 + beta21*df$m)/beta1)
# naive, uncorrected WTP with 100% response rate
res1 <- clogit(best ~ cost + x + strata(id),data=df)
WTP_1 <- -res1$coefficients[["x"]]/res1$coefficients[["cost"]]

# uncorrected WTP estimate with partial response
res2 <- clogit(best ~ cost + x + strata(id),data=df[which(df$respond==1),])
WTP_2 <- -res2$coefficients[["x"]]/res2$coefficients[["cost"]]

