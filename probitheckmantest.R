require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,
       doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo,haven,readstata13)

# heckman probit
set.seed(0)
library("sampleSelection")
library("mvtnorm")

n <- 1000

eps <- rmvnorm(n, c(0,0), matrix(c(1,-0.7,
                                     -0.7,1), 
                                   2, 2))
xs <- round(runif(n))
ys <- xs + eps[,1] > 0
xo <- runif(n)
yoX <- xo + eps[,2]
yo <- (yoX>0)*(ys > 0)
eps


summary( selection(ys~xs, factor(yo) ~1,method="ml"))


df <- data.frame(xs = xs,
                 ys = ys,
                 xo = xo,
                 yo = yo)

df$yo <- ifelse(!df$ys,NA,df$yo)

summary( selection(ys~factor(xs), factor(yo) ~ 1,data=df,method="ml"))
summary(naive <- probit(factor(yo) ~ 1,data=df))

pnorm(naive$estimate)


### CCO factor variables for selection equation ####

n <- 1000
iter <- 500

results <- data.frame(true=rep(NA,iter),wrong=NA,corrected=NA,se=NA,r_rate=NA)

for (i in 100:iter) {

eps <- rmvnorm(n*3, c(-0.5,0), matrix(c(1,0.7,
                                     0.7,1), 
                                   2, 2))

df <- data.frame(CCO = factor(x=rep(c("A","B","C"),n)),
                 ys = NA,
                 yo = NA,
                 eps_s = eps[,1],
                 eps_o = eps[,2])

df <- df %>% mutate(ys = (-0.2*(CCO=="A") + 0*(CCO=="B")+ 0.5*(CCO=="C") + eps_s)>0) %>% 
  mutate(yo = eps_o+0.4>0)

est1 <- mean(df$yo)
responserate <- paste0(round(100*mean(df$ys)))

nocorrect <- mean(df$yo[df$ys])

summary(good <- selection(ys~CCO, yo ~ 1,data=df,method="ml"))
good2 <- summary(good)
good3 <- good2$estimate %>% as.data.frame()
se <- good3$`Std. Error`[4]
corrected <- pnorm(good$estimate[4])
ub <- pnorm(corrected+1.96*se)
lb <- pnorm(corrected-1.96*se)

results[i,] <- c(est1,nocorrect,corrected,se,mean(df$ys))

}


ggplot(results) +
  geom_histogram(aes(x=wrong),fill="red",color="red", alpha=0.5,binwidth=0.01) +
  geom_histogram(aes(x=corrected),fill="forestgreen",color="forestgreen",alpha=0.5,binwidth=0.025) +
  geom_vline(aes(xintercept=mean(true,na.rm=T))) +
  labs(x="Estimates",y="Count")
  #geom_vline(aes(xintercept=mean(corrected,na.rm=T)),color="darkgreen") + 

meanA <- mean(df$yo[df$ys & df$CCO=="A"])
meanB <- mean(df$yo[df$ys & df$CCO=="B"])
meanC <- mean(df$yo[df$ys & df$CCO=="C"])
weighted <- mean(c(meanA, meanB, meanC))

cat("Full participation:    ",round(100*est1),"%",
    "\nNo correction:         ",round(100*nocorrect),
    "\nCorrected (weighting): ",round(weighted*100),
    "\nCorrected (Heckman):   ",round(corrected*100),
    "\nwith response rate of  ",responserate,"\n")

cat("CI for heckman:","(",floor(100*lb),",",floor(100*ub),")")
